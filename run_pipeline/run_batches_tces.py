"""Script created to run LEO Vetter in batch mode for SPOC TCEs.

Requires a CSV file containing TCEs with the following columns:
- target_id (TIC ID), int
- uid (unique identifier for the SPOC TCE. E.g., using format <tic_id>_<planet_candidate>_<sector_run>), string
- sector_run (identifier for the sector run, e.g. "1" for sector run S1, "1-6" for sector run S1-S6), string
- tce_time0bk (TCE epoch in BTJD), float
- tce_period (TCE period in days), float
- tce_duration (TCE duration in hours), float
- tce_plnt_num (TCE planet candidate number, e.g. 1 for the first planet candidate, 2 for the second, etc.), int
- sectors_observed (string indicating which sectors the TIC was observed in for the specific sector run, either as a binary string like "100100" where each digit represents a sector, or as an explicit list of sector numbers like "1_4_27"), string
"""

# Suppress expected RuntimeWarnings from edge cases (divide by zero, empty arrays, etc)
# This must be set before any leo_vetter imports to apply globally
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# imports
import argparse
import pandas as pd
from pathlib import Path
import lightkurve as lk
import numpy as np
from astroquery.mast import Catalogs
from multiprocessing import Pool
from tqdm import tqdm
import logging
import time

from requests.exceptions import RequestException
from urllib3.exceptions import HTTPError

from leo_vetter.stellar import quadratic_ldc
from leo_vetter.main import TCELightCurve
from leo_vetter.plots import plot_modshift, plot_summary
from leo_vetter.thresholds import check_thresholds


TIC_COLUMNS = ["rad", "mass", "rho", "Teff", "logg"]
TCE_COLUMNS = ["target_id", "uid", "sector_run", "tce_time0bk", "tce_period", "tce_duration", "tce_plnt_num", "sector_run", "sectors_observed"]
LC_SOURCE_OPTIONS = ("2min", "ffi")
REMOTE_RETRY_ATTEMPTS = 4
REMOTE_RETRY_BASE_DELAY_SECONDS = 2.0


def get_logger():
    """Get the logger for the pipeline."""
    return logging.getLogger("leo_vetter_pipeline")


def is_retryable_remote_error(error):
    """Determine if an error is a retryable remote error (e.g., network issues, timeouts).
    
    :param Exception error: the exception to check
    :return bool: True if the error is a retryable remote error, False otherwise
    """
    
    return isinstance(
        error,
        (
            RequestException,
            HTTPError,
            ConnectionError,
            TimeoutError,
            OSError,
        ),
    )


def retry_remote_call(operation, description, attempts=REMOTE_RETRY_ATTEMPTS):
    """Retry a remote operation with exponential backoff if it fails due to a retryable remote error.
    
    :param function operation: the remote operation to perform, as a function that takes no arguments
    :param str description: a description of the operation for logging purposes
    :param int attempts: number of retry attempts, defaults to REMOTE_RETRY_ATTEMPTS
    :return: result of the remote operation
    :raises Exception: the last exception raised if all retry attempts fail
    """
    
    logger = get_logger()

    for attempt in range(1, attempts + 1):
        try:
            return operation()
        except Exception as error:
            if not is_retryable_remote_error(error) or attempt == attempts:
                raise

            delay_seconds = REMOTE_RETRY_BASE_DELAY_SECONDS * attempt
            logger.warning(
                "%s failed on attempt %s/%s: %s. Retrying in %.1fs",
                description,
                attempt,
                attempts,
                error,
                delay_seconds,
            )
            time.sleep(delay_seconds)


def convert_sectors_observed_binary_string_to_int_string(sectors_observed_binary_string):
    """Convert `sectors_observed` from a binary string to a string with sectors numbers separated by "_".

    :param str sectors_observed_binary_string: sectors observed in binary format
    :raises ValueError: `sectors_observed` is NaN
    :raises ValueError: `sectors_observed` is empty
    :return str: `sectors_observed` in new format
    """

    if pd.isna(sectors_observed_binary_string):
        raise ValueError("sectors_observed_binary_string cannot be NaN")

    sectors_observed_binary_string = str(sectors_observed_binary_string).strip()
    if not sectors_observed_binary_string:
        raise ValueError("sectors_observed_binary_string cannot be empty")

    # Already in explicit format like "1_2_27".
    if '_' in sectors_observed_binary_string:
        return '_'.join(
            sector for sector in sectors_observed_binary_string.split('_') if sector
        )

    # Convert long binary mask (e.g., "100100...") to 1-indexed sector list.
    if set(sectors_observed_binary_string).issubset({'0', '1'}) and len(sectors_observed_binary_string) > 3:
        sectors = [
            str(sector_i)
            for sector_i, sector_bin_i in enumerate(sectors_observed_binary_string)
            if sector_bin_i == '1'
        ]
        return '_'.join(sectors)

    return sectors_observed_binary_string


def setup_logging(log_dir):
    """Configure logging to capture info and error messages to a file.
    
    :param Path log_dir: directory to save log files
    :return logging.Logger: configured logger instance
    """
    
    log_dir.mkdir(exist_ok=True, parents=True)
    log_file = log_dir / "pipeline.log"
    
    # Create logger
    logger = logging.getLogger("leo_vetter_pipeline")
    logger.setLevel(logging.DEBUG)
    
    # Clear any existing handlers
    logger.handlers.clear()
    
    # File handler (captures everything)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(file_formatter)
    logger.addHandler(file_handler)
    
    # Console handler (only INFO and above)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Redirect warnings to logger
    logging.captureWarnings(True)
    warnings_logger = logging.getLogger("py.warnings")
    warnings_logger.addHandler(file_handler)
    
    logger.info("Logging initialized - RuntimeWarnings suppressed")
    
    return logger


def get_cached_lc_files(tic, sector_numbers, save_lc_dir, lc_source):
    """Gets cached light curve FITS files for TIC ID `tic` in sectors `sector_numbers` found in `save_lc_dir`.

    :param int tic: TIC ID
    :param list sector_numbers: sector IDs
    :param str save_lc_dir: light curve directory
    :param str lc_source: either "2min" or "ffi" for SPOC 2-min/FFI light curves, respectively
    :raises ValueError: light curve source not supported/known
    :return set: sorted set of Paths for the cached light curves for the target in the given sectors
    """

    tic_pattern = f"{int(tic):016d}"
    sector_tokens = {f"s{sector:04d}" for sector in sector_numbers}

    if lc_source == "2min":
        local_lc_files = []
        for sector_token in sector_tokens:
            local_lc_files.extend(
                save_lc_dir.rglob(f"tess*-{sector_token}-{tic_pattern}-*_lc.fits")
            )
        return sorted(set(local_lc_files))

    if lc_source == "ffi":
        local_lc_files = []
        for sector_token in sector_tokens:
            local_lc_files.extend(
                save_lc_dir.rglob(f"hlsp_tess-spoc_tess_phot_{tic_pattern}-{sector_token}_tess_v1_tp.fits")
            )
        return sorted(set(local_lc_files))

    raise ValueError(f"Unsupported lc_source: {lc_source}")


def get_lc_data(tic, sectors_observed, save_lc_dir, lc_source):
    """Gets light curve data for TIC ID `tic` in sectors `sectors_observed`.
    
    :param int tic: TIC ID
    :param str sectors_observed: sectors observed, separated by "_"
    :param str save_lc_dir: light curve directory
    :param str lc_source: either "2min" or "ffi" for SPOC 2-min/FFI light curves, respectively
    :return tuple: light curve object and list of local light curve files
    """
    
    sectors_numbers = [int(sector) for sector in sectors_observed.split('_')]

    local_lc_files = get_cached_lc_files(tic, sectors_numbers, save_lc_dir, lc_source)

    if local_lc_files:
        lcs = lk.LightCurveCollection([lk.read(local_lc_file) for local_lc_file in local_lc_files])
    else:
        author = "SPOC" if lc_source == "2min" else "TESS-SPOC"
        search_result = retry_remote_call(
            lambda: lk.search_lightcurve(
                f"TIC {tic}",
                mission="TESS",
                author=author,
                sector=sectors_numbers,
            ),
            f"Light curve search for TIC {tic} sectors {sectors_observed}",
        )
        lcs = retry_remote_call(
            lambda: search_result.download_all(download_dir=str(save_lc_dir)),
            f"Light curve download for TIC {tic} sectors {sectors_observed}",
        )
        local_lc_files = get_cached_lc_files(tic, sectors_numbers, save_lc_dir, lc_source)

    if lcs is None or len(lcs) == 0:
        raise FileNotFoundError(
            f"No light curves available for TIC {tic} sectors {sectors_observed}"
        )
    
    # Stitch light curves together into a multi-sector light curve
    lc = lcs.stitch()

    # Remove NaNs and poor quality cadences
    lc = lc[~np.isnan(lc["flux"]) & (lc["quality"] == 0)]

    return lc, local_lc_files


def get_lc_data_for_tce(lc, epo, per, dur):
    """Extracts and flattens the light curve data for a given TCE.
    
    :param lk.LightCurve lc: light curve object
    :param float epo: epoch of the transit (BTJD)
    :param float per: period of the transit (day)
    :param float dur: duration of the transit (day)
    :return tuple: arrays of time, raw flux, flattened flux, and flux errors
    """
    
    # Flatten the light curve
    # Highly recommend that the transits are masked out when flattening the light curve
    transit_mask = lc.create_transit_mask(transit_time=epo, period=per, duration=dur)
    lc_flat = lc.flatten(mask=transit_mask)

    # Extract the relevant arrays
    time = np.asarray(lc_flat["time"].value)
    raw = np.asarray(lc["flux"].value)
    flux = np.asarray(lc_flat["flux"].value)
    flux_err = np.asarray(lc_flat["flux_err"].value)
    
    return time, raw, flux, flux_err

def query_tic_for_stellar_parameters(tic):
    """Queries the TIC catalog for stellar parameters for a given TIC ID and returns them in a dictionary, along with limb-darkening coefficients.
    
    :param int tic: TIC ID
    :return dict: dictionary containing stellar parameters and limb-darkening coefficients for the target TIC
    """

    def _to_finite_float(value, fallback):
        try:
            cast_value = float(value)
        except (TypeError, ValueError):
            return fallback
        return cast_value if np.isfinite(cast_value) else fallback

    result = retry_remote_call(
        lambda: Catalogs.query_criteria(catalog="TIC", ID=tic),
        f"TIC catalog query for TIC {tic}",
    )
    star = {"tic": tic}

    # Use solar-like defaults when TIC metadata are missing/non-finite.
    stellar_defaults = {
        "rad": 1.0,
        "mass": 1.0,
        "rho": 1.0,
        "Teff": 5777.0,
        "logg": 4.44,
    }

    for key in TIC_COLUMNS:
        star[key] = _to_finite_float(result[key], stellar_defaults[key])
        star["e_" + key] = _to_finite_float(result["e_" + key], 0.0)

    # Get limb-darkening parameters from sanitized Teff/logg values.
    star["u1"], star["u2"] = quadratic_ldc(star["Teff"], star["logg"])

    return star

def generate_tce_metrics(tic_id, per, epo, dur, lc_tic, tic_params, metrics_save_fp, planetno=1):
    """Generates metrics for the TCE.

    :param int tic_id: TIC ID
    :param float per: period (day)
    :param float epo: epoch (BTJD)
    :param float dur: transit duration (day)
    :param lk.Lightcurve lc_tic: target light curve object
    :param dict tic_params: target stellar parameters
    :param Path metrics_save_fp: filepath used to save metrics CSV file
    :param int planetno: SPOC planet number, defaults to 1
    :return TCELightCurve: TCE object
    """
    
    time, raw, flux, flux_err = get_lc_data_for_tce(lc_tic, epo, per, dur)
        
    tlc = TCELightCurve(tic_id, time, raw, flux, flux_err, per, epo, dur, planetno=planetno)
        
    tlc.compute_flux_metrics(tic_params, verbose=True)
    
    tlc.save_metrics(save_file=metrics_save_fp)
    
    return tlc

def check_thresholds_tce(tlc, decision_thresholds, tce_uid, verbose=False):
    """Checks the metrics for a given TCE against the decision thresholds to determine FA/FP classification and which tests were failed.
    
    :param TCELightCurve tlc: TCE object containing the metrics to check
    :param dict decision_thresholds: dictionary containing the decision thresholds for each test
    :param str tce_uid: unique identifier for the TCE
    :param bool verbose: whether to print verbose output, defaults to False
    :return pd.DataFrame: DataFrame containing FA/FP classification and failed tests
    """
    
    # FA is True if any tests failed; False otherwise
    FA, FA_failed_tests = check_thresholds(tlc.metrics, "FA", thresholds=decision_thresholds, verbose=True) 
    # FP is True if any tests failed; False otherwise
    FP, FP_failed_tests = check_thresholds(tlc.metrics, "FP", thresholds=decision_thresholds, verbose=True)
    
    failed_tests = '_'.join(FA_failed_tests + FP_failed_tests) if (FA_failed_tests or FP_failed_tests) else 'None'
    
    fa_fp_tests_df = pd.DataFrame({
        'FA': [FA],
        'FP': [FP],
        'Failed Tests': [failed_tests],
    }, index=[tce_uid])

    if not FA and not FP and verbose:
        print(f"TIC-{tlc.tic}.{tlc.planetno} is a planet candidate!")   

    return fa_fp_tests_df

def process_tic(tic_id, tic_data, decision_thresholds, save_lc_dir, lc_source, delete_lc_after_target=False, plot_modshift_flag=False, plot_summary_flag=False, plot_modshift_save_dir=None, plot_summary_save_dir=None, metrics_save_dir=None, fa_fp_tests_save_dir=None):
    """Processes a single TIC through the pipeline, including light curve retrieval, metric generation, FA/FP classification, and plotting.
    
    :param int tic_id: TIC ID
    :param pd.DataFrame tic_data: DataFrame containing TCE data for the TIC, indexed by TCE UID
    :param dict decision_thresholds: dictionary containing the decision thresholds for each test
    :param Path save_lc_dir: directory to save light curves
    :param str lc_source: source of the light curves. Either "2min" for SPOC 2-min light curves or "ffi" for SPOC FFI light curves.
    :param bool delete_lc_after_target: whether to delete light curves after processing, defaults to False
    :param bool plot_modshift_flag: whether to plot modshift, defaults to False
    :param bool plot_summary_flag: whether to plot summary, defaults to False
    :param Path plot_modshift_save_dir: directory to save modshift plots
    :param Path plot_summary_save_dir: directory to save summary plots
    :param Path metrics_save_dir: directory to save metrics CSV files
    :param Path fa_fp_tests_save_dir: directory to save FA/FP tests CSV files
    :return dict: dictionary containing processing results for the TIC
    """
    
    logger = get_logger()

    try:
        tic_params = query_tic_for_stellar_parameters(tic_id)
    except Exception as error:
        logger.exception("Failed to fetch stellar parameters for TIC %s", tic_id)
        return {
            "tic_id": tic_id,
            "status": "failed",
            "processed_tces": 0,
            "failed_sector_runs": tic_data["sector_run"].nunique(),
            "error": str(error),
        }

    fa_fp_tests_df_tic_lst = []
    lc_files_to_cleanup = set()
    processed_tces = 0
    failed_sector_runs = []
    for sector_run, tic_data_sector_run in tqdm(tic_data.groupby('sector_run'), desc=f'Processing TIC {tic_id}', unit='sector run', total=tic_data["sector_run"].nunique()):

        sectors_observed = tic_data_sector_run["sectors_observed"].iloc[0]

        try:
            lc_tic, lc_files_used = get_lc_data(tic_id, sectors_observed, save_lc_dir, lc_source)
        except Exception as error:
            logger.exception(
                "Skipping TIC %s sector run %s after light-curve retrieval failure",
                tic_id,
                sector_run,
            )
            failed_sector_runs.append(
                {
                    "sector_run": sector_run,
                    "sectors_observed": sectors_observed,
                    "error": str(error),
                }
            )
            continue

        lc_files_to_cleanup.update(lc_files_used)

        for tce_uid, tce_data in tqdm(tic_data_sector_run.iterrows(), desc=f'Processing TIC {tic_id} sector run {sector_run}', unit='TCE', total=tic_data_sector_run.shape[0]):

            epo = tce_data["tce_time0bk"]
            per = tce_data["tce_period"]
            dur = tce_data["tce_duration"]
            planetno = tce_data["tce_plnt_num"]

            tlc = generate_tce_metrics(tic_id, per, epo, dur, lc_tic, tic_params, metrics_save_dir / f"metrics_tic{tic_id}_tce{tce_uid}.csv", planetno=planetno)

            fa_fp_tests_df_tce = check_thresholds_tce(tlc, decision_thresholds, tce_uid, verbose=True)
            fa_fp_tests_df_tic_lst.append(fa_fp_tests_df_tce)
            processed_tces += 1

            if plot_modshift_flag:
                plot_modshift(tlc, save_fig=plot_modshift_flag, save_file=plot_modshift_save_dir / f"modshift_tic{tic_id}_tce{tce_uid}.png")
            if plot_summary_flag:
                plot_summary(tlc, tic_params, save_fig=plot_summary_flag, save_file=plot_summary_save_dir / f"summary_tic{tic_id}_tce{tce_uid}.png")

    fa_fp_tests_df_tic_lst = [df for df in fa_fp_tests_df_tic_lst if df is not None]
    if fa_fp_tests_df_tic_lst:
        fa_fp_tests_df_tic_df = pd.concat(fa_fp_tests_df_tic_lst, axis=0)
        fa_fp_tests_df_tic_df.to_csv(fa_fp_tests_save_dir / f"fa_fp_tests_tic{tic_id}.csv", index=True)

    if delete_lc_after_target:
        for lc_file in sorted(lc_files_to_cleanup):
            try:
                lc_file.unlink(missing_ok=True)
            except OSError:
                pass

    status = "success" if not failed_sector_runs else "partial"
    return {
        "tic_id": tic_id,
        "status": status,
        "processed_tces": processed_tces,
        "failed_sector_runs": len(failed_sector_runs),
        "error": " | ".join(
            f"{item['sector_run']}: {item['error']}" for item in failed_sector_runs
        ),
    }


def read_tce_table(tce_tbl_fp):
    """Reads TCE table and prepares it for the run.

    :param Path tce_tbl_fp: filepath to TCE table
    :return pd.DataFrame: loaded TCE table
    """
    
    tce_tbl = pd.read_csv(tce_tbl_fp, usecols=TCE_COLUMNS, on_bad_lines='skip', engine='python', dtype={'sectors_observed': str})
    
    tce_tbl = tce_tbl.rename(columns={"target_id": "tic"})
    tce_tbl['tce_duration'] = tce_tbl['tce_duration'] / 24. # convert duration from hours to days
    
    tce_tbl['sectors_observed'] = tce_tbl['sectors_observed'].apply(convert_sectors_observed_binary_string_to_int_string)
    
    return tce_tbl

def run_pipeline(tce_tbl_fp, decision_thresholds, save_lc_dir, res_dir, lc_source="2min", delete_lc_after_target=False, plot_modshift_flag=False, plot_summary_flag=False, num_processes=4, additional_metadata=None):
    """Runs the LEO-vetter pipeline for a batch of TCEs specified in a TCE table CSV file.
    
    :param Path tce_tbl_fp: filepath to TCE table CSV file
    :param dict decision_thresholds: dictionary containing the decision thresholds for each test
    :param Path save_lc_dir: directory to save light curves
    :param Path res_dir: directory to save results (metrics, FA/FP tests, plots, logs)
    :param str lc_source: source of the light curves. Either "2min" for SPOC 2-min light curves or "ffi" for SPOC FFI light curves, defaults to "2min"
    :param bool delete_lc_after_target: whether to delete light curves after processing each target, defaults to False
    :param bool plot_modshift_flag: whether to generate and save modshift plots, defaults to False
    :param bool plot_summary_flag: whether to generate and save summary plots, defaults to False
    :param int num_processes: number of parallel processes to use for processing TICs, defaults to 4
    :param dict additional_metadata: optional dictionary of additional metadata to include in the saved decision thresholds CSV file
    """
    
    res_dir.mkdir(exist_ok=True)
    logger = setup_logging(res_dir / "logs")
    
    metrics_save_dir = res_dir / "metrics"
    fa_fp_tests_save_dir = res_dir / "fa_fp_tests"
    plot_modshift_save_dir = res_dir / 'modshift_plots'
    plot_summary_save_dir = res_dir / "summary_plots"
    
    metrics_save_dir.mkdir(exist_ok=True)
    fa_fp_tests_save_dir.mkdir(exist_ok=True)
    plot_modshift_save_dir.mkdir(exist_ok=True)
    plot_summary_save_dir.mkdir(exist_ok=True)
    
    save_lc_dir.mkdir(exist_ok=True, parents=True)
    
    # save decision thresholds dictionary to a CSV in the results directory for record-keeping
    decision_thresholds_df = pd.DataFrame.from_dict(decision_thresholds, orient='index', columns=['threshold'])
    # add metadata to the dataframe
    decision_thresholds_df.attrs['description'] = "Decision thresholds used for FA/FP classification in the LEO-vetter pipeline. If a TCE's metric value exceeds the threshold for a given test, it fails that test. FA is True if any tests failed; FP is True if any tests failed. These thresholds are applied to the metrics computed for each TCE to determine its FA/FP classification."
    decision_thresholds_df.attrs['source'] = "Defined in run_batches_tces.py and saved here for record-keeping."
    decision_thresholds_df.attrs['notes'] = "These thresholds can be adjusted based on the desired balance between false positives and false negatives. They were chosen based on analysis of known planets and false positives in TESS data, but may be further refined with additional data and testing."
    decision_thresholds_df.attrs['created'] = pd.Timestamp.now().isoformat()
    if additional_metadata:
        for key, value in additional_metadata.items():
            decision_thresholds_df.attrs[key] = value
    with open(res_dir / "decision_thresholds.csv", "w") as f:
        for key, value in decision_thresholds_df.attrs.items():
            f.write(f"{key}: {value}\n")
        decision_thresholds_df.to_csv(f, index=True)
    
    tce_tbl = read_tce_table(tce_tbl_fp)
    
    # Check for previously processed TCEs in aggregate results and skip them to avoid redundant processing
    agg_metrics_fp = res_dir / "agg_metrics.csv"
    agg_fa_fp_tests_fp = res_dir / "agg_fa_fp_tests.csv"
    if agg_metrics_fp.exists():
        agg_metrics_df = pd.read_csv(agg_metrics_fp, index_col=0)
        processed_uids_metrics = set(agg_metrics_df.index)
    if agg_fa_fp_tests_fp.exists():
        agg_fa_fp_tests_df = pd.read_csv(agg_fa_fp_tests_fp, index_col=0)
        processed_uids_fa_fp = set(agg_fa_fp_tests_df.index)
    processed_uids = set()
    if agg_metrics_fp.exists() and agg_fa_fp_tests_fp.exists():
        processed_uids = processed_uids_metrics.intersection(processed_uids_fa_fp)
        tce_tbl = tce_tbl[~tce_tbl['uid'].isin(processed_uids)]
        logger.info(f"Found {len(processed_uids)} previously processed TCEs in aggregate results. These will be skipped in this run.")
        logger.info(f"{len(tce_tbl)} TCEs remain to be processed after skipping previously processed TCEs.")
    else:
        logger.info("No previously processed TCEs found in aggregate results. All TCEs will be processed in this run.")
    
    # split the TCE table by TIC and run each TIC through the pipeline separately to avoid memory issues with loading in all the light curves at once
    tic_jobs = [(tic_id, tic_data.set_index("uid")) for tic_id, tic_data in tce_tbl.groupby("tic")]
    pipeline_results = []
    
    logger.info(
        f"Starting pipeline with {len(tic_jobs)} TICs, {num_processes} processes, "
        f"lc_source={lc_source}, delete_lc_after_target={delete_lc_after_target}"
    )
    
    with Pool(processes=num_processes) as pool, tqdm(total=len(tic_jobs), desc='Processing TICs', unit='TIC') as pbar:
        
        async_results = []
        for tic_job in tic_jobs:
            async_result = pool.apply_async(
                process_tic, 
                args=(*tic_job, decision_thresholds, save_lc_dir, lc_source, delete_lc_after_target, plot_modshift_flag, plot_summary_flag, plot_modshift_save_dir, plot_summary_save_dir, metrics_save_dir, fa_fp_tests_save_dir),
                callback=lambda _: pbar.update(),
            )
            async_results.append(async_result)
        
        for async_result in async_results:
            pipeline_results.append(async_result.get())

    pipeline_results_df = pd.DataFrame(pipeline_results)
    if not pipeline_results_df.empty:
        pipeline_results_df.sort_values(["status", "tic_id"], inplace=True)
        pipeline_results_df.to_csv(res_dir / "pipeline_status.csv", index=False)

        failed_results = pipeline_results_df[pipeline_results_df["status"] != "success"]
        if failed_results.empty:
            logger.info("Pipeline processing completed successfully")
        else:
            logger.warning(
                "Pipeline completed with %s TICs requiring attention. See %s",
                len(failed_results),
                res_dir / "pipeline_status.csv",
            )
    else:
        logger.info("Pipeline completed with no TICs to process")

def aggregate_metrics(metrics_dir, output_fp):
    """Aggregates individual TCE metrics CSV files into a single CSV file for all TCEs.
    
    :param Path metrics_dir: directory containing individual TCE metrics CSV files
    :param Path output_fp: filepath to save the aggregated metrics CSV file
    """
    
    all_metrics = []
    for metrics_file in metrics_dir.glob("metrics_tic*_tce*.csv"):
        df = pd.read_csv(metrics_file)
        df['uid'] = metrics_file.stem.split("tce")[1] # extract the TCE uid from the filename and add it as a column to the metrics dataframe
        all_metrics.append(df)
    
    all_metrics_df = pd.concat(all_metrics, ignore_index=True)
    all_metrics_df.sort_values(["tic", "uid"], inplace=True)
    all_metrics_df.set_index("uid", inplace=True)
    all_metrics_df.to_csv(output_fp, index=True)
    
    print(f"Aggregated metrics for {len(all_metrics_df)} TCEs saved to {output_fp}")
    

def aggregate_fa_fp_tests(fa_fp_tests_dir, output_fp):
    """Aggregates individual TCE FA/FP tests CSV files into a single CSV file for all TCEs.
    
    :param Path fa_fp_tests_dir: directory containing individual TCE FA/FP tests CSV files
    :param Path output_fp: filepath to save the aggregated FA/FP tests CSV file
    """
    
    all_fa_fp_tests = []
    for fa_fp_tests_file in fa_fp_tests_dir.glob("fa_fp_tests_tic*.csv"):
        df = pd.read_csv(fa_fp_tests_file, index_col=0)
        df['uid'] = df.index # add the uid as a column to the fa_fp_tests dataframe
        all_fa_fp_tests.append(df)
    
    all_fa_fp_tests_df = pd.concat(all_fa_fp_tests, ignore_index=True)
    all_fa_fp_tests_df.sort_values(["uid"], inplace=True)
    all_fa_fp_tests_df.set_index("uid", inplace=True)
    all_fa_fp_tests_df.to_csv(output_fp, index=True)
    
    print(f"Aggregated FA/FP tests for {len(all_fa_fp_tests_df)} TCEs saved to {output_fp}")
    
           
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Run LEO-vetter on TESS TCE tables.")
    parser.add_argument(
        "--lc-source",
        choices=LC_SOURCE_OPTIONS,
        default="2min",
        help="Choose which SPOC light-curve source to use for cache matching and downloads.",
    )
    parser.add_argument(
        "--tce-table",
        type=Path,
        default=Path('/data/exoplnt_dl/ephemeris_tables/tess/tess_spoc_2min/tess-spoc-2min-tces-dv_s1-s94_s1s92_9-19-2025_1518_exofop-sg1-tois_9-22-2025.csv'),
        help="Path to the TCE table CSV.",
    )
    parser.add_argument(
        "--delete-lc-after-target",
        action="store_true",
        help="Delete cached light-curve FITS files after all TCEs for the target TIC are processed.",
    )
    args = parser.parse_args()
    
    additional_metadata = {
        'TCEs catalog': 'TESS SPOC 2-min TCEs S1-S94 (up to S14-S78 multisector run)'
    }

    tce_tbl_fp = args.tce_table
    num_processes = 6
    decision_thresholds = {
        "MES": 6.2,
        "N_transit": 3,
        "SHP": 0.6,
        "MS1": 0.2,
        "MS2": 0.8,
        "MS3": 0.8,
        "chases": 0.78,
        "DMM": 1.5,
        "max_SES_to_MES": 0.8,
        "AIC1": -60,
        "AIC2": -30,
        "SWEET": 15,
        "ASYM": 8,
        "CHI": 7.8,
        "frac_gap": 0.5,
        "V_shape": 1.5,
        "size": 22,
        "MS4": 0,
        "MS5": -1,
        "MS6": -1,
        "offset": 15,
    }
    
    plot_modshift_flag = False  # if True, modshift plots will be generated and saved for each TCE. If False, modshift plots will not be generated.
    plot_summary_flag = False  # if True, summary plots will be generated and saved for each TCE. If False, summary plots will not be generated.
    res_dir = Path('/data/LEO-vetter/results/')
    lc_data_dir = Path('/data/LEO-vetter/data/lcs')
    
    run_pipeline(tce_tbl_fp, decision_thresholds, save_lc_dir=lc_data_dir, res_dir=res_dir, lc_source=args.lc_source, delete_lc_after_target=args.delete_lc_after_target,
                 plot_modshift_flag=plot_modshift_flag, plot_summary_flag=plot_summary_flag, num_processes=num_processes, additional_metadata=additional_metadata)
    
