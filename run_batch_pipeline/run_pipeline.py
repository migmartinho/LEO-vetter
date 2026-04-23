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

Optional columns for stellar parameters:
- tce_srad_dv (tce_srad_err_dv): stellar radius (uncertainty) in R_s
- tce_smass_dv (tce_smass_err_dv): stellar mass (uncertainty) in M_s
- tce_sdens_dv (tce_sdens_err_dv): stellar density (uncertainty) in 
- tce_steff_dv (tce_steff_err_dv): stellar effective temperature (uncertainty) in K
- tce_slogg_dv (tce_slogg_err_dv): stellar surface gravity (uncertainty) in log
"""

# Suppress expected RuntimeWarnings from edge cases (divide by zero, empty arrays, etc)
# This must be set before any leo_vetter imports to apply globally
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# imports
import argparse
import itertools
import pandas as pd
from pathlib import Path
import lightkurve as lk
import numpy as np
from astroquery.mast import Catalogs, Observations
from astroquery.exceptions import RemoteServiceError
from multiprocessing import Pool
from tqdm import tqdm
import logging
import time
import yaml
from requests.exceptions import RequestException
from urllib3.exceptions import HTTPError
import csv
import shutil

from leo_vetter.stellar import quadratic_ldc
from leo_vetter.main import TCELightCurve
from leo_vetter.plots import plot_modshift, plot_summary
from leo_vetter.thresholds import check_thresholds

Observations.enable_cloud_dataset()


TIC_COLUMNS = ["rad", "mass", "rho", "Teff", "logg"]
TCE_COLUMNS = ["target_id", "uid", "sector_run", "tce_time0bk", "tce_period", "tce_duration", "tce_plnt_num", "sector_run", "sectors_observed"]
# TCE_STELLAR_COLUMNS = ["tce_smass_dv", "tce_smass_err_dv", "tce_srad_dv", "tce_srad_err_dv", "tce_sdens_dv", "tce_sdens_err_dv", "tce_steff_dv", "tce_steff_err_dv", "tce_slogg_dv", "tce_slogg_err_dv"]
TCE_STELLAR_COLUMNS = ["tic_smass", "tic_smass_err", "tic_sradius", "tic_sradius_err", "tic_sdens", "tic_sdens_err", "tic_steff", "tic_steff_err", "tic_slogg", "tic_slogg_err"]
LC_SOURCE_OPTIONS = ("2min", "ffi")
REMOTE_RETRY_ATTEMPTS = 4
REMOTE_RETRY_BASE_DELAY_SECONDS = 2.0

STELLAR_DEFAULTS = {
    "rad": 1.0,
    "mass": 1.0,
    "rho": 1.0,
    "Teff": 5777.0,
    "logg": 4.44,
}
MAP_STELLAR_TCE_TABLE_NAMES = {
    "tic_smass": 'mass',
    "tic_smass_err": 'e_mass',
    "tic_sradius": 'rad',
    "tic_sradius_err": 'e_rad',
    "tic_sdens": 'rho',
    "tic_sdens_err": 'e_rho',
    "tic_steff": 'Teff',
    "tic_steff_err": 'e_Teff',
    "tic_slogg": 'logg',
    "tic_slogg_err": 'e_logg',
}


def get_logger():
    """Get the logger for the pipeline."""
    return logging.getLogger("leo_vetter_pipeline")


def is_retryable_remote_error(error):
    """Determine if an error is a retryable remote error (e.g., network issues, timeouts).
    
    :param Exception error: the exception to check
    :return bool: True if the error is a retryable remote error, False otherwise
    """
    
    retryable_types = (
        RequestException,
        HTTPError,
        ConnectionError,
        TimeoutError,
        OSError,
        RemoteServiceError,
    )
    
    if not isinstance(error, retryable_types):
        return False
    
    # Check for transient MAST error indicators in message
    error_msg = str(error).lower()
    transient_indicators = [
        'timeout',
        'pool',
        'gateway',
        'service unavailable',
        'tempdb',
        'connection',
    ]
    
    return any(indicator in error_msg for indicator in transient_indicators)


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
    :return str|None: `sectors_observed` in new format, or None when missing/empty
    """

    if pd.isna(sectors_observed_binary_string):
        return None

    sectors_observed_binary_string = str(sectors_observed_binary_string).strip()
    if not sectors_observed_binary_string:
        return None

    if sectors_observed_binary_string.lower() in {"nan", "none"}:
        return None

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

    if lc_source == "2min":
        if sector_numbers is None:
            # Get all 2-min light curves for this TIC (no sector filtering)
            local_lc_files = list(save_lc_dir.rglob(f"tess*-s????-{tic_pattern}-*_lc.fits"))
        else:
            sector_tokens = {f"s{sector:04d}" for sector in sector_numbers}
            local_lc_files = []
            for sector_token in sector_tokens:
                local_lc_files.extend(
                    save_lc_dir.rglob(f"tess*-{sector_token}-{tic_pattern}-*_lc.fits")
                )
        return sorted(set(local_lc_files))

    if lc_source == "ffi":
        if sector_numbers is None:
            # Get all FFI light curves for this TIC (no sector filtering)
            local_lc_files = list(save_lc_dir.rglob(f"hlsp_tess-spoc_tess_phot_{tic_pattern}-s????_tess_v1_lc.fits"))
        else:
            sector_tokens = {f"s{sector:04d}" for sector in sector_numbers}
            local_lc_files = []
            for sector_token in sector_tokens:
                local_lc_files.extend(
                    save_lc_dir.rglob(f"hlsp_tess-spoc_tess_phot_{tic_pattern}-{sector_token}_tess_v1_lc.fits")
                )
        return sorted(set(local_lc_files))

    raise ValueError(f"Unsupported lc_source: {lc_source}")


def cleanup_cached_lc_files_for_tic(tic, save_lc_dir, lc_source, lc_files_to_cleanup=None):
    """Delete cached light-curve files and their parent directories for a TIC.

    MAST downloads create a directory per product with the FITS file inside.
    This deletes both the files and their parent directories.
    Uses both explicitly tracked files and a full TIC cache scan so cleanup still
    works when failures happen before files are tracked.
    """

    files_to_delete = set(lc_files_to_cleanup or [])

    try:
        files_to_delete.update(get_cached_lc_files(tic, None, save_lc_dir, lc_source))
    except Exception:
        # Cleanup should be best-effort and never break TIC processing.
        pass

    # Delete parent directories (MAST creates a dir per product with the FITS file inside)
    dirs_to_delete = set()
    for lc_file in files_to_delete:
        parent = lc_file.parent
        if parent.exists() and parent != save_lc_dir:
            dirs_to_delete.add(parent)

    for lc_dir in sorted(dirs_to_delete):
        try:
            shutil.rmtree(lc_dir, ignore_errors=True)
        except Exception:
            # Best-effort; don't fail TIC processing if cleanup fails
            pass


def get_lc_data(tic, sectors_observed, save_lc_dir, lc_source):
    """Gets light curve data for TIC ID `tic` in sectors `sectors_observed`.
    
    :param int tic: TIC ID
    :param str|None sectors_observed: sectors observed separated by "_", or None to use all sectors
    :param str save_lc_dir: light curve directory
    :param str lc_source: either "2min" or "ffi" for SPOC 2-min/FFI light curves, respectively
    :return tuple: light curve object and list of local light curve files
    """

    sectors_numbers = None if sectors_observed is None else [int(sector) for sector in sectors_observed.split('_')]

    local_lc_files = get_cached_lc_files(tic, sectors_numbers, save_lc_dir, lc_source)

    if len(local_lc_files) == 0:
        # TODO: temporarily disable MAST download in favor of manual download until we can implement robust retry logic and cleanup to handle MAST instability issues
        raise ValueError(f'No download permitted. TIC {tic}')
        # lcs = lk.LightCurveCollection([lk.read(local_lc_file) for local_lc_file in local_lc_files])
    # else:
        # author = "SPOC" if lc_source == "2min" else "TESS-SPOC"
        # search_result = retry_remote_call(
        #     lambda: lk.search_lightcurve(
        #         f"TIC {tic}",
        #         mission="TESS",
        #         author=author,
        #         sector=sectors_numbers,
        #     ),
        #     f"Light curve search for TIC {tic} sectors {sectors_observed}",
        # )
        # lcs = retry_remote_call(
        #     lambda: search_result.download_all(download_dir=str(save_lc_dir)),
        #     f"Light curve download for TIC {tic} sectors {sectors_observed}",
        # )
        # Build query parameters conditionally.
        # If sectors_observed is None, we skip sector filtering in product filtering.
        # FFI SPOC data is in HLSP collection with provenance_name='TESS-SPOC'
        # 2-min SPOC data is in TESS collection
        if lc_source == 'ffi':
            query_kwargs = {
                'target_name': tic,
                'obs_collection': 'HLSP',
                'provenance_name': 'TESS-SPOC',
            }
        else:  # 2min
            query_kwargs = {
                'target_name': tic,
                'obs_collection': 'TESS',
            }
        obs_table = retry_remote_call(
            lambda: Observations.query_criteria(**query_kwargs),
            f"MAST observation query for TIC {tic} ({lc_source})",
        )
        if len(obs_table) == 0:
            raise ValueError(f'No observations found for TIC {tic} at the MAST for {lc_source} data. Skipping.')
        # get table with all available products for queried observations
        products = retry_remote_call(
            lambda: Observations.get_product_list(obs_table),
            f"MAST product list query for TIC {tic} ({lc_source})",
        )

        if len(products) == 0:
            raise ValueError(f'No products found for TIC {tic} at the MAST for {lc_source} data. Skipping.')

        # Apply sector filtering only if sectors_numbers is not None
        product_filenames = products['productFilename']
        if sectors_numbers is not None:
            sector_strings = [f'-s{str(sector).zfill(4)}' for sector in sectors_numbers]
            mask = [
                fn.endswith('lc.fits') and
                'fast-lc' not in fn and
                any(sector_str in fn for sector_str in sector_strings)
                for fn in product_filenames
            ]
        else:
            # No sector filtering - get all light curves
            mask = [
                fn.endswith('lc.fits') and
                'fast-lc' not in fn
                for fn in product_filenames
            ]

        lc_products = products[mask]
        if len(lc_products) == 0:
            raise ValueError(f'No TESS light curve files found for TIC {tic} in {lc_source} data. Skipping.')

        _ = Observations.download_products(lc_products, download_dir=str(save_lc_dir), mrp_only=False)

        local_lc_files = get_cached_lc_files(tic, sectors_numbers, save_lc_dir, lc_source)

    lcs = lk.LightCurveCollection([lk.read(local_lc_file) for local_lc_file in local_lc_files])
    
    if lcs is None or len(lcs) == 0:
        raise FileNotFoundError(
            f"No TESS light curves available for TIC {tic} in sectors {sectors_observed} for {lc_source} data either locally or at the MAST. Skipping."
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

def get_stellar_parameters_for_tic(tic, tic_data=None, query_tic=False):
    """Queries the TIC catalog for stellar parameters for a given TIC ID and returns them in a dictionary, along with limb-darkening coefficients.
    
    :param int tic: TIC ID
    :param pd.DataFrame tic_data: DataFrame containing TCE data for the TIC, indexed by TCE UID
    :param bool query_tic: if True, TIC is queried to get stellar paramters. Otherwise, they are grabbed from the TIC TCE data
    :return dict: dictionary containing stellar parameters and limb-darkening coefficients for the target TIC
    """
    
    if tic_data is None and not query_tic:
        raise ValueError(f'No stellar parameters were provided for TIC {tic} and querying TIC catalog was disabled.')

    def _to_finite_float(value, fallback):
        try:
            cast_value = float(value)
        except (TypeError, ValueError):
            return fallback
        return cast_value if np.isfinite(cast_value) else fallback

    if query_tic:
        result = retry_remote_call(
            lambda: Catalogs.query_criteria(catalog="TIC", ID=tic),
            f"TIC catalog query for TIC {tic}",
        )
    else:            
        result = {tic_col: tic_data[tic_col].iloc[0] for tic_col in TIC_COLUMNS}
        result.update({f'e_{tic_col}': tic_data[f'e_{tic_col}'].iloc[0] for tic_col in TIC_COLUMNS})
        
    star = {"tic": tic}

    # Use solar-like defaults when TIC metadata are missing/non-finite.
    for key in TIC_COLUMNS:
        star[key] = _to_finite_float(result[key], STELLAR_DEFAULTS[key])
        star["e_" + key] = _to_finite_float(result["e_" + key], 0.0)

    # Get limb-darkening parameters from sanitized Teff/logg values.
    star["u1"], star["u2"] = quadratic_ldc(star["Teff"], star["logg"])

    return star

def generate_tce_metrics(tic_id, per, epo, dur, lc_tic, tic_params, planetno=1, verbose=False):
    """Generates metrics for the TCE.

    :param int tic_id: TIC ID
    :param float per: period (day)
    :param float epo: epoch (BTJD)
    :param float dur: transit duration (day)
    :param lk.Lightcurve lc_tic: target light curve object
    :param dict tic_params: target stellar parameters
    # :param Path metrics_save_fp: filepath used to save metrics CSV file
    :param int planetno: SPOC planet number, defaults to 1
    :param bool verbose: whether to print verbose output, defaults to False
    :return TCELightCurve: TCE object
    """
    
    time, raw, flux, flux_err = get_lc_data_for_tce(lc_tic, epo, per, dur)
        
    tlc = TCELightCurve(tic_id, time, raw, flux, flux_err, per, epo, dur, planetno=planetno)
        
    tlc.compute_flux_metrics(tic_params, verbose=verbose)
    
    # tlc.save_metrics(save_file=metrics_save_fp)
    
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
    FA, FA_failed_tests = check_thresholds(tlc.metrics, "FA", thresholds=decision_thresholds, verbose=verbose) 
    # FP is True if any tests failed; False otherwise
    FP, FP_failed_tests = check_thresholds(tlc.metrics, "FP", thresholds=decision_thresholds, verbose=verbose)
    
    failed_tests = '_'.join(FA_failed_tests + FP_failed_tests) if (FA_failed_tests or FP_failed_tests) else 'None'
    
    fa_fp_tests_df = pd.DataFrame({
        'uid': [tce_uid],
        'FA': [FA],
        'FP': [FP],
        'Failed Tests': [failed_tests],
    })

    if not FA and not FP and verbose:
        print(f"TCE {tce_uid} is a planet candidate!")   

    return fa_fp_tests_df

def process_tic(tic_id, tic_data, decision_thresholds, save_lc_dir, lc_source, delete_lc_after_target=False, plot_modshift_flag=False, plot_summary_flag=False, plot_modshift_save_dir=None, plot_summary_save_dir=None, metrics_save_dir=None, fa_fp_tests_save_dir=None, query_tic=False, verbose=False):
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
    :param bool query_tic: if True, TIC catalog is queried for stellar parameters for each target
    :param bool verbose: whether to print verbose output during processing, defaults to False
    :return dict: dictionary containing processing results for the TIC
    """
    
    logger = get_logger()

    try:
        tic_params = get_stellar_parameters_for_tic(tic_id, tic_data, query_tic)
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
    tlc_tcelst = []
    
    # Use all-sectors mode when all rows for this TIC have empty sectors_observed.
    use_all_sectors_mode = tic_data["sectors_observed"].isna().all()
    
    # If using all sectors, fetch light curves once before the sector_run loop
    lc_tic = None
    if use_all_sectors_mode:
        try:
            lc_tic, lc_files_used = get_lc_data(tic_id, None, save_lc_dir, lc_source)
            lc_files_to_cleanup.update(lc_files_used)
        except Exception as error:
            logger.exception(
                "Skipping TIC %s after light-curve retrieval failure (all sectors mode)",
                tic_id,
            )
            if delete_lc_after_target:
                cleanup_cached_lc_files_for_tic(tic_id, save_lc_dir, lc_source, lc_files_to_cleanup)
            return {
                "tic_id": tic_id,
                "status": "failed",
                "processed_tces": 0,
                "failed_sector_runs": tic_data["sector_run"].nunique(),
                "error": str(error),
            }
    
    for sector_run, tic_data_sector_run in tqdm(tic_data.groupby('sector_run'), desc=f'Processing TIC {tic_id}', unit='sector run', total=tic_data["sector_run"].nunique()):

        sectors_observed = tic_data_sector_run["sectors_observed"].iloc[0]

        # Allow mixed tables: empty sectors_observed in a sector_run means all sectors for that run.
        if pd.isna(sectors_observed):
            sectors_observed = None

        # Skip light curve retrieval if already fetched in all-sectors mode
        if not use_all_sectors_mode:
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

            tlc_tce = generate_tce_metrics(tic_id, per, epo, dur, lc_tic, tic_params, planetno=planetno, verbose=verbose)
            tlc_tce.metrics['uid'] = tce_uid

            tlc_tcelst.append(tlc_tce)

            fa_fp_tests_df_tce = check_thresholds_tce(tlc_tce, decision_thresholds, tce_uid, verbose=verbose)
            fa_fp_tests_df_tic_lst.append(fa_fp_tests_df_tce)
            processed_tces += 1

            if plot_modshift_flag:
                plot_modshift(tlc_tce, save_fig=plot_modshift_flag, save_file=plot_modshift_save_dir / f"modshift_tic{tic_id}_tce{tce_uid}.png")
            if plot_summary_flag:
                plot_summary(tlc_tce, tic_params, save_fig=plot_summary_flag, save_file=plot_summary_save_dir / f"summary_tic{tic_id}_tce{tce_uid}.png")

    # if len(tlc_tcelst) == 0:
    #     logger.warning(f"No TCEs were successfully processed for TIC {tic_id} in any sector run. TIC will be marked as failed.")
    #     no_tces = True
    # if len(fa_fp_tests_df_tic_lst) == 0:
    #     logger.warning(f"No FA/FP tests were successfully generated for TIC {tic_id} in any sector run. TIC will be marked as partial failure.")
    #     no_fa_fp_tests = True  
    # if not no_fa_fp_tests and not no_tces:
    try:
        saved_metrics_tic, saved_fa_fp_tests_tic = False, False
        
        # save metrics for target's TCEs
        if len(tlc_tcelst) > 0:
            metrics_tic_df = metrics_save_dir / f"metrics_tic{tic_id}.csv"
            with open(metrics_tic_df, "w") as f:
                writer = csv.writer(f, delimiter=",")
                for tce_i, tlc_tce in enumerate(tlc_tcelst):
                    if tce_i == 0:
                        writer.writerow(tlc_tce.metrics.keys())
                    writer.writerow(tlc_tce.metrics.values())
            
            saved_metrics_tic = True
                
        # save FA/FP threshold tests for target's TCEs
        fa_fp_tests_df_tic_lst = [df for df in fa_fp_tests_df_tic_lst if df is not None]
        if fa_fp_tests_df_tic_lst:
            fa_fp_tests_df_tic_df = pd.concat(fa_fp_tests_df_tic_lst, axis=0)
            fa_fp_tests_df_tic_df.to_csv(fa_fp_tests_save_dir / f"fa_fp_tests_tic{tic_id}.csv", index=False)
        
            saved_fa_fp_tests_tic = True
        
    except Exception as error:
        if not saved_metrics_tic:  # if metrics failed to save, we consider the whole TIC a failure and delete any partial outputs
            logger.exception(
                "Failed to save metrics for TIC %s", tic_id
            )
            metrics_tic_df.unlink(missing_ok=True)
            (fa_fp_tests_save_dir / f"fa_fp_tests_tic{tic_id}.csv").unlink(missing_ok=True)
        if not saved_fa_fp_tests_tic:  # if FA/FP tests failed to save, we consider it a partial failure and only delete the FA/FP tests output
            logger.exception(
                "Failed to save FA/FP tests for TIC %s", tic_id
            )
            (fa_fp_tests_save_dir / f"fa_fp_tests_tic{tic_id}.csv").unlink(missing_ok=True)

    if delete_lc_after_target:
        cleanup_cached_lc_files_for_tic(tic_id, save_lc_dir, lc_source, lc_files_to_cleanup)

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
    

def read_tce_table(tce_tbl_fp, get_stellar_parameters_tic_from_table=False):
    """Reads TCE table and prepares it for the run.

    :param str tce_tbl_fp: filepath to TCE table
    :param bool get_stellar_parameters_tic_from_table: if False, stellar parameters are fetched from the TCE table
    :return pd.DataFrame: loaded TCE table
    """
    
    if get_stellar_parameters_tic_from_table:
        tce_tbl = pd.read_csv(tce_tbl_fp, usecols=TCE_COLUMNS + TCE_STELLAR_COLUMNS, on_bad_lines='skip', engine='python', dtype={'sectors_observed': str})
        tce_tbl = tce_tbl.rename(columns=MAP_STELLAR_TCE_TABLE_NAMES)
    else:
        tce_tbl = pd.read_csv(tce_tbl_fp, usecols=TCE_COLUMNS, on_bad_lines='skip', engine='python', dtype={'sectors_observed': str})

    tce_tbl = tce_tbl.loc[tce_tbl['uid'] == '229685063-1-S14-41']
    
    tce_tbl = tce_tbl.rename(columns={"target_id": "tic"})
    tce_tbl['tce_duration'] = tce_tbl['tce_duration'] / 24. # convert duration from hours to days
    
    tce_tbl['sectors_observed'] = tce_tbl['sectors_observed'].apply(convert_sectors_observed_binary_string_to_int_string)
    
    return tce_tbl

def run_pipeline(tce_tbl_fp, decision_thresholds, save_lc_dir, res_dir, lc_source="2min", delete_lc_after_target=False, plot_modshift_flag=False, plot_summary_flag=False, num_processes=4, 
                 additional_metadata=None, query_tic=False, aggregate_checkpoint_tces=0, verbose=False, tic_timeout=600, use_all_observed_sectors=False):
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
    :param bool query_tic: if True, TIC catalog is queried for stellar parameters for each target
    :param int aggregate_checkpoint_tces: if >0, aggregate results and delete individual CSV files whenever this many new TCEs are processed
    :param bool verbose: whether to print verbose output during processing, defaults to False
    :param int tic_timeout: seconds to wait for a single TIC worker before recording it as timed-out and moving on; 0 means no timeout, defaults to 7200
    :param bool use_all_observed_sectors: if True, ignore table sector filters and process each TIC using all available sectors from the light-curve query (no pre-query pass), defaults to False
    """
    
    res_dir = Path(res_dir)
    save_lc_dir = Path(save_lc_dir)
    
    res_dir.mkdir(exist_ok=True)
    logger = setup_logging(res_dir / "logs")
    
    metrics_save_dir = res_dir / "metrics"
    fa_fp_tests_save_dir = res_dir / "fa_fp_tests"
    plot_modshift_save_dir = res_dir / 'modshift_plots'
    plot_summary_save_dir = res_dir / "summary_plots"
    
    metrics_save_dir.mkdir(exist_ok=True)
    fa_fp_tests_save_dir.mkdir(exist_ok=True)
    if plot_modshift_flag:
        plot_modshift_save_dir.mkdir(exist_ok=True)
    if plot_summary_flag:
        plot_summary_save_dir.mkdir(exist_ok=True)
    
    save_lc_dir.mkdir(exist_ok=True, parents=True)
    
    # save decision thresholds dictionary to a CSV in the results directory for record-keeping
    decision_thresholds_df = pd.DataFrame.from_dict(decision_thresholds, orient='index', columns=['threshold'])
    # add metadata to the dataframe
    decision_thresholds_df.attrs['description'] = "Decision thresholds used for FA/FP classification in the LEO-vetter pipeline. If a TCE's metric value exceeds the threshold for a given test, it fails that test. FA is True if any tests failed; FP is True if any tests failed. These thresholds are applied to the metrics computed for each TCE to determine its FA/FP classification."
    decision_thresholds_df.attrs['source'] = "Defined in run_pipeline.py and saved here for record-keeping."
    decision_thresholds_df.attrs['notes'] = "These thresholds can be adjusted based on the desired balance between false positives and false negatives. They were chosen based on analysis of known planets and false positives in TESS data, but may be further refined with additional data and testing."
    decision_thresholds_df.attrs['created'] = pd.Timestamp.now().isoformat()
    if additional_metadata:
        for key, value in additional_metadata.items():
            decision_thresholds_df.attrs[key] = value
    with open(res_dir / "decision_thresholds.csv", "w") as f:
        for key, value in decision_thresholds_df.attrs.items():
            f.write(f"{key}: {value}\n")
        decision_thresholds_df.to_csv(f, index=True)
    
    if query_tic:
        get_stellar_params_from_tce_table = False
    else:
        get_stellar_params_from_tce_table = True
    tce_tbl = read_tce_table(tce_tbl_fp, get_stellar_parameters_tic_from_table=get_stellar_params_from_tce_table)
    
    # Optional override: ignore table-provided sector filters.
    # This does NOT do a pre-query; get_lc_data will fetch without sector filtering.
    if use_all_observed_sectors:
        logger.info("use_all_observed_sectors=True: ignoring sectors_observed from table (single-query mode)")
        tce_tbl['sectors_observed'] = None
    
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
        f"lc_source={lc_source}, delete_lc_after_target={delete_lc_after_target}, "
        f"use_all_observed_sectors={use_all_observed_sectors}"
    )
    
    completed_tces = 0
    last_checkpoint_tces = 0
    with Pool(processes=num_processes) as pool, tqdm(total=len(tic_jobs), desc='Processing TICs', unit='TIC') as pbar:

        pending_results = {}
        tic_jobs_iter = iter(tic_jobs)

        def _submit_next():
            """Submit the next job from the queue if one is available."""
            tic_job = next(tic_jobs_iter, None)
            if tic_job is None:
                return
            ar = pool.apply_async(
                process_tic,
                args=(*tic_job, decision_thresholds, save_lc_dir, lc_source, delete_lc_after_target, plot_modshift_flag, plot_summary_flag, plot_modshift_save_dir, plot_summary_save_dir, metrics_save_dir, fa_fp_tests_save_dir, query_tic, verbose),
            )
            pending_results[ar] = (time.time(), tic_job[0])  # (start_time, tic_id)

        # Prime the pool with up to num_processes jobs so start_time reflects
        # actual execution start rather than queue-submission time.
        for _ in range(min(num_processes, len(tic_jobs))):
            _submit_next()

        while pending_results:
            finished_results = []

            for async_result, (start_time, tic_id_for_result) in list(pending_results.items()):
                result = None

                if async_result.ready():
                    try:
                        result = async_result.get(timeout=0)
                    except Exception as error:
                        logger.exception(
                            "A TIC worker failed with an exception and will be recorded as failed: %s",
                            error,
                        )
                        result = {
                            "tic_id": tic_id_for_result,
                            "status": "failed",
                            "processed_tces": 0,
                            "failed_sector_runs": 0,
                            "error": str(error),
                        }
                # elif tic_timeout and (time.time() - start_time) > tic_timeout:
                #     logger.error(
                #         f"A TIC worker timed out after {tic_timeout} seconds and will be recorded as failed: start {start_time}, end {time.time()}."
                #     )
                #     result = {
                #         "tic_id": tic_id_for_result,
                #         "status": "timeout",
                #         "processed_tces": 0,
                #         "failed_sector_runs": 0,
                #         "error": f"timed out after {tic_timeout}s",
                #     }

                if result is None:
                    continue

                # A slot freed up — submit the next queued job immediately so
                # its start_time is recorded close to when it will execute.
                _submit_next()

                finished_results.append(async_result)
                pipeline_results.append(result)
                completed_tces += int(result.get("processed_tces", 0) or 0)
                pbar.update(1)

                if aggregate_checkpoint_tces > 0 and (completed_tces - last_checkpoint_tces) >= aggregate_checkpoint_tces:
                    logger.info(
                        "Reached checkpoint at %s processed TCEs. Aggregating and cleaning up individual CSV files.",
                        completed_tces,
                    )
                    aggregate_and_cleanup_results(res_dir, logger=logger)
                    last_checkpoint_tces = completed_tces

            for async_result in finished_results:
                pending_results.pop(async_result, None)

            if pending_results:
                time.sleep(1)

    if aggregate_checkpoint_tces > 0:
        logger.info("Final aggregate/cleanup pass after pipeline completion")
        aggregate_and_cleanup_results(res_dir, logger=logger)

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
    :return int: number of target metrics CSV files found for the run.
    """
    
    metrics_files = sorted(metrics_dir.glob("metrics_tic*.csv"))
    if not metrics_files:
        return 0
    
    all_metrics = []
    skipped_files = 0
    for metrics_file in metrics_files:
        try:
            df = pd.read_csv(metrics_file)
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
            print(f"Skipping unreadable metrics file {metrics_file.name}: {error}")
            skipped_files += 1
            continue
        all_metrics.append(df)

    if not all_metrics:
        print(
            "No valid per-target metrics files were found to aggregate "
            f"(skipped {skipped_files} files)."
        )
        return 0
    
    all_metrics_df = pd.concat(all_metrics, ignore_index=True)
    
    output_fp = Path(output_fp)
    if output_fp.exists():
        prev_df = pd.read_csv(output_fp)
        all_metrics_df = pd.concat([prev_df, all_metrics_df], ignore_index=True)
    all_metrics_df.drop_duplicates(subset=["uid"], keep="last", inplace=True)
    
    all_metrics_df.sort_values(["tic", "uid"], inplace=True)
    all_metrics_df.set_index("uid", inplace=True)
    all_metrics_df.to_csv(output_fp, index=True)
    
    print(f"Aggregated metrics for {len(all_metrics_df)} TCEs saved to {output_fp}")
    return len(metrics_files)
    

def aggregate_fa_fp_tests(fa_fp_tests_dir, output_fp):
    """Aggregates individual TCE FA/FP tests CSV files into a single CSV file for all TCEs.
    
    :param Path fa_fp_tests_dir: directory containing individual TCE FA/FP tests CSV files
    :param Path output_fp: filepath to save the aggregated FA/FP tests CSV file
    :return int: number of target FA/FP test CSV files found for the run.
    """
    
    fa_fp_files = sorted(fa_fp_tests_dir.glob("fa_fp_tests_tic*.csv"))
    if not fa_fp_files:
        return 0
    
    all_fa_fp_tests = []
    skipped_files = 0
    for fa_fp_tests_file in fa_fp_files:
        try:
            df = pd.read_csv(fa_fp_tests_file)
        except (pd.errors.EmptyDataError, pd.errors.ParserError) as error:
            print(f"Skipping unreadable FA/FP file {fa_fp_tests_file.name}: {error}")
            skipped_files += 1
            continue
        all_fa_fp_tests.append(df)

    if not all_fa_fp_tests:
        print(
            "No valid per-target FA/FP files were found to aggregate "
            f"(skipped {skipped_files} files)."
        )
        return 0
    
    all_fa_fp_tests_df = pd.concat(all_fa_fp_tests, ignore_index=True)
    
    output_fp = Path(output_fp)
    if output_fp.exists():
        prev_df = pd.read_csv(output_fp)
        all_fa_fp_tests_df = pd.concat([prev_df, all_fa_fp_tests_df], ignore_index=True)
    all_fa_fp_tests_df.drop_duplicates(subset=["uid"], keep="last", inplace=True)
    
    all_fa_fp_tests_df.sort_values(["uid"], inplace=True)
    all_fa_fp_tests_df.set_index("uid", inplace=True)
    all_fa_fp_tests_df.to_csv(output_fp, index=True)
    
    print(f"Aggregated FA/FP tests for {len(all_fa_fp_tests_df)} TCEs saved to {output_fp}")
    return len(fa_fp_files)
    
    
def aggregate_and_cleanup_results(res_dir, logger=None):
    """Aggregate per-target CSV outputs and remove individual files after successful writes.
    
    :param str res_dir: path to results directory
    :param Logger logger: logger object
    """

    res_dir = Path(res_dir)
    metrics_dir = res_dir / "metrics"
    fa_fp_tests_dir = res_dir / "fa_fp_tests"
    agg_metrics_fp = res_dir / "agg_metrics.csv"
    agg_fa_fp_tests_fp = res_dir / "agg_fa_fp_tests.csv"

    aggregated_metrics_count = 0
    aggregated_fa_fp_count = 0

    try:
        aggregated_metrics_count = aggregate_metrics(metrics_dir, agg_metrics_fp)
    except Exception as error:
        if logger:
            logger.exception("Failed to aggregate metrics: %s", error)
        else:
            print(f"Failed to aggregate metrics: {error}")

    try:
        aggregated_fa_fp_count = aggregate_fa_fp_tests(fa_fp_tests_dir, agg_fa_fp_tests_fp)
    except Exception as error:
        if logger:
            logger.exception("Failed to aggregate FA/FP tests: %s", error)
        else:
            print(f"Failed to aggregate FA/FP tests: {error}")

    # Only clean up per-target files when both aggregations succeeded, to keep
    # agg_metrics.csv and agg_fa_fp_tests.csv in sync.
    if aggregated_metrics_count > 0 and aggregated_fa_fp_count > 0:
        for metrics_file in metrics_dir.glob("metrics_tic*.csv"):
            try:
                metrics_file.unlink(missing_ok=True)
            except OSError:
                pass
        for fa_fp_file in fa_fp_tests_dir.glob("fa_fp_tests_tic*.csv"):
            try:
                fa_fp_file.unlink(missing_ok=True)
            except OSError:
                pass

    if logger:
        logger.info(
            "Aggregate/cleanup summary: metrics_files=%s, fa_fp_files=%s",
            aggregated_metrics_count,
            aggregated_fa_fp_count,
        )
        
           
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Run LEO-vetter on TESS TCE tables.")
    parser.add_argument(
        "--run_dir",
        type=str,
        required=True,
        help="Path to run directory.",
    )
    parser.add_argument(
        "--lc_dir",
        type=str,
        required=True,
        help="Path to root directory to target light curve files.",
    )
    parser.add_argument(
        "--run_config",
        type=str,
        required=True,
        help="Configuration YAML file for the run.",
    )
    parser.add_argument(
        "--lc_source",
        type=str,
        choices=LC_SOURCE_OPTIONS,
        required=True,
        help="Choose which SPOC light-curve source to use for cache matching and downloads.",
    )
    parser.add_argument(
        "--tce_table",
        required=True,
        type=str,
        help="Path to the TCE table CSV.",
    )
    
    parser.add_argument(
        "--num_processes",
        type=int,
        default=1,
        help="Number of processes used for parallelizing the run.",
    )
    
    parser.add_argument(
        "--query_tic_catalog",
        action="store_true",
        help="Query TIC catalog for stellar parameters for each target.",
    )
    
    parser.add_argument(
        "--delete_lc_after_target",
        action="store_true",
        help="Delete cached light-curve FITS files after all TCEs for the target TIC are processed.",
    )
    
    parser.add_argument(
        "--aggregate_checkpoint_tces",
        type=int,
        default=0,
        help="If >0, periodically aggregate and delete individual CSV outputs every N processed TCEs.",
    )
    
    parser.add_argument(
        "--plot_modshift_flag",
        action="store_true",
        help="Plot modshift plot for each TCE.",
    )
    
    parser.add_argument(
        "--plot_summary_flag",
        action="store_true",
        help="Plot summary plot for each TCE.",
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output during processing.",
    )

    parser.add_argument(
        "--tic_timeout",
        type=int,
        default=600,
        help="Seconds to wait for a single TIC worker before recording it as timed-out and moving on. 0 means no timeout. Default: 600.",
    )

    parser.add_argument(
        "--use_all_observed_sectors",
        action="store_true",
        help="Ignore sectors_observed from the input table and process using all sectors returned by the light-curve query (no pre-query pass).",
    )

    args = parser.parse_args()
    
    with open(args.run_config, 'r') as file:
        run_config = yaml.load(file, Loader=yaml.SafeLoader)
        
        decision_thresholds = run_config['decision_thresholds']
        additional_metadata = run_config['additional_metadata']
    
    run_pipeline(args.tce_table, decision_thresholds, save_lc_dir=args.lc_dir, res_dir=args.run_dir, lc_source=args.lc_source, delete_lc_after_target=args.delete_lc_after_target,
                 plot_modshift_flag=args.plot_modshift_flag, plot_summary_flag=args.plot_summary_flag, num_processes=args.num_processes, additional_metadata=additional_metadata, 
                 query_tic=args.query_tic_catalog, aggregate_checkpoint_tces=args.aggregate_checkpoint_tces, verbose=args.verbose, tic_timeout=args.tic_timeout,
                 use_all_observed_sectors=args.use_all_observed_sectors)
    
