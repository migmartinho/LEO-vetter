# LEO-Vetter Batch Pipeline

This folder contains scripts to run LEO-Vetter over a batch (table) of TCEs, classify each TCE, and write per-target and aggregated outputs (metrics and FA/FP threshold-based tests).

Current implementation:
- Designed for both TESS SPOC 2-min and FFI light curves
- Flux-level vetting
	- When metric cannot be computed, metric value is returned as NaN
	- When test cannot be evaluated, the TCE does not fail the test
- **Pixel-level vetting has not been included in the pipeline** (but can be easily added; requirements: install [transit-diffImage](https://github.com/stevepur/transit-diffImage) which requires target pixel files; see [main README.md](/README.md#installation))

## Input Format

The main input is a CSV TCE table passed with `--tce_table`.

Required columns:

1. `target_id` (int): TIC ID
2. `uid` (str): unique ID per TCE (recommended format: `<tic>-<planetno>_<sector_run>`)
3. `sector_run` (str): sector run label (e.g., 1-6)
4. `tce_time0bk` (float): epoch in BTJD
5. `tce_period` (float): period in days
6. `tce_duration` (float): duration in hours
7. `tce_plnt_num` (int): candidate number for the target
8. `sectors_observed` (str): sectors in which the target was observed for this run. It can be a mix of the following options for different TCEs!):
	 - explicit sectors joined by `_` (example: `1_4_27`), or
	 - binary-like mask string (auto-converted by the pipeline to explicit sectors joined by `_`) (example: `0100001` converts to `1_6`)
	 - `None` or `np.nan`; if observed sectors are not provided, then all available sectors for the target will be used

Optional stellar columns (used when not querying TIC for the stellar parameters):

1. `tic_smass`, `tic_smass_err`: stellar mass ($M_{\odot}$)
2. `tic_sradius`, `tic_sradius_err`: stellar radius ($R_{\odot}$)
3. `tic_sdens`, `tic_sdens_err`: stellar density ($\rho_{\odot}$)
4. `tic_steff`, `tic_steff_err`: stellar effective temperature ($K$)
5. `tic_slogg`, `tic_slogg_err`: stellar surface gravity ($log_{10}(cm/s^2)$)

Notes:

1. `tce_duration` is converted from hours to days internally.
2. If `--query_tic_catalog` is used, stellar parameters are pulled from TIC and table stellar columns are not required.
3. You can choose which sectors to use for each TCE. Keep in mind that some TCE ephemerides, especially those derived from early single-sector TESS runs, may be stale.

## Example TCE Table

Minimal CSV (required columns only):

- without observed sectors
```csv
target_id,uid,sector_run,tce_time0bk,tce_period,tce_duration,tce_plnt_num,sectors_observed
1003831,1003831-1_S8,,1518.203536,1.651142,0.758184,1,8
```

- with observed sectors
```csv
target_id,uid,sector_run,tce_time0bk,tce_period,tce_duration,tce_plnt_num,sectors_observed
1003831,1003831-1_S8,8,1518.203536,1.651142,0.758184,1,8
```

CSV including optional stellar columns (used when you do not pass `--query_tic_catalog`):

```csv
target_id,uid,sector_run,tce_time0bk,tce_period,tce_duration,tce_plnt_num,sectors_observed,tic_smass,tic_smass_err,tic_sradius,tic_sradius_err,tic_sdens,tic_sdens_err,tic_steff,tic_steff_err,tic_slogg,tic_slogg_err
1003831,1003831-1_S8,8,1518.203536,1.651142,0.758184,1,8,0.977,0,1.12196,0,0.691766,0,5550,0,4.32801,0
```

## Environment And Installation

Run these from the repository root (the folder that contains [setup.py](/setup.py)).

Important:

1. Installing only from [requirements.txt](requirements.txt) may not include everything needed for this pipeline workflow and CLI entry usage. Option 1 is recommended.
2. For full functionality, install the package itself with `pip install -e .` after dependencies.

Option 1: [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)

```bash
micromamba env create /path/to/new/env -f leo-vetter_env.yaml
micromamba activate leo-vetter
pip install -e .
```

Option 2: venv

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install -e .
```

Quick import check:

```bash
python -c "import leo_vetter; print('leo_vetter import ok')"
```

## Pre-Run Checklist

1. Environment is activated and `pip install -e .` completed.
2. `--tce_table` exists and contains all required columns with the right formats/data types.
3. `--run_config` exists and includes `decision_thresholds` and `additional_metadata`.
4. `--lc_dir` points to your cache/download root and matches `--lc_source` (`2min` or `ffi`).
5. Confirm your sector list corresponds to sectors with available TESS SPOC 2-min/HLSP products for TICs of interest.

## Run With Shell Script

Template script: [run_pipeline/run_pipeline.sh](/run_pipeline/run_pipeline.sh)

Typical steps:

1. Edit path variables in [run_pipeline/run_pipeline.sh](/run_pipeline/run_pipeline.sh):
	 - `RUN_PIPELINE_PY`
	 - `RUN_DIR`
	 - `LC_DIR`
	 - `LC_SOURCE`
	 - `RUN_CONFIG`
	 - `TCE_TABLE`
	 - You can also adjust other variables that have default values, see both shell and Python scripts.
2. Run:

```bash
bash run_pipeline/run_pipeline.sh
```

## Run With Python Script

Main runner: [run_pipeline/run_pipeline.py](/run_pipeline/run_pipeline.py)

Example:

```bash
python run_pipeline/run_pipeline.py \
	--run_dir /path/to/results/run_YYYYMMDD \
	--lc_dir /path/to/lc_cache \
	--run_config run_pipeline/run_config.yaml \
	--lc_source 2min \
	--tce_table /path/to/tces.csv \
	--num_processes 6 \
	--aggregate_checkpoint_tces 500 \
	--query_tic_catalog
```

Common flags:

1. `--lc_source {2min,ffi}`
2. `--delete_lc_after_target`
3. `--plot_modshift_flag`
4. `--plot_summary_flag`
5. `--aggregate_checkpoint_tces N` (periodically aggregates and removes per-target CSVs)
6. `--use_all_observed_sectors` (overwrites `sectors_observed` value in the input table and uses all available sectors for the target)

For full CLI options, run:

```bash
python run_pipeline/run_pipeline.py --help
```

Configuration YAML example is in [run_pipeline/run_config.yaml](/run_pipeline/run_config.yaml) and must define:

1. `decision_thresholds`
2. `additional_metadata`

## Expected Outputs

All outputs are written under `--run_dir`. See output of example run in [run_pipeline/example](/run_pipeline/example/).

Always created:

1. `decision_thresholds.csv`: thresholds + metadata snapshot used for the run
2. `pipeline_status.csv`: per-TIC run status summary
3. `logs/pipeline.log`: detailed logs from the pipeline logger (stdout/err output can also be seen from `run_output.txt`)
4. `metrics/metrics_tic*.csv`: per-target metrics files
5. `fa_fp_tests/fa_fp_tests_tic*.csv`: per-target FA/FP outcomes

Optional plots (if flags enabled):

1. `modshift_plots/modshift_tic*_tce*.png`
2. `summary_plots/summary_tic*_tce*.png`

Aggregated outputs:

1. `agg_metrics.csv`
2. `agg_fa_fp_tests.csv`

These are produced by periodic checkpoint aggregation (`--aggregate_checkpoint_tces > 0`) and/or by running [run_pipeline/run_aggregate_results.py](/run_pipeline/run_aggregate_results.py).

## Aggregate Results Separately

You can aggregate an existing run directory manually:

```bash
python run_pipeline/run_aggregate_results.py \
	--input_dir /path/to/results/run_YYYYMMDD \
	--output_metrics_file /path/to/results/run_YYYYMMDD/agg_metrics.csv \
	--output_fa_fp_file /path/to/results/run_YYYYMMDD/agg_fa_fp_tests.csv
```

## Troubleshooting

1. Problem: `ModuleNotFoundError: No module named 'leo_vetter'`
	Fix: activate the environment and run `pip install -e .` from the repository root.

2. Problem: No light curves found / download failures from MAST
	Fix: verify `--lc_dir`, `--lc_source`, sector list format in `sectors_observed`, and network access. The pipeline retries transient remote errors, but persistent failures still need path/source/network fixes. Ensure that the TICs of interest were observed for the requested sectors in the data collection mode (i.e., TESS SPOC 2min or FFI) you set.

3. Problem: Files were downloaded, but pipeline still says no light curves were found
	Fix: check for `*_lc.fits` files under your `--lc_dir` tree, and verify file permissions/readability. If needed, clear stale partial downloads and rerun for the TIC.

4. Problem: Aggregates missing (`agg_metrics.csv`, `agg_fa_fp_tests.csv`)
	Fix: either run with `--aggregate_checkpoint_tces > 0` or run [run_pipeline/run_aggregate_results.py](/run_pipeline/run_aggregate_results.py) after the pipeline run.

5. Problem: Run appears to skip many TCEs unexpectedly
	Fix: check for existing aggregate files in your run directory. The pipeline skips previously processed TCEs found in both aggregate files.

6. Problem: Plots are not generated
	Fix: pass `--plot_modshift_flag` and/or `--plot_summary_flag`.

7. Problem: MAST is down (you see no light curve files being downloaded and run seems to hang/timeout per TIC is reached)
	Fix: wait until it comes back up or use local target light curve FITS files and add your own stellar parameters to the TCE input table.
