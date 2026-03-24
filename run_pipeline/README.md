# LEO Batch Pipeline

This folder contains scripts to run LEO-Vetter over a table of TCEs, classify each TCE, and write per-target and aggregated outputs (metrics and FA/FP threshold-based tests).

## Input Format

The main input is a CSV TCE table passed with `--tce_table`.

Required columns:

1. `target_id` (int): TIC ID
2. `uid` (str): unique ID per TCE (recommended format: `<tic>_<planetno>_<sector_run>`)
3. `sector_run` (str): sector-run label
4. `tce_time0bk` (float): epoch in BTJD
5. `tce_period` (float): period in days
6. `tce_duration` (float): duration in hours
7. `tce_plnt_num` (int): candidate number for the target
8. `sectors_observed` (str): either:
	 - explicit sectors joined by `_` (example: `1_4_27`), or
	 - binary-like mask string (auto-converted by the pipeline)

Optional stellar columns (used when not querying TIC):

1. `tic_smass`, `tic_smass_err`
2. `tic_sradius`, `tic_sradius_err`
3. `tic_sdens`, `tic_sdens_err`
4. `tic_steff`, `tic_steff_err`
5. `tic_slogg`, `tic_slogg_err`

Notes:

1. `tce_duration` is converted from hours to days internally.
2. If `--query_tic_catalog` is used, stellar parameters are pulled from TIC and table stellar columns are not required.

## Example TCE Table

Minimal CSV (required columns only):

```csv
target_id,uid,sector_run,tce_time0bk,tce_period,tce_duration,tce_plnt_num,sectors_observed
1003831,1003831_1_s1-s3,s1-s3,1355.1234,5.678901,2.40,1,1_2_3
```

CSV including optional stellar columns (used when you do not pass `--query_tic_catalog`):

```csv
target_id,uid,sector_run,tce_time0bk,tce_period,tce_duration,tce_plnt_num,sectors_observed,tic_smass,tic_smass_err,tic_sradius,tic_sradius_err,tic_sdens,tic_sdens_err,tic_steff,tic_steff_err,tic_slogg,tic_slogg_err
1003831,1003831_1_s1-3,1-3,1355.1234,5.678901,2.40,1,1_2_3,0.98,0.05,1.02,0.04,1.10,0.15,5750,80,4.42,0.08
```

## Environment And Installation

Run these from the repository root (the folder that contains [setup.py](setup.py)).

Important:

1. Installing only from [requirements.txt](requirements.txt) may not include everything needed for this pipeline workflow and CLI entry usage.
2. For full functionality, install the package itself with `pip install -e .` after dependencies.

Option 1: venv

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
pip install -e .
```

Option 2: conda

```bash
conda env create -f leo-vetter_env.yaml
conda activate leo-vetter
pip install -e .
```

Quick import check:

```bash
python -c "import leo_vetter; print('leo_vetter import ok')"
```

## Pre-Run Checklist

1. Environment is activated and `pip install -e .` completed.
2. `--tce_table` exists and contains all required columns.
3. `--run_config` exists and includes `decision_thresholds` and `additional_metadata`.
4. `--lc_dir` points to your cache/download root and matches `--lc_source` (`2min` or `ffi`).
5. You are using `--num_processes` (plural) and have chosen an appropriate process count.

## Run With Shell Script

Template script: [run_pipeline/run_pipeline.sh](run_pipeline/run_pipeline.sh)

Typical steps:

1. Edit path variables in [run_pipeline/run_pipeline.sh](run_pipeline/run_pipeline.sh):
	 - `RUN_PIPELINE_PY`
	 - `RUN_DIR`
	 - `LC_DIR`
	 - `RUN_CONFIG`
	 - `TCE_TABLE`
2. Ensure argument name is `--num_processes` (plural) when running [run_pipeline/run_batches_tces.py](run_pipeline/run_batches_tces.py).
3. Run:

```bash
bash run_pipeline/run_pipeline.sh
```

## Run With Python Script

Main runner: [run_pipeline/run_batches_tces.py](run_pipeline/run_batches_tces.py)

Example:

```bash
python run_pipeline/run_batches_tces.py \
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

Configuration YAML example is in [run_pipeline/run_config.yaml](run_pipeline/run_config.yaml) and must define:

1. `decision_thresholds`
2. `additional_metadata`

## Expected Outputs

All outputs are written under `--run_dir`.

Always created:

1. `decision_thresholds.csv`: thresholds + metadata snapshot used for the run
2. `pipeline_status.csv`: per-TIC run status summary
3. `logs/pipeline.log`: detailed logs
4. `metrics/metrics_tic*.csv`: per-target metrics files
5. `fa_fp_tests/fa_fp_tests_tic*.csv`: per-target FA/FP outcomes

Optional plots (if flags enabled):

1. `modshift_plots/modshift_tic*_tce*.png`
2. `summary_plots/summary_tic*_tce*.png`

Aggregated outputs:

1. `agg_metrics.csv`
2. `agg_fa_fp_tests.csv`

These are produced by periodic checkpoint aggregation (`--aggregate_checkpoint_tces > 0`) and/or by running [run_pipeline/run_aggregate_results.py](run_pipeline/run_aggregate_results.py).

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

2. Problem: Run script fails with unrecognized argument `--num_process`
	Fix: use `--num_processes` (plural) with [run_pipeline/run_batches_tces.py](run_pipeline/run_batches_tces.py).

3. Problem: No light curves found / download failures from MAST
	Fix: verify `--lc_dir`, `--lc_source`, sector list format in `sectors_observed`, and network access. The pipeline retries transient remote errors, but persistent failures still need path/source/network fixes.

4. Problem: Aggregates missing (`agg_metrics.csv`, `agg_fa_fp_tests.csv`)
	Fix: either run with `--aggregate_checkpoint_tces > 0` or run [run_pipeline/run_aggregate_results.py](run_pipeline/run_aggregate_results.py) after the batch.

5. Problem: Run appears to skip many TCEs unexpectedly
	Fix: check for existing aggregate files in your run directory. The pipeline skips previously processed UIDs found in both aggregate files.

6. Problem: Plots are not generated
	Fix: pass `--plot_modshift_flag` and/or `--plot_summary_flag`. Without flags, plot folders may exist but remain empty.
