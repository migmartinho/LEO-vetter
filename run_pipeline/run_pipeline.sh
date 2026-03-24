
RUN_PIPELINE_PY=/home6/msaragoc/work_dir/LEO-vetter/run_pipeline/run_batches_tces.py
RUN_DIR=/home6/msaragoc/work_dir/LEO-vetter/results/test-run_3-23-2026_1322
LC_DIR=/nobackup/msaragoc/work_dir/Kepler-TESS_exoplanet/data/FITS_files/TESS/spoc_2min/lc/sectors/
RUN_CONFIG=/home6/msaragoc/work_dir/LEO-vetter/run_pipeline/run_config.yaml
LC_SOURCE=2min
TCE_TABLE=/home6/msaragoc/work_dir/LEO-vetter/results/test-run_3-23-2026_1322/test_tce_tbl.csv
NUM_PROCESSES=1
AGG_CHECKPT_TCES=2
LOG_FP=$RUN_DIR/run_output.txt
export PYTHONPATH=/home6/msaragoc/work_dir/LEO-vetter/

mkdir -p $RUN_DIR

python $RUN_PIPELINE_PY --run_dir=$RUN_DIR --lc_dir=$LC_DIR --run_config=$RUN_CONFIG --lc_source=$LC_SOURCE --tce_table=$TCE_TABLE --num_process=$NUM_PROCESSES --aggregate_checkpoint_tces=$AGG_CHECKPT_TCES &> $LOG_FP
