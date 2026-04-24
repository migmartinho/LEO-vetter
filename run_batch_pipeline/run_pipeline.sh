
RUN_PIPELINE_PY=/path/to/leo-vetter/codebase/run_batch_pipeline/run_pipeline.py
RUN_DIR=/path/to/save/pipeline/results
LC_DIR=/path/to/light/curve/files
RUN_CONFIG=/path/to/leo-vetter/codebase/run_batch_pipeline/run_config.yaml
LC_SOURCE=2min  # 2min, ffi
TCE_TABLE=/path/to/tce_table.csv
NUM_PROCESSES=1
AGG_CHECKPT_TCES=2
LOG_FP=$RUN_DIR/run_output_$(date +%Y-%m-%d_%H%M%S).txt
export PYTHONPATH=/path/to/leo-vetter/codebase

mkdir -p $RUN_DIR

python $RUN_PIPELINE_PY --run_dir=$RUN_DIR --lc_dir=$LC_DIR --run_config=$RUN_CONFIG --lc_source=$LC_SOURCE --tce_table=$TCE_TABLE --num_processes=$NUM_PROCESSES --aggregate_checkpoint_tces=$AGG_CHECKPT_TCES &> $LOG_FP
