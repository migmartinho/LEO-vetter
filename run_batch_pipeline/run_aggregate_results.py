"""Aggregate results from the LEO-vetter pipeline."""

import argparse
import glob
from pathlib import Path

from run_batch_pipeline.run_pipeline import aggregate_metrics, aggregate_fa_fp_tests

def main():
    
    parser = argparse.ArgumentParser(
        description="Aggregate results from the LEO-vetter pipeline."
    )
    # Support both legacy positional args and explicit --flag style args.
    parser.add_argument(
        "input_dir",
        nargs="?",
        type=str,
        help="Directory containing the output files from the LEO-vetter pipeline.",
    )
    parser.add_argument(
        "output_metrics_file",
        nargs="?",
        type=str,
        help="Path to the output file where the aggregated metrics will be saved.",
    )
    parser.add_argument(
        "output_fa_fp_file",
        nargs="?",
        type=str,
        help="Path to the output file where the aggregated FA/FP tests will be saved.",
    )
    parser.add_argument(
        "--input_dir",
        dest="input_dir_opt",
        type=str,
        help="Directory containing the output files from the LEO-vetter pipeline.",
    )
    parser.add_argument(
        "--output_metrics_file",
        dest="output_metrics_file_opt",
        type=str,
        help="Path to the output file where the aggregated metrics will be saved.",
    )
    parser.add_argument(
        "--output_fa_fp_file",
        dest="output_fa_fp_file_opt",
        type=str,
        help="Path to the output file where the aggregated FA/FP tests will be saved.",
    )
    args = parser.parse_args()

    input_dir = args.input_dir_opt or args.input_dir
    output_metrics_file = args.output_metrics_file_opt or args.output_metrics_file
    output_fa_fp_file = args.output_fa_fp_file_opt or args.output_fa_fp_file

    if not input_dir or not output_metrics_file or not output_fa_fp_file:
        parser.error(
            "the following arguments are required: input_dir, output_metrics_file, output_fa_fp_file"
        )

    
    metrics_dir = Path(input_dir) / 'metrics'
    # Find all output files in the metrics directory
    input_files = glob.glob(f"{metrics_dir}/*.csv")
    if not input_files:
        print(f"No CSV files found in {metrics_dir}.")
        return
    aggregate_metrics(metrics_dir, output_metrics_file)
    fa_fp_dir = Path(input_dir) / 'fa_fp_tests'
    # Find all output files in the FA/FP tests directory
    input_files = glob.glob(f"{fa_fp_dir}/*.csv")
    if not input_files:
        print(f"No CSV files found in {fa_fp_dir}.")
        return
    aggregate_fa_fp_tests(fa_fp_dir, output_fa_fp_file)
    
if __name__ == "__main__":
    
    main()