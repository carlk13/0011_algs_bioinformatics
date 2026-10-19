import argparse
import csv
import json
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from preprocessor.exact_matcher import ExactAdapterTrimmer
from preprocessor.io_utils import stream_lines


def main():
    parser = argparse.ArgumentParser(
        description="Fast KMP Adapter Trimmer for Sequencing Reads"
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/s_3_sequence_1M.txt"),
        help="Path to input .txt or .gz file",
    )
    parser.add_argument(
        "--adapter",
        type=str,
        default="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        help="Adapter sequence to trim",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=1,
        help="Minimum overlap length required for trimming",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/day1"),
        help="Directory to save output stats and reports",
    )

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading input file: {args.input}")
    print(f"Adapter: {args.adapter}")
    print(f"Min Overlap: {args.min_overlap}")
    print(f"Output Directory: {args.output_dir}")

    # Initialize Streamer and Trimmer
    read_stream = stream_lines(args.input, show_progress=True)
    trimmer = ExactAdapterTrimmer(
        adapter_sequence=args.adapter, min_overlap=args.min_overlap
    )

    print("Starting stream processing...")
    (
        total_reads,
        trimmed_reads_count,
        overlap_dist,
        length_dist,
        elapsed_time,
        throughput,
    ) = trimmer.process_stream(read_stream)

    # Derived metrics
    untrimmed_reads_count = total_reads - trimmed_reads_count
    trimmed_pct = (trimmed_reads_count / total_reads * 100) if total_reads > 0 else 0
    untrimmed_pct = (
        (untrimmed_reads_count / total_reads * 100) if total_reads > 0 else 0
    )

    print("\n" + "=" * 40)
    print("       TASK 1 PROCESSING SUMMARY       ")
    print("=" * 40)
    print(f"Total Reads Processed : {total_reads:,}")
    print(f"Trimmed Reads         : {trimmed_reads_count:,} ({trimmed_pct:.2f}%)")
    print(f"Untrimmed Reads       : {untrimmed_reads_count:,} ({untrimmed_pct:.2f}%)")
    print(f"Elapsed Time          : {elapsed_time:.2f} seconds")
    print(f"Throughput            : {throughput:,.2f} reads/sec")

    print("\nTop 5 Remaining Lengths:")
    for length, count in length_dist.most_common(5):
        pct = (count / total_reads * 100) if total_reads > 0 else 0
        print(f"  - Length {length:2d}: {count:,} reads ({pct:.2f}%)")
    print("=" * 40)

    summary_data = {
        "metadata": {
            "input_file": str(args.input),
            "adapter_sequence": args.adapter,
            "min_overlap": args.min_overlap,
        },
        "metrics": {
            "total_reads": total_reads,
            "trimmed_reads": trimmed_reads_count,
            "untrimmed_reads": untrimmed_reads_count,
            "trimmed_percentage": round(trimmed_pct, 4),
            "elapsed_time_seconds": round(elapsed_time, 4),
            "throughput_reads_per_sec": round(throughput, 2),
        },
    }
    json_path = args.output_dir / "task1_summary.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=4)
    print(f"Saved JSON report to: {json_path}")

    rem_csv_path = args.output_dir / "remaining_lengths.csv"
    with open(rem_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["remaining_length", "count", "percentage"])
        for length in sorted(length_dist.keys()):
            count = length_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"Saved remaining lengths CSV to: {rem_csv_path}")

    overlap_csv_path = args.output_dir / "overlap_lengths.csv"
    with open(overlap_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["overlap_length", "count", "percentage"])
        for length in sorted(overlap_dist.keys()):
            count = overlap_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"Saved overlap lengths CSV to: {overlap_csv_path}")


if __name__ == "__main__":
    main()
