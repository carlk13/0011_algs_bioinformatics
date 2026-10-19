import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

# Add project root and src to sys.path
sys.path.append(str(Path(__file__).resolve().parent.parent))

from preprocessor.approx_matcher import HammingAdapterTrimmer, SemiGlobalDPTrimmer
from preprocessor.error_analyzer import ErrorProfiler
from preprocessor.io_utils import stream_lines


def save_remaining_lengths_csv(
    path: Path, length_dist: Counter, total_reads: int
) -> None:
    """Exports remaining read length distribution to CSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["remaining_length", "count", "percentage"])
        for length in sorted(length_dist.keys()):
            count = length_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"Saved remaining lengths CSV to: {path}")


def save_overlap_lengths_csv(
    path: Path, overlap_dist: Counter, total_reads: int
) -> None:
    """Exports adapter overlap length distribution to CSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["overlap_length", "count", "percentage"])
        for length in sorted(overlap_dist.keys()):
            count = overlap_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"Saved overlap lengths CSV to: {path}")


def save_positional_distribution_csv(path: Path, cycle_dist: list[dict]) -> None:
    """Exports positional error rates per sequencing cycle to CSV."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["cycle_position", "errors", "observations", "error_rate"])
        for entry in cycle_dist:
            writer.writerow(
                [
                    entry["cycle"],
                    entry["errors"],
                    entry["observations"],
                    f"{entry['error_rate']:.6f}",
                ]
            )
    print(f"Saved positional error distribution CSV to: {path}")


def load_task1_baseline(day1_dir: Path, report_cutoff: int):
    """Loads Task 1 (Exact Matching) baseline metrics for comparative summary if available."""
    json_path = day1_dir / "task1_summary.json"
    csv_path = day1_dir / "overlap_lengths.csv"
    if not json_path.exists():
        return None

    try:
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        metrics = data.get("metrics", {})
        total_reads = metrics.get("total_reads", 0)
        trimmed_1 = metrics.get("trimmed_reads", 0)
        elapsed = metrics.get("elapsed_time_seconds", 0.0)

        # Compute cutoff (L >= report_cutoff) from overlap CSV
        trimmed_cutoff = 0
        if csv_path.exists():
            with open(csv_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if int(row["overlap_length"]) >= report_cutoff:
                        trimmed_cutoff += int(row["count"])
        else:
            trimmed_cutoff = trimmed_1

        return {
            "total_reads": total_reads,
            "trimmed_min_overlap": trimmed_1,
            "trimmed_report_cutoff": trimmed_cutoff,
            "elapsed_time_seconds": elapsed,
        }
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Preprocessor Task 2 (Approximate Matching) and Task 3 (Sequencing Error Profiling)"
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/s_3_sequence_1M.txt"),
        help="Path to input .txt or .gz file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/day2"),
        help="Directory to save output stats and reports",
    )
    parser.add_argument(
        "--adapter",
        type=str,
        default="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        help="Adapter sequence to trim",
    )
    parser.add_argument(
        "--min-overlap",
        "--min_overlap",
        type=int,
        default=1,
        dest="min_overlap",
        help="Minimum overlap length required for trimming (default: 1)",
    )
    parser.add_argument(
        "--report-overlap-cutoff",
        type=int,
        default=10,
        help="Minimum overlap length threshold for reliable biological reporting (default: 10)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Limit number of reads processed for quick testing / dry-runs (default: None for full file)",
    )
    parser.add_argument(
        "--dp-sample-size",
        type=int,
        default=20000,
        help="Number of reads to evaluate with Semi-Global DP for indel profiling (0 for full run, default: 20000)",
    )
    parser.add_argument(
        "--no-dp",
        action="store_true",
        help="Skip Semi-Global DP indel evaluation",
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 65)
    print("       TASK 2 & TASK 3 EXECUTION PIPELINE       ")
    print("=" * 65)
    print(f"Input file            : {args.input}")
    print(f"Adapter               : {args.adapter}")
    print(f"Min Overlap           : {args.min_overlap}")
    print(f"Report Overlap Cutoff : {args.report_overlap_cutoff}")
    print(
        f"Read Limit            : {args.limit if args.limit else 'None (Full Dataset)'}"
    )
    print(f"Output Directory      : {args.output_dir}")
    print("=" * 65)

    # -------------------------------------------------------------
    # 1. Hamming Trimming with 10% Mismatch Threshold & Error Profiling
    # -------------------------------------------------------------
    print("\n[Step 1/3] Running Hamming Trimmer (10% mismatch threshold)...")
    read_stream_10 = stream_lines(args.input, limit=args.limit, show_progress=True)
    error_profiler = ErrorProfiler(
        max_cycles=50, min_overlap_cutoff=args.report_overlap_cutoff
    )
    hamming_trimmer_10 = HammingAdapterTrimmer(
        adapter_sequence=args.adapter,
        max_error_rate=0.10,
        min_overlap=args.min_overlap,
    )

    (
        total_reads_10,
        trimmed_reads_count_10,
        overlap_dist_10,
        length_dist_10,
        mismatch_dist_10,
        elapsed_time_10,
        throughput_10,
    ) = hamming_trimmer_10.process_stream(read_stream_10, error_profiler=error_profiler)

    trimmed_cutoff_10 = sum(
        c for l, c in overlap_dist_10.items() if l >= args.report_overlap_cutoff
    )
    pct_cutoff_10 = (
        (trimmed_cutoff_10 / total_reads_10 * 100) if total_reads_10 > 0 else 0
    )

    print(
        f"10% Trimming Complete in {elapsed_time_10:.2f}s ({throughput_10:,.1f} reads/s)"
    )
    print(
        f"  - Matches (L >= {args.min_overlap}): {trimmed_reads_count_10:,} ({trimmed_reads_count_10 / total_reads_10 * 100:.2f}%)"
    )
    print(
        f"  - Confident Matches (L >= {args.report_overlap_cutoff}): {trimmed_cutoff_10:,} ({pct_cutoff_10:.2f}%)"
    )

    # -------------------------------------------------------------
    # 2. Hamming Trimming with 25% Mismatch Threshold
    # -------------------------------------------------------------
    print("\n[Step 2/3] Running Hamming Trimmer (25% mismatch threshold)...")
    read_stream_25 = stream_lines(args.input, limit=args.limit, show_progress=True)
    hamming_trimmer_25 = HammingAdapterTrimmer(
        adapter_sequence=args.adapter,
        max_error_rate=0.25,
        min_overlap=args.min_overlap,
    )

    (
        total_reads_25,
        trimmed_reads_count_25,
        overlap_dist_25,
        length_dist_25,
        mismatch_dist_25,
        elapsed_time_25,
        throughput_25,
    ) = hamming_trimmer_25.process_stream(read_stream_25)

    trimmed_cutoff_25 = sum(
        c for l, c in overlap_dist_25.items() if l >= args.report_overlap_cutoff
    )
    pct_cutoff_25 = (
        (trimmed_cutoff_25 / total_reads_25 * 100) if total_reads_25 > 0 else 0
    )

    print(
        f"25% Trimming Complete in {elapsed_time_25:.2f}s ({throughput_25:,.1f} reads/s)"
    )
    print(
        f"  - Matches (L >= {args.min_overlap}): {trimmed_reads_count_25:,} ({trimmed_reads_count_25 / total_reads_25 * 100:.2f}%)"
    )
    print(
        f"  - Confident Matches (L >= {args.report_overlap_cutoff}): {trimmed_cutoff_25:,} ({pct_cutoff_25:.2f}%)"
    )

    # -------------------------------------------------------------
    # 3. Semi-Global DP Trimming with Indels
    # -------------------------------------------------------------
    total_reads_dp = 0
    trimmed_reads_dp = 0
    trimmed_cutoff_dp = 0
    overlap_dist_dp = Counter()
    length_dist_dp = Counter()
    total_subs_dp = 0
    total_ins_dp = 0
    total_dels_dp = 0
    elapsed_time_dp = 0.0
    throughput_dp = 0.0

    if not args.no_dp:
        # Determine sample limit for DP
        if args.limit is not None and (
            args.dp_sample_size <= 0 or args.limit < args.dp_sample_size
        ):
            dp_limit = args.limit
        else:
            dp_limit = None if args.dp_sample_size <= 0 else args.dp_sample_size

        limit_desc = f"{dp_limit:,} reads sample" if dp_limit else "full dataset"
        print(
            f"\n[Step 3/3] Running Semi-Global DP Trimmer with Indels ({limit_desc})..."
        )

        read_stream_dp = stream_lines(args.input, limit=dp_limit, show_progress=True)
        semi_global_trimmer = SemiGlobalDPTrimmer(
            adapter_sequence=args.adapter,
            min_score_threshold=args.report_overlap_cutoff,
        )

        (
            total_reads_dp,
            trimmed_reads_dp,
            overlap_dist_dp,
            length_dist_dp,
            total_subs_dp,
            total_ins_dp,
            total_dels_dp,
            elapsed_time_dp,
            throughput_dp,
        ) = semi_global_trimmer.process_stream(read_stream_dp, limit=dp_limit)

        trimmed_cutoff_dp = sum(
            c for l, c in overlap_dist_dp.items() if l >= args.report_overlap_cutoff
        )
        pct_cutoff_dp = (
            (trimmed_cutoff_dp / total_reads_dp * 100) if total_reads_dp > 0 else 0
        )

        print(
            f"DP Trimming Complete in {elapsed_time_dp:.2f}s ({throughput_dp:,.1f} reads/s)"
        )
        print(
            f"  - Confident Matches (L >= {args.report_overlap_cutoff}): {trimmed_cutoff_dp:,} ({pct_cutoff_dp:.2f}%)"
        )
        print(
            f"  - Traceback Errors: Subs={total_subs_dp:,}, Ins={total_ins_dp:,}, Dels={total_dels_dp:,}"
        )
    else:
        print("\n[Step 3/3] Skipping Semi-Global DP as requested (--no-dp).")

    # -------------------------------------------------------------
    # 4. Task 3 Sequencing Error Profile Extraction
    # -------------------------------------------------------------
    task3_metrics = error_profiler.get_summary()

    # Incorporate empirical indel counts from DP traceback if available
    total_dp_errors = total_subs_dp + total_ins_dp + total_dels_dp
    if total_dp_errors > 0:
        task3_metrics["error_type_breakdown"] = {
            "substitutions": {
                "count": total_subs_dp,
                "percentage": round(total_subs_dp / total_dp_errors * 100, 2),
            },
            "insertions": {
                "count": total_ins_dp,
                "percentage": round(total_ins_dp / total_dp_errors * 100, 2),
            },
            "deletions": {
                "count": total_dels_dp,
                "percentage": round(total_dels_dp / total_dp_errors * 100, 2),
            },
        }

    # -------------------------------------------------------------
    # 5. Print Summary Tables to Console
    # -------------------------------------------------------------
    task1_baseline = load_task1_baseline(
        args.output_dir.parent / "day1", args.report_overlap_cutoff
    )

    print("\n" + "=" * 65)
    print("                TASK 2 RESULTS SUMMARY                  ")
    print("=" * 65)
    print(
        f"{'Condition':<23} | {'Matches (L>=1)':<16} | {'Matches (L>=10)':<16} | {'Runtime':<10}"
    )
    print("-" * 65)
    if task1_baseline:
        t1_total = task1_baseline["total_reads"]
        t1_min = task1_baseline["trimmed_min_overlap"]
        t1_cut = task1_baseline["trimmed_report_cutoff"]
        t1_time = task1_baseline["elapsed_time_seconds"]
        print(
            f"{'Exact (Task 1)':<23} | {t1_min:>7,} ({t1_min / t1_total * 100:>5.2f}%) | {t1_cut:>7,} ({t1_cut / t1_total * 100:>5.2f}%) | {t1_time:>7.2f}s"
        )
    print(
        f"{'10% Mismatch':<23} | {trimmed_reads_count_10:>7,} ({trimmed_reads_count_10 / total_reads_10 * 100:>5.2f}%) | {trimmed_cutoff_10:>7,} ({pct_cutoff_10:>5.2f}%) | {elapsed_time_10:>7.2f}s"
    )
    print(
        f"{'25% Mismatch':<23} | {trimmed_reads_count_25:>7,} ({trimmed_reads_count_25 / total_reads_25 * 100:>5.2f}%) | {trimmed_cutoff_25:>7,} ({pct_cutoff_25:>5.2f}%) | {elapsed_time_25:>7.2f}s"
    )
    if not args.no_dp:
        print(
            f"{'Semi-Global DP (Indels)':<23} | {trimmed_reads_dp:>7,} ({trimmed_reads_dp / total_reads_dp * 100:>5.2f}%) | {trimmed_cutoff_dp:>7,} ({pct_cutoff_dp:>5.2f}%) | {elapsed_time_dp:>7.2f}s"
        )
    print("=" * 65)

    print("\nTop 5 Remaining Lengths (10% Mismatch):")
    for length, count in length_dist_10.most_common(5):
        pct = (count / total_reads_10 * 100) if total_reads_10 > 0 else 0
        print(f"  - Length {length:2d}: {count:,} reads ({pct:.2f}%)")

    print("\nTop 5 Remaining Lengths (25% Mismatch):")
    for length, count in length_dist_25.most_common(5):
        pct = (count / total_reads_25 * 100) if total_reads_25 > 0 else 0
        print(f"  - Length {length:2d}: {count:,} reads ({pct:.2f}%)")

    print("\n" + "=" * 65)
    print("          TASK 3 SEQUENCING ERROR PROFILE SUMMARY       ")
    print("=" * 65)
    print(
        f"Adapter-Containing Reads (L >= {args.report_overlap_cutoff}) : {task3_metrics['total_adapter_reads']:,}"
    )
    print(
        f"Reads with >= 1 Error                        : {task3_metrics['reads_with_errors']:,}"
    )
    print(
        f"Total Aligned Adapter Bases                  : {task3_metrics['total_aligned_bases']:,}"
    )
    print(
        f"Per-Nucleotide Error Rate                    : {task3_metrics['per_nucleotide_error_rate'] * 100:.3f}% ({task3_metrics['per_nucleotide_error_rate']:.6f})"
    )
    print(
        f"Per-Sequence Error Rate                      : {task3_metrics['per_sequence_error_rate'] * 100:.2f}%"
    )

    print("\nError Type Breakdown (Empirical Traceback):")
    for err_type, stats in task3_metrics["error_type_breakdown"].items():
        print(
            f"  - {err_type.capitalize():<14}: {stats['count']:>6,} ({stats['percentage']:>5.2f}%)"
        )

    cycle_data = task3_metrics["cycle_distribution"]
    print("\nPositional Error Rates Along Read (Sample Cycles):")
    sample_cycles = [0, 10, 20, 30, 40, 45, 49]
    for c_idx in sample_cycles:
        if c_idx < len(cycle_data):
            c_info = cycle_data[c_idx]
            print(
                f"  - Cycle {c_info['cycle']:2d} : {c_info['error_rate'] * 100:.3f}% ({c_info['errors']:,} errors / {c_info['observations']:,} obs)"
            )
    print("=" * 65)

    # -------------------------------------------------------------
    # 6. Save JSON & CSV Deliverables
    # -------------------------------------------------------------
    task2_summary_data = {
        "metadata": {
            "input_file": str(args.input),
            "adapter_sequence": args.adapter,
            "min_overlap": args.min_overlap,
            "report_overlap_cutoff": args.report_overlap_cutoff,
            "read_limit": args.limit,
        },
        "exact_task1_baseline": task1_baseline,
        "mismatch_10pct": {
            "total_reads": total_reads_10,
            "trimmed_reads_min_overlap": trimmed_reads_count_10,
            "trimmed_reads_report_cutoff": trimmed_cutoff_10,
            "trimmed_pct_min_overlap": round(
                trimmed_reads_count_10 / total_reads_10 * 100, 4
            )
            if total_reads_10 > 0
            else 0,
            "trimmed_pct_report_cutoff": round(pct_cutoff_10, 4),
            "elapsed_time_seconds": round(elapsed_time_10, 4),
            "throughput_reads_per_sec": round(throughput_10, 2),
        },
        "mismatch_25pct": {
            "total_reads": total_reads_25,
            "trimmed_reads_min_overlap": trimmed_reads_count_25,
            "trimmed_reads_report_cutoff": trimmed_cutoff_25,
            "trimmed_pct_min_overlap": round(
                trimmed_reads_count_25 / total_reads_25 * 100, 4
            )
            if total_reads_25 > 0
            else 0,
            "trimmed_pct_report_cutoff": round(pct_cutoff_25, 4),
            "elapsed_time_seconds": round(elapsed_time_25, 4),
            "throughput_reads_per_sec": round(throughput_25, 2),
        },
        "indels_dp": {
            "evaluated_reads": total_reads_dp,
            "trimmed_reads_min_overlap": trimmed_reads_dp,
            "trimmed_reads_report_cutoff": trimmed_cutoff_dp,
            "trimmed_pct_report_cutoff": round(pct_cutoff_dp, 4)
            if total_reads_dp > 0
            else 0,
            "substitutions": total_subs_dp,
            "insertions": total_ins_dp,
            "deletions": total_dels_dp,
            "elapsed_time_seconds": round(elapsed_time_dp, 4),
            "throughput_reads_per_sec": round(throughput_dp, 2),
        },
    }

    task2_json_path = args.output_dir / "task2_summary.json"
    with open(task2_json_path, "w", encoding="utf-8") as f:
        json.dump(task2_summary_data, f, indent=4)
    print(f"Saved Task 2 JSON summary to: {task2_json_path}")

    # Remaining lengths CSVs
    save_remaining_lengths_csv(
        args.output_dir / "remaining_lengths_10pct.csv",
        length_dist_10,
        total_reads_10,
    )
    save_remaining_lengths_csv(
        args.output_dir / "remaining_lengths_25pct.csv",
        length_dist_25,
        total_reads_25,
    )
    if not args.no_dp:
        save_remaining_lengths_csv(
            args.output_dir / "remaining_lengths_indels.csv",
            length_dist_dp,
            total_reads_dp,
        )

    # Overlap lengths CSVs (for Task 2 comparative distribution plots)
    save_overlap_lengths_csv(
        args.output_dir / "overlap_lengths_10pct.csv",
        overlap_dist_10,
        total_reads_10,
    )
    save_overlap_lengths_csv(
        args.output_dir / "overlap_lengths_25pct.csv",
        overlap_dist_25,
        total_reads_25,
    )

    # Positional error distribution CSV
    save_positional_distribution_csv(
        args.output_dir / "positional_error_distribution.csv",
        task3_metrics["cycle_distribution"],
    )

    # Task 3 Error Summary JSON
    task3_json_path = args.output_dir / "error_summary_task3.json"
    task3_output_data = {
        "metadata": {
            "input_file": str(args.input),
            "adapter_sequence": args.adapter,
            "report_overlap_cutoff": args.report_overlap_cutoff,
            "dp_sample_evaluated": total_reads_dp,
            "read_limit": args.limit,
        },
        "metrics": {
            "total_adapter_reads": task3_metrics["total_adapter_reads"],
            "reads_with_errors": task3_metrics["reads_with_errors"],
            "total_aligned_bases": task3_metrics["total_aligned_bases"],
            "per_nucleotide_error_rate": task3_metrics["per_nucleotide_error_rate"],
            "per_sequence_error_rate": task3_metrics["per_sequence_error_rate"],
            "error_type_breakdown": task3_metrics["error_type_breakdown"],
        },
    }
    with open(task3_json_path, "w", encoding="utf-8") as f:
        json.dump(task3_output_data, f, indent=4)
    print(f"Saved Task 3 JSON error summary to: {task3_json_path}")


if __name__ == "__main__":
    main()
