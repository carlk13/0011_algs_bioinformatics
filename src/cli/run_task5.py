import argparse
import csv
import json
import sys
import time
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from preprocessor.demultiplexer import BarcodeIdentifier, LibraryDemultiplexer
from preprocessor.io_utils import stream_lines


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Task 5: Barcode Identification & Library Demultiplexing"
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/MultiplexedSamples"),
        help="Path to pooled multiplexed sequencing reads (.txt or .gz)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/day3/task5"),
        help="Directory to save demultiplexing stats and per-sample reports",
    )
    parser.add_argument(
        "--adapter",
        type=str,
        default="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
        help="3' adapter sequence to trim after barcode excision",
    )
    parser.add_argument(
        "--barcode-length",
        type=int,
        default=4,
        help="Expected length of sample barcodes in bp (default: 4)",
    )
    parser.add_argument(
        "--barcodes",
        nargs="*",
        default=None,
        help="Optional explicit list of known barcodes (if omitted, auto-discovered)",
    )
    parser.add_argument(
        "--max-mismatches",
        type=int,
        default=1,
        help="Maximum allowed mismatches for barcode assignment (default: 1)",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=10,
        help="Minimum overlap length required for 3' adapter trimming (default: 10)",
    )
    parser.add_argument(
        "--sample-limit",
        type=int,
        default=100000,
        help="Sample size for discovering candidate barcodes (default: 100000)",
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("      TASK 5: BARCODE IDENTIFICATION & DEMULTIPLEXING")
    print("=" * 60)
    print(f"Input File     : {args.input}")
    print(f"Barcode Length : {args.barcode_length} bp")
    print(f"Max Mismatches : {args.max_mismatches}")
    print(f"Min Overlap    : {args.min_overlap} bp")
    print(f"Output Dir     : {args.output_dir}")
    print("-" * 60)

    # discover sample barcodes
    if args.barcodes:
        barcodes = args.barcodes
        print(f"\n[Phase 1] Using specified barcodes: {barcodes}")
    else:
        print("\n[Phase 1] Auto-discovering barcodes from read prefix spectrum...")
        identifier = BarcodeIdentifier(
            barcode_length=args.barcode_length, min_distance=3
        )
        sample_stream = stream_lines(
            args.input, limit=args.sample_limit, show_progress=False
        )
        barcodes = identifier.discover_barcodes(
            read_stream=sample_stream,
            sample_limit=args.sample_limit,
            min_frequency_ratio=0.05,
        )
        print(f"  Discovered {len(barcodes)} barcodes (min distance >= 3): {barcodes}")

    if not barcodes:
        print("Error: No valid barcodes could be identified.", file=sys.stderr)
        sys.exit(1)

    # Resolve active 3' adapter sequence
    adapter_to_use = args.adapter
    if "TCGTATGCCGTCTTCTGCTTG" in args.adapter:
        idx = args.adapter.find("TCGTATGCCGTCTTCTGCTTG")
        adapter_to_use = args.adapter[idx:] + "AAAAAA"

    # demultiplexing
    print(f"\n[Phase 2] Demultiplexing library across {len(barcodes)} samples...")
    demux = LibraryDemultiplexer(
        known_barcodes=barcodes,
        adapter_sequence=adapter_to_use,
        max_mismatches=args.max_mismatches,
        min_overlap=args.min_overlap,
    )

    total_reads = 0
    start_time = time.time()
    for read in stream_lines(args.input, show_progress=True):
        total_reads += 1
        demux.process_read(read)

    elapsed_time = time.time() - start_time
    throughput = (total_reads / elapsed_time) if elapsed_time > 0 else 0.0

    print("\n[Phase 3] Exporting demultiplexing statistics and sample profiles...")

    assigned_reads = sum(s.total_reads for s in demux.samples.values())
    unassigned_reads = demux.unassigned_count
    assigned_pct = (assigned_reads / total_reads * 100) if total_reads > 0 else 0.0
    unassigned_pct = (unassigned_reads / total_reads * 100) if total_reads > 0 else 0.0

    samples_summary = {}
    for bc, stats in demux.samples.items():
        sample_pct = (stats.total_reads / total_reads * 100) if total_reads > 0 else 0.0
        samples_summary[bc] = {
            "total_reads": stats.total_reads,
            "percentage_of_library": round(sample_pct, 4),
            "unique_sequences": len(stats.sequence_frequencies),
        }

    summary_data = {
        "metadata": {
            "input_file": str(args.input),
            "barcode_length": args.barcode_length,
            "max_mismatches": args.max_mismatches,
            "adapter_sequence": adapter_to_use,
            "min_overlap": args.min_overlap,
        },
        "metrics": {
            "total_reads": total_reads,
            "assigned_reads": assigned_reads,
            "assigned_percentage": round(assigned_pct, 4),
            "unassigned_reads": unassigned_reads,
            "unassigned_percentage": round(unassigned_pct, 4),
            "elapsed_time_seconds": round(elapsed_time, 4),
            "throughput_reads_per_sec": round(throughput, 2),
        },
        "samples": samples_summary,
    }
    summary_json_path = args.output_dir / "demultiplex_summary.json"
    with open(summary_json_path, "w", encoding="utf-8") as f:
        json.dump(summary_data, f, indent=4)
    print(f"  Saved summary JSON to: {summary_json_path}")

    for bc, stats in demux.samples.items():
        sample_csv_path = args.output_dir / f"sample_lengths_{bc}.csv"
        with open(sample_csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["remaining_length", "count", "percentage"])
            for length in sorted(stats.remaining_length_distribution.keys()):
                count = stats.remaining_length_distribution[length]
                pct = (
                    (count / stats.total_reads * 100) if stats.total_reads > 0 else 0.0
                )
                writer.writerow([length, count, f"{pct:.4f}"])
        print(f"  Saved length distribution for {bc} to: {sample_csv_path}")

    top_seqs_csv_path = args.output_dir / "top_sequences_per_sample.csv"
    with open(top_seqs_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["barcode", "rank", "sequence", "length", "count", "percentage"]
        )
        for bc, stats in demux.samples.items():
            for rank, (seq, count) in enumerate(
                stats.sequence_frequencies.most_common(10), start=1
            ):
                pct = (
                    (count / stats.total_reads * 100) if stats.total_reads > 0 else 0.0
                )
                writer.writerow([bc, rank, seq, len(seq), count, f"{pct:.4f}"])
    print(f"  Saved top sequences CSV to: {top_seqs_csv_path}")

    print("\n" + "=" * 60)
    print("                TASK 5 PROCESSING SUMMARY")
    print("=" * 60)
    print(f"Total Reads Processed  : {total_reads:,}")
    print(f"Assigned Reads         : {assigned_reads:,} ({assigned_pct:.2f}%)")
    print(f"Unassigned Reads       : {unassigned_reads:,} ({unassigned_pct:.2f}%)")
    print(f"Elapsed Time           : {elapsed_time:.2f} seconds")
    print(f"Throughput             : {throughput:,.2f} reads/sec")
    print("\nPer-Sample Read Allocation:")
    for bc, stats in demux.samples.items():
        pct = (stats.total_reads / total_reads * 100) if total_reads > 0 else 0.0
        top_seq, top_cnt = (
            stats.sequence_frequencies.most_common(1)[0]
            if stats.sequence_frequencies
            else ("N/A", 0)
        )
        print(
            f"  - Sample [{bc}]: {stats.total_reads:,} reads ({pct:.2f}%) | Top: {top_seq[:18]}... ({top_cnt:,})"
        )
    print("=" * 60)


if __name__ == "__main__":
    main()
