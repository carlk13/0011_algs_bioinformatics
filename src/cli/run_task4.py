import argparse
import csv
import json
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from preprocessor.adapter_finder import UnknownAdapterFinder
from preprocessor.exact_matcher import ExactAdapterTrimmer
from preprocessor.io_utils import stream_lines


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Task 4: Unsupervised 3' Adapter Discovery and Trimming"
    )

    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/tdt4287-unknown-adapter.txt"),
        help="Path to input .txt or .gz file containing unknown adapter reads",
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=10,
        help="Minimum overlap length required for trimming",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/day3/task4"),
        help="Directory to save output stats and reports",
    )

    parser.add_argument(
        "--sample-limit",
        type=int,
        default=100000,
        help="Sample size for seed extraction and consensus assembly (default: 100000)",
    )

    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("      TASK 4: UNKNOWN ADAPTER INFERENCE & TRIMMING")
    print("=" * 60)
    print(f"Input File     : {args.input}")
    print(f"Min Overlap    : {args.min_overlap}")
    print(f"Sample Limit   : {args.sample_limit:,} reads")
    print(f"Output Dir     : {args.output_dir}")
    print("-" * 60)

    # Find frequent suffices from a sample
    print("\n[Phase 1] Sampling reads for terminal suffix frequency spectrum...")
    finder = UnknownAdapterFinder(seed_length=12)
    sample_stream = stream_lines(
        args.input, limit=args.sample_limit, show_progress=False
    )
    frequent_suffixes = finder.find_frequent_suffixes(
        sample_stream, sample_limit=args.sample_limit
    )

    top_seed, seed_count = frequent_suffixes[0]
    seed_pct = (seed_count / args.sample_limit * 100) if args.sample_limit > 0 else 0
    print(f"  Top terminal 12-mer: {top_seed} ({seed_count:,} reads, {seed_pct:.2f}%)")

    # Save frequent_suffixes.csv
    freq_suffixes_csv_path = args.output_dir / "frequent_suffixes.csv"
    with open(freq_suffixes_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "suffix", "count", "percentage"])
        for rank, (suffix, count) in enumerate(frequent_suffixes[:50], start=1):
            pct = (count / args.sample_limit * 100) if args.sample_limit > 0 else 0
            writer.writerow([rank, suffix, count, f"{pct:.4f}"])
    print(f"  Saved frequent suffixes to: {freq_suffixes_csv_path}")

    # Assemble consensus adapter
    print("\n[Phase 2] Assembling full adapter consensus sequence...")
    assembly_stream = stream_lines(
        args.input, limit=args.sample_limit, show_progress=False
    )
    inferred_adapter = finder.assemble_adapter_consensus(
        seed=top_seed,
        read_stream=assembly_stream,
        max_length=65,
        sample_limit=args.sample_limit,
    )

    print(f"  Inferred Adapter: {inferred_adapter}")
    print(f"  Length          : {len(inferred_adapter)} bp")

    # Save inferred_adapter.json
    inferred_adapter_data = {
        "metadata": {
            "input_file": str(args.input),
            "sample_limit": args.sample_limit,
            "min_overlap": args.min_overlap,
        },
        "inference": {
            "top_seed_kmer": top_seed,
            "seed_observations": seed_count,
            "seed_percentage_in_sample": round(seed_pct, 4),
            "inferred_adapter_sequence": inferred_adapter,
            "inferred_adapter_length": len(inferred_adapter),
        },
    }
    inferred_json_path = args.output_dir / "inferred_adapter.json"
    with open(inferred_json_path, "w", encoding="utf-8") as f:
        json.dump(inferred_adapter_data, f, indent=4)
    print(f"  Saved inferred adapter JSON to: {inferred_json_path}")

    # full libraray trimming with inferred adapter
    print(
        f"\n[Phase 3] Trimming full library using inferred adapter (min_overlap={args.min_overlap})..."
    )
    trimmer = ExactAdapterTrimmer(
        adapter_sequence=inferred_adapter, min_overlap=args.min_overlap
    )

    total_reads = 0
    trimmed_reads_count = 0
    overlap_dist: Counter[int] = Counter()
    length_dist: Counter[int] = Counter()
    sequence_freqs: Counter[str] = Counter()

    start_time = time.time()
    for read in stream_lines(args.input, show_progress=True):
        total_reads += 1
        trimmed_seq, L, is_trimmed = trimmer.trim_read(read)

        if is_trimmed:
            trimmed_reads_count += 1
            sequence_freqs[trimmed_seq] += 1

        overlap_dist[L] += 1
        length_dist[len(trimmed_seq)] += 1

    elapsed_time = time.time() - start_time
    throughput = (total_reads / elapsed_time) if elapsed_time > 0 else 0.0

    untrimmed_reads_count = total_reads - trimmed_reads_count
    trimmed_pct = (trimmed_reads_count / total_reads * 100) if total_reads > 0 else 0.0
    untrimmed_pct = (
        (untrimmed_reads_count / total_reads * 100) if total_reads > 0 else 0.0
    )

    rem_csv_path = args.output_dir / "remaining_lengths_task4.csv"
    with open(rem_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["remaining_length", "count", "percentage"])
        for length in sorted(length_dist.keys()):
            count = length_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0.0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"\n  Saved remaining lengths CSV to: {rem_csv_path}")

    # 2. top_sequences_task4.csv
    top_seqs_csv_path = args.output_dir / "top_sequences_task4.csv"
    with open(top_seqs_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["rank", "sequence", "length", "count", "percentage_of_trimmed"]
        )
        for rank, (seq, count) in enumerate(sequence_freqs.most_common(50), start=1):
            pct = (
                (count / trimmed_reads_count * 100) if trimmed_reads_count > 0 else 0.0
            )
            writer.writerow([rank, seq, len(seq), count, f"{pct:.4f}"])
    print(f"  Saved top sequences CSV to: {top_seqs_csv_path}")

    # 3. overlap_lengths_task4.csv
    overlap_csv_path = args.output_dir / "overlap_lengths_task4.csv"
    with open(overlap_csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["overlap_length", "count", "percentage"])
        for length in sorted(overlap_dist.keys()):
            count = overlap_dist[length]
            pct = (count / total_reads * 100) if total_reads > 0 else 0.0
            writer.writerow([length, count, f"{pct:.4f}"])
    print(f"  Saved overlap lengths CSV to: {overlap_csv_path}")

    print("\n" + "=" * 60)
    print("                TASK 4 PROCESSING SUMMARY")
    print("=" * 60)
    print(f"Total Reads Processed  : {total_reads:,}")
    print(f"Trimmed Reads          : {trimmed_reads_count:,} ({trimmed_pct:.2f}%)")
    print(f"Untrimmed Reads        : {untrimmed_reads_count:,} ({untrimmed_pct:.2f}%)")
    print(f"Elapsed Time           : {elapsed_time:.2f} seconds")
    print(f"Throughput             : {throughput:,.2f} reads/sec")
    print(f"Inferred Adapter       : {inferred_adapter}")
    print("\nTop 5 Remaining Sequence Lengths:")
    for length, count in length_dist.most_common(5):
        pct = (count / total_reads * 100) if total_reads > 0 else 0.0
        print(f"  - Length {length:2d} bp: {count:,} reads ({pct:.2f}%)")
    print("\nTop 5 Biological Sequences (microRNAs):")
    for rank, (seq, count) in enumerate(sequence_freqs.most_common(5), start=1):
        pct = (count / trimmed_reads_count * 100) if trimmed_reads_count > 0 else 0.0
        print(f"  {rank}. {seq} ({len(seq)} bp) : {count:,} reads ({pct:.2f}%)")
    print("=" * 60)


if __name__ == "__main__":
    main()
