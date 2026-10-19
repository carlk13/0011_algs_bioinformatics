import time
from collections import Counter
from collections.abc import Generator

from preprocessor.kmp import compute_kmp_pi, kmp_suffix_prefix


class ExactAdapterTrimmer:
    def __init__(self, adapter_sequence: str, min_overlap: int = 1):
        self.adapter = adapter_sequence
        self.n = len(adapter_sequence)
        self.min_overlap = min_overlap

        # // Precompute the KMP failure table (pi table) for the adapter once
        self.pi = compute_kmp_pi(self.adapter)

    def trim_read(self, read: str) -> tuple[str, int, bool]:
        m = len(read)

        # find the longest suffix-prefix overlap length L using KMP
        L = kmp_suffix_prefix(read, self.adapter, self.pi)

        # check if the overlap meets the minimum threshold
        if L >= self.min_overlap:
            trimmed_seq = read[0 : m - L]
            return (trimmed_seq, L, True)

        else:
            return (read, 0, False)

    def process_stream(self, read_stream: Generator[str]) -> tuple:
        total_reads = 0
        trimmed_reads_count = 0

        longest_overlap_distribution: Counter[int] = Counter()
        remaining_length_distribution: Counter[int] = Counter()

        # Start timer
        start_time = time.time()

        for read in read_stream:
            total_reads += 1
            m = len(read)

            # Trim the read
            trimmed_seq, L, is_trimmed = self.trim_read(read)

            if is_trimmed:
                trimmed_reads_count += 1

            # Update distributions
            longest_overlap_distribution[L] += 1
            remaining_length_distribution[len(trimmed_seq)] += 1

        # Stop timer
        end_time = time.time()
        elapsed_time = end_time - start_time

        # Calculate throughput
        throughput = total_reads / elapsed_time if elapsed_time > 0 else 0

        return (
            total_reads,
            trimmed_reads_count,
            longest_overlap_distribution,
            remaining_length_distribution,
            elapsed_time,
            throughput,
        )
