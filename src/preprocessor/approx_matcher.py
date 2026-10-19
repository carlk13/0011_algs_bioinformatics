import time
from collections import Counter
from collections.abc import Generator
from math import floor

import numpy as np

from preprocessor.error_analyzer import AlignmentRecord, ErrorProfiler


class HammingAdapterTrimmer:
    """
    Trims adapters allowing for substitution errors (Hamming distance)
    up to a specified maximum error rate (e.g., 10% or 25%).
    """

    def __init__(
        self, adapter_sequence: str, max_error_rate: float = 0.10, min_overlap: int = 1
    ):
        self.adapter_seq = adapter_sequence
        self.n = len(adapter_sequence)
        self.max_error_rate = max_error_rate
        self.min_overlap = min_overlap

    def trim_read(self, read: str) -> tuple[str, int, int, bool]:
        """
        Check potential suffix-prefix overlaps of the read against the adapter.
        Returns: (trimmed_seq, overlap_length, mismatch_count, is_trimmed)
        """
        m = len(read)

        # Determine the maximum possible overlap
        max_possible_L = min(m, self.n)

        for L in range(max_possible_L, self.min_overlap - 1, -1):
            # max allowed mismatches for this specific overlap length L
            allowed_mismatches = floor(L * self.max_error_rate)
            mismatch_count = 0
            is_valid = True

            # Compare read suffix against adapter prefix
            for i in range(L):
                read_char = read[m - L + i]
                adapter_char = self.adapter_seq[i]

                if read_char != adapter_char:
                    mismatch_count = mismatch_count + 1

                    if mismatch_count > allowed_mismatches:
                        is_valid = False
                        break

            if is_valid:
                trimmed_seq = read[0 : m - L]
                return (trimmed_seq, L, mismatch_count, True)

        # If no overlap length met the criteria
        return (read, 0, 0, False)

    def process_stream(
        self,
        read_stream: Generator[str],
        error_profiler: ErrorProfiler | None = None,
    ) -> tuple:
        total_reads = 0
        trimmed_reads_count = 0

        longest_overlap_distribution = Counter()
        remaining_length_distribution = Counter()
        mismatch_count_distribution = Counter()

        # Start timer
        start_time = time.time()

        for read in read_stream:
            total_reads += 1
            m = len(read)

            trimmed_seq, L, mismatch_count, is_trimmed = self.trim_read(read)

            if is_trimmed:
                trimmed_reads_count += 1
                if error_profiler is not None:
                    # Map mismatch positions to read cycle coordinates (m - L + i)
                    error_positions = [
                        m - L + i
                        for i in range(L)
                        if read[m - L + i] != self.adapter_seq[i]
                    ]
                    record = AlignmentRecord(
                        read_length=m,
                        overlap_length=L,
                        num_substitutions=mismatch_count,
                        num_insertions=0,
                        num_deletions=0,
                        error_positions=error_positions,
                        is_trimmed=True,
                    )
                    error_profiler.update(record)

            # Update distributions
            longest_overlap_distribution[L] += 1
            remaining_length_distribution[len(trimmed_seq)] += 1
            if is_trimmed:
                mismatch_count_distribution[mismatch_count] += 1

        # Stop timer
        end_time = time.time()
        elapsed_time = end_time - start_time

        throughput = total_reads / elapsed_time if elapsed_time > 0 else 0

        return (
            total_reads,
            trimmed_reads_count,
            longest_overlap_distribution,
            remaining_length_distribution,
            mismatch_count_distribution,
            elapsed_time,
            throughput,
        )


class SemiGlobalDPTrimmer:
    """
    Match reward: +1
    Mismatch penalty: -2
    Gap penalty: -3


    If your penalties are too harsh, valid mutated adapters will be missed;
    if they are too lenient, random genomic sequence will accidentally match the adapter.

    """

    def __init__(
        self, adapter_sequence, match=1, mismatch=-2, gap=-3, min_score_threshold=10
    ):
        self.adapter = adapter_sequence
        self.n = len(adapter_sequence)
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.min_threshold = min_score_threshold

    def trim_read(self, read: str, detailed: bool = False) -> tuple:
        m = len(read)

        # DP table of dimensions (m + 1) x (n + 1)
        DP = np.zeros((m + 1, self.n + 1))

        # INITIALIZATION FOR SEMI-GLOBAL (Fitting alignment):
        # Row 0 is all 0s -> Allows the adapter to start aligning anywhere
        # at the tail of the read without penalty.

        # FILL THE MATRIX:
        for i in range(1, m + 1):
            for j in range(1, self.n + 1):
                # Calculate score for substitution/match
                score = (
                    self.match if read[i - 1] == self.adapter[j - 1] else self.mismatch
                )
                match_or_sub = DP[i - 1, j - 1] + score

                # Calculate score for deletion in read / insertion in adapter
                del_in_read = DP[i - 1, j] + self.gap

                # Calculate score for insertion in read / deletion in adapter
                ins_in_read = DP[i, j - 1] + self.gap

                # Take the maximum score (maximizing score rather than minimizing distance)
                DP[i, j] = max(0, match_or_sub, del_in_read, ins_in_read)

        # FIND THE BEST MATCH IN THE LAST ROW:
        best_j = int(np.argmax(DP[m, :]))
        max_score = float(DP[m, best_j])

        # Check if the score meets our quality threshold
        if max_score >= self.min_threshold and best_j > 0:
            i_start, subs, ins, dels, err_pos = self.find_alignment_start(
                DP, read, m, best_j
            )
            trimmed_seq = read[:i_start]
            overlap_length = m - i_start
            if detailed:
                return (
                    trimmed_seq,
                    overlap_length,
                    max_score,
                    subs,
                    ins,
                    dels,
                    err_pos,
                    True,
                )
            return (trimmed_seq, overlap_length, max_score, True)
        else:
            if detailed:
                return (read, 0, 0.0, 0, 0, 0, [], False)
            return (read, 0, 0.0, False)

    def find_alignment_start(self, DP: np.ndarray, read: str, m: int, best_j: int):
        i_start = m
        j = best_j
        num_subs = 0
        num_ins = 0
        num_dels = 0
        error_positions = []

        # Walk backward from the bottom-right endpoint until we hit the top of the adapter (j = 0)
        while j > 0 and i_start > 0:
            current_score = DP[i_start, j]
            is_match = read[i_start - 1] == self.adapter[j - 1]
            score_match = self.match if is_match else self.mismatch

            # Check which neighboring cell generated the current score
            if (
                j > 0
                and i_start > 0
                and current_score == DP[i_start - 1, j - 1] + score_match
            ):
                if not is_match:
                    num_subs += 1
                    error_positions.append(i_start - 1)
                i_start -= 1
                j -= 1

            elif i_start > 0 and current_score == DP[i_start - 1, j] + self.gap:
                num_ins += 1
                error_positions.append(i_start - 1)
                i_start -= 1

            elif j > 0 and current_score == DP[i_start, j - 1] + self.gap:
                num_dels += 1
                error_positions.append(i_start)
                j -= 1

            else:
                # Safety break if we hit a boundary or score of 0
                break

        return i_start, num_subs, num_ins, num_dels, error_positions

    def process_stream(
        self,
        read_stream: Generator[str],
        limit: int | None = None,
        error_profiler: ErrorProfiler | None = None,
    ) -> tuple:
        total_reads = 0
        trimmed_reads_count = 0
        longest_overlap_distribution = Counter()
        remaining_length_distribution = Counter()
        total_subs = 0
        total_ins = 0
        total_dels = 0

        start_time = time.time()
        for read in read_stream:
            total_reads += 1
            (
                trimmed_seq,
                L,
                score,
                subs,
                ins,
                dels,
                err_pos,
                is_trimmed,
            ) = self.trim_read(read, detailed=True)

            if is_trimmed:
                trimmed_reads_count += 1
                total_subs += subs
                total_ins += ins
                total_dels += dels
                if error_profiler is not None:
                    error_profiler.update(
                        AlignmentRecord(
                            read_length=len(read),
                            overlap_length=L,
                            num_substitutions=subs,
                            num_insertions=ins,
                            num_deletions=dels,
                            error_positions=err_pos,
                            is_trimmed=True,
                        )
                    )

            longest_overlap_distribution[L] += 1
            remaining_length_distribution[len(trimmed_seq)] += 1

            if limit is not None and total_reads >= limit:
                break

        end_time = time.time()
        elapsed_time = end_time - start_time
        throughput = total_reads / elapsed_time if elapsed_time > 0 else 0

        return (
            total_reads,
            trimmed_reads_count,
            longest_overlap_distribution,
            remaining_length_distribution,
            total_subs,
            total_ins,
            total_dels,
            elapsed_time,
            throughput,
        )
