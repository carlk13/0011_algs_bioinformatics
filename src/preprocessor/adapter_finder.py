from collections import Counter
from collections.abc import Generator


class UnknownAdapterFinder:
    """
    Infers unknown 3' adapter sequences from high-throughput reads using
    suffix-frequency spectrum analysis and greedy consensus extension.
    """

    def __init__(self, seed_length: int = 12, min_frequency_threshold: int = 1000):
        self.k = seed_length
        self.min_freq = min_frequency_threshold

    def find_frequent_suffixes(
        self, read_stream: Generator[str], sample_limit: int | None = 200000
    ) -> list[tuple[str, int]]:
        """Extracts top frequent terminal suffix k-mers."""
        suffix_counts: Counter[str] = Counter()
        for count, read in enumerate(read_stream):
            if sample_limit and count >= sample_limit:
                break
            if len(read) >= self.k:
                suffix_counts[read[-self.k :]] += 1
        return suffix_counts.most_common(20)

    def assemble_adapter_consensus(
        self,
        seed: str,
        read_stream: Generator[str],
        max_length: int = 65,
        sample_limit: int | None = 200000,
    ) -> str:
        """
        Extends the seed rightward by analyzing positional successor frequencies
        in reads containing the current consensus prefix.
        """
        # Load sample into memory for multi-round extension
        reads: list[str] = []
        for count, r in enumerate(read_stream):
            if sample_limit and count >= sample_limit:
                break
            reads.append(r)

        consensus = list(seed)

        # Leftward extension
        while len(consensus) < max_length:
            query = "".join(consensus[:12])
            prev_base_counter: Counter[str] = Counter()
            for r in reads:
                idx = r.find(query)
                if idx > 0:
                    prev_base_counter[r[idx - 1]] += 1

            if not prev_base_counter:
                break

            total_obs = sum(prev_base_counter.values())
            best_prev_base, count = prev_base_counter.most_common(1)[0]
            # Stop when consensus drops
            if count / total_obs < 0.70 or count < 50:
                break

            consensus.insert(0, best_prev_base)

        #  Rightward extension (towards 3' end)
        while len(consensus) < max_length:
            curr_str = "".join(consensus)
            query = curr_str[-12:]
            next_base_counter: Counter[str] = Counter()

            for r in reads:
                idx = r.find(query)
                if idx != -1 and idx + len(query) < len(r):
                    next_base_counter[r[idx + len(query)]] += 1

            if not next_base_counter:
                break

            total_obs = sum(next_base_counter.values())
            best_next_base, count = next_base_counter.most_common(1)[0]
            if count / total_obs < 0.70 or count < 50:
                break

            consensus.append(best_next_base)

        return "".join(consensus)
