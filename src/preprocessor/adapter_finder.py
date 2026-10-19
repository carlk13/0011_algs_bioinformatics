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
        suffix_counts = Counter()
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
        reads = []
        for count, r in enumerate(read_stream):
            if sample_limit and count >= sample_limit:
                break
            reads.append(r)

        consensus = list(seed)
        while len(consensus) < max_length:
            curr_str = "".join(consensus)
            # Use the last 12-15 bp of current consensus as search query
            query = curr_str[-12:]
            next_base_counter = Counter()

            for r in reads:
                idx = r.find(query)
                if idx != -1 and idx + len(query) < len(r):
                    next_base_counter[r[idx + len(query)]] += 1

            if not next_base_counter or next_base_counter.most_common(1)[0][1] < 50:
                break

            best_next_base, _ = next_base_counter.most_common(1)[0]
            consensus.append(best_next_base)

        return "".join(consensus)
