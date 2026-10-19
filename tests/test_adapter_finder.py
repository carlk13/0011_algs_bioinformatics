import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

import random

from preprocessor.adapter_finder import UnknownAdapterFinder


def test_find_frequent_suffixes_synthetic() -> None:
    """find_frequent_suffixes identifies the dominant terminal k-mer against random noise."""
    finder = UnknownAdapterFinder(seed_length=12)
    target_seed = "TGGAATTCTCGA"

    # Generate 150 reads ending with target_seed and 50 ending with random k-mers
    random.seed(42)
    bases = ["A", "C", "G", "T"]
    reads = []
    for _ in range(150):
        prefix = "".join(random.choices(bases, k=25))
        reads.append(prefix + target_seed)
    for _ in range(50):
        reads.append("".join(random.choices(bases, k=37)))

    suffixes = finder.find_frequent_suffixes(iter(reads), sample_limit=200)
    assert len(suffixes) > 0
    top_kmer, count = suffixes[0]
    assert top_kmer == target_seed
    assert count == 150


def test_assemble_adapter_consensus_recovery() -> None:
    """assemble_adapter_consensus reconstructs the full adapter from a terminal seed."""
    finder = UnknownAdapterFinder(seed_length=12)
    known_adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACACTGATCTCGTATGCCGTCTTCTGCTTGAA"
    seed = known_adapter[-12:]  # terminal 12-mer

    random.seed(42)
    bases = ["A", "C", "G", "T"]
    reads = []
    # Reads have diverse biological prefixes (20 bp) followed by the known adapter
    for _ in range(100):
        bio_insert = "".join(random.choices(bases, k=20))
        reads.append(bio_insert + known_adapter)

    consensus = finder.assemble_adapter_consensus(
        seed=seed,
        read_stream=iter(reads),
        max_length=len(known_adapter),
        sample_limit=100,
    )

    assert consensus == known_adapter


def test_assemble_adapter_consensus_stops_at_junction() -> None:
    """Assembler halts leftward extension when predecessor bases become non-consensus."""
    finder = UnknownAdapterFinder(seed_length=10)
    known_adapter = "TCGTATGCCGTCTTCTGCTTGAA"
    seed = known_adapter[-10:]

    random.seed(123)
    bases = ["A", "C", "G", "T"]
    reads = []
    # Diverse random inserts preceding the adapter
    for _ in range(80):
        bio_insert = "".join(random.choices(bases, k=22))
        reads.append(bio_insert + known_adapter)

    consensus = finder.assemble_adapter_consensus(
        seed=seed,
        read_stream=iter(reads),
        max_length=50,
        sample_limit=80,
    )

    # Consensus must start with the adapter, not absorbing random insert
    assert consensus.startswith("TCGTATGCCG")
    assert consensus == known_adapter


def test_sample_limit_respected() -> None:
    """Processing terminates once sample_limit is reached."""
    finder = UnknownAdapterFinder(seed_length=8)

    def infinite_stream():
        while True:
            yield "ACGTACGTACGTACGT"

    suffixes = finder.find_frequent_suffixes(infinite_stream(), sample_limit=50)
    assert len(suffixes) == 1
    assert suffixes[0][1] == 50


def test_empty_stream() -> None:
    """Empty generator returns empty results without error."""
    finder = UnknownAdapterFinder(seed_length=10)
    assert finder.find_frequent_suffixes(iter([])) == []
    assert finder.assemble_adapter_consensus("ACGT", iter([])) == "ACGT"
