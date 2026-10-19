import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from preprocessor.demultiplexer import BarcodeIdentifier, LibraryDemultiplexer

TEST_BARCODES = ["TAGC", "TGAG", "TCTT", "TTCA"]
TEST_ADAPTER = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"


def test_barcode_exact_match() -> None:
    """Read starting with exact barcode assigns to the correct sample."""
    demux = LibraryDemultiplexer(
        known_barcodes=TEST_BARCODES,
        adapter_sequence=TEST_ADAPTER,
        max_mismatches=1,
    )
    read = "TAGC" + "ACGTACGT" + TEST_ADAPTER[:15]
    bc, trimmed, L = demux.process_read(read)
    assert bc == "TAGC"
    assert trimmed == "ACGTACGT"
    assert L == 15
    assert demux.samples["TAGC"].total_reads == 1
    assert demux.unassigned_count == 0


def test_barcode_1_mismatch_correction() -> None:
    """Read with 1 mismatch in barcode (e.g. AAGC -> TAGC) correctly assigns."""
    demux = LibraryDemultiplexer(
        known_barcodes=TEST_BARCODES,
        adapter_sequence=TEST_ADAPTER,
        max_mismatches=1,
    )
    # Mutate position 0: 'T' -> 'A' (AAGC is distance 1 from TAGC, >= 2 from others)
    read = "AAGC" + "ACGTACGT" + TEST_ADAPTER[:10]
    bc, trimmed, L = demux.process_read(read)
    assert bc == "TAGC"
    assert trimmed == "ACGTACGT"
    assert demux.samples["TAGC"].total_reads == 1
    assert demux.unassigned_count == 0


def test_barcode_2_mismatches_rejection() -> None:
    """Read with 2 mismatches from all known barcodes is marked unassigned."""
    demux = LibraryDemultiplexer(
        known_barcodes=TEST_BARCODES,
        adapter_sequence=TEST_ADAPTER,
        max_mismatches=1,
    )
    # AAAA has Hamming distance 3 from TAGC (AAGC), 3 from TGAG, 4 from TCTT, 3 from TTCA
    read = "AAAA" + "ACGTACGT"
    bc, trimmed, L = demux.process_read(read)
    assert bc is None
    assert trimmed == read
    assert demux.unassigned_count == 1
    assert all(s.total_reads == 0 for s in demux.samples.values())


def test_double_trimming_clean_bounds() -> None:
    """Barcode excised from 5' end and adapter from 3' end without off-by-one."""
    demux = LibraryDemultiplexer(
        known_barcodes=TEST_BARCODES,
        adapter_sequence=TEST_ADAPTER,
        max_mismatches=1,
        min_overlap=10,
    )
    bio_insert = "GATTACAGATTACA"  # 14 bp
    read = "TGAG" + bio_insert + TEST_ADAPTER[:20]

    bc, trimmed, L = demux.process_read(read)
    assert bc == "TGAG"
    assert trimmed == bio_insert
    assert len(trimmed) == 14
    assert L == 20

    stats = demux.samples["TGAG"]
    assert stats.total_reads == 1
    assert stats.remaining_length_distribution[14] == 1
    assert stats.sequence_frequencies[bio_insert] == 1


def test_read_without_adapter_preserves_insert() -> None:
    """Read matching barcode but having no adapter keeps full insert length."""
    demux = LibraryDemultiplexer(
        known_barcodes=TEST_BARCODES,
        adapter_sequence=TEST_ADAPTER,
        max_mismatches=1,
        min_overlap=10,
    )
    bio_insert = "CCCCCCCCCCCCCCCC"  # 16 bp
    read = "TCTT" + bio_insert

    bc, trimmed, L = demux.process_read(read)
    assert bc == "TCTT"
    assert trimmed == bio_insert
    assert L == 0
    assert len(trimmed) == 16


def test_barcode_identifier_discovery() -> None:
    """BarcodeIdentifier extracts candidate barcodes separated by distance >= 3."""
    identifier = BarcodeIdentifier(barcode_length=4, min_distance=3)

    reads = []
    # 4 dominant barcodes, 20 reads each
    for bc in TEST_BARCODES:
        for _ in range(20):
            reads.append(bc + "ACGTACGTACGT")
    # Low frequency noise barcode
    for _ in range(2):
        reads.append("AAAAACGTACGTACGT")

    discovered = identifier.discover_barcodes(
        read_stream=iter(reads),
        sample_limit=100,
        min_frequency_ratio=0.05,
    )

    assert set(discovered) == set(TEST_BARCODES)
