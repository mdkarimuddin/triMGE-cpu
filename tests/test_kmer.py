from trimgecpu.fasta import clean_sequence, parse_fasta
from trimgecpu.kmer import compute_kmer_frequencies, kmer_vocabulary


def test_clean_sequence():
    assert clean_sequence("acgtnxyz") == "ACGT"


def test_parse_fasta_filters_short_sequences(tmp_path):
    fasta = tmp_path / "tiny.fasta"
    fasta.write_text(">short\nACGT\n>long\n" + ("ACGT" * 30) + "\n")
    records = parse_fasta(fasta)
    assert len(records) == 1
    assert records[0].sequence_id == "long"


def test_kmer_dimension():
    vec = compute_kmer_frequencies("ACGTACGTACGT", k=5)
    assert vec.shape == (len(kmer_vocabulary(5)[0]),)
    assert abs(vec.sum() - 1.0) < 1e-9
