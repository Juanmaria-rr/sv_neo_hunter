from Bio.Seq import Seq


def set_aa_seq(svfusion):
    """
    :function: translate the nucelotide sequence in svfusion to amino acid sequence
    :param: svfusion: a SVFusion class
    :return: the corresponding amino acid sequence
    :NOTE: if there is a stop codon, the process will stop before translating all nucleotides
    """
    dna_seq = Seq(trim_to_3x(svfusion.nt_sequence))
    mrna_seq = dna_seq.transcribe()
    aa_seq = str(mrna_seq.translate(to_stop=True))
    # PATCHED (see vendor/patches/004-non-aug-initiation.md).
    # A transcript annotated with a non-AUG start (CUG, GUG, ...) is still
    # initiated with methionine; Ensembl's protein for it begins with M. Read
    # literally, the first codon gives L/V/..., so the N-terminal windows differ
    # from wild type by one residue and pass the wild-type subtraction as
    # spurious neopeptides. When the fusion begins at the 5' transcript's own
    # start codon, its first residue is methionine.
    if aa_seq and svfusion.starts_at_native_start and not aa_seq.startswith('M'):
        aa_seq = 'M' + aa_seq[1:]
    return aa_seq


def set_nt_seq(svfusion):
    """
    :function: given a svfusion, paste the nucleotides in cc_1, cc_2 and 3utr
    :param svfusion: a SVFusion class
    :return: the nucleotide sequence from start codon to the end of 3utr
    """
    return svfusion.nt_sequence_cds + svfusion.nt_sequence_3utr


def trim_to_3x(nt_sequence):
    """
    :function: trim a sequence to be divisible by 3
    :return: trimmed sequence
    """
    remainder = len(nt_sequence) % 3
    if remainder:
        return nt_sequence[: -remainder]
    else:
        return nt_sequence


def generate_neoepitopes(svfusion, window_range):
    """
    :function: cut the WT and MUT protein sequence by window_range, then get the MUT specific peptides
    :param svfusion: a SVFusion class
    :param window_range: a list, specifying the range of window size, e.g. [8,9,10,11]
    :return: a list of neopeptides for svfusion
    """
    mut_peptides = []
    for window in window_range:
        mut_peptides = mut_peptides + cut_sequence(svfusion.aa_sequence, window)
    wt_peptides = []
    for window in window_range:
        wt_peptides = wt_peptides + cut_sequence(svfusion.cc_1.transcript.protein_sequence, window)
        wt_peptides = wt_peptides + cut_sequence(svfusion.cc_2.transcript.protein_sequence, window)
    return list(set(mut_peptides)-set(wt_peptides))


def cut_sequence(aa_sequence, window):
    """
    :param aa_sequence: amino acid sequence
    :param window: the size of window, window should be smaller than len(aa_sequence)
    :return: all possible slices with length = window in this sequence
    """
    if len(aa_sequence) < window:
        return []
    else:
        return [aa_sequence[i: i+window] for i in range(len(aa_sequence)-window+1)]
