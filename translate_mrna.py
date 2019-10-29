"""Translate the DNA or RNA sequence to amino acid sequence provided by the
user.

The user provide the DNA and RNA sequence with its sequence type, then the
sequence is translated at three reading frames and the molecular weight of
the protein is calculated after the the protein sequence is defined.

If no seuquence is provided the default sequence will be Prkcb transcript
variant 2 mRNA sequence which contains 671 amino acids and reading frame 2.

Author: Jia Yi Terri Shen
Date: October 2019
"""

import sys

def codon_lookup(codon):
    """Return the associating amino acid in the single-letter format.

    Args:
        codon(str): A 3-nt sequence.
    Returns:
        string: 1-nt letter.
    """

    codon_table = {
        'TTT': "F", 'TTC': "F", 'TTA': "L", 'TTG': "L", 'CTT': "L",
        'CTA': "L", 'CTG': "L", 'ATT': "I", 'ATC': "I", 'ATA': "I",
        'GTT': "V", 'GTC': "V", 'GTA': "V", 'GTG': "V", 'TCT': "S",
        'TCA': "S", 'TCG': "S", 'CCT': "P", 'CCC': "P", 'CCA': "P",
        'ACT': "T", 'ACC': "T", 'ACA': "T", 'ACG': "T", 'GCT': "A",
        'GCA': "A", 'GCG': "A", 'TAT': "Y", 'TAC': "Y", 'TAA': "-",
        'CAT': "H", 'CAC': "H", 'CAA': "Q", 'CAG': "Q", 'AAT': "N",
        'AAA': "K", 'AAG': "K", 'GAT': "D", 'GAC': "D", 'GAA': "E",
        'TGT': "C", 'TGC': "C", 'TGA': "-", 'TGG': "W", 'CGT': "R",
        'CGA': "R", 'CGG': "R", 'AGT': "S", 'AGC': "S", 'AGA': "R",
        'GGT': "G", 'GGC': "G", 'GGA': "G", 'GGG': "G", 'CTC': "L",
        'ATG': "M", 'TCC': "S", 'CCG': "P", 'GCC': "A", 'TAG': "-",
        'AAC': "N", 'GAG': "E", 'CGC': "R", 'AGG': "R",
    }

    return codon_table.get(codon, "x")

def translate_sequence(sequence, start):
    """Find the amino acid of the mRNA sequence from the starting point.

    Args:
        sequence(string): A string containing the mRNA sequence
        start(interger): The starting point of the reading frame
    Returns:
        string: A proper translated protein seq starting with M and end with -
    """

    # Declare empty string to store amino acid.
    protein = ""

    # At the reading frame, split the sequence into codons.
    for i in range(start, len(sequence), 3):

        # For each codon, find the associating amino acid.
        codon = sequence[i : i+3]
        amino_acid = codon_lookup(codon)

        # If the amino acid is "M", store it in the empty protein string.
        if not protein and amino_acid == 'M':
            protein += amino_acid
        # Started to translate when the starting amino acid is found.
        elif protein and amino_acid == '-':
            return protein
        elif protein:
            protein += amino_acid

    return


def calculate_molecular_weight(sequence):
    """Find the average molecular weight of the DNA sequence.

    Args:
        sequence(string): A string containing the protein sequence

    Returns:
        numeric: A float real number with two decimal pointsthat represent the
        average molecualar weight in kilodaltons of the protein
    """

    # Define the average weight of one amino acid in kilos
    avg_weight_kilo = 110 / 1000

    # Calcualte the sequence weight
    protein_avg_weight = float(len(sequence) * avg_weight_kilo)

    return protein_avg_weight

def main():
    """The main function of the script."""

    # Define a variable to store the Prkcb.transcript vairant2 mRNA sequence
    mrna_seq = """
GCGGCCCTGCGGTCCCCGGGCGGCAGCAGCGGCCGCCTAGTCCCGCGCCTCTCCGGGCTTACAGCCCCGC
GGTCCCGCCGCCCCGGGGCCGCCACCTCTCGGGGCTCCCCCCAGTCCCCGCGCGCGCAAGATGGCTGACC
CGGCTGCGGGGCCGCCGCCGAGCGAGGGCGAGGAGAGCACGGTGCGCTTCGCCCGCAAAGGCGCCCTCCG
GCAGAAGAACGTGCACGAGGTGAAGAACCACAAATTCACCGCCCGCTTCTTCAAGCAGCCCACCTTCTGC
AGCCACTGCACCGACTTCATTTGGGGCTTAGGTTTGCAGGGATTCCAGAGTCAGGTCTGCTGCTTTGTTG
TACACAAGCGCTGCCATGAATTCGTCACGTTCTCCTGCCCTGGTGCAGACAAGGGCCCGGCCTCTGATGA
CCCACGGAGCAAACACAAGTTTAAGATCCACACCTACTCCAGCCCTACCTTCTGTGACCACTGTGGATCA
CTGCTGTATGGGCTCATCCACCAGGGGATGAAATGCGACACCTGTATGATGAATGTCCACAAGCGCTGCG
TGATGAACGTCCCCAGCCTCTGTGGCACCGACCACACAGAACGCCGTGGCCGCATCTACATCCAGGCCCA
CATCGACAGGGAGGTCCTCATCGTTGTTGTAAGAGATGCTAAAAATCTGGTACCTATGGACCCCAACGGC
TTGTCAGATCCCTACGTAAAACTGAAACTGATCCCTGATCCCAAAAGTGAGAGCAAGCAGAAGACCAAGA
CTATCAAATGCTCCCTCAACCCGGAGTGGAACGAAACCTTCAGATTTCAGCTGAAGGAATCAGACAAAGA
CAGAAGACTGTCCGTAGAGATCTGGGATTGGGACCTGACCAGCAGGAATGACTTCATGGGATCTCTGTCG
TTTGGGATTTCAGAACTACAGAAAGCCGGAGTGGATGGCTGGTTCAAGTTACTAAGCCAGGAAGAAGGCG
AGTACTTTAATGTGCCGGTGCCGCCGGAAGGAAGCGAGGGCAATGAAGAGCTGCGGCAGAAGTTTGAGAG
AGCCAAGATTGGCCAAGGTACCAAGGCTCCAGAAGAAAAGACAGCGAACACTATATCCAAATTTGACAAC
AATGGCAACAGGGACCGGATGAAACTGACCGATTTTAACTTCCTGATGGTGCTGGGGAAAGGCAGCTTTG
GCAAGGTCATGCTCTCAGAGCGGAAGGGTACAGATGAACTCTATGCCGTGAAGATCCTGAAGAAAGATGT
GGTGATCCAAGATGACGATGTGGAGTGCACAATGGTGGAGAAGAGGGTGCTGGCCCTGCCTGGGAAGCCC
CCATTCCTGACTCAGCTCCATTCCTGCTTCCAGACCATGGACCGCCTCTACTTTGTGATGGAGTATGTGA
ACGGGGGCGACCTCATGTACCACATCCAACAAGTTGGCCGTTTCAAGGAGCCCCATGCTGTATTTTACGC
TGCAGAGATTGCCATCGGTCTTTTCTTCTTGCAGAGCAAGGGCATCATTTACCGTGACCTGAAACTTGAC
AACGTGATGCTGGATTCCGAGGGGCACATCAAAATCGCTGACTTTGGCATGTGTAAAGAGAATATCTGGG
ATGGGGTGACAACCAAGACATTCTGTGGCACTCCAGACTACATTGCCCCAGAGATCATTGCTTATCAGCC
CTACGGGAAGTCTGTGGACTGGTGGGCGTTTGGAGTCCTGCTGTATGAAATGTTGGCTGGCCAGGCACCT
TTTGAAGGGGAGGATGAGGATGAACTCTTCCAGTCAATCATGGAGCACAACGTGGCGTATCCCAAGTCCA
TGTCTAAGGAAGCTGTGGCAATCTGCAAAGGGCTAATGACCAAACACCCAGGCAAGCGCCTGGGTTGTGG
GCCTGAAGGGGAACGAGACATTAAGGAGCATGCATTTTTCCGGTATATCGACTGGGAGAAACTCGAACGC
AAGGAGATTCAGCCACCTTATAAACCAAAAGCTAGAGACAAGCGAGACACCTCCAACTTCGACAAAGAGT
TCACCAGGCAGCCTGTGGAACTGACTCCCACTGACAAACTCTTCATCATGAACTTGGACCAAAATGAATT
TGCTGGCTTCTCGTATACTAACCCAGAGTTTGTCATTAATGTGTAGGTGAATGCAGATTCCATCGCTGAG
CCTGTGTGTAAGGCTGCAGGCTGAATGTCTATTATCAATTCCAGTCTTCCAGGATTCATGGTGCCTCTGT
TGGCATCCGTCATGTGGAGAGCTTGTCTTAGAGGGCTTTTCTTTGTATGTATAGCTTGCTAGTTTGTTTT
CTACATTTCAAAATGTTTAGTTTAGAATAAGTGCATTGCCCACTGATAGAGGTACAATTTTCCAGACTTC
CAGAAACTCATCCAATGAACCAACAGTGTCAAAACTTAACTGTGTCCGATACCAAAATGCTTCAGTATTT
GTAATTTTTAAAGTCAGATGCTGATGTTCCTGGTCAAAGTTTTTACAGTTACTCTCGAATATCTCCTTTG
AATGCTACCTAAGCATGACCGGTATTTTTAAAAGTTGTGAGTAAGCTTTGCAGTTACTGTGAACTCTTGT
CTCTTGGAGGAAACTTTTTGTTTAAGAATTGGTATGATTAAACTGAATTCATATATGCAAAAAAAAAAAA
AAAAAA """

    # Remove any new line or spaces in between.
    mrna_seq = mrna_seq.replace("\n", "").replace("\r", "")

    mrna = "mrna"
    seq_type = "sequence_type"
    valid = ['DNA', 'RNA']

    # Set the default seq and seq type when there is no user input.
    if len(sys.argv) < 2:
        mrna = mrna_seq
        seq_type = "DNA"
    # Quit the function if the seq is provided but seq type is not.
    elif len(sys.argv) == 2:
        sys.exit("Please indicate if the sequence is 'DNA' or 'RNA'.")
    elif len(sys.argv) == 3:
        mrna = sys.argv[1]
        seq_type = sys.argv[2]
        # Quit the function if the sequence type is not RNA and DNA.
        if seq_type not in valid:
            sys.exit("Sorry, " + seq_type + " is not recognized. " +
                     "Please indicate if the sequence is 'DNA' or 'RNA'.")
        # Translate the DNA sequence if it is RNA
        elif seq_type == "RNA":
            mrna = mrna.replace('U', 'T')
        else:
            pass

    # Printing a summary report to the user.
    for frame in range(3):

        # Define output variables and stores it into a list
        protein = translate_sequence(mrna, frame)
        weight = calculate_molecular_weight(protein) if protein else 'N/A'

        frame = f"Reading frame: {frame + 1}"
        sequence = f"Sequence: {protein if protein else 'None'}"
        length = f"Length: {len(protein) if protein else 'N/A'}"
        molecular_weight = f"Weight (kilodaltons): {weight}"
        dash = "-" * 4
        output = [protein, frame, sequence, length, molecular_weight, dash]
        print(*output[1:], sep="\n")

if __name__ == "__main__":
    main()
