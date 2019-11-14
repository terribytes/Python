"""Parse trinity output of differentially expressed genes with its Swissport
ID and the list of associated GO terms and corresponding GO terms descriptions.

This script is the refactor code of diff_exp_annotations-messy.py, which
makes annotation report.tsv file for the trinity output of differentially
expressed genes under 4 stress conditions with their association protein name
, GO IDs, and GO description.

This file contains below functions:
    * transcript_protein_dict - accept blast .outfmt6 file return a dictionary
                                that maps transcript ID to SwissPort ID.
    * gene_go_dict - accept GO annotation .gaf file return a dictionary that
                     maps SwissPort ID to its lists of unique GO IDs.
    * go_name_dict - accept GO terms .obo file return a dictionary that maps
                     GO ID to its description.

Author: Jia Yi Terri Shen
Date: November 2019
"""

import re


def transcript_protein_dict(blast_filename):
    """Load transcript IDs of query sequence (qseqid) and SwissProt ID
    of subject sequence (sseqid) without version number to the dictionary.

    Args:
        filename (str): File path to the blast output .outfmt6 file.

    Returns:
        dictionary: A dictionary matching the transctipt(key) to sp_id
                    (value) that is the first corresponding BLAST hit with
                    identity over 99.
    """

    transcript_to_protein = {}
    if not blast_filename.endswith(".outfmt6"):
        print("Please provide valid BLAST .outfmt6 file.")

    else:
        blast_file = open(blast_filename)
        for line in blast_file:
            qseqid, sseqid, pident, *others = line.split("\t")

            # Only screen for identity over 99
            if float(pident) > 99:
                transcript = re.search(r"(\S+)\|\S+", qseqid).group(1)
                sp_id = re.search(r".+sp\|(.+)\..+", sseqid).group(1)

                # Mapped the transcript with the first corresponding BLAST hit
                if transcript not in transcript_to_protein:
                    transcript_to_protein[transcript] = sp_id

        blast_file.close()

    return transcript_to_protein


def gene_go_dict(gene_to_go_filename):
    """Load protein IDs and corresponding unique GO terms to the dictionary.

    Args:
        filename (str): File path to the gene association .gaf file.
    Returns:
        dictionary: A dictionary matching the Swissport ID(key) to unique GO
                    IDs (value).
    """

    gene_to_go = {}
    if not gene_to_go_filename.endswith(".gaf"):
        print("Please provide valid gene association .gaf file.")

    else:
        gene_to_go_file = open(gene_to_go_filename)
        for line in gene_to_go_file:
            column_info = re.split(r"\t", line)
            object_id = column_info[1]
            go_id = column_info[4]

            # Combine unique GO IDs into existing keys
            if object_id in gene_to_go:
                gene_to_go[object_id].add(go_id)
            else:
                gene_to_go[object_id] = {go_id}

        gene_to_go_file.close()

    return gene_to_go


def go_name_dict(go_terms_filename):
    """Load GO IDs and their names to the dictionary.

    Args:
        filename (str): A GO terms file.

    Returns:
        dictionary: A dictionary matching the GO ID(key) to its description
        (value).
    """

    go_to_desc = {}
    if not go_terms_filename.endswith(".obo"):
        print("Please provide valid GO terms .obo file.")

    else:
        go_terms_file = open(go_terms_filename)
        terms = go_terms_file.read()
        terms = re.findall(r"\[Term]\n(.*?)\n\n", terms, re.DOTALL)

        for term in terms:
            go_id = re.search(r"^id:\s+(GO:\d+?)\n", term).group(1)
            go_name = re.search(r"^name:\s+(.+?)\n", term, re.M).group(1)

            # Check if both ID and name have a value before adding
            if go_id and go_name:
                go_to_desc[go_id] = go_name

        go_terms_file.close()

    return go_to_desc


def wrtie_annotations(diff_exp_filename, transcript_to_protein, gene_to_go,
                      go_to_desc, report_filename):
    """Loop through the differential expression file and dictionaries then
    write annotation to output .tsv file.

    Args:
        diff_exp_filename (str): A differential expression file with four
                                 stressed conditions.
        transcript_to_protein (dictionary): A dict mapping transcript to its
                                            SwissPort ID.
        gene_to_go (dictionary): A dict mapping SwissPort ID to its GO terms.
        go_to_desc (dictionary): A dict mapping GO term to its name.
        report_file (str): An output .tsv file.

    Return:
        file: A formated file.
    """

    report_file = open(report_filename, "w")
    if not diff_exp_filename.endswith(".matrix"):
        print("Please provide valid differential expression .matrix file.")

    else:
        # Loop through differential expression file
        diff_exp_file = open(diff_exp_filename)
        diff_exp_file.readline()  # skip header

        for transcript_info in diff_exp_file:

            # Lookup the protein ID
            transcript, *conditions = transcript_info.rstrip().split("\t")
            protein = transcript_to_protein.get(transcript, "NA")

            # Lookup GO IDs and GO description
            count = 0
            go_ids = gene_to_go.get(protein, ["NA"])
            for go_id in sorted(go_ids):
                go_desc = go_to_desc.get(go_id, "NA")

                # Print results to REPORT output
                if count == 0:
                    op_1 = [transcript, protein, *conditions, go_id, go_desc]
                    report_file.write("\t".join(op_1) + "\n")
                    count += 1
                else:
                    op_2 = [5*"\t", go_id, go_desc]
                    report_file.write("\t".join(op_2) + "\n")

        diff_exp_file.close()
        report_file.close()

    return report_file


def main():
    "The main function of the script."

    # Define file name for reading and writing.
    blast_filename = "blastp.outfmt6"
    gene_to_go_filename = "gene_association_subset.gaf"
    diff_exp_filename = "diffExpr.P1e-3_C2.matrix"
    go_terms_filename = "go-basic.obo"
    report_filename = "report.tsv"

    # Make dictionaries
    transcript_to_protein = transcript_protein_dict(blast_filename)
    gene_to_go = gene_go_dict(gene_to_go_filename)
    go_to_desc = go_name_dict(go_terms_filename)

    # Produce output file
    wrtie_annotations(diff_exp_filename, transcript_to_protein, gene_to_go,
                      go_to_desc, report_filename)


if __name__ == "__main__":
    main()
