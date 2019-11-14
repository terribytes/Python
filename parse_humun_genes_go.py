#!/usr/bin/env python3
"""Parse human genes and their list of assoicated GO terms

This script create an annotation report of human genes and their list of
associating GO terms, including the parent terms.

This file contains below function to parse the gene information:
    * split_terms - Function accepts .obo file and a list that contains
                    individual GO temrms.
    * map_protein_to_go - Function accept .gaf file and return a dictionary
                          that maps the protein ID with its GO terms.
    * parse_go_term - Function accept a GO term and return a list of
                      original GO term with it's direct parent GO terms.
    * find_parent_terms - Fiction accept a GO term and an empty dictionary
                          which recursively look for all parent GO terms.

Author: Jia Yi Terri Shen
Date: November 2019
"""

import sys
import re


def split_terms(filename):
    """Open the GO term file and split the file into individual GO terms.
    Args:
        filename(string): A file path to a GO terms file.
    Return:
        list: Sepatated terms as a list or an empty list if the file is not
              found.
    """

    try:
        with open(filename) as go_term_file:
            go_term_info = go_term_file.read()
            split_go_term_list = re.split(r"\[Term\]", go_term_info)
        return split_go_term_list

    except FileNotFoundError:
        return []


def map_protein_to_go(filename):
    """Open the GO annotation file (GAF) and build the mapping relationship
    between the protein ID and its list of associating GO terms.

    Arg:
        filename(string): A file path to a GO annotation file (GAF).
    Return:
        dictionary: A protein and GO terms mapping dictionary or an empty
        dictionary if the file is not found.
    """

    try:
        with open(filename) as go_association_file:
            go_association = go_association_file.read()
            split_go_association = re.split(r"\n", go_association)

        # Ignore the general file information, which is the line starting
        # with "!"".
        go_association_info = []
        for line in split_go_association:
            if line and not line.startswith("!"):
                go_association_info.append(line)

        # Declare the tuple to parse the protein and go term as a pair and
        # store it in the set to avoid duplicate situation
        go_protein_dict = {}
        for column in go_association_info:
            column_info = re.split(r"\t", column)
            protein_id = column_info[1]
            go_term = column_info[4]

            if protein_id in go_protein_dict:
                go_protein_dict[protein_id].add(go_term)
            else:
                go_protein_dict[protein_id] = {go_term}
        return go_protein_dict

    except FileNotFoundError:
        return {}


def parse_go_term(term):
    """Parsing the ID and is_a in the GO term file and return a collection.

    Arg:
        term(string): A single GO term.
    Return:
        tuple: A collection that correspond the ID in it's corresponding
        is_a values.
    """

    id_pattern = re.compile(r"^id:\s+(GO:[0-9]+)\s+", re.M)
    is_a_pattern = re.compile(r"^is_a:\s+(.*?) !", re.M)

    id_parse = re.findall(id_pattern, term)
    is_a_parse = re.findall(is_a_pattern, term)

    return id_parse, is_a_parse


def find_parent_terms(go_id, go_dict):
    """Search for all the parent terms of the starting GO term and return
    them as a collection.

    Args:
        go_id: A single GO term.
        go_dict: A dictionary of GO terms.
    Return:
        tuple: A collection that correspond the parent GO term to it's
        starting GO term.
    """

    go_set = set()
    values = go_dict[go_id]

    for value in values:
        go_set.add(value)
        more_values = find_parent_terms(value, go_dict)
        for more_value in more_values:
            go_set.add(more_value)

    return go_set


def main():
    """The main function of the script that will validate the user inputs and
    output some warining messages or export an annotation file that contains
    three columns which entitle protein ID, GO terms, and it's parent terms.

    Args:
        input_string(string): A mandatory genbank file that needs to be
                              provided by the users
        filename(string): A non-mandatory file name that might be provided by
                              the users

    Returns:(One of the below)
        fasta_outputs(string): An output file that contains the gene
                               information will be provided with defined or
                               default filename if the genbank file is valid
        sys.exit(string): Two different error messages will be provided as
                          warnings to the user depending on their inputs
    """

    # Accept up to three command-line arguments
    input_terms = "<input_GO_terms_file>"
    input_annotations = "<input_gene_associations_file>"
    output_filename = "<output_filename>"


    # The first two arguments are required GO terms file ending with .obo
    # and gene association GAF file ending with .gaf
    if len(sys.argv) < 3:
        sys.exit("Please provide required GO terms .obo file and gene " +
                 "assocatiion .gaf file.")
    elif not sys.argv[1].endswith(".obo"):
        sys.exit("Please provide a GO terms .obo file.")
    elif not sys.argv[2].endswith(".gaf"):
        sys.exit("Please provide a gene association .gaf file.")
    else:
        input_terms = sys.argv[1]
        input_annotations = sys.argv[2]


    # Check if the provided import .obo or .gaf files exist
    if not input_terms:
        sys.exit(input_terms + " not found. Check the file path and try again.")
    elif not input_annotations:
        sys.exit(input_annotations + " not found. Check the file path and try again.")
    elif len(sys.argv) == 3:
        output_filename = "results.tsv"
        sys.stdout = open("results.tsv", "w")
    elif len(sys.argv) == 4:
        output_filename = sys.argv[3] + ".tsv"
        sys.stdout = open(output_filename, "w")


    # parse id and is_valeus and make a go_dict
    split_input_terms = split_terms(input_terms)
    go_dict = {}
    for record in split_input_terms:
        (go_id, is_a) = parse_go_term(record)
        key_go_dict = "".join(go_id)
        go_dict[key_go_dict] = is_a


    # Export an annotation gene information to tsv format into the output file
    gene_association_map = map_protein_to_go(input_annotations)
    for protein, go_ids in sorted(gene_association_map.items()):
        print(protein, end="")

        for go_id in sorted(go_ids):
            parent_go_ids = find_parent_terms(go_id, go_dict)

            count = 0
            for parent_go_id in sorted(parent_go_ids):

                if count == 0:
                    print("\t", go_id, "\t", parent_go_id)
                    count += 1
                else:
                    print("\t", parent_go_id, sep="\t")

    sys.stdout.close()


if __name__ == "__main__":
    main()
