#!/usr/bin/env python3
"""Parse genbank file into FASTA format

The script allows user to parse genbank file into FASTA format.
The genebank file can be imported per user's input and fasta file will be
exported after parsing.

This file contains below function to parse the gene information:
    * get_header - returns the FASTA header that contains info about VERSION
                   and DEFINITION
    * get_sequence - returns the FASTA format sequence under ORIGIN
    * split_records - returns a list of all possible gene record in a genbank
                      file

Author: Jia Yi Terri Shen
Date: October 2019
"""

import sys
import re
import textwrap


def get_header(record):
    """Parse for the VERSION and DEFINITION values within the file and stitch
    them together to form and return a >FASTA header.

    Args:
        record(string): A multi-line GenBank record
    Returns:
        string: A signle-line that contains the VERSION and DEFINITION values
        in fasta format
    """

    version_parse = "".join(re.findall(r"^VERSION\s*(\S+)\s*", record, re.M))
    definition_parse = "".join(re.findall(r"^DEFINITION\s*(.*\s[\n\s]*.*)\.",
                                          record, re.M))

    return ">" + version_parse + " " + definition_parse


def get_sequence(record):
    """Parse for the ORIGIN value within record and convert the sequence into
    a FASTA-formatted sequence.

    Args:
        record(string): A multi-line GenBank record
    Returns:
        string: A multi-line sequence that is all CAPs, wrap around 70
        characters, and no whitespaces in between.
    """

    origin_parse = re.findall(r"ORIGIN.+((?:\n.+)+)", record, re.M)
    aaseq_parse = "".join(re.findall(r"([a-z]+)", origin_parse[0], re.M))

    return "\n".join(textwrap.wrap(aaseq_parse.upper(), width=70))


def split_records(filename):
    """Open the input file, read it, and split the file into individual
    GenBank records.

    Args:
        filename(string): A file path to the input GenBank file
    Returns:
        lists: Seperated entries that contains gene info as a list or an
               empty list if the file is not found
    """

    try:
        with open(filename) as gb_file:
            gb_info = gb_file.read()
            split_list = re.split(r"//\n", gb_info)
        record_list = []
        for record in split_list:
            if "LOCUS" in record:
                record_list.append(record)
        return record_list
    except FileNotFoundError:
        return []

    gb_file.close()


def main():
    """The main function of the script that will validate the user inputs and
    output some warining messages or export a file that contains gene
    information in FASTA format.

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

    # Accept up to two command-line arguments
    input_file = "<input_file>"
    output_filename = "<output_filename>"

    # The first argument is required and should be a genbank file
    if len(sys.argv) < 2 or not sys.argv[1].endswith(".gb"):
        sys.exit("Provide a GenBank file to convert to FASTA.")
    else:
        input_file = sys.argv[1]
        gb_records = split_records(input_file)

    # Check if the provided import file exist
    if not gb_records:
        sys.exit(input_file + " not found. Check the file path and try again.")
    elif len(sys.argv) == 2:
        output_filename = "sequences.fasta"
        sys.stdout = open("sequences.fasta", "w")
    elif len(sys.argv) == 3:
        output_filename = sys.argv[2] + ".fasta"
        sys.stdout = open(output_filename, "w")

    # Export gene information to FASTA format into the output file
    fasta_outputs = []
    for record in gb_records:
        fasta_record = "\n".join([get_header(record), get_sequence(record)])
        fasta_outputs.append(fasta_record)
    print("\n\n".join(fasta_outputs), end="\n\n", flush=False)

    sys.stdout.close()


if __name__ == "__main__":
    main()
