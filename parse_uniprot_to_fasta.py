"""Convert Uniport record into FASTA format.

This script parse the Uniport record (uniprot-neurofibromas.txt) and create
FASTA format file.
Eache Swiss-Prot record will be parse and converted to a FASTA-formatted
header and sequence.

This file imported an output file called "output.fasta" and contains the
following class and functions:
    * UniprotEntry - A class defines the collection of UniprotEntry object
                     and following method:
                        * is_reviewed: To crean if the record is reviewed.
                        * to_fasta: To construct FASTA header and sequence.
    * main - The main function of the script to create FASTA output file.

Authout: Jia Yi Terri Shen
Date: November 2019
"""

import re
import textwrap

class UniprotEntry():
    """Collection of UniprotEntry object per Uniprot record.

    Attributes:
        db(str): Indicate the Uniprot records are from Swiss-Prot database.
        entry_name(str): Entry name at line ID.
        status(str): Review status of the entry at line ID.
        accession(str): Accession ID at line AC.
        gene_name(str): Gene name, or (OrderedLousName/ORFNames) at line GN.
        protein_name(str): Recommended protein name or (Fragment) at line DE.
        protein_evi(int): Index for proteion evidence at line PE.
        organism(str): Formal organism name at line OS.
        taxon_id(int): Taxon_id at line OX.
        seq_version(int): Index for sequence version at line DT.
        seq(str): Sequence at line block follow by line SQ.
    Methods:
        is_reviewed: Returns TRUE if the record is reviewed.
        to_fasta: Returns header and sequence in fasta format.
    """

    db = "sp"

    def __init__(self, record):
        # Find the status and entry name - ID
        entry_name_pattern = re.compile(r"ID\s+(\S+)\s+(\S+);.+")
        self.status = re.search(entry_name_pattern, record).group(2)
        self.entry_name = re.search(entry_name_pattern, record).group(1)

        # Find the accession ID - AC
        self.accession = re.search(r"^AC\s+(.+?);", record, re.M).group(1)

        # Find the gene name - GN
        gene_name_pattern = re.compile(r"GN\s+\w+=(.+?)\s*({.+})*;", re.M)
        if re.search(gene_name_pattern, record):
            self.gene_name = re.search(gene_name_pattern, record).group(1)
        else:
            self.gene_name = ""

        # Find the protein name or (Fragment) - DE
        protein_pattern = re.compile(r"DE\s+RecName:\s+Full=(.*?)\s*({.+})*;" \
                                    , re.M)
        flags_pattern = re.compile(r"^DE\s+Flags:\s+(Fragment);", re.M)

        if re.search(protein_pattern, record):
            protein_name = re.search(protein_pattern, record).group(1)
        else:
            protein_name = None

        if re.search(flags_pattern, record):
            self.protein_name = f"{protein_name} " \
                                f"({re.search(flags_pattern,record).group(1)})"
        else:
            self.protein_name = f"{protein_name}"

        # Find protein eviende - PE
        self.protein_evi = re.search(r"^PE\s+(\d+):", record, re.M).group(1)

        # Find organism name - OS
        self.organism = re.search(r"^OS\s+(.+)\s\(.+\).", record, re.M). \
                                 group(1)

        # Find taxon_id - OX
        self.taxon_id = re.search(r"^OX\s+NCBI_TaxID=(\d+)\s*({.+})*;", \
                                 record, re.M).group(1)

        # Find seq_version - DT
        self.seq_version = re.search(r"^DT\s+.+sequence\sversion.(\d).", \
                                    record, re.M).group(1)

        # Find seq - SQ
        seq_block = re.search(r"SQ.+;\n(.+)", record, re.S).group(1)
        seq_combine = "".join(re.findall(r"(\w+)", seq_block, re.M))
        self.seq = "\n".join(textwrap.wrap(seq_combine, width=60))


    def is_reviewed(self):
        """Returns if the record is reviewed or not.

        Returns:
            boolean: True if the record is reviewed, False otherwise.
        """
        return self.status == "Reviewed"


    def to_fasta(self):
        """Parse into the FASTA format.

        Returns:
            String: The relevant attributes in UniProtKB FASTA format.
        """
        header = f">{UniprotEntry.db}|{self.accession}|{self.entry_name}" \
                 f" {self.protein_name} OS={self.organism}" \
                 f" OX={self.taxon_id} GN={self.gene_name}" \
                 f" PE={self.protein_evi} SV={self.seq_version}"
        sequence = self.seq
        return  header + "\n" + sequence + "\n"


def main():
    """ The main function of the script."""

    with open("uniprot-neurofibromas.txt") as uniport_file:
        records = uniport_file.read()
        records = re.split(r"//\n", records)

        with open("output.fasta", "w") as output_fasta:
            for record in records:
                if record:
                    uniprot_entry = UniprotEntry(record)
                    if uniprot_entry.is_reviewed():
                        output_fasta.write(uniprot_entry.to_fasta())


if __name__ == "__main__":
    main()
