"""File for defining the BlastHit and Blast classes.

Author: Jia Yi Terri Shen
Date: December 2019
"""

import re

class BlastHit:
    """Collection of Blast Hit record from one Blast Hit.

    Args:
        blast_hit(string): One blast hit.

    Attributes:
        transcript_id(str): The transcript ID within the query sequence ID.
        sp_id(str): The SwissPort ID within the subject sequence ID.
        pident(float): The percent of identical match.
        mismatch(int): The number of mismatches.

    Methods:
        __lt__: Returns if the hit has less mismatch.
        hit_good_match: Return True if the record is a good really match.
    """

    def __init__(self, blast_hit):
        self.blast_hit = blast_hit
        column = blast_hit.split("\t")
        self.transcript_id = re.search(r"(\S+)\|\S+", column[0]).group(1)
        self.sp_id = re.search(r".+sp\|(.+)\..+", column[1]).group(1)
        self.pident = float(column[2])
        self.mismatch = int(column[4])

    def __repr__(self):
        return f"BlastHit({self.blast_hit})"

    def __lt__(self, other):
        """Return boolean True if the hit has less mismatch."""

        if isinstance(self.mismatch) != isinstance(other.mismatch):
            raise Exception("Mismatch object must be the same type.")
        return self.mismatch < other.mismatch

    def hit_good_match(self):
        """Return boolean True if the record is a good really match (>95%)"""

        return self.pident > 95


class Blast:
    """Collections of Blast objects.

    Arg:
        blast_filename(string): One blast file name.

    Attribute:
        hits(list): A list of Blast objects from .outfmt6 file.

    Method:
        __iter__: Return iterator of the blast input.
    """

    def __init__(self, blast_filename):
        self.blast_filename = blast_filename
        with open(blast_filename) as blast_file:
            blast_hits = blast_file.read()
            blast_hits = blast_hits.rstrip("\n").split("\n")
            self.blast_hit_list = [BlastHit(hit) for hit in blast_hits]

    def __repr__(self):
        return f"Blast({self.blast_filename})"

    def __iter__(self):
        """Return iterator of the blast input."""

        return iter(self.blast_hit_list)
