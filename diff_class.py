"""File for defining the DiffExp and Matirx classes.

Author: Jia Yi Terri Shen
Date: December 2019
"""

class DiffExp:
    """Collection of Differential Expression info from .matrix file.

    Args:
        transcript_info(string): One differential expression informtation.

    Attributes:
        transcript(str):
        sp_ds(float): diauxic shift
        sp_hs(float): heat shock
        sp_log(float): logarithmic growth
        sp_plat(float): plateau phase

    Methods:
        data_attributes: Return the data sample attributes as a tuple.
    """

    def __init__(self, transcript_info):
        self.transcript_info = transcript_info
        self.transcript, self.sp_ds, self.sp_hs, self.sp_log, self.sp_plat \
            = transcript_info.rstrip().split("\t")

    def __repr__(self):
        return f"DiffExp({self.transcript_info})"

    def data_attributes(self):
        """Return tuple that contains data sample attributes."""

        return (self.sp_ds, self.sp_hs, self.sp_log, self.sp_plat)


class Matrix:
    """Collections of differential expression objects.

    Arg:
        diff_exp_filename: One differential expresssion matrix file name.

    Attribute:
        expression(list): A list of DiffExp objects.

    Method:
        __iter__: Return iterator of the diff_exp input.
    """

    def __init__(self, diff_exp_filename):
        self.diff_exp_filename = diff_exp_filename
        with open(diff_exp_filename) as diff_exp_file:
            self.expressions = [DiffExp(info) for info in diff_exp_file]

    def __repr__(self):
        return f"Matrix({self.diff_exp_filename})"

    def __iter__(self):
        """Return iterator of the diff_exp input."""

        return iter(self.expressions)
