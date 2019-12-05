"""Create annotation file from Blast Hit file and Diff Exp matrix.

This script finds the transcript ID and SwissProt ID that are a good match
(pident >95%) in the Blast Hit file and map them with the transcript ID in
the Diff Exp matrix and produce an output file.

This file can be imported and contains the following functions:
    * tuple_to_stirng: Accept a tuple and retrun it as a tab-separated string.

Author: Jia Yi Terri Shen
Date: December 2019
"""

from diff_class import Matrix
from blast_class import Blast

def tuple_to_string(transcript_info):
    """Accept a tuple and retrun it as a tab-separated string.

    Arg:
        diff_exp_tuple(tuple): A differntial expression tuple.

    Return:
        str: A tab separated string.
    """

    return "\t".join(transcript_info.data_attributes())


def main():
    """ The main function of script"""

    #Open the BLAST and differential expressions files and create its objects
    blast_filename = "blastp.outfmt6"
    diff_exp_filename = "diffExpr.P1e-3_C2.matrix"
    blast = Blast(blast_filename)
    matrix = Matrix(diff_exp_filename)

    #Load transcript_id and sp_id within the good BlastHit into dictionary
    blast_dict = {blast.transcript_id: blast.sp_id \
                 for blast in blast.blast_hit_list if blast.hit_good_match()}

    #Perform lookup and produce an output annotation file
    with open("output.txt", "w") as output:
        for info in matrix.expressions:
            matrix_info = blast_dict.get(info.transcript, info.transcript) \
                                             + "\t" + tuple_to_string(info)
            output.write(matrix_info + "\n")


if __name__ == "__main__":
    main()
