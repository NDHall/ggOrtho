
"""
code taken from
https://biopython.org/wiki/AlignIO a
and modified slightly.


"""





from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def convertFasta2Phylip(instring,outstring):

    input_handle = open(instring, "rU")
    output_handle = open(outstring, "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip")

    output_handle.close()
    input_handle.close()


def selectFasta2Phylip(instring,outstring,seq_list ):
    with open(seq_list) as f:
        _seq_list = f.read().rstrip('\n').split('\n')

    input_handle = open(instring, "rU")


    alignments = AlignIO.read(input_handle, "fasta")
    rec_add = []
    for rec in alignments:
        if str(rec.id) in _seq_list  and len(rec.seq) > 0:
            rec_add += [str(rec.id), rec]
    if set(rec_add[0::2]) == set(_seq_list):
        outAlign = MultipleSeqAlignment(rec_add[1::2])
        output_handle = open(outstring, "w")
        AlignIO.write(outAlign, output_handle, "phylip")
        output_handle.close()

    input_handle.close()










if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description="""
    This is a fasta to phylip script it takes a list of fasta file and converts them to phylip files. """)
    parser.add_argument('-i', '--infile',
                        help='list of fasta files to convert to phylip files. ',
                        required=True,
                        type=str,
                        nargs='+',
                        dest='in_files')
    parser.add_argument('-o', '--out',
                        help='directory to save results to, default is current working directory',
                        required=False,
                        default='out',
                        type=str,
                        dest='out')
    parser.add_argument('-l', '--list',
                        help='directory to save results to, default is current working directory',
                        required=False,
                        default=None,
                        type=str,
                        dest='list')


    argv = parser.parse_args()

    for in_file in argv.in_files:
        outstring = '{f}_{o}.phy'.format(f=in_file, o=argv.out)
        if argv.list is None:
            convertFasta2Phylip(in_file,outstring)
        else :
            selectFasta2Phylip(in_file, outstring, argv.list)


