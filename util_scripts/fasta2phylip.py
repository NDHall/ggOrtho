
"""
code taken from
https://biopython.org/wiki/AlignIO a
and modified slightly.


"""





from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def convertFasta2Phylip(instring,outstring):
    """

    :param instring: in fasta
    :param outstring: out phylip handle
    :return: basic phylip
    """


    input_handle = open(instring, "rU")
    output_handle = open(outstring, "w")

    alignments = AlignIO.parse(input_handle, "fasta")
    AlignIO.write(alignments, output_handle, "phylip")

    output_handle.close()
    input_handle.close()


def selectFasta2Phylip(instring,outstring,seq_list):
    """

    :param instring: in fasta
    :param outstring: out phylip name
    :param seq_list: list of required sequences
    :return:  phylip
    """

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

        output_handle.close()

    input_handle.close()

def writePaml(rec_add, outstring, column_width=None):
    """

    :param rec_add: list that contains [rec.id, rec,...]
    :param outstring: output handle
    :param column_width: either None type or pre-calculated width based on seq.list ids. If None type, it is calculated
           within the function.
    :return:  writes out PAML style phylip.
    """

    if column_width is None:
        rec_title_len = []
        for rec in rec_add[0::2]:
            rec_title_len.append(len(rec))
        column_width = max(rec_title_len) + 2
    with open(outstring, 'w') as fout:
        aln_len = len(str(rec_add[1].seq).replace(" ", ""))
        num_alns = len(rec_add[0::2])
        fout.write("  {n} {l}\n".format(n=num_alns, l=aln_len))

        for record in rec_add[1::2]:
            space_to_add = "".join([ " " for i in range((column_width-len(record.id)))])
            print("{r}{s}RIGHT".format(r=record.id,s=space_to_add))
            fout.write('{rec}{s}{seq}\n'.format(rec=str(record.id),
                                              s=space_to_add,
                                              seq=str(record.seq).replace(" ", "")))


def fasta2PamlPhylip(instring, outstring, seq_list=None):
    """
    This function will convert fasta to paml style phylip. It can take a list of required taxa, but the list is not
    required.

    :param instring: target fasta
    :param outstring: out phylip handle
    :param seq_list: list of sequences to return.
    :return: paml style phylip.
    """
    input_handle = open(instring, "rU")


    alignments = AlignIO.read(input_handle, "fasta")
    rec_add = []
    _seq_list = None
    for rec in alignments:
        if seq_list is not None:
            with open(seq_list) as f:
                _seq_list = f.read().rstrip('\n').split('\n')

            input_handle = open(instring, "rU")
            rec_title_len = []
            for seq in _seq_list:
                rec_title_len.append(len(seq))
            column_width = max(rec_title_len) + 2
            if str(rec.id) in _seq_list  and len(rec.seq) > 0:
                rec_add += [str(rec.id), rec]
        else :
            rec_add += [str(rec.id), rec]

    if _seq_list is not None:
        if set(rec_add[0::2]) == set(_seq_list):
            writePaml(rec_add,outstring,column_width)
    else:
        writePaml(rec_add, outstring)


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
    parser.add_argument( '--paml',
                        help='convert to paml out',
                        required=False,
                        default=False,
                        action='store_true',
                        dest='paml')


    argv = parser.parse_args()

    for in_file in argv.in_files:
        outstring = '{f}_{o}.phy'.format(f=in_file, o=argv.out)
        if argv.list is None and argv.paml is None:
            convertFasta2Phylip(in_file,outstring)
        elif argv.list is not None and argv.paml is False :
            selectFasta2Phylip(in_file, outstring, argv.list)
        elif  argv.paml is True:
            fasta2PamlPhylip(in_file, outstring, argv.list)




