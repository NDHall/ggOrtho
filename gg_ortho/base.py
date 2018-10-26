import numpy as np
from Bio import SearchIO

class Syn:
    """
    class for Syntetic blocks
    """


    def __init__(self,block_num, score, org_a, org_b, orient, num_gene_pairs ):
        self.block_num = int(block_num)
        self.score = float(score)
        self.org_a = org_a
        self.org_b = org_b
        self.orient = orient
        self.num_gene_pairs = int(num_gene_pairs)
        self.genes = []
        self.avg_per_id = 'Null'

class Gene:
    """
    for genes in syntentic blocks
    """


    def __init__(self,org_a_chrom, org_a_region, org_a_start, org_a_stop,
            org_b_chrom, org_b_region, org_b_start, org_b_stop,
            evalue, accumulative_score ):
        self.org_a_chrom = org_a_chrom
        self.org_a_region = org_a_region
        self.org_a_start = int(org_a_start)
        self.org_a_stop = int(org_a_stop)
        self.org_b_chrom = org_b_chrom
        self.org_b_region = org_b_region
        self.org_b_start = int(org_b_start)
        self.org_b_stop = int(org_b_stop)
        self.evalue = evalue
        self.accumulative_scorer = accumulative_score



class Region:
    """
    for lastz hit reports in syntentic blocks.
    """

    def __init__(self, chrom, start,stop,gene,orientation,seq_type,score1,score2,percent_id):
        self.chrom = chrom
        self.start = int(start)
        self.stop = int(stop)
        self.gene = gene
        self.orientation = int(orientation)
        self.seq_type = seq_type
        self.score1 = int(score1)
        self.score2 = float(score2)
        self.percent_id = float(percent_id)


class KsSyn(Syn):

    def process_float(self,val):
        indicator = False
        try :
            val = float(val)
            if type(val) == float:
                indicator = True
        except:
            val = str(val)
        return [val, indicator]



    def __init__(self,block_num, score, org_a, org_b, orient, num_gene_pairs,ks,ka ):
        super().__init__(block_num, score, org_a, org_b, orient, num_gene_pairs)
        self.ks, self.ks_is_float = self.process_float(ks)
        self.ka, self.ka_is_float = self.process_float(ka)


class KsGene(Gene):

    def process_float(self,val):
        indicator = False
        try :
            val = float(val)
            if type(val) == float:
                indicator = True
        except:
            val = str(val)
        return [val, indicator]


    def __init__(self, ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,
            org_b_chrom, org_b_region, org_b_start, org_b_stop,
            evalue, accumulative_score ):
        super().__init__(org_a_chrom, org_a_region, org_a_start, org_a_stop,
                     org_b_chrom, org_b_region, org_b_start, org_b_stop,
                     evalue, accumulative_score)
        self.ks, self.ks_is_float = self.process_float(ks)
        self.ka, self.ka_is_float = self.process_float(ka)


class GeneInBag():
    def __init__(self, call, pgene, total_cregion,total_pregion, perid,meanperid,
                 ks, meanks,ka,meanka,
                 cgene, cgene_region,
                 pgene_region,syn_len,
                 pgene_orient, cgene_orient):
        self.call = call
        self.pgene = pgene
        self.ctregion = total_cregion
        self.ptregion = total_pregion
        self.perid = perid
        self.meanperid = meanperid
        self.ks = ks
        self.meanks = meanks
        self.cgene = cgene
        self.cgene_region = cgene_region
        self.pgene_region = pgene_region
        self.ka = ka
        self.meanka = meanka
        self.perid_len = None
        self.ka_len = None
        self.ks_len = None
        self.syn_len = syn_len
        self.cgene_orient = cgene_orient
        self.pgene_orient = pgene_orient
    def __str__(self):
        return '{p}\t{c}'.format(p=self.pgene,c=self.cgene)

"""
Functions to build classes. 
"""

def parse_ks(infile,pid_cutoff,ks_cutoff,syn_len_cutoff,
             call,parent='b',strict_ks=True, qac=True):
    """

    :param infile: This the handle for opening DAGChainer output.
    :param pid_cutoff: Minimum percent identitty to accepth
    :param ks_cutoff:  Maximum ks to accept.
    :param syn_len_cutoff: Minimum number of genes required to include genes from a syntentic block.
    :param out_file: This is a prefix that will be used for writitng out files.
    :param call: String this is the default input if no call is made for a gene. It can be used to label
                 species when calls are not being made.
    :param parent: This is the col of the DAGChainer file to use as the reference 'a'==0  and 'b'==1
    :param strict_ks: Only report calls that pass ks cutoffs. This drastically reduces call numbers.
    :param call_ab: Bool, call ab if the run is being used to separate A and B syntentic blocks.
    :param qac:  Was quota align used in the process if so True. This will directly affect how file is parsed.
    :param in_len: Default == None. If not none it is the handle for a tab delimited file that has chrom lengths.
                   These are turned into a dict of values so that we can spoof length of sequence.
                   structure ==
                   <length>\t<chrom name>\n
                   <length>\t<chrom name>\n
                   ...

    :return: output 5 files
            1. a set of _ab_gene.tsv calls.
            2. gff of A genes as children of syntenic regions.
            3. gff of B genes as children of syntenic regions.
            4. gff of all genes that were eligble to be called.
            5. _bag_of_genes.tsv. contains all the information for each gene stored in gene in bag variable.
    """
    bag_of_genes_dict = {}
    f = open(infile, 'r')
    ks_file = f.read()
    if qac is False:
        ks_file = ks_file.rstrip("\n").split('#')[2::]
    else:
        ks_file = ks_file.rstrip("\n").replace('\n###','\n#')
        ks_file = ks_file.split("#")[2::]
    for block, body in zip( ks_file[0::2],ks_file[1::2]):
        # we are not using block for this, right now. May want it later.
        body = body.rstrip("\n").split("\n")[1:] # we are dropping the comment line with column names.
        gene_list = ks_body_parser(body)
        bag_of_genes_dict = parse_gene_list(gene_list,pid_cutoff,ks_cutoff,syn_len_cutoff, bag_of_genes_dict,strict_ks,
                                            call,
                                            ref=parent)
    return bag_of_genes_dict

def ks_body_parser(body):
    """
    Here we are parsing the body and changing into a list of Ks_gene objects.
    I want to keep this function so that it is portable. So I am willing to
    run through the list a second time to parse out desired values.
    :param body:
    :return:
    """
    return_list = []  #
    for gene in body:
        ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,\
        org_b_chrom, org_b_region, org_b_start, org_b_stop,\
        evalue, accumulative_score = gene.split("\t")
        gene_obj = KsGene(ks, ka, org_a_chrom, org_a_region, org_a_start, org_a_stop,
                                   org_b_chrom, org_b_region, org_b_start, org_b_stop,
                          evalue, accumulative_score)
        #print(gene_obj.ka_is_float,gene_obj.ka,gene_obj.ks,gene_obj.org_a_stop)
        return_list.append(gene_obj)
    return return_list
def parse_gene_list(gene_list,pid_cutoff,ks_cutoff,syn_len_cutoff,bag_of_genes_dict, strict_ks,
                    call ,ref='b',cutoff=3.0):
    """

    :param gene_list: When the DAGChainer file is parsed it is split into a block(the header) and a gene_list(the body)
                      We can do this by reading the whole DAGChainer file into memory. then splitting based on `#` then
                      splitting resulting blocks by line breaks. in the case of ks files. gene list == body[1:] because
                      the body comes with its own header that needs discarded.
    :param pid_cutoff: The percent id cutoff for the mean pid of a syntenic block. I set the default in argparse object.
                       default == 90.0 for this program. This is because we are hunting a recent AB divergence.
    :param ks_cutoff:  maximum acceptable mean ks value for a syntenic block for a gene to be considered. This
                       is really a bit redundant after quota align. But it can be used as a filter. Though it mostly
                       succeeds in over-filtering the results, since ks and ka reports are not output for every gene.
                       Even if it is calculated for every gene.
    :param syn_len_cutoff: Minimum number of genes used a syntenic region for it to be used in calls.
    :param bag_of_genes_dict: It is  dict of dict {gene:{'pass':[Geneinbag,...j], 'fail':[Geneinbag,...j]}  genes that
                              pass all cutoff filters are added 'pass' else: gene is added in 'fail'.
    :param strict_ks:   Bool if true ks_cutoff is used to filter genes that are added bag_of_genes_dict{}.
                        default==False. I have found that calls get really sparse with ks
    :param call:        Default value of call variabe. Useful in case call_ab is true.
    :param ref:         This takes the parent variable. The defaut is b which corresponds the the b genome in
                        DAGChainer output or column 2 of the DAGChainer report. Some may prefer to run SynMap with
                        the AB ref in col 'a'. This can specified here.
    :param cutoff:     codeml and by extension CoGe report ks values 0-inf, even though they are based on the
                       Nei and Gojobori 1986 approach. Though this is not the exact model used proxmial to the output.
                       Adjustments are made for models of nucleotide substitions. I have limited the max ks to 3.0 which
                       is pretty high since without accounting for models of substition this would be max. We don't want
                       highly different sites for this comparison.
    :return:
    """

    baglet_of_genes = []
    syn_len = len(gene_list)
    if ref == 'a':
        parent_start = region_parser(gene_list[0].org_a_region)
        child_start = region_parser(gene_list[0].org_b_region)
        parent_stop = region_parser(gene_list[-1].org_a_region)
        child_stop = region_parser(gene_list[-1].org_b_region)

    else:
        child_start = region_parser(gene_list[0].org_a_region)
        parent_start = region_parser(gene_list[0].org_b_region)
        child_stop = region_parser(gene_list[-1].org_a_region)
        parent_stop = region_parser(gene_list[-1].org_b_region)
    if parent_start.start  > parent_stop.stop :
        pstart = parent_stop.stop -1
        pstop = parent_start.start
    else :
        pstart = parent_start.start -1
        pstop = parent_stop.stop
    if child_start.start > child_stop.stop :
        cstart = child_stop.stop -1
        cstop = child_start.start
    else :
        cstart = child_start.start -1
        cstop = child_stop.stop
    ks_vals = []
    ka_vals = []
    pid_vals = []
    for gene in gene_list:


        if ref == 'a' :
            parent_region = region_parser(gene.org_a_region)
            child_region = region_parser(gene.org_b_region)

        else:
            parent_region = region_parser(gene.org_b_region)
            child_region = region_parser(gene.org_a_region)


        pgene = parent_region.gene
        cgene_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=child_region.chrom,
                  start=child_region.start -1 ,
                  stop=child_region.stop)
        pgene_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=parent_region.chrom,
                  start=parent_region.start -1 ,
                  stop=parent_region.stop)
        perid = float(child_region.percent_id)
        if type(perid) is float:
            pid_vals.append(perid)
        ks = gene.ks
        ka = gene.ka
        #print(gene.ks, type(gene.ks))
        meanka = None
        if type(gene.ks) is float and \
            gene.ks <= cutoff:
            ks_vals.append(gene.ks)
        if type(gene.ka) is float and \
            gene.ka <= cutoff:
            ka_vals.append(gene.ka)
        cgene = child_region.gene
        ctot_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=child_region.chrom,
                  start=cstart ,
                  stop=cstop)
        ptot_region = '{chrom}\t{start}\t{stop}'.format(
                  chrom=parent_region.chrom,
                  start=pstart,
                  stop=pstop)
        classy_gene = GeneInBag(
            call=call,
            pgene=pgene,
            total_cregion=ctot_region,
            total_pregion=ptot_region,
            perid=perid,
            ks=ks,
            ka=ka,
            meanka=None,
            meanks=None,
            cgene=cgene,
            cgene_region=cgene_region,
            pgene_region=pgene_region,
            syn_len=syn_len,
            meanperid=0,
            cgene_orient=int(child_region.orientation),
            pgene_orient=int(parent_region.orientation))
        baglet_of_genes.append(classy_gene)

    ka_len = None
    perid_len = None
    #print('ks vals:',ks_vals)
    if len(ks_vals) > 0:
        meanks = np.mean(ks_vals)
        ks_len = len(ks_vals)
    else:
        meanks = None
        ks_len = None
    if len(ka_vals) > 0:
        meanka = np.mean(ka_vals)
        ka_len = len(ka_vals)
    if len(pid_vals) >0:
        meanperid = np.mean(pid_vals)
        perid_len = len(pid_vals)
    # update all values.
    # add to dictionary in dictionary classify as a high quality or low quality based on ks mean for block and
    # or mean percent id per gene or per block.
    for classy_gene in baglet_of_genes:
        #print(meanperid)
        classy_gene.meanka = meanka
        classy_gene.meanks = meanks
        classy_gene.meanperid = meanperid
        classy_gene.perid_len = perid_len
        classy_gene.ka_len = ka_len
        classy_gene.ks_len = ks_len
        # this is where we filter for gene comparisons.
        #print(classy_gene.meanperid)
        if classy_gene.pgene not in bag_of_genes_dict:
            bag_of_genes_dict[classy_gene.pgene] = {'pass':[], 'fail':[]}
        if strict_ks == True:
            if classy_gene.meanks is not None and \
               classy_gene.meanks <= ks_cutoff and \
               classy_gene.meanperid >= pid_cutoff and \
               classy_gene.syn_len >= syn_len_cutoff:
                    bag_of_genes_dict[classy_gene.pgene]['pass'].append(classy_gene)
            else:
                    bag_of_genes_dict[classy_gene.pgene]['fail'].append(classy_gene)
        else:
            # now we just call based on perid cutoff and syntenic block length.
            if classy_gene.meanperid >= pid_cutoff and \
               classy_gene.syn_len >= syn_len_cutoff:
                    bag_of_genes_dict[classy_gene.pgene]['pass'].append(classy_gene)
            else:
                    bag_of_genes_dict[classy_gene.pgene]['fail'].append(classy_gene)

    return bag_of_genes_dict



def region_parser( region):
    """

    :param region: takes region out of diag chain file. Essentially
                   these are blast hits
    :return: returns class region.
    """
    #print(region)
    chrom, start, stop, gene, orientation, seq_type, score1, score2, percent_id = region.split('||')
    region = Region(chrom, start, stop, gene, orientation, seq_type, score1, score2, percent_id)
    assert region.start <region.stop, "blast hits out of order in region %s"%(region)
    return region


"""
Functions for writing seqs out
"""

import os
from Bio import SeqIO
"""
Many thanks to keithweaver 
code takenfrom and modified 
https://gist.github.com/keithweaver/562d3caa8650eefe7f84fa074e9ca949
"""


def createFolder(dir):
    try:
        if not os.path.exists(dir):
            os.makedirs(dir)
    except :
        pass



def createSeqFailTable(pickle_handle):
    _target_dir = '/'.join(pickle_handle.split('/')[:-1])
    exclude_seqs = '{t}/exclude.tsv'.format(t=_target_dir)
    if os.path.isfile(exclude_seqs):
        f = open(exclude_seqs,'a' )
    else:
        f = open(exclude_seqs,'w' )
    return f

# Creates a folder in the current directory called data

def write_parent_fastas(pc_dict,parent_seqs,target_dir):
    target_dir = target_dir.rstrip('/') + '/'
    createFolder(target_dir)
    with open(parent_seqs,'r') as seq_handle:
        seqs= SeqIO.parse(seq_handle,'fasta')
        for record in seqs :
            """seqs fasta cds names must perfectly match the fasta used here.
               If we started something like re matching. The problem would be
               that g1 would could find g1,g10,g100 without specific knowledge
               of the naming conventions and rules which are varied.  
            """
            key = str(record.id)
            outhandle = '{g}.fa'.format(g=key)
            if key in pc_dict:
                classy_gene = pc_dict[key]
                outstring = '>{pg}\n{pseq}\n'.format(
                    pg=str(record.id),
                    pseq=str(record.seq)
                )
                block_folder = '{t}/{b}'.format(t=target_dir,b=classy_gene.ptregion.replace("\t",'_'))
                createFolder(block_folder)
                full_out_path = '{b}/{g}'.format(b=block_folder,g=outhandle)
                f = open(full_out_path,'w')
                f.write(outstring)
                f.close()

def append_seqs_to_parent_fastas(pc_dict,cp_dict, child_seqs,target_dir):
    target_dir = target_dir.rstrip('/') + '/'
    with open(child_seqs,'r') as seq_handle:
        seqs= SeqIO.parse(seq_handle,'fasta')
        for record in seqs :
            #print(str(record.id))
            if str(record.id) in cp_dict:
                #print('found')
                append_seq(cp_dict,pc_dict,record,target_dir)



def append_blast_seqs_to_parent_fastas(target_dir,blast_seqs,tg_dict,pc_dict,pickle_handle,label=None):
    target_dir = target_dir.rstrip('/') + '/'
    f = createSeqFailTable(pickle_handle)
    with open(blast_seqs,'r') as seq_handle:
        seqs= SeqIO.parse(seq_handle,'fasta')
        for record in seqs :
            #print(str(record.id))
            if str(record.id) in tg_dict and \
                   len(tg_dict[record.id]) == 1 :
                print('added {r}'.format(r=record.id))
                tg_dict[record.id] = tg_dict[record.id][0]
                append_seq(tg_dict,pc_dict,record,target_dir,label=label)
            elif str(record.id) in tg_dict and \
                   len(tg_dict[record.id]) > 1 :
                record_len = len(tg_dict[record.id])
                for parent_gene in  tg_dict[record.id]:
                    if label is None:
                        label = 'None'
                    outstring = '{pg}\t{l}\t{tg}\t{lab}\n'.format(
                        pg=parent_gene,
                        l=record_len,
                        tg=record.id,
                        lab=label
                    )
                    f.write(outstring)
    f.close()
import glob
import shutil



def get_fasta_paths(target_dir):
    fa_path = glob.glob('{td}/*/*.fa'.format(td=target_dir.rstrip('/')))
    return fa_path

def get_exclude_names(tsv):
    ret_exclude = []
    with open(tsv) as f :
        exclude = f.read().rstrip('\n').split('\n')
        for line in exclude:
            gene_name = line.split('\t')[0]
            ret_exclude.append(gene_name)
    return set(ret_exclude)



def filter_fasta_paths(fasta_paths,exclude_names,acceptable_minimum=2):
    for fasta in fasta_paths:
        entries = []
        fail = False
        wd = '/'.join( fasta.split('/')[:-1])
        fa = '/'.join( fasta.split('/')[-1])
        gene_name = '.'.join( fa.split('.')[:-1])
        fail_path = '{wd}/fail'.format(wd=wd)
        pass_path = '{wd}/pass'.format(wd=wd)
        createFolder(pass_path)
        createFolder(fail_path)
        if gene_name in exclude_names:
            shutil.move(fasta, fail_path)
        else:
            with open(fasta) as f :
                records = f.read().rstrip('\n').split('>')[1:]
                for record in records :
                    id = record.split('_')[0]
                    if id in entries:
                        fail = True
                    entries.append(id)
            if fail is True or len(entries) < acceptable_minimum:
               shutil.move(fasta, fail_path)
            else:
                shutil.move(fasta, pass_path)


def append_seq(cp_dict,pc_dict, record,target_dir,label=None):
    key = cp_dict[str(record.id)]
    outhandle = '{g}.fa'.format(g=key)
    if key in pc_dict :
        print('found_ke in pc_dict key == {k}'.format(k=key))
        classy_gene = pc_dict[key]
        if label is None:
            outstring = '>{pg}\n{pseq}\n'.format(
                pg=str(record.id),
                pseq=str(record.seq)
            )
        elif label is not None:
            outstring = '>{l}_{pg}\n{pseq}\n'.format(
                l=label,
                pg=str(record.id),
                pseq=str(record.seq)
            )
        block_folder = '{t}{b}'.format(t=target_dir, b=classy_gene.ptregion.replace("\t", '_'))
        full_out_path = '{b}/{g}'.format(b=block_folder, g=outhandle)
        f = open(full_out_path, 'a')
        f.write(outstring)
        print('wrote{o} \nto {p}'.format(o=outstring, p=full_out_path))
        f.close()

from Bio.Blast import NCBIXML

def filter_blast_results(blastp_xml_handle,e_cutoff, score_cutoff, per_len_cutoff):
    _xml_handle = open(blastp_xml_handle)
    blast_res = NCBIXML.parse(_xml_handle)
    transcript_to_genome_dict = {}
    for hit in blast_res:
        for aln in hit.alignments:
            for hsp in aln.hsps:
                if hsp.expect <= e_cutoff and \
                    hsp.score >= score_cutoff  and \
                        (hsp.identities/hit.query_length) >= per_len_cutoff:
                    """
                    Notes: Once a hit passes we can add it to new child to parent library.
                           Here also, names must be identical, this is not really an issue 
                           but it bears noting since sometimes, fasta names are altered to 
                           avoid characters such as spaces.
                           
                    Percent ID Proxy:
                    I don't believe we are using a true precent id here, but it comes close.
                    Basically we want to know how much of the total query length can be accounted
                    for by a postional match to the subject. If it passes we will add it.
                    """
                    print(aln.hit_def, hit.query, hit.query_length, aln.length)
                    clean_hit_name = hit.query.split(' ')[0]
                    if clean_hit_name not in transcript_to_genome_dict:
                        transcript_to_genome_dict[clean_hit_name] = [aln.hit_def]
                    else:
                        transcript_to_genome_dict[clean_hit_name].append(aln.hit_def)


    _xml_handle.close()
    return transcript_to_genome_dict


