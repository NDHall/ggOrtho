import gg_ortho.base as gg
import networkx as nx



def get_bog(dcfile, pid_cutoff=0, ks_cutoff=20, syn_len_cutoff=5,
            call='', parent='b', qac=False, strict_ks=False):
    '''

    :param dcfile: dagchainer file for output
    :param pid_cutoff: set to zero because we are just accepting whatever is being output by synmap here
    :param ks_cutoff: set high, but not used
    :param syn_len_cutoff: 5 genes because we again are accepting synmap output here basically unchanged.
    :param call: '' no labels.
    :param parent: b, but it does not really matter in this case.
    :param qac:  False we are not filtering results here
    :param strict_ks: False, we want all calls. here.
    :return:  bag of gene dict.
    '''
    bog = gg.parse_ks(dcfile,
                      pid_cutoff,
                      ks_cutoff=ks_cutoff,
                      syn_len_cutoff=syn_len_cutoff,
                      call=call,
                      parent=parent,
                      qac=qac,
                      strict_ks=strict_ks)
    return bog


def parse_bog(bog):
    passing_genes_cp = {}
    passing_genes_pc = {}
    graph = nx.Graph()
    for gene in bog:
        for classy_gene in bog[gene]['pass']:  # this can be 1 but that should be the lower limit.
            print(classy_gene.pgene, classy_gene.cgene)
            if graph.has_edge(u=str(classy_gene.pgene) , v=str(classy_gene.cgene) ) is False: # this way we will only add each gene relationship once.
                graph.add_edge(str(classy_gene.pgene),str(classy_gene.cgene))
                graph.add_edge(str(classy_gene.cgene),str(classy_gene.pgene))
                passing_genes_cp[classy_gene.cgene] = classy_gene.pgene
                passing_genes_pc[classy_gene.pgene] = classy_gene  # this is important for writing files, out.
                print(classy_gene.cgene)
    return [passing_genes_cp, passing_genes_pc]

if __name__ =='__main__':

     dcfile = '/home/ndh0004/Downloads/tmp_passthrough_files/' \
            '51576_51576.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.gcoords.ks'
     parent_seqs = '/home/ndh0004/Documents/gg_ortho_indica/Ecor_PR202_maker_predicted_transcripts2.fa'
     child_seqs = parent_seqs
     target_dir = '/home/ndh0004/Downloads/ksSeqs'
     bog = get_bog(dcfile=dcfile,
            pid_cutoff=0,
            ks_cutoff=20,
            syn_len_cutoff=5,
            strict_ks=False,
            qac=False)
     cp_dict, pc_dict = parse_bog(bog)
     gg.write_parent_fastas(pc_dict,parent_seqs,target_dir)
     gg.append_seqs_to_parent_fastas(pc_dict,cp_dict,child_seqs,target_dir)




