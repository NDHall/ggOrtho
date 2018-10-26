import gg_ortho.base as gg
import logging
import pickle
import yaml

logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


def create_bag_of_genes(dcfile,parent, ks, qac=True, pid_cutoff=70, syn_len_cutoff=5,ks_cutoff=3, strict_ks=False,call=False):
    if ks is True :
        bog = gg.parse_ks(dcfile,pid_cutoff,ks_cutoff=ks_cutoff, syn_len_cutoff=syn_len_cutoff, call=call,parent=parent,
                          qac=qac,strict_ks=strict_ks)
        return bog
    else :
        logging.ERROR('non-ks DAGChainer files not currently supported.')
        return None



def parse_bog(bog,prop):
    passing_genes_cp = {}
    passing_genes_pc = {}
    for gene in bog:
        if len(bog[gene]['pass']) == prop:
            for classy_gene in bog[gene]['pass']:  # this will typically be one.
                passing_genes_cp[classy_gene.cgene] = classy_gene.pgene
                passing_genes_pc[classy_gene.pgene] = classy_gene  # this is important for writing files, out.
                print(classy_gene.cgene)
    return [passing_genes_cp, passing_genes_pc]










if __name__ == '__main__':
    #   dcfile = '/home/ndh0004/Documents/gg_ortho_indica/' \
    #            '51527_52024.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A5.aligncoords.qac1.1.50.gcoords.ks'
    #   prop = 1
    #   parent = 'b'
    #   ks = True
    yaml_config = '/home/ndh0004/Documents/gg_ortho_indica/config_4.yaml'

    yobj = yaml.load(open(yaml_config, 'r'))
    prop = yobj['prop']
    parent = yobj['parent']
    ks = yobj['ks']
    target_dir = yobj['target_dir']
    dcfile = yobj['dc_chain_and_seqs']['dcfile']
    parent_seqs = yobj['dc_chain_and_seqs']['parent_seqs']
    child_seqs = yobj['dc_chain_and_seqs']['child_seqs']


    bog = create_bag_of_genes(dcfile,parent, ks)
    info_dir = '{t}/{i}'.format(t=target_dir.rstrip("/"), i=yobj['info_dir'])
    gg.createFolder(info_dir)



    if bog is not None:
       cp_dict, pc_dict = parse_bog(bog,prop)
       gg.write_parent_fastas(pc_dict,parent_seqs,target_dir)
       gg.append_seqs_to_parent_fastas(pc_dict,cp_dict,child_seqs,target_dir)
       #write dictionary to json file.


    pickle.dump(pc_dict, open('{i}/{p}'.format(i=info_dir,p=yobj['pickle_handle']), 'wb'))






