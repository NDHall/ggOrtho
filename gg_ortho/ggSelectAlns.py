import gg_ortho.base as gg
import yaml
import glob




def select_alns(yaml_config):
    yobj = yaml.load(open(yaml_config,'r'))

    fa_path = gg.get_fasta_paths(yobj['target_dir'])

    info_dir = '{t}/{i}'.format(t=yobj['target_dir'].rstrip("/"), i=yobj['info_dir'])
    exclude_seqs = '{t}/exclude.tsv'.format(t=info_dir)
    exclude_names = gg.get_exclude_names(exclude_seqs)
    gg.filter_fasta_paths(fa_path,exclude_names,acceptable_minimum=yobj['acceptable_min'])






if __name__ == '__main__':

    yaml_config = '/home/ndh0004/Documents/gg_ortho_indica/config_6.yaml'

    select_alns(yaml_config)