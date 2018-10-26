import gg_ortho.base as gg
import logging
import pickle
import pyaml
import yaml


logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)





def add_new_seqs(add_seq_blast_handle,e_cutoff,score_cutoff, per_len_cutoff):
    tg_dict = gg.filter_blast_results(add_seq_blast_handle, e_cutoff,score_cutoff,per_len_cutoff)
    return tg_dict


def add_blast_seqs(target_dir, p_dict, tg_dict, add_seq_handle,pickle_handle,label=None):
    gg.append_blast_seqs_to_parent_fastas(target_dir,add_seq_handle,tg_dict,p_dict,pickle_handle,label)


def add_seq(seq_blast_handle,seq_handle,target_dir,pickle_handle,e_cutoff,score_cutoff, per_len_cutoff,
            yaml_seq_handle=None,label=''):

    if yaml_seq_handle is None and seq_blast_handle is not None:
        tg_dict = add_new_seqs(seq_blast_handle, e_cutoff, score_cutoff, per_len_cutoff)

    elif yaml_seq_handle is not None and seq_blast_handle is None:
        tg_dict = yaml.load(open(yaml_seq_handle,'r'))

    p_dict = pickle.load(open(pickle_handle,'rb'))
    add_blast_seqs(target_dir, p_dict, tg_dict, seq_handle, pickle_handle,label)


if __name__ == '__main__':


    yaml_config = '/home/ndh0004/Documents/gg_ortho_indica/config_5.yaml'
    yobj = yaml.load(open(yaml_config, 'r'))

    if 'blast_and_seq' in yobj:
        for key in yobj['blast_and_seq']:
            if 'e_cutoff' in yobj['blast_and_seq'][key]:
                e_cutoff = float(yobj['blast_and_seq'][key]['e_cutoff'])
            else:
                e_cutoff = float(yobj['e_cutoff'])

            if 'score_cutoff' in yobj['blast_and_seq'][key]:
                score_cutoff = float(yobj['blast_and_seq'][key]['score_cutoff'])
            else:
                score_cutoff = float(yobj['score_cutoff'])
            if 'per_len_cutoff' in yobj['blast_and_seq'][key]:
                per_len_cutoff = int(yobj['blast_and_seq'][key]['per_len_cutoff'])
            else:
                per_len_cutoff = int(yobj['per_len_cutoff'])
            if 'label' in yobj['blast_and_seq'][key]:
                label = yobj['blast_and_seq'][key]['label']
            else:
                label = None

            target_dir = yobj['target_dir']
            blast = yobj['blast_and_seq'][key]['blast']
            seqs = yobj['blast_and_seq'][key]['seq']
            info_dir = '{t}/{i}'.format(t=target_dir.rstrip("/"), i=yobj['info_dir'])
            pickle_handle='{i}/{p}'.format(i=info_dir, p=yobj['pickle_handle'])
            add_seq(blast,seqs,target_dir, pickle_handle, e_cutoff,score_cutoff, per_len_cutoff,label=label)

    if 'yaml_and_seq' in yobj:
        for key in yobj['yaml_and_seq']:
            if 'e_cutoff' in yobj['yaml_and_seq'][key]:
                e_cutoff = float(yobj['yaml_and_seq'][key]['e_cutoff'])
            else:
                e_cutoff = float(yobj['e_cutoff'])

            if 'score_cutoff' in yobj['yaml_and_seq'][key]:
                score_cutoff = float(yobj['yaml_and_seq'][key]['score_cutoff'])
            else:
                score_cutoff = float(yobj['score_cutoff'])
            if 'per_len_cutoff' in yobj['yaml_and_seq'][key]:
                per_len_cutoff = int(yobj['yaml_and_seq'][key]['per_len_cutoff'])
            else:
                per_len_cutoff = int(yobj['per_len_cutoff'])
            if 'label' in yobj['yaml_and_seq'][key]:
                label = yobj['yaml_and_seq'][key]['label']
            else:
                label = None
            target_dir = yobj['target_dir']
            seqs = yobj['yaml_and_seq'][key]['seq']
            info_dir = '{t}/{i}'.format(t=target_dir.rstrip("/"), i=yobj['info_dir'])
            pickle_handle = '{i}/{p}'.format(i=info_dir, p=yobj['pickle_handle'])
            yaml_handle = yobj['yaml_and_seq'][key]['yaml_handle']
            add_seq(None,seqs,target_dir, pickle_handle, e_cutoff,score_cutoff,
                    per_len_cutoff,yaml_seq_handle=yaml_handle,
            label=label)



