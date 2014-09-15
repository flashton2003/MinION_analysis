__author__ = 'flashton'

try:
    import cPickle as pickle
except:
    import pickle

from minion_core import BlastRes, BlastHit, MinionRead

contigs = '/Users/flashton/Dropbox/H58_from_iMac/H58/data/refs/ST1_pAKU_contigs.fa'

with open('/Users/flashton/projects/H58/res_dict.pick', 'rb') as inhandle:
    minion_mapping = pickle.load(inhandle)

#print_res_dict(minion_mapping)



def read_walk(res_dict):

    for read in res_dict:
        minion_read = MinionRead()
        minion_read.read_name = res_dict[read].read_name

        if len(res_dict[read].hits) > 1:
            print res_dict[read].read_name, len(res_dict[read].hits)


read_walk(minion_mapping)
