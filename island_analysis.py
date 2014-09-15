__author__ = 'flashton'

try:
    import cPickle as pickle
except:
    import pickle

from minion_core import print_res_dict, find_best_hits

############################## variables ##############################

island_contigs = ['NODE_26_length_65708_cov_26.6671_ID_2073446', 'NODE_72_length_768_cov_81.9565_ID_2067122', 'NODE_61_length_2869_cov_9.97726_ID_2076452', 'NODE_71_length_820_cov_103.867_ID_2070328', 'NODE_64_length_2029_cov_7.17072_ID_2023902', 'NODE_25_length_66162_cov_25.497_ID_2076656', 'NODE_50_length_9081_cov_17.9961_ID_2075706', 'NODE_48_length_14772_cov_23.1447_ID_2076292', 'NODE_32_length_48020_cov_25.0245_ID_1886496', 'NODE_65_length_1933_cov_29.8807_ID_2073682', 'NODE_154_length_59_cov_109.75_ID_545390', 'NODE_56_length_4160_cov_15.6957_ID_2074936']

############################## functions ##############################

def reads_that_map_to_multiple_contigs(res_dict):
    out_dict = {}
    for read in res_dict:
        #print read, res_dict[read]
        #for hit in res_dict[read].hits:
        if len(set(res_dict[read].hits)) > 1:
            out_dict[read] = res_dict[read]
            #print read
    return out_dict

def reads_that_map_to_island(res_dict, island_contigs):
    out_dict = {}
    for read in res_dict:
        hit_sbjct = []
        for hit in res_dict[read].hits:
            hit_sbjct.append(hit.sbjct)
        hit_sbjct = set(hit_sbjct)
        if len(hit_sbjct.intersection(island_contigs)) > 1:
            out_dict[read] = res_dict[read]
            #print read + '\t' + str(len(hit_sbjct))
    return out_dict

############################## main ##############################

res_dict = pickle.load(open('/Users/flashton/Dropbox/H58_from_iMac/H58/res_dict.pick'))

res_dict = find_best_hits(res_dict)

multi_contig_reads = reads_that_map_to_multiple_contigs(res_dict)

multi_contig_island_reads = reads_that_map_to_island(multi_contig_reads, island_contigs)

print_res_dict(multi_contig_island_reads)