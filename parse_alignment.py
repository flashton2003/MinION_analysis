from __future__ import division

__author__ = 'flashton'

import re
from __init__ import MinionRead, ReadContigMatch
try:
    import cPickle as pickle
except:
    import pickle

try:
    import numpy as np
except ImportError:
    print 'No Numpy, mapping_stats wont work but is not accessible via the command line interface anyway'


def parse_blast_text(blast_file):
    res_dict = {}
    read_name = ''
    with open(blast_file) as fi:
        for line in fi.readlines():
            line = line.strip()

            if line.startswith('Query='):
                ## want to refresh sbjct every time there is a new query
                sbjct = ''

                if read_name != '':
                    #read_contig_match.print_res(blast_res.read_name, blast_res.read_len)
                    read_contig_match.query_start = min(read_contig_match.query_coordinates)
                    read_contig_match.query_stop = max(read_contig_match.query_coordinates)
                    read_contig_match.sbjct_start = min(read_contig_match.sbjct_coordinates)
                    read_contig_match.sbjct_stop = max(read_contig_match.sbjct_coordinates)
                    minion_read.hits.append(read_contig_match)

                read_name = line.split(' ')[-1]
                minion_read = MinionRead()
                #print dir(blast_res)

                res_dict[read_name] = minion_read
                minion_read.read_name = read_name
                #print blast_res.read_name

            if line.startswith('('):
                read_len = int(line.replace('(', ' ').split(' ')[1])
                minion_read.read_len = read_len

            if line.startswith('>NODE'):
                if sbjct != '':

                    #read_contig_match.print_res(blast_res.read_name, blast_res.read_len)
                    read_contig_match.query_start = min(read_contig_match.query_coordinates)
                    read_contig_match.query_stop = max(read_contig_match.query_coordinates)
                    read_contig_match.sbjct_start = min(read_contig_match.sbjct_coordinates)
                    read_contig_match.sbjct_stop = max(read_contig_match.sbjct_coordinates)
                    minion_read.hits.append(read_contig_match)

                sbjct = line[1:]
                #print sbjct
                read_contig_match = ReadContigMatch()
                read_contig_match.sbjct = sbjct

            if line.startswith('Scor'):
                score = line.split(' ')[-1]
                read_contig_match.score = int(score)

            if line.startswith('Identi'):
                split_line = line.replace('/', ' ').split(' ')
                read_contig_match.match_len = int(split_line[3])
                read_contig_match.match_pos = int(split_line[2])
                try:
                    read_contig_match.match_gap = int(split_line[7])
                except IndexError:
                    read_contig_match.match_gap = 0

            if line.startswith('Query:'):
                split_line = line.split()

                read_contig_match.query_coordinates.append(int(split_line[1]))
                read_contig_match.query_coordinates.append(int(split_line[3]))

            if line.startswith('Sbjct:'):
                split_line = line.split()
                read_contig_match.sbjct_coordinates.append(int(split_line[1]))
                read_contig_match.sbjct_coordinates.append(int(split_line[3]))

            if line.startswith('Strand'):
                split_line = line.split()
                to_join = [split_line[2], split_line[4]]
                ori = '/'.join(to_join)
                read_contig_match.orientation = ori



        #read_contig_match.print_res(blast_res.read_name, blast_res.read_len)
        read_contig_match.query_start = min(read_contig_match.query_coordinates)
        read_contig_match.query_stop = max(read_contig_match.query_coordinates)
        read_contig_match.sbjct_start = min(read_contig_match.sbjct_coordinates)
        read_contig_match.sbjct_stop = max(read_contig_match.sbjct_coordinates)
        minion_read.hits.append(read_contig_match)

    return res_dict

def find_best_hits(res_dict):

    for read in res_dict:
        hit_dict = {}
        #print res_dict[read].read_name
        for hit in res_dict[read].hits:
            #print (hit.query_start, hit.query_stop)
            hit_dict[hit] = range(hit.query_start, hit.query_stop)
        #print hit_dict
        overlap_dict = {}
        for hit1 in hit_dict:
            overlap_dict[hit1] = {}
            for hit2 in hit_dict:
                if hit1 != hit2:
                    if hit2 in overlap_dict:
                        if hit1 not in overlap_dict[hit2]:
                            overlap_dict[hit1][hit2] = len(set(hit_dict[hit1]).intersection(hit_dict[hit2]))
        #print overlap_dict
        to_delete = []
        if len(overlap_dict) > 1:
            for hit1 in overlap_dict:
                for hit2 in overlap_dict[hit1]:
                    if hit1 != hit2:
                        #print hit1, hit1.match_len, hit2, hit2.match_len, overlap_dict[hit1][hit2]
                        hits = [hit1, hit2]
                        for h in hits:
                        #    #print h.match_len, overlap_dict[hit1][hit2] * 0.75
                            #print h.match_len
                            if h.match_len * 0.75 < overlap_dict[hit1][hit2]:
                                if hit1.score < hit2.score:
                                    to_delete.append(hit1)
                                if hit2.score < hit1.score:
                                    to_delete.append(hit2)
                                if hit1.score == hit2.score:
                                    to_delete = []
                                    to_delete.append(hit2)
        #print res_dict[read].hits
        #print 'to delete', to_delete
        for h in to_delete:
            for i, hit in enumerate(res_dict[read].hits):
                if hit == h:
                    #print len(res_dict[read].hits)
                    del res_dict[read].hits[i]
                    #print 'left', res_dict[read].read_name, res_dict[read].hits
    return res_dict



    #print(hit_dict)

def mapping_stats(res_dict):

    ## outdict{read_type:[accuracy, number reads aligned, gaps, length aligned], ...}

    #out_dict = {'template':[[], 0, [], 0], 'complement':[[], 0, [], 0], 'twodirections':[[], 0, [], 0]}

    out_dict = {'template':{'lAccuracy':[], 'iNumAligned':0, 'lGaps':[], 'iLengthAligned':0}, 'complement':{'lAccuracy':[], 'iNumAligned':0, 'lGaps':[], 'iLengthAligned':0}, 'twodirections':{'lAccuracy':[], 'iNumAligned':0, 'lGaps':[], 'iLengthAligned':0}}


    for read in res_dict:
        for type in out_dict:
            if re.search(type, read):
                out_dict[type]['iNumAligned'] += 1
                #print out_dict[type]
                for hit in res_dict[read].hits:
                    out_dict[type]['lGaps'].append(hit.match_gap / hit.match_len)
                    out_dict[type]['lAccuracy'].append(hit.match_pos / hit.match_len)
                    out_dict[type]['iLengthAligned'] += hit.match_len



    ## this code gives the total results
    total_acc = []
    total_gap = []
    total_len_aligned = 0
    total_number_reads = 0
    for each in out_dict:
        total_number_reads += out_dict[each]['iNumAligned']
        total_acc += out_dict[each]['lAccuracy']
        total_gap += out_dict[each]['lGaps']
        total_len_aligned += out_dict[each]['iLengthAligned']

    number_alignments = len(total_gap)

    print 'total', '\t', total_number_reads, '\t', number_alignments, '\t', total_len_aligned, '\t', np.mean(np.array(total_acc)), '\t', np.mean(np.array(total_gap))

    ## this code breaks down the results by type
    for each in out_dict:
        ## this gives the number of aligned reads, total number of alignments, mean accuracy of reads broken down by read type
        print each, '\t', out_dict[each]['iNumAligned'], '\t', len(out_dict[each]['lAccuracy']), '\t', out_dict[each]['iLengthAligned'], '\t', np.mean(np.array(out_dict[each]['lAccuracy'])), '\t', np.mean(np.array(out_dict[each]['lGaps']))
        ## this gives the number of aligned reads of each type and the total number of alignments of those reads
        #print each, out_dict[each][1], len(out_dict[each][0])

























