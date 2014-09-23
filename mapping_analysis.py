from __future__ import division

__author__ = 'flashton'

from minion_core import BlastRes, BlastHit, find_best_hits

try:
    import cPickle as pickle
except:
    import pickle

import re
import numpy as np
from Bio import SeqIO
############################## variables ##############################

raw = '/Users/flashton/Dropbox/H58_from_iMac/H58/data/H566_30min/2014.08.19.H566_30min.mapped.fastq'

blast_format = '/Users/flashton/Dropbox/ulf_minion/TXT_last_map_against_spades_contigs_unmapped-EcWSU1.blast'

pileup = '/Users/flashton/Dropbox/minion/results/minion_vs_ct18.pileup'
reference = '/Users/flashton/Dropbox/minion/data/ct18.fa'


############################## functions ##############################

def parse_blast_text(blast_file):
    """

    :rtype : dictionary
    """
    res_dict = {}
    read_name = ''
    with open(blast_file) as fi:
        for line in fi.readlines():
            line = line.strip()

            if line.startswith('Query='):
                ## want to refresh sbjct every time there is a new query
                sbjct = ''

                if read_name != '':
                    #blast_hit.print_res(blast_res.read_name, blast_res.read_len)
                    blast_hit.query_start = min(blast_hit.query_coordinates)
                    blast_hit.query_stop = max(blast_hit.query_coordinates)
                    blast_hit.sbjct_start = min(blast_hit.sbjct_coordinates)
                    blast_hit.sbjct_stop = max(blast_hit.sbjct_coordinates)
                    blast_res.hits.append(blast_hit)

                read_name = line.split(' ')[-1]
                blast_res = BlastRes()
                #print dir(blast_res)

                res_dict[read_name] = blast_res
                blast_res.read_name = read_name
                #print blast_res.read_name

            if line.startswith('('):
                read_len = int(line.replace('(', ' ').split(' ')[1])
                blast_res.read_len = read_len

            if line.startswith('>NODE'):
                if sbjct != '':

                    #blast_hit.print_res(blast_res.read_name, blast_res.read_len)
                    blast_hit.query_start = min(blast_hit.query_coordinates)
                    blast_hit.query_stop = max(blast_hit.query_coordinates)
                    blast_hit.sbjct_start = min(blast_hit.sbjct_coordinates)
                    blast_hit.sbjct_stop = max(blast_hit.sbjct_coordinates)
                    blast_res.hits.append(blast_hit)

                sbjct = line[1:]
                #print sbjct
                blast_hit = BlastHit()
                blast_hit.sbjct = sbjct

            if line.startswith('Scor'):
                score = line.split(' ')[-1]
                blast_hit.score = int(score)

            if line.startswith('Identi'):
                split_line = line.replace('/', ' ').split(' ')
                blast_hit.match_len = int(split_line[3])
                blast_hit.match_pos = int(split_line[2])
                try:
                    blast_hit.match_gap = int(split_line[7])
                except IndexError:
                    blast_hit.match_gap = 0

            if line.startswith('Query:'):
                split_line = line.split()

                blast_hit.query_coordinates.append(int(split_line[1]))
                blast_hit.query_coordinates.append(int(split_line[3]))

            if line.startswith('Sbjct:'):
                split_line = line.split()
                blast_hit.sbjct_coordinates.append(int(split_line[1]))
                blast_hit.sbjct_coordinates.append(int(split_line[3]))

            if line.startswith('Strand'):
                split_line = line.split()
                to_join = [split_line[2], split_line[4]]
                ori = '/'.join(to_join)
                blast_hit.orientation = ori



        #blast_hit.print_res(blast_res.read_name, blast_res.read_len)
        blast_hit.query_start = min(blast_hit.query_coordinates)
        blast_hit.query_stop = max(blast_hit.query_coordinates)
        blast_hit.sbjct_start = min(blast_hit.sbjct_coordinates)
        blast_hit.sbjct_stop = max(blast_hit.sbjct_coordinates)
        blast_res.hits.append(blast_hit)

    return res_dict

def print_res_dict(res_dict):
    """

    :param res_dict:
    """
    #outhandle = open('/Users/flashton/Dropbox/H58_from_iMac/H58/blast_txt', 'w')
    print 'query    read len	subject	orientation	score	match len	match pos	match gap	q start	q stop	s start	s stop'
    for read in res_dict:
    #print each, res_dict[each]
    #print res_dict[each].hits
        for every in res_dict[read].hits:
            print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (res_dict[read].read_name, res_dict[read].read_len,
                                                                    every.sbjct, every.orientation, every.score, every.match_len, every.match_pos, every.match_gap, every.query_start, every.query_stop, every.sbjct_start, every.sbjct_stop))
    #outhandle.close()

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

def pull_out_mapped_reads(raw, res_dict):
    outhandle_temp = '/Users/flashton/Dropbox/H58_from_iMac/H58/data/H566_30min/2014.08.19.H566_30min.mapped.template.fastq'
    outhandle_comp = '/Users/flashton/Dropbox/H58_from_iMac/H58/data/H566_30min/2014.08.19.H566_30min.mapped.complement.fastq'
    outhandle_twod = '/Users/flashton/Dropbox/H58_from_iMac/H58/data/H566_30min/2014.08.19.H566_30min.mapped.2D.fastq'

    out_dict = {'complement':[[], outhandle_comp], 'template':[[], outhandle_temp], 'twodirections':[[], outhandle_twod]}

    i = 0
    out_list = []
    for rec in SeqIO.parse(raw, 'fastq'):
        read_name = rec.id.split(':')[0]

        if read_name in res_dict:
            for each in out_dict:
                if re.search(each, read_name):

                    out_dict[each][0].append(rec)
            #SeqIO.append([rec], outhandle, 'fastq')
            i += 1

    for each in out_dict:
        SeqIO.write(out_dict[each][0], out_dict[each][1], 'fastq')

    print i


############################## functions ##############################

res_dict = parse_blast_text(blast_format)
res_dict = find_best_hits(res_dict)
print_res_dict(res_dict)

#with open('/Users/flashton/Dropbox/H58_from_iMac/H58/res_dict.pick', 'wb') as outhandle:
#    pickle.dump(res_dict, outhandle)

#res_dict = pickle.load(open('/Users/flashton/Dropbox/H58_from_iMac/H58/res_dict.pick'))

#mapping_stats(res_dict)

#pull_out_mapped_reads(raw, res_dict)























