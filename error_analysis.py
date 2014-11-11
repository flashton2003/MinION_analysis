from __future__ import division

__author__ = 'flashton'

from Bio import SeqIO
import re
import numpy as np
import os


def find_indels(pileup):
    deleted_sections = []
    inserted_sections = []
    with open(pileup) as fi:
        for line in fi:
            split_line = line.split('\t')
            read_bases = split_line[4]
            m = re.findall('-([0-9]+)+([ACGTNacgtn]+)+', read_bases)
            if m:
                for each in m:
                    ## if there is a substitution after the indel, this regular expression will include it. this if statement
                    ## checks that the length of the nucleotide string matches the integer.
                    if len(each[1]) != int(each[0]):
                        deleted_sections.append(each[1][:int(len(each[0]))])
                    else:
                        deleted_sections.append(each[1])
            ## same for insertions
            n = re.findall('\+([0-9]+)+([ACGTNacgtn]+)+', read_bases)
            if n:
                for each in n:
                    if len(each[1]) != int(each[0]):
                        inserted_sections.append(each[1][:int(len(each[0]))])
                    else:
                        inserted_sections.append(each[1])

    ## commented out the histogram printing as not very interesting, but wanted to leave it in just in case
    #histogram_of_error_length(deleted_sections, 'deleted')
    #histogram_of_error_length(inserted_sections, 'inserted')


    for i, each in enumerate(deleted_sections):
        deleted_sections[i] = each.upper()
    for i, each in enumerate(inserted_sections):
        inserted_sections[i] = each.upper()

    ## basically just default dict on deletion/insertion lists.
    sorted_deletion = common_error_profiles(deleted_sections)
    sorted_insertion = common_error_profiles(inserted_sections)
    return sorted_deletion, sorted_insertion


def common_error_profiles(error_list):
    res_dict = {}
    for each in error_list:
        if each in res_dict:
            res_dict[each] += 1
        else:
            res_dict[each] = 1
    #sorted_res = sorted(res_dict.iteritems(), key=operator.itemgetter(1))
    return res_dict
    # for each in sorted_res:
    #     print each[0], '\t', each[1]


def find_kmer_freq(reference):
    '''
    This is quick but scales badly with increasing k mer range (runs out of ram)
    '''
    kmers = range(1, 7)
    res_dict = {}
    for kmer in kmers:
        for each_contig in SeqIO.parse(reference, 'fasta'):
            i = 0
            each_contig = str(each_contig.seq)
            for x in range(len(each_contig) - kmer)[kmer:]:
                seq = each_contig[i:x]
                if seq in res_dict:
                    res_dict[seq] += 1
                else:
                    res_dict[seq] = 1
                i += 1
    return res_dict


def slow_find_kmer_freq(reference, error_kmers):
    '''
    This is slow but has scale
    '''
    res_dict = {}
    for kmer in error_kmers:
        for each_contig in SeqIO.parse(reference, 'fasta'):
            each_contig = str(each_contig.seq)
            m = re.findall(kmer, each_contig)
            if kmer in res_dict:
                res_dict[kmer] += len(m)
            else:
                res_dict[kmer] = len(m)
    return res_dict


def characterise_deletions(deleted_kmers, typhi_kmers, output_dir):
    # iupac = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    ## normalise the deleted kmers by the occurence of that kmer in the reference genome
    error_proportion_by_kmer = {}
    for each in typhi_kmers:
        if each in deleted_kmers:
            error_proportion = deleted_kmers[each] / typhi_kmers[each]
            error_proportion_by_kmer[each] = error_proportion

    ## group the normalised deletion occurences by kmer length
    error_proportion_by_k_len = {}
    for kmer in error_proportion_by_kmer:
        if len(kmer) in error_proportion_by_k_len:
            error_proportion_by_k_len[len(kmer)].append(error_proportion_by_kmer[kmer])
        else:
            error_proportion_by_k_len[len(kmer)] = []
            error_proportion_by_k_len[len(kmer)].append(error_proportion_by_kmer[kmer])

    ## calculate the per-k-length mean and stdev of the deletion proportion,
    ## for comparison against each k-mer to allow calc
    ## of z-score
    mean_stdev_error_proportion_by_k_len = {}
    for each in error_proportion_by_k_len:
        a = np.array(error_proportion_by_k_len[each])
        mean_stdev_error_proportion_by_k_len[each] = []
        mean_stdev_error_proportion_by_k_len[each].append(np.mean(a))
        mean_stdev_error_proportion_by_k_len[each].append(np.std(a))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open('%s/deleted_kmers.txt' % (output_dir), 'w') as fo:
        fo.write('kmer\tkmer_length\tnumber_in_reference\tnumber_of_deletions\tdeletions_normalised_by_reference\tz-score\n')
        for each in error_proportion_by_kmer:
            fo.write(each + '\t' + str(len(each)) + '\t' + str(typhi_kmers[each]) + '\t'
                     + str(deleted_kmers[each]) + '\t' + str(error_proportion_by_kmer[each]) + '\t'
                     + str((error_proportion_by_kmer[each] - mean_stdev_error_proportion_by_k_len[len(each)][0]) /
                                                               mean_stdev_error_proportion_by_k_len[len(each)][1]) + '\n')



def analyse_insertions(inserted_kmers, output_dir):

    sorted_by_k_len = {}
    for kmer in inserted_kmers:
        if len(kmer) in sorted_by_k_len:
            sorted_by_k_len[len(kmer)].append(inserted_kmers[kmer])
        else:
            sorted_by_k_len[len(kmer)] = []
            sorted_by_k_len[len(kmer)].append(inserted_kmers[kmer])

    mean_stdev_error_proportion_by_k_len = {}

    for each in sorted_by_k_len:
        a = np.array(sorted_by_k_len[each])
        mean_stdev_error_proportion_by_k_len[each] = []
        mean_stdev_error_proportion_by_k_len[each].append(np.mean(a))
        mean_stdev_error_proportion_by_k_len[each].append(np.std(a))

    print mean_stdev_error_proportion_by_k_len

    print inserted_kmers

    with open('%s/inserted_kmers.txt' % (output_dir), 'w') as fo:
        for each in inserted_kmers:
            ## if the stdev (i.e. mean_stdev_error_proportion_by_k_len[len(each)][1]) is 0 then the z-score cannot be calculated
            ## this if statement stops the RuntimeError being raised
            if mean_stdev_error_proportion_by_k_len[len(each)][1] != 0:
                fo.write(str(each) + '\t' + str(len(each)) + '\t' + str(inserted_kmers[each]) + '\t' + str((inserted_kmers[each] -
                                                                                                        mean_stdev_error_proportion_by_k_len[len(each)][0]) / mean_stdev_error_proportion_by_k_len[len(each)][1]) + '\n')



def total_len_error(erroneous_kmers, error_type):
    total = 0
    l = []
    for each in erroneous_kmers:
        l.append(len(each))
        total += len(each) * erroneous_kmers[each]
    print error_type + '\t' + str(total)


def histogram_of_error_length(erroneous_kmers, name):
    error_lengths = []
    for x in erroneous_kmers:
        error_lengths.append(len(x))
    hist, bins = np.histogram(error_lengths, bins = range(0, max(error_lengths)))
    print '=' * 20
    #for x, y in zip(hist, bins):
    #    print x, y
    a = np.array(error_lengths)
    print np.median(a)
    print np.mean(a)


def find_substitutions(pileup):
    res_dict = {'A':{'T':0, 'C':0, 'G':0}, 'T':{'A':0, 'C':0, 'G':0}, 'C':{'A':0, 'T':0, 'G':0}, 'G':{'A':0, 'T':0, 'C':0}}
    with open(pileup) as fi:
        for line in fi:
            split_line = line.split('\t')
            ref_base = split_line[2]
            read_bases = split_line[4]
            ## substitutions are all the characters (ATCG) in the 4th column of the pileup that aren't indels.
            ## first, check if there is an indel in the 4th column
            m = re.findall('([-+][0-9]+)+([ACGTNacgtn]+)+', read_bases)
            if m:
                for each in m:
                    ## if there is, remove that indel from the 4th column
                    if len(each[1]) != int(each[0][1:]):
                        ## however, need to check that there isn't a substitution next to the indel, so count the nucleotide
                        ## string and check is equivalent to the integer following the +/-
                        cor_str =  each[1][:int(each[0][1:])]
                        e = ''.join([each[0], cor_str])
                        read_bases = read_bases.replace(e, '')
                    else:
                        e = ''.join([each[0], each[1]])
                        read_bases = read_bases.replace(e, '')
                read_bases = read_bases.upper()
                for b in read_bases:
                    if b in ['A', 'C', 'T', 'G']:
                        res_dict[ref_base][b] += 1
            else:
                ## if there is no indel, can do a simple count of nucleotide characters
                read_bases = read_bases.upper()
                for b in read_bases:
                    if b in ['A', 'C', 'T', 'G']:
                        res_dict[ref_base][b] += 1
    print 'Substitutions {ref_base:{changed_to:frequency}}'
    print res_dict