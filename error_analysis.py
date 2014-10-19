from __future__ import division

__author__ = 'flashton'

from Bio import SeqIO
import re
import os

def make_pileup(reference, input_file):
    print 'Running mpileup'
    path_and_filename = os.path.splitext(os.path.splitext(input_file)[0])[0]
    os.system('samtools mpileup -BQ0 -f %s %s > %s.pileup' % (reference, input_file, path_and_filename))

def error_profile(input_file):
    pileup = '{0}.pileup'.format(os.path.splitext(os.path.splitext(input_file)[0])[0])
    deleted_sections = []
    inserted_sections = []
    with open(pileup) as fi:
        for line in fi.readlines():
            split_line = line.split('\t')
            read_bases = split_line[4]
            #m = re.search('\-', read_bases)
            m = re.search('-([0-9]+)+([ACGTNacgtn]+)+', read_bases)
            if m:
                if len(m.groups()[1]) != int(m.groups()[0]):
                    deleted_sections.append(m.groups()[1][:int(len(m.groups()[0]))])
                else:
                    #print line
                    #print m.groups()
                    #print m.groups()[1][:int(len(m.groups()[0]))]
                    deleted_sections.append(m.groups()[1])

            # n = re.search('\+', read_bases)
            # if n:
            #     split_read_bases = list(read_bases)
            #     indicies = [i for i, x in enumerate(split_read_bases) if x == '+']
            #     for i in indicies:
            #         ins_len = int(split_read_bases[i + 1])
            #         #print del_len
            #         start = i + 2
            #         stop = i + 2 + ins_len
            #         inserted_section = split_read_bases[start:stop]
            #         inserted_sections.append(''.join(inserted_section))
    for i, each in enumerate(deleted_sections):
        deleted_sections[i] = each.upper()
    #print len(deleted_sections)
    sorted_deletion = common_error_profiles(deleted_sections)
    #print len(inserted_sections)
    # sorted_insertion = common_error_profiles(inserted_sections)
    return sorted_deletion

def common_error_profiles(deleted_sections):
    res_dict = {}
    for each in deleted_sections:
        if each in res_dict:
            res_dict[each] += 1
        else:
            res_dict[each] = 1
    #sorted_res = sorted(res_dict.iteritems(), key=operator.itemgetter(1))
    return res_dict
    # for each in sorted_res:
    #     print each[0], '\t', each[1]

def find_kmer_freq(reference):
    kmers = range(1, 21)
    res_dict = {}
    for kmer in kmers:
        for each_contig in SeqIO.parse(reference, 'fasta'):
            i = 0
            each_contig = str(each_contig.seq)
            #print range(len(each_contig) - kmer)[kmer:]
            for x in range(len(each_contig) - kmer)[kmer:]:
                #print i, x
                seq = each_contig[i:x]
                #print seq
                if seq in res_dict:
                    res_dict[seq] += 1
                else:
                    res_dict[seq] = 1
                i += 1
    return res_dict
        # for each in res_dict:
        #     print each, res_dict[each]

def compare_kmers(error_kmers, typhi_kmers):
    res_dict = {}
    iupac = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    for each in typhi_kmers:
        if each in error_kmers:
            error_proportion = error_kmers[each] / typhi_kmers[each]
            res_dict[each] = error_proportion
    total = 0
    with open('/Users/flashton/missing_20_mers.txt', 'w') as fo:
        for each in res_dict:
        # if any(x in each for x in iupac):
        #     pass
        # else:
            total += len(each) * error_kmers[each]
            fo.write(each + '\t' + str(len(each)) + '\t' + str(typhi_kmers[each]) + '\t' + str(error_kmers[each]) + \
                            '\t' + str(res_dict[each]) + '\n')
    print total
