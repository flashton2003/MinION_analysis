__author__ = 'flashton'

from Bio import SeqIO
import re
import glob
import numpy as np
import os
import operator

######################## variables etc ########################

script_dir = os.path.dirname(os.path.realpath(__file__))

raw = '/Users/flashton/Dropbox/H58_minion/data/2014.07.26.Rabsch/2014.07.26.Rabsch.fasta'
root_dir = '/Users/flashton/Dropbox/H58_minion/data/2014.07.26.Rabsch_R7'
fast5_folder = '/Users/flashton/projects/H58/H566_30min_inc'

types = ['template', 'complement', 'twodirection']
pt_types = ['all', 'fwd', 'rev', '2D']

fastas = ['2014.08.11.H566_30_min_inc.complement.fasta', '2014.08.11.H566_30_min_inc.fasta', '2014.08.11.H566_30_min_inc.template.fasta', '2014.08.11.H566_30_min_inc.twod.fasta']

######################## functions ########################

def mapping_stats(blast_text):
    '''
    blast format will include multiple hits for each read
    '''
    gaps = []
    identities = []
    lens = []
    res_dict = {}

    with open(blast_text) as fi:
        for line in fi.readlines():
            line = line.strip()
            m = re.search('Iden', line)
            if m:
                try:
                    #print line

                    identity = int(line.split(' ')[3].replace('%', '(').split('(')[1])
                    identities.append(identity)
                    gap = int(line.split(' ')[7].replace('%', '(').split('(')[1])
                    gaps.append(gap)
                    match_len = int(line.split(' ')[2].split('/')[1])
                    lens.append(match_len)
                except:
                    #print line
                    pass


    i = np.array(identities)
    g = np.array(gaps)
    l = np.array(lens)

    print np.median(i)
    print np.median(g)
    print np.median(l)
    print np.max(l)
    print np.sum(l)
    print len(i)

def compare_raw_mapped_reads(mapped):
    read_len = []
    res_dict = {}

    raw_names = []

    for each_contig in SeqIO.parse(raw, 'fasta'):
        raw_names.append(each_contig.id)
        res_dict[each_contig.id] = {}
        res_dict[each_contig.id]['raw'] = []
        res_dict[each_contig.id]['raw'].append(len(each_contig.seq))
        #read_len.append(len(each_contig.seq))
        #print len(each_contig.seq)

    mapped_names = []

    j = 0
    for each_contig in SeqIO.parse(mapped, 'fasta'):
        mapped_names.append(each_contig.id)

        if 'mapped' not in res_dict[each_contig.id]:
            res_dict[each_contig.id]['mapped'] = []
            res_dict[each_contig.id]['mapped'].append(len(each_contig.seq))

        else:
            res_dict[each_contig.id]['mapped'].append(len(each_contig.seq))

    #print len(mapped_names)
    #print len(set(mapped_names))



    temp_dict = {}

    for each in mapped_names:
        if each in temp_dict:
            temp_dict[each] += 1
        else:
            temp_dict[each] = 1

    for each in temp_dict:
        if temp_dict[each] == 2:
            print each

    '''
    for each in res_dict:
        if 'mapped' in res_dict[each]:
            if len(res_dict[each]['mapped']) > 1:
                print each, res_dict[each]


    raw_len = []
    mapped_len = []

    for each in res_dict:
        if len(res_dict[each]) > 1:
            #print each, res_dict[each]
            raw_len.append(res_dict[each]['raw'])
            mapped_len.append(res_dict[each]['mapped'])
            #print 'match'

    i = 0
    for each in res_dict:
        if 'mapped' in res_dict[each]:
            i += 1

    a = np.array(raw_len)
    b = np.array(mapped_len)

    print len(raw_len), len(mapped_len)
    print np.median(a), np.median(b)
    #print np.max(a)
    #print np.sum(a)
    '''

def error_profile(pileup):
    deleted_sections = []
    with open(pileup) as fi:
        for line in fi.readlines():
            split_line = line.split('\t')
            read_bases = split_line[4]
            m = re.search('-', read_bases)
            if m:
                split_read_bases = list(read_bases)
                indicies = [i for i, x in enumerate(split_read_bases) if x == '-']


                #print split_read_bases
                #print indicies
                for i in indicies:
                    del_len = int(split_read_bases[i + 1])
                    #print del_len
                    start = i + 2
                    stop = i + 2 + del_len
                    deleted_section = split_read_bases[start:stop]
                    deleted_sections.append(''.join(deleted_section))

    for i, each in enumerate(deleted_sections):
        deleted_sections[i] = each.upper()

    return deleted_sections

def common_error_profiles(deleted_sections):
    res_dict = {}

    for each in deleted_sections:
        if each in res_dict:
            res_dict[each] += 1
        else:
            res_dict[each] = 1

    sorted_res = sorted(res_dict.iteritems(), key=operator.itemgetter(1))

    return res_dict

    #for each in sorted_res:
    #	print each[0], '\t', each[1]

def find_nucleotide_freq(reference):
    for each_contig in SeqIO.parse(reference, 'fasta'):
        #print len(each_contig)
        each_contig = str(each_contig.seq)
        #print each_contig
        #typhi = Seq(each_contig)
        print each_contig.count('A')
        print each_contig.count('T')
        print each_contig.count('C')
        print each_contig.count('G')

def find_kmer_freq(reference):
    kmers = range(1, 6)
    res_dict = {}
    for kmer in kmers:
        print kmer
        i = 0

        for each_contig in SeqIO.parse(reference, 'fasta'):
            each_contig = str(each_contig.seq)

            for x in range(len(each_contig) - kmer)[kmer:]:
                #print x
                seq = each_contig[i:x]

                #print seq
                if seq in res_dict:
                    res_dict[seq] += 1
                else:
                    res_dict[seq] = 1

                i += 1
    return res_dict
        #for each in res_dict:
        #	print each, res_dict[each]

def compare_kmers(error_kmers, typhi_kmers):
    res_dict = {}
    iupac = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N']
    for each in typhi_kmers:
        try:
            if each in error_kmers:
                error_proportion = error_kmers[each] / typhi_kmers[each]
                res_dict[each] = error_proportion
        except:
            pass
    for each in res_dict:
        try:
            if any(x in each for x in iupac):
                pass
            else:
                print each, typhi_kmers[each], error_kmers[each], res_dict[each]
        except:
            print

def how_many_2d_reads(raw):
    res_dict = {}
    with open(raw) as fi:
        for line in fi.readlines():
            if line.startswith('>'):
                read_type = line.split(' ')[0].split('_')[-1]
                #read_type = split_line[-1]
                #print read_type

                if read_type in res_dict:
                    res_dict[read_type] += 1
                else:
                    res_dict[read_type]= 1

    for each in res_dict:
        print each, res_dict[each]

def divide_reads_by_type(raw, types):
    root_name = '.'.join(raw.split('.')[:-1])
    print root_name
    template_out = open('%s.template.fasta' % root_name, 'w')
    complement_out = open('%s.complement.fasta' % root_name, 'w')
    twod_out = open('%s.twod.fasta' % root_name, 'w')

    out_dict = {'template':template_out, 'complement':complement_out, 'twodirection':twod_out}

    for each_read in SeqIO.parse(raw, 'fasta'):
        for type in types:
            m = re.search(type, each_read.id)
            if m:
                SeqIO.write(each_read, out_dict[type], 'fasta')
                #print type, each_read.id

def read_length(raw):

    out_dict = {'template':[], 'complement':[], 'twodirections':[]}

    for each_read in SeqIO.parse(raw, 'fasta'):
        for each in out_dict:
            if re.search(each, each_read.id):
                out_dict[each].append(len(each_read.seq))


    total = []
    for type in out_dict:
        total += out_dict[type]
        a = np.array(out_dict[type])

        print type, '\t', len(a), '\t', np.median(a), '\t', np.max(a), '\t', np.sum(a)

    a = np.array(total)
    print 'total', '\t', len(a), '\t', np.median(a), '\t', np.max(a), '\t', np.sum(a)

    #for each in out_dict:
    #    for every in out_dict[each]:
    #        print each, every

def get_read_qual(root_dir):

    res = glob.glob('%s/*fastq' % root_dir)

    for each in res:

        lPhred = []
        lAccuracy = []

        #outhandle = '.'.join(each.split('.')[:-1])
        sample = each.split('/')[-1]

        fFile = open(each)
        while True:
            sLine = fFile.readline()
            sLine = fFile.readline()
            sLine = fFile.readline()
            sLine = fFile.readline().strip()
            if not sLine:
                break
            for c in sLine:

                phred = float(ord(c)-33)
                lPhred.append(phred)

                accuracy = pow(10.0, -1.0*(phred/10.0))
                lAccuracy.append(accuracy)

        aPhred = np.array(lPhred)
        aAccuracy = np.array(lAccuracy)

        print sample, np.median(aPhred), np.median(aAccuracy)




######################### main ########################




#read_length(raw)

get_read_qual(root_dir)

