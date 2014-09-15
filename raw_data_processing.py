__author__ = 'flashton'

import os
import sys
import glob
import StringIO
import h5py
from Bio import SeqIO



################################# variables #################################

root_dir = '/Users/flashton/Dropbox/H58_minion/data/2014.07.26.Rabsch'

fast5_folder = '/Users/flashton/projects/H58/data/raw_3_Rabsch/downloads'



pt_types = ['all', 'fwd', 'rev', '2D']

keys = {'template' : '/Analyses/Basecall_2D_000/BaseCalled_template/Fastq',
        'complement' : '/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq',
        'twodirections' : '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'}

################################# functions #################################

def extract_fasta_nl():
    '''
    Use this one!
    '''

    for each_channel in os.listdir(root_dir):
        #outhandle = open('/Users/flashton/data/H58/raw/%s/%s.nl.minion.fasta' % (each_channel, each_channel), 'w')
        #if each_channel == 'ch99':
        fast5s = glob.glob('%s/%s/*.fast5' % (root_dir, each_channel))
        for each_file in fast5s:
            sn = each_file.split('/')[-1]
            fn = '%s' % each_file
            #print fn
            try:
                hdf = h5py.File(fn, 'r')
            except Exception, e:
                #print >>sys.stderr, "Error opening %s: %s" % (fn, e)
                continue

            for id, key in keys.iteritems():
                try:
                    fq = hdf[key][()]
                    #print fq
                    rec = SeqIO.read(StringIO(fq), "fastq")
                    rec.id += "_" + id
                    rec.description = sn
                    SeqIO.write([rec], sys.stdout, "fasta")
                except Exception, e:
                #	print >>sys.stderr, e
                    continue
            hdf.close()

def sort_data_out(root_dir):
    for each in os.listdir(root_dir):
        if not os.path.isdir('%s/%s' % (root_dir, each)):

            split_name = each.split('_')
            channel = split_name[-3]
            #print each
            #print split_name
            print channel

            if not os.path.exists('%s/%s' % (root_dir, channel)):
                os.system('mkdir %s/%s' % (root_dir, channel))

            os.system('mv %s/*%s_* %s/%s' % (root_dir, channel, root_dir, channel))

def single(inhandle):
    hdf = h5py.File(inhandle, 'r')
    print hdf.filename
    print hdf.keys()
    print hdf['Analyses'].keys()
    print hdf['Analyses']['Basecall_2D_000']['BaseCalled_template']['Fastq'][()]
    #fq = hdf[datasetname][()]
    #print fq
    hdf.close()

def extract_fastqs_pt(root_dir, fast5_folder):
    for type in pt_types:
        os.system('poretools fastq --type %s %s > %s/2014.07.26.Rabsch.%s.fastq' % (type, fast5_folder, root_dir, type))

################################# main #################################

extract_fastqs_pt(root_dir, fast5_folder)