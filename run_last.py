__author__ = 'flashton'

import os

reference = '/Users/flashton/Dropbox/H58_minion/data/refs/ST2_contigs'
query = '/Users/flashton/Dropbox/H58_minion/data/H566_ON/2014.08.11.H566_ON_both.fasta'
outhandle = '/Users/flashton/Dropbox/H58_minion/results/H566_ON/2014.08.11.H566_ON_both_vs_ST2_contigs'

def run_last(reference, query, outhandle):
    os.system('lastdb -Q 0 %s.lastindex %s.fasta' % (reference, reference))
    os.system('lastal -s 2 -T 0 -Q 0 -a 1 %s.lastindex %s > %s.last.txt' % (reference, query, outhandle))
    os.system('maf-convert.py sam %s.last.txt > %s.sam' % (outhandle, outhandle))
    os.system('samtools view -T %s.fa -bS %s.sam | samtools sort - %s.sorted' % (reference, outhandle, outhandle))
    os.system('samtools index %s.sorted.bam' % (outhandle))
    os.system('maf-convert.py blast %s.last.txt > %s.blast' % (outhandle, outhandle))



if __name__ == '__main__':
    run_last(reference, query, outhandle)