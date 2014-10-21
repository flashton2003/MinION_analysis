__author__ = 'flashton'

import os


def run_last(reference, query, outdir):
    reference_path = os.path.splitext(reference)[0]
    reference_name = os.path.basename(reference_path)
    query_name = os.path.splitext(os.path.basename(query))[0]
    print '### indexing reference for last alignment ###'
    os.system('lastdb -Q 0 %s.lastindex %s' % (reference_path, reference))
    print '### LAST is aligning the query against the reference ###'
    os.system('lastal -s 2 -T 0 -Q 0 -a 1 %s.lastindex %s > %s/%s_vs_%s.last.txt' % (reference_path, query, outdir,
    query_name, reference_name))
    print '### converting the last format to blast format using maf-convert.py ###'
    os.system('maf-convert.py blast %s/%s_vs_%s.last.txt > %s/%s_vs_%s.blast.txt' % (outdir, query_name, reference_name, outdir,
                                                                                   query_name, reference_name))
    os.system('maf-convert.py sam %s/%s_vs_%s.last.txt > %s/%s_vs_%s.sam' % (outdir, query_name, reference_name, outdir, query_name, reference_name))
    os.system('samtools view -T %s -bS %s/%s_vs_%s.sam | samtools sort - %s/%s_vs_%s.sorted' % (reference, outdir, query_name, reference_name, outdir, query_name, reference_name))
    os.system('samtools index %s/%s_vs_%s.sorted.bam' % (outdir, query_name, reference_name))
    return '%s/%s_vs_%s.sorted.bam' % (outdir, query_name, reference_name)

def make_pileup(reference, input_file):
    print 'Running mpileup'
    path_and_filename = os.path.splitext(os.path.splitext(input_file)[0])[0]
    os.system('samtools mpileup -BQ0 -f %s %s > %s.pileup' % (reference, input_file, path_and_filename))