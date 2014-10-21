MinION_analysis
===============

Useage: python minion_analysis_main.py -h

There are two commands
  run_last - aligns a set of Nanopore reads in a fasta file against a reference genome (typically, an illumina assembly of the same isolate), converts the last output into BLAST format using maf-convert.py (bundled with LAST)
  parse_last_output - parses the LAST alignment which has been converted to BLAST format file produced by maf-convert.py and outputs a tab delimited file
  
The example data presented here (specifically the Illumina assembly) was produced using SPAdes v.2.5.1 and is not the final assembly in the associated paper (need to insert reference to paper).
  

Output of parse_last_output
==============================
The output of parse_last_output is intended as a flexible output to be useful for downstream analysis. In our workflow, it is loaded into Excel and sorted so that matches to the same MinION read are grouped together. Within the matches for each MinION read, they are sorted by their query start position. Additionally, it is useful to have the reads sorted by the number of different contigs matched, in our case, we are most interested in reads that match more than 2 contigs.

This setup allows the analyser to focus on the reads that map to multiple contigs. Here, the analytical process for one read from our example dataset will be gone through in detail, hopefully allowing you to replicate our analysis. The read we choose is 'channel\_473\_read\_80\_twodirections', which will henceforth be refered to as the read of interest (ROI).

The ROI maps to 4 different contigs in our Illumina assembly (nodes 25, 72, 61 and 50). As we have the matches ordered by the query start position, you should should be able to see that the ROI matches node 25 between position 31 (on the ROI) and 4541 (also on ROI). This match is to positions 61406 to 66162 on the contig (the end of the contig). The next match is between ROI position 4499 and 5232 and node 72 position 1-757 (entire length of node 72). In the same region of the ROI, there is another match to node 61 which is much shorter (44 bp compared with 785 bp), this match will be ignored. The next match is between ROI position 5181 and 12051 to node 50. The match in node 50 is position 2043-9080 (which is the end of node 50). 

Therefore, this ROI provides evidence that node 25 and node 50 are neighbours in the genome and that the short read genome assembly has been interupted by node 72 (an insertion sequence).


Requirements
===================

The LAST aligner needs to be installed (available via homebrew on Mac OS X) and in the path. Similarly, the maf-convert.py script bundled with LAST needs to be in the path.
Python >= 2.7.6
Numpy >= 1.6.2
