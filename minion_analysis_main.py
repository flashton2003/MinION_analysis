__author__ = 'flashton'

'''
This section of the code is heavily influenced by the design of Aaron Quinlan and Nick Loman's poretools package - check it out!
https://github.com/arq5x/poretools
'''
import argparse
import os
from __init__ import __version__, print_res_dict
import run_last
import parse_alignment
import error_analysis


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr", action="store_true", dest="quiet")


def run_command(args):
    if args.command == 'run_last':
        if os.path.exists(args.output_dir):
            pass
        else:
            os.makedirs(args.output_dir)
        run_last.run_last(args.reference, args.minion_reads, args.output_dir)

    if args.command == 'parse_last_output':
        res_dict = parse_alignment.parse_blast_text(args.input_file)
        res_dict = parse_alignment.find_best_hits(res_dict)
        parse_alignment.mapping_stats(res_dict)
        print_res_dict(res_dict)

    if args.command == 'error_profile':
        error_analysis.make_pileup(args.reference, args.input_file)
        deleted_kmers = error_analysis.error_profile(args.input_file)
        ref_kmers = error_analysis.find_kmer_freq(args.reference)
        error_analysis.compare_kmers(deleted_kmers, ref_kmers)


def main():
    parser = argparse.ArgumentParser(prog='minion_analysis', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed minion_analysis version", action="version",
                        version="%(prog)s " + str(__version__))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    parser_run_last = subparsers.add_parser('run_last', help='Takes MinION reads and maps then to a reference using LAST')
    parser_run_last.add_argument('reference', help='Assembly of Illumina data from the same organism as the MinION reads')
    parser_run_last.add_argument('minion_reads', help='MinION reads in FASTA format')
    parser_run_last.add_argument('output_dir', help='Where the output of the LAST alignment and subsequent processing '
                                                    'will be stored')

    parser_parse_last_output = subparsers.add_parser('parse_last_output', help='Parses LAST output (in BLAST text format '
                                                                               'produced by maf-convert) to produce a format '
                                                                               'useful for scaffolding')
    parser_parse_last_output.add_argument('input_file', help='The BLAST text (*.blast.txt) formatted LAST output (will be in '
                                                             'the output directory you passed to run_last')

    parser_error_profile = subparsers.add_parser('error_profile', help='Converts LAST output to pileup, identifies different '
                                                                       'error types (indel, miscalled bases)')
    parser_error_profile.add_argument('input_file', help='sorted bam of the LAST alignment')
    parser_error_profile.add_argument('reference', help='Reference genome that was aligned to')

    args = parser.parse_args()
    run_command(args)


if __name__ == "__main__":
    main()