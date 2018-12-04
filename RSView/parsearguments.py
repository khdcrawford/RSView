"""Argument parser for `RSView`"""

import sys
import argparse
import RSView


class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""
    def error(self, message):
        """Prints error message, then help."""
        sys.stderr.write('error: {0}\n\n'.format(message))
        self.print_help()
        sys.exit(2)


def mapParser():
    """Returns argument parser for `map_rsv.py`"""
    parser = argparse.ArgumentParser(description="Plot global distribution of RSV")
    parser.add_argument(
        'level', type=str, choices=['subtype', 'genotype'],
        help="Specify whether the subtype or genotype of RSV sequences should be plotted")
    parser.add_argument(
        '--genotype-level', type=str, choices=['collapse', 'all'], default='collapse',
        help="Specify whether to plot all genotypes of RSV or collapse them into major clades")
    parser.add_argument(
        '--years', default=[1990,2018],
        help="Specify a range of years to plot")
    return parser

def seqParser():
    """Returns argument parser for `seq_download.py`."""
    parser = ArgumentParserNoArgHelp(
            description='Downloads RSV G protein sequences & metadata from '\
            'Genbank. This program is part of {0} (version {1}) written by '\
            '{2}.'.format(RSView.__name__, RSView.__version__, 
                    RSView.__author__) + 'For documentation, see '
                    '{0}'.format(RSView.__url__), 
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--email', required=True, type=str, help='User email'\
            ' for GenBank download.')
    parser.add_argument('--query', required=True, type=str, help='Search '\
            'term(s) for downloading sequences from GenBank.')
    parser.add_argument('--outdir', default='./data', type=str,
            help='Directory for downloaded data.')
    parser.add_argument('--outprefix', default='RSVgb_metadata', type=str, 
            help='Beginning of file name for output `.csv`. Suffix will '\
            'specify number of sequences downloaded.')
    parser.add_argument('--db', default='nuccore', type=str, help='Entrez'\
            ' database to search.')
    parser.add_argument('--firstseq', default=0, type=int, help='Index of '\
            'first sequence to download.')
    parser.add_argument('--filesize', default=5000, type=int, help='Number of '\
            'seqs to download into one file. Default of 5000 balances '\
            'download time and minimizing the number of separate files.')
    parser.add_argument('--maxseqs', default=20000, type=int, help='Maximum '\
            'number of sequences to download across all output files.')
    parser.add_argument('--batchsize', default=100, type=int, help='Number '\
            'seqs to download in one retrieval. If much larger than 100, '\
            'download will be quite slow.')
    parser.add_argument('--filetype', default='gb', type=str, help='File '\
            'type of sequences downloaded from GenBank.')
    parser.add_argument('--outmode', default='xml', type=str, help='File '\
            'type for results of GenBank query.')

    return parser

#def plotParser():

