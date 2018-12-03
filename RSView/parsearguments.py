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
    parser.add_argument('--email', required=True, type=str, help='User email for '\
            'GenBank download.')
    parser.add_argument('--outprefix', default='./data/RSVgb_metadata_', 
            type=string, help='Path and beginning of file name for output '\
            '`.csv`. Suffix will specify number of seqs downloaded.')
    parser.add_argument('--db', default='nuccore', type=str, help='Entrez'\
            ' database to search.')
    parser.add_argument('--retnum', default=5000, type=int, help='')
    parser.add_argument('--batchsize', )
    parser.add_argument('--query', required=True, type=str, help='')
    return parser

#def plotParser():

