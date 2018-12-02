import argparse

def get_arguments():
    PARSER = argparse.ArgumentParser(description="Plot global distribution of RSV")
    PARSER.add_argument(
        'level', type=str, choices=['subtype', 'genotype'],
        help="Specify whether the subtype or genotype of RSV sequences should be plotted")
    PARSER.add_argument(
        '--genotype-level', type=str, choices=['collapse', 'all'], default='collapse',
        help="Specify whether to plot all genotypes of RSV or collapse them into major clades")
    PARSER.add_argument(
        '--years', default=[1990,2018],
        help="Specify a range of years to plot")
    ARGS = PARSER.parse_args()

    return ARGS
