import os
from map_argparser import get_arguments
import seq_download
import map_rsv

def main():
    args = get_arguments()
    for filename in map_rsv.DATAFILES:
        if os.path.isfile(filename):
            pass
        else:
            seq_download.main()
    map_rsv.main(args.level, genotype_level=args.genotype_level, years=args.years)

if __name__ == '__main__':
	main()
