"""Script for assigning genotypes"""
import pandas as pd
import os
import glob
import RSView.parsearguments


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    assert len(s1) == len(s2), "Undefined for sequences of unequal length"
    return sum(site1 != site2 for site1, site2 in zip(s1, s2))


def main():

    parser = RSView.parsearguments.genotypeParser()
    args = vars(parser.parse_args())
    prog = parser.prog

    print("\nExecuting {0} ({1}) in {2} at {3}.\n".format(
            prog, RSView.__version__, os.getcwd(), time.asctime()))

    files = [filename for filename in glob.glob(args['infile'])]

    print(files)

    file_frames = []
    for file in files:
        if not os.path.isfile(file):
            raise ValueError('Sequence data not downloaded. Run '\
                    '`seq_download.py`.')
        file_df=pd.read_csv(file)
        file_frames.append(fild_df)

    rsv_df = pd.concat(file_frames, ignore_index=True, sort=False)

    print(rsv_df)


if __name__ == '__main__':
    main()