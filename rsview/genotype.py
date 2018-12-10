"""Assign RSV G sequences genotypes.

This program relies on some seqs already being genotyped. It will
only add sequences to a genotype that already has a high quality
(< 60 gaps in the alignment) sequence in the input data.
"""

import time
import os
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
import rsview.parsearguments


#GT_A and GT_B lists based on genotypes already known to be present.
GT_A_LIST = ['ON1', 'GA1', 'GA2', 'GA3', 'GA5', 'GA6', 'GA7', 'NA1', 'NA2',
             'NA3', 'SAA1', 'SAA2']
GT_B_LIST = ['BA', 'BA2', 'BA4', 'BA5', 'BA7', 'BA8', 'BA9', 'BA10', 'BA11',
             'BA12', 'BA14', 'THB', 'SAB1', 'SAB3', 'SAB4', 'GB1', 'GB2',
             'GB3', 'GB4', 'GB13', 'GB12']

def hamming_distance(seq1, seq2):
    """Return the Hamming distance between equal-length sequences."""
    assert len(seq1) == len(seq2), "Undefined for sequences of unequal length"
    return sum(site1 != site2 for site1, site2 in zip(seq1, seq2))


def merge_csvs(csv_files):
    """Merge given *csv_files* into one dataframe.

    Args:
        `csv_files` (list)
            list of .csv files to merge into dataframe
            
    Returns:
        `full_df` (dataframe)
            pandas dataframe containing data concatenated from csv_files
    """
    file_frames = []
    for file in csv_files:
        if not os.path.isfile(file):
            raise ValueError('Sequence data not downloaded. Run '\
                    '`seq_download.py`.')
        file_df = pd.read_csv(file)
        file_frames.append(file_df)

    full_df = pd.concat(file_frames, ignore_index=True, sort=False)
    full_df.drop('Unnamed: 0', axis=1, inplace=True)

    return full_df


def seqstofastas(seqs_df, outfiles):
    """
    Takes a dataframe containing sequence, subtype, and genotype info and
    outputs specified `.fasta` files.

    Args:
        `outfiles` (list)
            List of `.fasta` files to output. Expects 3 files:
                1. Long G seqs already genotyped
                2. Long G seqs not yet genotyped
                3. All short G seqs

        `seqs_df` (dataframe)
            pandas dataframe containing sequence, subtype, and genotype info
    """
    assert len(outfiles) == 3, 'Unexpected number of files to output.'

    lt_fasta = outfiles[0]
    l_fasta = outfiles[1]
    s_fasta = outfiles[2]

    with open(lt_fasta, 'w') as longtyped, \
         open(l_fasta, 'w') as long_nogt, \
         open(s_fasta, 'w') as short:
        for i in range(len(seqs_df)):
            stype = seqs_df.at[i, 'subtype']
            gtype = seqs_df.at[i, 'genotype']
            header = '>{0} {1} {2}\n'.format(i, stype, gtype)
            seq = seqs_df.at[i, 'G_seq']
            if seq != 'NaN':
                if gtype != 'NaN':
                    if len(seq) > 290:
                        longtyped.write('{0}{1}\n'.format(header, seq))
                    else:
                        short.write('{0}{1}\n'.format(header, seq))
                else:
                    if len(seq) > 290:
                        long_nogt.write('{0}{1}\n'.format(header, seq))
                    else:
                        short.write('{0}{1}\n'.format(header, seq))


def align_seqs(infiles, outfiles):
    """Use mafft to align RSV G sequences.

    Due to significant disparities in length, use a three step approach.
        1. Align long G sequences with known genotypes
        2. Add in all long G sequences using --add with --keeplength
        3. Add in short G sequences using --addfragments with --keeplength

    Expects to align 1st infile and ouput to first outfile, then add 2nd
    infile and output as 2nd outfile, and finally add 3rd infile and output
    as 3rd outfile.
    """

    assert len(infiles) == len(outfiles) == 3, 'Incorrect number of files '\
            'for i/o.'

    subprocess.check_call('mafft --auto {0} > {1}'.format(infiles[0],
            outfiles[0]), shell=True)

    subprocess.check_call('mafft --add {0} --reorder --keeplength {1} > {2}'\
            .format(infiles[1], outfiles[0], outfiles[1]), shell=True)

    subprocess.check_call('mafft --addfragments {0} --reorder --6merpair '\
            '--thread -1 --keeplength {1} > {2}'.format(infiles[2],
            outfiles[1], outfiles[2]), shell=True)


def getrefseqs(ltyped_alignment, full_alignment):
    """From alignments, finds the longest sequences with genotype info.
    
    Args:
        `ltyped_alignment` (.fasta alignment)
            Alignment of long sequences with genotypes.
            Try to assign reference genotype sequences from these seqs first.
        
        `full_alignment` (.fasta alignment)
            Alignment of all sequences of interest.
            Look through these sequences for genotypes not found in the
            `ltyped_alignment` and assign them a reference sequence if there
            is a sequence of sufficient length (< 60 gaps compared to full).
    
    Returns:
        `gt_seqs` (dict)
            Dictionary of a reference sequence for each genotype already
            called in the data.

    """
    gt_seqs = {}
    # Set reference seqs with long seqs first
    for record in SeqIO.parse(ltyped_alignment, 'fasta'):
        genotype_info = record.description.split(' ')
        if len(genotype_info) == 3:
            genotype = genotype_info[-1]
        elif len(genotype_info) == 4:
            genotype = genotype_info[2] + genotype_info[3]
            if genotype == 'BAIV':
                genotype = 'BA4'
        if genotype not in gt_seqs.keys():
            if 'X' not in str(record.seq): # No ambiguous amino acids
                gt_seqs[genotype] = str(record.seq)

    # Add ref seqs for genotypes with no 'long' seq
    added_gts = {}
    for record in SeqIO.parse(full_alignment, 'fasta'):
        genotype_info = record.description.split(' ')
        if len(genotype_info) == 3:
            if genotype_info[-1] != 'NaN':
                genotype = genotype_info[-1]
        elif len(genotype_info) == 4:
            genotype = genotype_info[2] + genotype_info[3]
            if genotype == 'BAIV':
                genotype = 'BA4'
        if genotype not in gt_seqs.keys(): #no full seq for gt
            if 'X' not in str(record.seq):
                if genotype not in added_gts.keys():
                    added_gts[genotype] = [str(record.seq)]
                else:
                    added_gts[genotype].append(str(record.seq))

    # Pick added reference seq with least number of gaps. Must have < 60.
    for added_gt in added_gts:
        possible_refs = added_gts[added_gt]
        gt_refseq = min(possible_refs, key=lambda seq: seq.count('-'))
        if gt_refseq.count('-') < 60:
            gt_seqs[added_gt] = gt_refseq

    return gt_seqs

def assign_gt(alignment, gt_refdict, hd_threshold):
    """
    Assign a genotype to a non-genotyped sequence as long as the match
    meets a certain hamming distance threshold.
    
    Args:
        `alignment` (.fasta alignment)
            Alignment of all sequences (genotyped and non) to analyze.

        `gt_refdict` (dict)
            Dictionary of genotypes and reference sequences

        `hd_threshold` (int)
            Hamming distance threshold that sets the maximum number of 
            mismatches allowed between a sequence and the most similar 
            genotype reference sequence. If the minimum hamming distance
            between a sequence and its most similar reference sequence is 
            greater than `hd_threshold`, the genotype will remain 'NaN'.

    Returns:
        `updated_gts` (list of tuples)
            List of df index and new genotype for sequences with new gt

        `mistyped` (int)
            Number of sequences that were mistyped using this genotyping 
            method. Being mistyped means the method assigned them a 
            genotype that did not agree with the downloaded subtype. 
            These genotypes are reset to 'NaN'.

    """

    updated_gts = []

    #Keep track of number seqs mistyped based on subtype/genotype disagreement
    mistyped = 0

    for record in SeqIO.parse(alignment, 'fasta'):
        genotype = record.description.split(' ')[2]
        subtype = record.description.split(' ')[1]
        # only add genotypes for seqs with subtypes so can check concordance
        if genotype == 'NaN' and subtype != 'NaN':
            gt_hds = {}
            for gtype in gt_refdict:
                gt_ref = gt_refdict[gtype]
                gt_hds[gtype] = hamming_distance(gt_ref, str(record.seq))
            # At least *hd_threshold* sites must match to call genotype.
            if min(gt_hds.values()) < (len(record.seq) - hd_threshold):
                new_gt = min(gt_hds, key=gt_hds.get)
                if record.description.split(' ')[1] == 'B':
                    if new_gt not in GT_B_LIST:
                        mistyped += 1
                        print("\nSeq: {0}. Genotype {1} and subtype {2} "\
                              "discordant. Reset genotype to 'NaN'.".format(
                              record.name, new_gt, 'B'))
                        new_gt = 'NaN'

                elif record.description.split(' ')[1] == 'A':
                    if new_gt not in GT_A_LIST:
                        mistyped += 1
                        print("\nSeq: {0}. Genotype {1} and subtype {2} "\
                              "discordant. Reset genotype to 'NaN'.".format(
                              record.name, new_gt, 'A'))
                        new_gt = 'NaN'

                updated_gts.append((record.name, new_gt))

    return [updated_gts, mistyped]


def main():
    """Align downloaded sequences, call genotypes, and return final df"""

    parser = rsview.parsearguments.genotype_parser()
    args = vars(parser.parse_args())
    prog = parser.prog

    print("\nExecuting {0} ({1}) in {2} at {3}.\n".format(
            prog, rsview.__version__, os.getcwd(), time.asctime()))

    files = [filename for filename in glob.glob('{0}*.csv'.format(
            args['inprefix']))]

    if not os.path.isdir('{0}'.format(args['seqsdir'])):
        os.makedirs('{0}'.format(args['seqsdir']))

    if not os.path.isdir('{0}'.format(args['outdir'])):
        os.makedirs('{0}'.format(args['outdir']))

    outfile = '{0}/RSVG_all_genotyped.csv'.format(args['outdir'])

    hd_threshold = args['threshold']

    print('Input files: {0}'.format(files))

    file_frames = []
    for file in files:
        if not os.path.isfile(file):
            raise ValueError('Sequence data not downloaded. Run '\
                    '`seq_download.py`.')
        file_df = pd.read_csv(file)
        file_frames.append(file_df)

    rsv_df = merge_csvs(files)

    seqs = rsv_df[['G_seq', 'genotype', 'subtype']]

    already_genotyped = sum(seqs['genotype'].value_counts())

    print('\nStarting with {0} of {1} seqs genotyped.\n'.format(
            already_genotyped, len(seqs)))

    seqs = seqs.fillna(value='NaN') #easily callable placeholder

    #Establish files for seqs. Make 3 files for iterative alignment.
    longtyped_fasta = '{0}/G_seqs_longtyped.fasta'.format(args['seqsdir'])
    long_fasta = '{0}/G_seqs_long_nogt.fasta'.format(args['seqsdir'])
    short_fasta = '{0}/G_seqs_short.fasta'.format(args['seqsdir'])

    seqs_files = [longtyped_fasta, long_fasta, short_fasta]

    seqstofastas(seqs, seqs_files)

    # Establish files for alignments
    aligned_ltyped = '{0}/G_longtyped_aligned.fasta'.format(args['seqsdir'])
    aligned_long = '{0}/G_long_all_aligned.fasta'.format(args['seqsdir'])
    aligned_all = '{0}/G_all_aligned.fasta'.format(args['seqsdir'])

    alignment_files = [aligned_ltyped, aligned_long, aligned_all]

    # Make alignments
    align_seqs(seqs_files, alignment_files)

    # Set reference seqs for each genotype
    gt_refs = getrefseqs(aligned_ltyped, aligned_all)

    # Assign new genotypes and calculate number mistyped
    new_gt_info = assign_gt(aligned_all, gt_refs, hd_threshold)
    new_gts = new_gt_info[0]
    num_mistyped = new_gt_info[1]

    print("\n{0} genotypes mistyped and reset to 'NaN'.".format(num_mistyped))
    print('{0} genotypes added.'.format(len(new_gts) - num_mistyped))
    print("{0} seqs now genotyped.".format(already_genotyped + len(new_gts)
            - num_mistyped))

    # Assign genotypes back to full dataframe
    for new_gt in new_gts:
        rsv_df.loc[int(new_gt[0]), 'genotype'] = new_gt[1]

    # Make a clean df with relevant columns and save as `.csv`
    clean_df = rsv_df[['collection_date', 'country', 'subtype', 'genotype', \
                       'G_seq']]
    clean_df.to_csv(outfile)


if __name__ == '__main__':
    START_TIME = time.time()
    main()
    END_TIME = time.time()
    print('Finished at {0}. Took {1:.3f} minutes to run.'.format(
            time.asctime(), (END_TIME - START_TIME)/60))
