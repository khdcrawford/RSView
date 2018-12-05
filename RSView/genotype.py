"""Assign RSV G sequences genotypes.

This program relies on some seqs already being genotyped. It will
not add sequences to a genotype that does not already have a high
quality (< 60 gaps in the alignment) sequence in the input data.
"""
import pandas as pd
import numpy as np
import time
import os
import glob
import subprocess
import RSView.parsearguments

from Bio import SeqIO

#GT_A and GT_B lists based on genotypes already known to be present.
GT_A_LIST = ['ON1', 'GA1', 'GA2', 'GA3', 'GA5', 'GA6', 'GA7', 'NA1', 'NA2',
             'NA3', 'SAA1', 'SAA2']
GT_B_LIST = ['BA', 'BA2', 'BA4', 'BA5', 'BA7', 'BA8', 'BA9', 'BA10', 'BA11',
             'BA12', 'BA14', 'THB', 'SAB1', 'SAB3', 'SAB4', 'GB1', 'GB2', 
             'GB3', 'GB4','GB13', 'GB12']

def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    assert len(s1) == len(s2), "Undefined for sequences of unequal length"
    return sum(site1 != site2 for site1, site2 in zip(s1, s2))

def merge_csvs(csv_files):
    """Merge given csv *files* into one dataframe"""
    file_frames = []
    for file in csv_files:
        if not os.path.isfile(file):
            raise ValueError('Sequence data not downloaded. Run '\
                    '`seq_download.py`.')
        file_df=pd.read_csv(file)
        file_frames.append(file_df)

    full_df = pd.concat(file_frames, ignore_index=True, sort=False)
    full_df.drop('Unnamed: 0', axis=1, inplace=True)

    return full_df

def main():
    """Align downloaded sequences, call genotypes, and return final df"""
    parser = RSView.parsearguments.genotypeParser()
    args = vars(parser.parse_args())
    prog = parser.prog

    print("\nExecuting {0} ({1}) in {2} at {3}.\n".format(
            prog, RSView.__version__, os.getcwd(), time.asctime()))

    files = [filename for filename in glob.glob('{0}*.csv'.format(
            args['inprefix']))]
    
    if not os.path.isdir('{0}'.format(args['seqsdir'])):
        os.makedirs('{0}'.format(args['seqsdir']))
    
    seqsdir = args['seqsdir']
    outfile = args['outfile']

    print('Input files: {0}'.format(files))

    file_frames = []
    for file in files:
        if not os.path.isfile(file):
            raise ValueError('Sequence data not downloaded. Run '\
                    '`seq_download.py`.')
        file_df=pd.read_csv(file)
        file_frames.append(file_df)

    rsv_df = merge_csvs(files)

    seqs=rsv_df[['G_seq', 'genotype', 'subtype']]
    seqs = seqs.fillna(value='NaN') #easily callable placeholder

    #Establish files for seqs. Make 3 files for iterative alignment.
    longtyped_fasta = '{0}/G_seqs_longtyped.fasta'.format(seqsdir)
    long_fasta = '{0}/G_seqs_long_nogt.fasta'.format(seqsdir)
    short_fasta = '{0}/G_seqs_short.fasta'.format(seqsdir)

    #Keep track of sequences with genotypes downloaded from GenBank
    already_genotyped = 0

    with open(longtyped_fasta, 'w') as longtyped, \
         open(long_fasta, 'w') as long_nogt, \
         open(short_fasta, 'w') as short:
        for i in range(len(seqs)):
            st = seqs.at[i, 'subtype']
            gt = seqs.at[i, 'genotype']
            header = '>{0} {1} {2}\n'.format(i, st, gt)
            seq = seqs.at[i, 'G_seq']
            if seq != 'NaN':
                if gt != 'NaN':
                    already_genotyped += 1
                    if len(seq) > 290:
                        longtyped.write('{0}{1}\n'.format(header, seq))
                    else:
                        short.write('{0}{1}\n'.format(header, seq))
                else:
                    if len(seq) > 290:
                        long_nogt.write('{0}{1}\n'.format(header, seq))
                    else:
                        short.write('{0}{1}\n'.format(header, seq))
    
    # Establish files for alignments
    aligned_ltyped = '{0}/G_longtyped_aligned.fasta'.format(seqsdir)
    aligned_long = '{0}/G_long_all_aligned.fasta'.format(seqsdir)
    aligned_all = '{0}/G_all_aligned.fasta'.format(seqsdir)

    # Make alignments 
    # subprocess.check_call('mafft --auto {0} > {1}'.format(longtyped_fasta, 
    #         aligned_ltyped), shell=True)

    # subprocess.check_call('mafft --add {0} --reorder --keeplength {1} > {2}'\
    #         .format(long_fasta, aligned_ltyped, aligned_long), shell=True)

    # subprocess.check_call('mafft --addfragments {0} --reorder --6merpair '\
    #         '--thread -1 --keeplength {1} > {2}'.format(short_fasta, 
    #         aligned_long, aligned_all), shell=True)

    # Set reference seqs for each genotype
    gt_seqs = {}
    # Set reference seqs with long seqs first
    for record in SeqIO.parse(aligned_ltyped, 'fasta'):
        genotype_info = record.description.split(' ')
        if len(genotype_info) == 3:
            genotype = genotype_info[-1]
        elif len(genotype_info) == 4:
            genotype = genotype_info[2] + genotype_info[3]
            if genotype == 'BAIV':
                genotype = 'BA4'
        if genotype not in gt_seqs.keys():
            if 'X' not in str(record.seq):
                gt_seqs[genotype] = str(record.seq)

    # Add ref seqs for genotypes with no 'long' seq
    added_gts = {}
    for record in SeqIO.parse(aligned_all, 'fasta'):
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
    
    # Pick new reference seq with least number of gaps
    for added_gt in added_gts:
        possible_refs = added_gts[added_gt]
        gt_refseq = min(possible_refs, key=lambda seq: seq.count('-'))
        if gt_refseq.count('-') < 60:
            gt_seqs[added_gt] = gt_refseq

    # Keep track of number of new genotypes assigned and number of sequences
    # mistyped based on subtype/genotype disagreement
    gt_assigned = 0
    mistyped = 0
    for record in SeqIO.parse(aligned_all, 'fasta'):
        genotype = record.description.split(' ')[2]
        subtype = record.description.split(' ')[1]
        # only add genotypes for seqs with subtypes so can check concordance
        if genotype == 'NaN' and subtype != 'NaN':
            gt_hds = {}
            for gt in gt_seqs:
                gt_seq = gt_seqs[gt]
                gt_hds[gt] = hamming_distance(gt_seq, str(record.seq))
            # At most 60 '-' in gt_refseq, so must match >= 80 addn'l sites
            if min(gt_hds.values()) < (len(record.seq) - 140):
                new_gt = min(gt_hds, key=gt_hds.get)
                gt_assigned += 1
                if record.description.split(' ')[1] == 'B':
                    if new_gt not in GT_B_LIST:
                        mistyped+=1
                        print("\nSeq: {0}. Genotype {1} and subtype {2} "\
                              "discordant. Reset genotype to 'NaN'.".format(
                              record.name, new_gt, 'B'))
                    new_gt = 'NaN'

                elif record.description.split(' ')[1] == 'A':
                    if new_gt not in GT_A_LIST:
                        mistyped+=1
                        print("\nSeq: {0}. Genotype {1} and subtype {2} "\
                              "discordant. Reset genotype to 'NaN'.".format(
                              record.name, new_gt, 'A'))
                        new_gt = 'NaN'

    print('\n{0} seqs previously genotyped.'.format(already_genotyped))
    print("{0} genotypes mistyped and reset to 'NaN'.".format(mistyped))
    print('{0} genotypes added.'.format(gt_assigned))
    print("{0} seqs now genotyped.".format(already_genotyped + gt_assigned))


if __name__ == '__main__':
    main()