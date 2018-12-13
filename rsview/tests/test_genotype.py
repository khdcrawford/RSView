"""Test rsview.genotype.py

Must be run from within `tests` due to file paths.
"""
import unittest
import os
import rsview.genotype as genotype


STR = 'MQGHQNCASSG'
CSVS = ['../data/RSVG_gb_metadata_0-5000.csv',
        '../data/RSVG_gb_metadata_5000-10000.csv']
COLUMNS = ['G_seq', 'genotype', 'subtype', 'country', 'collection_date']

TESTDIR = './testdata_genotype'

SEQSDIR = '../data/seqs'

TESTSEQS = [TESTDIR+'/seqs1.fasta', TESTDIR+'/seqs2.fasta',
             TESTDIR+'/seqs3.fasta']

ALIGNMENTS = [SEQSDIR+'/G_longtyped_aligned.fasta',
              SEQSDIR+'/G_long_all_aligned.fasta',
              SEQSDIR+'/G_all_aligned.fasta']

GT_A_LIST = ['ON1', 'GA1', 'GA2', 'GA3', 'GA5', 'GA6', 'GA7', 'NA1', 'NA2',
             'NA3', 'SAA1', 'SAA2']
GT_B_LIST = ['BA', 'BA2', 'BA4', 'BA5', 'BA7', 'BA8', 'BA9', 'BA10', 'BA11',
             'BA12', 'BA14', 'THB', 'SAB1', 'SAB3', 'SAB4', 'GB1', 'GB2',
             'GB3', 'GB4', 'GB13', 'GB12']

class TestGenotyping(unittest.TestCase):
    """
    Tests genotype.py
    """

    def test_hamming_distance(self):
        """Test get correct output and error from genotype.hamming_distance"""
        str1 = STR
        str2 = str1.replace('Q', 'T')
        self.assertEqual(genotype.hamming_distance(str1, str2),
                         str1.count('Q'))

        with self.assertRaises(AssertionError):
            genotype.hamming_distance(str1, str1[2:])


    def test_merge_csvs(self):
        """Test genotpe.merge_csvs gives right output df length"""
        merged_df = genotype.merge_csvs(CSVS)
        self.assertEqual(len(merged_df), 10000)

        for column in COLUMNS:
            self.assertTrue(column in list(merged_df))


    def test_seqstofastas(self):
        """Test creation of fasta files.

        Test that function can be called.

        Test that files are made.

        Test those files have expected `.fasta` headings.
        """
        # Assumes test_merge_csvs passed

        test_df = genotype.merge_csvs(CSVS)
        # Select columns and rid of 'nan', so all G_seqs are strs
        testseqsdf = test_df[['G_seq', 'genotype', 'subtype']].fillna('NaN')

        # Ensure some seqs will be "full length", so get seqs in all files
        testlength = testseqsdf.G_seq.map(len).max() - 20

        with self.assertRaises(AssertionError):
            genotype.seqstofastas(testseqsdf, TESTSEQS[:2], testlength)

        # Start with empty TESTDIR then test making expected files
        if os.path.isdir(TESTDIR):
            for seqfile in TESTSEQS:
                if os.path.isfile(seqfile):
                    os.remove(seqfile)
        else:
            os.makedirs(TESTDIR)

        for seqfasta in TESTSEQS:
            self.assertFalse(os.path.isfile(seqfasta))

        genotype.seqstofastas(testseqsdf, TESTSEQS, testlength)

        for seqfasta in TESTSEQS:
            self.assertTrue(os.path.isfile(seqfasta))
            with open(seqfasta) as fasta:
                for line in fasta:
                    if '>' in line:
                        self.assertTrue(line[0] == '>')
                        self.assertTrue(3 <= len(line.split(' ')) <= 4)
                        self.assertTrue(line.split(' ')[1] in ['A', 'B', 'NaN'])


    def test_alignseqs(self):
        """Test that alignseqs is called and outputs error if improper input.

        Not testing alignment by mafft as it is slow and mafft is supported
        elsewhere.

        Therefore, I will not test toy data, but will test processed data
        in `../data/seqs/` directory for appropriate headings.
        """

        with self.assertRaises(AssertionError):
            genotype.alignseqs(TESTSEQS, ALIGNMENTS[:2])
            genotype.alignseqs(TESTSEQS[0], ALIGNMENTS[0])

        align_seq = ''
        lenprevseq = 0

        # Only run tests on output if such files exist.
        for alignment in ALIGNMENTS:
            if os.path.isfile(alignment):
                with open(alignment) as alignfile:
                    for line in alignfile:
                        if '>' in line:
                            # Test all aligned seqs have appropriate headers
                            self.assertTrue(line[0] == '>')
                            self.assertTrue(3 <= len(line.split(' ')) <= 4)
                            self.assertTrue(line.split(' ')[1] in ['A', 'B',
                                    'NaN'])
                            # Test all aligned seqs have same length
                            if lenprevseq != 0:
                                self.assertEqual(len(align_seq), lenprevseq)
                            if align_seq:
                                lenprevseq = len(align_seq)
                                align_seq = ''
                        else:
                            align_seq += line


    def test_getrefseqs(self):
        """Test getting reference sequences yields expected output.

        Test on actual alignments from data rather than toy data.
        """

        # Only run tests on output if such files exist.
        for alignment in ALIGNMENTS:
            self.assertTrue(os.path.isfile(alignment))

        gt_list = GT_A_LIST + GT_B_LIST

        test_refs = genotype.getrefseqs(ALIGNMENTS[0], ALIGNMENTS[2])

        for key in test_refs:
            self.assertTrue(key in gt_list)
            # Empty strings are false. Ensure genotypes have seq.
            self.assertTrue(test_refs[key])



    def test_assign_gt(self):
        """Test assigning gneotypes with actual data.

        I have not created toy `.fasta` files for testing, so test on real
        data.
        There is a rough check for mistyping in the main code, so just make
        sure output is as expected. Overall, we only check if our genotypes
        are logical and make no claim that they are fully accurate.

        Check that if threshold is changed, number of new genotypes and 
        mistyped genotypes changes as expected.

        This test is slow and takes ~25 seconds to run.
        """

        self.assertTrue(os.path.isfile(ALIGNMENTS[2]))

        alignall = ALIGNMENTS[2]

        # Assumes test_getrefseqs passes
        test_refgts = genotype.getrefseqs(ALIGNMENTS[0], ALIGNMENTS[2])

        threshold_norm = 150
        threshold_lax = 50
        threshold_strict = 200

        assigngt_norm = genotype.assign_gt(alignall, test_refgts,
                threshold_norm)
        self.assertTrue(len(assigngt_norm[2]) == assigngt_norm[1])

        assigngt_lax = genotype.assign_gt(alignall, test_refgts,
                threshold_lax)
        self.assertTrue(len(assigngt_lax[2]) == assigngt_lax[1])

        assigngt_strict = genotype.assign_gt(alignall, test_refgts,
                threshold_strict)
        self.assertTrue(len(assigngt_strict[2]) == assigngt_strict[1])

        self.assertTrue(assigngt_strict[1] < assigngt_norm[1] \
                < assigngt_lax[1])

        self.assertTrue(len(assigngt_strict[0]) < len(assigngt_norm[0]) \
                < len(assigngt_lax[0]))


if __name__ == '__main__':
    unittest.main()
