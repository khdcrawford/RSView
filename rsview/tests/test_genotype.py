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
        """Test getting reference sequences from toy data"""


    def test_assigngt(self):
        """Test assigning gneotypes with toy data"""








if __name__ == '__main__':
    unittest.main()
