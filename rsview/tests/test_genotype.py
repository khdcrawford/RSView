"""Test rsview.genotype.py"""
import unittest
import os
from Bio import SeqIO
import rsview.genotype as genotype


STR = 'MQGHQNCASSG'
CSVS = ['../data/RSVG_gb_metadata_0-5000.csv',
        '../data/RSVG_gb_metadata_5000-10000.csv']      
COLUMNS = ['G_seq', 'genotype', 'subtype', 'country', 'collection_date']
TESTDIR = './testdata'


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
        """Test creation of fasta files"""
        # Assumes test_merge_csvs passed
        merged_df = genotype.merge_csvs(CSVS)
        test_seqsdf = merged_df[['G_seq', 'genotype', 'subtype']]
        # Ensure some seqs will be considered full length
        test_length = 100
        outfiles = ['{0}/test1.fasta'.format(TESTDIR), 
                    '{0}/test2.fasta'.format(TESTDIR), 
                    '{0}/test3.fasta'.format(TESTDIR)]
        
        outfiles_short = outfiles[:2]

        with self.assertRaises(AssertionError):
            genotype.seqstofastas(test_seqsdf, outfiles_short, test_length)




if __name__ == '__main__':
    unittest.main()
