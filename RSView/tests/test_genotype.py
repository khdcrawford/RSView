"""Test RSView.genotype.py"""
import unittest
import RSView.genotype as genotype

STR = 'MQGHQNCASSG'
CSVS = ['../data/RSVG_gb_metadata_0-5000.csv', 
        '../data/RSVG_gb_metadata_5000-10000.csv']

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
        self.assertEqual(len(genotype.merge_csvs(CSVS)), 10000)

    



if __name__ == '__main__':
    unittest.main()



