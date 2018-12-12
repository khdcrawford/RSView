"""Test rsview.seq_download.py

These tests rely on Entrez.efetch functioning properly and 
they do not test that function. 
"""
import unittest
from Bio import Entrez
import rsview.seq_download as seq_download

# constants for testing
EMAIL = 'test@uw.edu'
DB = 'nuccore'
RETMAX = 5
TERM = 'human respiratory syncytial virus G'

#GT_A and GT_B lists based on genotypes already known to be present.
GT_A_LIST = ['ON1', 'GA1', 'GA2', 'GA3', 'GA5', 'GA6', 'GA7', 'NA1', 'NA2',
             'NA3', 'SAA1', 'SAA2']
GT_B_LIST = ['BA', 'BA2', 'BA4', 'BA5', 'BA7', 'BA8', 'BA9', 'BA10', 'BA11',
             'BA12', 'BA14', 'THB', 'SAB1', 'SAB3', 'SAB4', 'GB1', 'GB2',
             'GB3', 'GB4', 'GB13', 'GB12']

TEST_DICTS = [{'organism': 'Human respiratory syncytial virus A'},
              {'organism': 'Human orthopneumovirus', 'strain': 'RSVB'},
              {'organism': 'Human paramyxovirus', 'note': 'subtype: A, '\
               'genotype: NA1'}, {'organism': 'Human orthopneumovirus'}]

TEST_SUBTYPES = ['A', 'B', 'A', '']
TEST_GENOTYPES = ['', '', 'NA1', '']

COLUMNS = ['G_seq', 'genotype', 'subtype', 'collection_date', 'country']

class TestSeqDownload(unittest.TestCase):
    """
    Tests seq_download.py
    """

    def test_getids(self):
        """Make sure Entrez query returns expected number of seq IDs."""
        Entrez.email = EMAIL
        self.assertTrue(len(seq_download.getids(DB, RETMAX, TERM)) == RETMAX)


    def test_find_subtype(self):
        """Make sure get correct subtypes from toy data."""
        for i in range(len(TEST_SUBTYPES)):
            self.assertTrue(seq_download.find_subtype(
                    TEST_DICTS[i]) == TEST_SUBTYPES[i])

    def test_find_genotype(self):
        """Makes sure genotype added to toy dicts.
        Assumes `test_find_subtype` passes.
        Make sure assigns correct genotypes basd on toy data.
        """
        for i in range(len(TEST_DICTS)):
            TEST_DICTS[i]['subtype'] = TEST_SUBTYPES[i] 
            genotyped_dict = seq_download.find_genotype(
                    TEST_DICTS[i], GT_A_LIST, GT_B_LIST)
            self.assertTrue('genotype' in genotyped_dict.keys())
            self.assertTrue(genotyped_dict['genotype'] == TEST_GENOTYPES[i])

    def test_makedf(self):
        """Test dataframe creation from small download.

        This does not test `seq_download.gethandle()` as that function is
        essentially just an implementation of `Entrez.efetch`, which is not
        maintained by this project.
        """
        # This assumes test_getids passes
        ids = seq_download.getids(DB, RETMAX, TERM)
        firstseq = 0
        rettype = 'gb'
        retmode = 'xml'
        handle = seq_download.gethandle(DB, ids, firstseq, RETMAX, 
                rettype, retmode)

        test_df = seq_download.makedf(handle)

        self.assertTrue(len(ids) == len(test_df))

        for column in COLUMNS:
            self.assertTrue(column in list(test_df))





if __name__ == '__main__':
    unittest.main()
