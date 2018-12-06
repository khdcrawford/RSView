"""Test RSView.seq_download.py"""
import unittest
from Bio import Entrez
import RSView.seq_download as seq_download

# constants for testing
EMAIL = 'test@uw.edu'
DB = 'nuccore'
RETMAX = 5
TERM = 'human respiratory syncytial virus G'

TEST_DICTS = [{'organism': 'Human respiratory syncytial virus A'}, 
              {'organism': 'Human orthopneumovirus', 'strain': 'RSVB'},
              {'organism': 'Human paramyxovirus', 'note': 'subtype: A, '\
               'genotype: NA1'}]

TEST_SUBTYPES = ['A', 'B', 'A']


class TestSeqDownload(unittest.TestCase):
    """
    Tests seq_download.py
    """

    def test_getIDs(self):
        """Make sure Entrez query returns expected number of seq IDs."""
        Entrez.email = EMAIL
        self.assertTrue(len(seq_download.getIDs(DB, RETMAX, TERM)) == RETMAX)


    def test_find_subtype(self):
        """Make sure get correct subtypes from toy data."""
        for i in range(len(TEST_SUBTYPES)):
            self.assertTrue(seq_download.find_subtype(TEST_DICTS[i]) == TEST_SUBTYPES[i])


if __name__ == '__main__':
    unittest.main()
