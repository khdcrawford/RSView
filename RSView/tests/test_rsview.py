import unittest

import map_rsv

class TestMapRsv(unittest.TestCase):
    def test_countries(self):
        assertTrue(len(organized_df[organized_df['adj_lon'].isnull()])==0)
