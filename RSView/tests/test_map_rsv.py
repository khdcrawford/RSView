import unittest
import os
import map_rsv

class TestMapRsv(unittest.TestCase):
    """
    Tests map_rsv.py
    """

    def test_df_structure(self):
        rsv_df = map_rsv.organize_data(map_rsv.DATAFILES, map_rsv.GENOTYPE_DICT)

        self.assertEqual(list(rsv_df.columns), ['collection_date', 'country', 'genotype',
                         'subtype', 'year', 'genotype_group'])



    # def test_countries(self):
    #     assertTrue(len(organized_df[organized_df['adj_lon'].isnull()])==0)

if __name__ == '__main__':
    unittest.main()
