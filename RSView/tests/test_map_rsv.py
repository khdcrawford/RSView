import unittest
import os
import RSView.map_rsv as map_rsv

class TestMapRsv(unittest.TestCase):
    """
    Tests map_rsv.py
    """

    def run_organize_data(self):
        rsv_df = map_rsv.organize_data(map_rsv.DATAFILES, map_rsv.GENOTYPE_DICT)
        return rsv_df

    def run_count_types(self, level, genotype_level='collapse'):
        organized_df = map_rsv.count_types(self.run_organize_data(), map_rsv.JITTER_DICT, level, genotype_level)
        return organized_df

    def run_map_rsv(self, level, genotype_level='collapse', years=[1980,2017]):
        fig = map_rsv.map_rsv(self.run_count_types(level, genotype_level), level,
                              genotype_level, years)
        return fig

    def test_organize_data(self):
        rsv_df = self.run_organize_data()

        self.assertEqual(list(rsv_df.columns), ['collection_date', 'country', 'genotype',
                         'subtype', 'year', 'genotype_group'])
        self.assertTrue(len(rsv_df[col].notnull()) == len(rsv_df) for col in rsv_df.columns)

    def test_count_types(self):
        #Test level='subtype'
        organized_df_subtype = self.run_count_types('subtype')

        self.assertEqual(list(organized_df_subtype.columns), ['country', 'subtype', 'year', 'count',
                         'country_code', 'Longitude', 'Latitude', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df_subtype[col].notnull()) == len(organized_df_subtype) for
                        col in organized_df_subtype.columns)

        #Test level='genotype'
        organized_df_genotype = self.run_count_types('genotype')

        self.assertEqual(list(organized_df_genotype.columns), ['country', 'subtype',
                         'genotype_group', 'year', 'count', 'country_code', 'Longitude',
                         'Latitude', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df_genotype[col].notnull()) == len(organized_df_genotype) for
                        col in organized_df_genotype.columns)

        #Test level='genotype', genotype_level='all'
        organized_df_genotypegroup = self.run_count_types('genotype', genotype_level='all')

        self.assertEqual(list(organized_df_genotypegroup.columns), ['country', 'subtype',
                         'genotype', 'year', 'count', 'country_code', 'Longitude',
                         'Latitude', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df_genotypegroup[col].notnull()) ==
                        len(organized_df_genotypegroup) for col in
                        organized_df_genotypegroup.columns)

    def test_map_rsv(self):
        #Test level='subtype'
        self.run_map_rsv('subtype')

        #Test level='subtype', years = 'all'
        self.run_map_rsv('subtype', years='all')

        #Test level='genotype'
        self.run_map_rsv('genotype')

        #Test level='genotype', genotype_level='all'
        self.run_map_rsv('genotype', genotype_level='all')

        #Test level='genotype', genotype_level='all', years='all'
        self.run_map_rsv('genotype', genotype_level='all', years='all')

if __name__ == '__main__':
    unittest.main()
