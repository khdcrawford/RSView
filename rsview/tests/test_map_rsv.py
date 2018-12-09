import unittest
import os
from unittest.mock import patch

import rsview.map_rsv as map_rsv

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
        organized_df = self.run_count_types('subtype')

        self.assertEqual(list(organized_df.columns), ['country', 'subtype', 'year', 'count',
                         'country_code', 'Longitude', 'Latitude', 'under_five_deaths', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df[col].notnull()) == len(organized_df) for
                        col in organized_df.columns)

        #Test level='genotype'
        organized_df = self.run_count_types('genotype')

        self.assertEqual(list(organized_df.columns), ['country', 'subtype',
                         'genotype_group', 'year', 'count', 'country_code', 'Longitude',
                         'Latitude', 'under_five_deaths', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df[col].notnull()) == len(organized_df) for
                        col in organized_df.columns)

        #Test level='genotype', genotype_level='all'
        organized_df = self.run_count_types('genotype', genotype_level='all')

        self.assertEqual(list(organized_df.columns), ['country', 'subtype',
                         'genotype', 'year', 'count', 'country_code', 'Longitude',
                         'Latitude', 'under_five_deaths', 'adj_lon', 'adj_lat'])
        self.assertTrue(len(organized_df[col].notnull()) ==
                        len(organized_df) for col in organized_df.columns)

    #Don't actually produce plot, just test function components
    @patch("rsview.map_rsv.py.plot")
    def test_map_rsv(self, mock_show):
        #Test level='subtype'
        fig = self.run_map_rsv('subtype')

        self.assertEqual(len(fig['data']), len(self.run_count_types('subtype'))+2)
        self.assertEqual(len(fig['layout']['sliders'][0]['steps']), int(2018-1980))
        self.assertTrue('subtype' in fig['data'][0]['hovertext'])

        #Test level='subtype', years = 'all'
        fig = self.run_map_rsv('subtype', years='all')
        year_range = [yr for yr in range(int(self.run_count_types('subtype').year.min()),
                      int(self.run_count_types('subtype').year.max()))]

        self.assertEqual(len(fig['layout']['sliders'][0]['steps']), len(year_range))

        #Test level='genotype'
        fig = self.run_map_rsv('genotype')
        organized_df = self.run_count_types('genotype')
        a_groups = list(set(organized_df[organized_df['subtype'] == 'A']
                               ['genotype_group'].tolist()))
        b_groups = list(set(organized_df[organized_df['subtype'] == 'B']
                               ['genotype_group'].tolist()))

        self.assertEqual(len(fig['data']), len(organized_df) +
                         len(a_groups+b_groups))
        self.assertTrue('genotype_group' in fig['data'][0]['hovertext'])

        #Test level='genotype', genotype_level='all'
        fig = self.run_map_rsv('genotype', genotype_level='all')
        organized_df = self.run_count_types('genotype', genotype_level='all')
        a_genotypes = list(set(organized_df[organized_df['subtype'] == 'A']
                               ['genotype'].tolist()))
        b_genotypes = list(set(organized_df[organized_df['subtype'] == 'B']
                               ['genotype'].tolist()))

        self.assertEqual(len(fig['data']), len(organized_df) +
                         len(a_genotypes+b_genotypes))
        self.assertTrue('genotype' in fig['data'][0]['hovertext'])


if __name__ == '__main__':
    unittest.main()
