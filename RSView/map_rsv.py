""" Map the global distribution of Respiratory Syncytial Virus (RSV) by collection date, location,
and viral subtype or genotype """

import os
import glob
import pandas as pd
import numpy as np
import plotly.plotly as py
import matplotlib.pyplot as plt
from matplotlib import colors

import RSView.parsearguments


JITTER_DICT = {'A':1.0, 'B':-1.0}
DATAFILES = [filename for filename in glob.glob('./data/RSVG_gb_metadata*.csv')]
GENOTYPE_DICT = {'GA2':'GA', 'GA5':'GA', 'GB12':'GB', 'GB13':'GB', 'GA3':'GA', 'NA1':'NA',
                 'NA2':'NA', 'ON1':'ON', 'SAA1':'SAA', 'BA':'BA', 'BA10':'BA', 'BA9':'BA',
                 'GB3':'GB', 'SAB1':'SAB', 'SAB3':'SAB', 'SAB4':'SAB', 'NA3':'NA', 'GB2':'GB',
                 'BA11':'BA', 'BA7':'BA', 'GA7':'GA', 'BA 10':'BA', 'BA 12':'BA', 'BA 14':'BA',
                 'BA 2':'BA', 'BA 8':'BA', 'BA 9':'BA', 'SAA2':'SAA', 'BA8':'BA', 'GA1':'GA',
                 'BA4':'BA', 'BA5':'BA', 'BA12':'BA', 'GB1':'GB', 'BA08':'BA', 'BA09':'BA',
                 'BA IV':'BA', 'THB':'TH'}

def organize_data(datafiles, genotype_dict):
    """
    Load .csv files containing RSV sequence data and extract relevant columns. Ensure country
    names to conform to standard names used to retrieve latitude and longitude data. Return a
    DataFrame where each sequence is row containing country, genotype, subtype, and collection date
    information.
    """

    #Return error if data files are not present
    for filename in datafiles:
        if os.path.isfile(filename):
            pass
        else:
            raise Error('Sequence data has not been downloaded yet. Run seq_download.py')

    #Append relevant columns from all data files to DataFrame
    rsv_df = pd.DataFrame()

    for datafile in datafiles:
        temp_df = pd.read_csv(datafile, usecols=['collection_date', 'country', 'subtype',
                                                 'genotype'], parse_dates=['collection_date'])
        rsv_df = rsv_df.append(temp_df, ignore_index=True)

    rsv_df['year'] = rsv_df['collection_date'].apply(lambda x: x.year)

    #Add column to group genotypes by clade
    rsv_df = rsv_df.assign(genotype_group=rsv_df['genotype'].map(genotype_dict))

    #Fix specific country names where city is given
    countries_with_cities = ['Brazil', 'China', 'Russia', 'New Zealand', 'Spain', 'Kenya',
                             'Germany', 'Egypt', 'India', 'Japan', 'Canada', 'Italy',
                             'Malaysia', 'Jordan', 'Saudi Arabia', 'Myanmar', 'Netherlands',
                             'France', 'Peru']
    for con in countries_with_cities:
        rsv_df['country'] = np.where(rsv_df['country'].str.contains(con), con, rsv_df['country'])

    #Fix specific country names where lat/lon table uses alternate country name
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('USA'), 'United States',
                                 rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('South Korea'), 'Korea',
                                 rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('Viet Nam'), 'Vietnam',
                                 rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('Laos'), 'Lao PDR',
                                 rsv_df['country'])

    #Only keep sequences with an assigned subtype
    rsv_df = rsv_df[rsv_df['subtype'].notnull()]

    return rsv_df


def count_types(rsv_df, jitter_dict, level, genotype_level='collapse'):
    """
    Restructure the DataFrame so that each row indicates the total number of RSV sequences
    found in each country, each year, for each subtype or genotype (specified by the level argument)
    """

    #use lat and long so datapoints can be jittered to show multiple subtypes
    #lat and long data from https://worldmap.harvard.edu/data/geonode:country_centroids_az8

    lat_lon = pd.read_csv('./data/country_centroids.csv',
                          usecols=['name', 'brk_a3', 'Longitude', 'Latitude']
                          ).rename(columns={'name':'country', 'brk_a3': 'country_code'})

    #Level specified by required argument
    if level == 'subtype':
        #count number of rows(seqs) from each country that are each subtype
        df_group = pd.DataFrame({'count' : rsv_df.groupby(['country', 'subtype',
                                                           'year']).size()}).reset_index()

        #compile country-specific subtype count data with lat and long for plotting
        organized_df = df_group.merge(lat_lon, how='left', left_on='country', right_on='country')

    elif level == 'genotype':
        #count number of rows(seqs) from each country that are each subtype
        genotype_subset = rsv_df[rsv_df['genotype'].notnull()]

        #genotype_level can be specified by an optional argument
        if genotype_level == 'collapse':
            df_group = pd.DataFrame(
                {'count' : genotype_subset.groupby(['country', 'subtype', 'genotype_group',
                                                    'year']).size()}).reset_index()
        else:
            df_group = pd.DataFrame(
                {'count' : genotype_subset.groupby(['country', 'subtype', 'genotype',
                                                    'year']).size()}).reset_index()

        #compile country-specific subtype count data with lat and long for plotting
        organized_df = df_group.merge(lat_lon, how='left', left_on='country', right_on='country')


    #Jitter points for countries that have multiple subtypes, so markers on map don't overlap
    country_group = organized_df.groupby('country').size()


    #With data separated by year
    organized_df['adj_lon'] = np.where(country_group[organized_df['country']] > 1,
                                       (organized_df['Longitude']+organized_df.subtype.map(
                                           lambda x: jitter_dict[x])), organized_df['Longitude'])

    organized_df['adj_lat'] = np.where(country_group[organized_df['country']] > 1,
                                       (organized_df['Latitude']+organized_df.subtype.map(
                                           lambda x: jitter_dict[x])), organized_df['Latitude'])

    #Find any country names that don't match between sequence DF and lat/lon database
    if len(organized_df[organized_df['adj_lon'].isnull()]) != 0:
        print('Warning: the following country names do not match between sequence DataFrames and\
               "country_centroids.csv"' +
              str(organized_df[organized_df['adj_lon'].isnull()]['country']))

    return organized_df


def map_rsv(organized_df, level, genotype_level='collapse', years=[1990,2017]):
    """
    Use ploy.ly to map RSV sequences onto a global map, with bubbles indicating the virus
    collection location. Bubbles are colored according to subtype or genotype (indicated by the
    level argument) and their size is proportional to the number of sequences collected from the
    given location. Maps are separated temporally by collection date and the slider at the bottom of
    the plot allows the user to scroll through time to see how RSV distribution changes over time.
    """

    #years can specified by an optional argument
    if years == 'all':
        year_range = [yr for yr in range(int(organized_df.year.min()),
                      int(organized_df.year.max()))]
    else:
        year_range = [yr for yr in range(years[0],years[1]+1)]

    if level == 'subtype':
        type_list = ['A', 'B']
        cmap = {'A':'royalblue', 'B':'salmon'}

    elif level == 'genotype':
        #genotype_level can be specified by an optional argument
        if genotype_level == 'collapse':
            a_genotypes = list(set(organized_df[organized_df['subtype'] == 'A']
                                   ['genotype_group'].tolist()))
            b_genotypes = list(set(organized_df[organized_df['subtype'] == 'B']
                                   ['genotype_group'].tolist()))
        else:
            a_genotypes = list(set(organized_df[organized_df['subtype'] == 'A']
                                   ['genotype'].tolist()))
            b_genotypes = list(set(organized_df[organized_df['subtype'] == 'B']
                                   ['genotype'].tolist()))
        type_list = a_genotypes + b_genotypes

        cmap = {}
        blues = plt.get_cmap('GnBu')
        reds = plt.get_cmap('OrRd')

        for a_genotype in a_genotypes:
            cmap[a_genotype] = colors.to_hex(
                blues((a_genotypes.index(a_genotype)+1.0)/len(a_genotypes)))
        for b_genotype in b_genotypes:
            cmap[b_genotype] = colors.to_hex(
                reds((b_genotypes.index(b_genotype)+1.0)/len(b_genotypes)))

    scale_markers = 1
    map_list = []

    #Reassign level for collapsed genotypes so 'genotype_group' will be referenced during plotting
    if level == 'genotype':
        if genotype_level == 'collapse':
            level = 'genotype_group'

    #Make dictionaries for each point to be plotted on a plotly map
    for i in range(len(organized_df)):

        map_country = dict(
            type='scattergeo',
            lat=[organized_df.loc[i, 'adj_lat']],
            lon=[organized_df.loc[i, 'adj_lon']],
            marker=dict(
                size=[np.min([organized_df.loc[i, 'count']*scale_markers, 75])],
                sizemin=5,
                color=cmap[organized_df.loc[i, level]],
                line=dict(width=0.5, color='rgb(40,40,40)'),
                opacity=0.75,
                sizemode='diameter'),
            hovertext=(organized_df.loc[i, 'country'] + ', ' + str(level) + ' ' +
                       organized_df.loc[i, level] + ' : ' + str(organized_df.loc[i, 'count'])+
                       'sequences'),
            name=organized_df.loc[i, 'country']+' '+organized_df.loc[i, level],
            legendgroup=organized_df.loc[i, level],
            showlegend=False,
            hoverinfo='text'
        )
        map_list.append(map_country)

    #Work around for showing legend
    for subtype in type_list:
        subtype_legend = dict(
            type='scattergeo',
            lat=[180.0],
            lon=[180.0],
            marker=dict(
                size=scale_markers*10,
                color=cmap[subtype],
                opacity=0.5,
                sizemode='area'),
            legendgroup=subtype,
            name=str(level)+ ' ' + subtype,
            showlegend=True,
            hovertext=None,
            )
        map_list.append(subtype_legend)

    steps = []
    for year in year_range:
        step = dict(
            method='restyle',
            label=year,
            args=['visible', [False] * (len(organized_df)+len(type_list))])
        for i in range(len(organized_df)):
            if organized_df.loc[i, 'year'] == year:
                step['args'][1][i] = True # Toggle i'th year to "visible"
        for subtype in type_list:
            step['args'][1][(len(organized_df)+type_list.index(subtype))] = True
        steps.append(step)

    layout = dict(
        title='Global distribution of RSV',
        sliders=[dict(
            steps=steps, y=0.185)],
        geo=dict(
            scope='world',
            showland=True,
            landcolor='rgb(217, 217, 217)',
            countrywidth=1,
            ),
        legend=dict(x=1.02, y=0.5))

    fig = dict(data=map_list, layout=layout)
    py.plot(fig)


def main(level, genotype_level, years):
    """
    Run organize_data, count_types, map_rsv
    """
    rsv_df = organize_data(DATAFILES, GENOTYPE_DICT)
    organized_df = count_types(rsv_df, JITTER_DICT, level, genotype_level=genotype_level)
    map_rsv(organized_df, level, genotype_level=genotype_level, years=years)

if __name__ == "__main__":

    parser = RSView.parsearguments.mapParser()
    args = parser.parse_args()

    main(args.level, genotype_level=args.genotype_level, years=args.years)
