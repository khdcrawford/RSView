import time
import pandas as pd
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go

CMAP= {'A':'royalblue','B':'salmon'}
JITTER_DICT= {'A':1.0, 'B':-1.0}
SUBTYPE_LIST=['A','B']
DATAFILES = ['../data/RSVG_gb_metadata.csv','../data/RSVG_gb_metadata_5000-10000.csv']

def organize_data(datafiles):

    rsv_df = pd.DataFrame()

    for datafile in datafiles:

        temp_df = pd.read_csv(datafile, usecols=['collection_date', 'country','subtype'],
                             parse_dates=['collection_date'])
        rsv_df = rsv_df.append(temp_df, ignore_index=True)

    rsv_df['year'] = rsv_df['collection_date'].apply(lambda x: x.year)
    
    #Fix specific country names where city is given
    countries_with_cities = ['Brazil','China','Russia','New Zealand','Spain',
                             'Germany','Egypt', 'India','Japan',
                             'Malaysia','Jordan', 'Saudi Arabia','Myanmar', 'Netherlands']
    for c in countries_with_cities:
        rsv_df['country'] = np.where(rsv_df['country'].str.contains(c),c,rsv_df['country'])

    #Fix specific country names where lat/lon table uses alternate country name
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('USA'),'United States',rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('South Korea'),'Korea',rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('Viet Nam'),'Vietnam',rsv_df['country'])
    rsv_df['country'] = np.where(rsv_df['country'].str.contains('Laos'),'Lao PDR',rsv_df['country'])


    #use lat and long so datapoints can be jittered to show multiple subtypes
    #lat and long data from https://worldmap.harvard.edu/data/geonode:country_centroids_az8

    lat_lon = pd.read_csv('../data/country_centroids_az8.csv', \
                          usecols=['name','brk_a3','Longitude','Latitude']\
                         ).rename(columns={'name':'country', 'brk_a3':'country_code'})

    #count number of rows(seqs) from each country that are each subtype
    df_count_time = pd.DataFrame({'count' : rsv_df.groupby(['country', 'subtype', 'year']).size()}).reset_index()

    #compile country-specific subtype count data with lat and long for plotting
    df_countries_time = df_count_time.merge(lat_lon, how='left', left_on='country', right_on='country')

    return df_countries_time


def jitter_locs(organized_df, jitter_dict):
    #Jitter points for countries that have multiple subtypes, so markers on map don't overlap
    country_group = organized_df.groupby('country').size()

    #With data separated by year
    organized_df['adj_lon'] = np.where(country_group[organized_df['country']]>1,
                                       (organized_df['Longitude']+organized_df.subtype.map(lambda x: jitter_dict[x])
                                       ), organized_df['Longitude'])

    organized_df['adj_lat'] = np.where(country_group[organized_df['country']]>1,
                                       (organized_df['Latitude']+organized_df.subtype.map(lambda x: jitter_dict[x])
                                       ), organized_df['Latitude'])

    return organized_df

def map_rsv(organized_df, cmap, year_range, subtype_list):
    scale_markers = 4

    map_list = []
    for i in range(len(organized_df)):

        map_country = dict(
            type = 'scattergeo',
    #         locationmode = 'country names',
    #         locations = [df_countries.loc[i,'country']],
            lat = [organized_df.loc[i,'adj_lat']],
            lon = [organized_df.loc[i,'adj_lon']],
            marker = dict(
                size = np.min([organized_df.loc[i,'count']*scale_markers, 75]), #Threshold max size of marker
                color = cmap[organized_df.loc[i,'subtype']],
                line = dict(width=0.5, color='rgb(40,40,40)'),
                opacity=0.5,
                sizemode = 'diameter'),
            hovertext = (organized_df.loc[i,'country']+', subtype '+organized_df.loc[i,'subtype']+
                         ' : '+str(organized_df.loc[i,'count'])+' sequences'),
            name = organized_df.loc[i,'country']+' '+organized_df.loc[i,'subtype'],
            legendgroup= organized_df.loc[i,'subtype'],
            showlegend=False,
            hoverinfo = 'text'
        )
        map_list.append(map_country)

    #Work around for showing legend
    for subtype in subtype_list:
        subtype_legend = dict(
                type = 'scattergeo',
                lat = [180.0],
                lon= [180.0],
                marker = dict(
                    size = scale_markers*10,
                    color = cmap[subtype],
                    opacity=0.5,
                    sizemode = 'area'),
                legendgroup = subtype,
                name = 'Subtype ' + subtype,
                showlegend=True,
                hovertext=None,
            )
        map_list.append(subtype_legend)

    steps = []
    for year in year_range:
        step = dict(
        method = 'restyle',
        label = year,
        args = ['visible', [False] * (len(organized_df)+len(subtype_list))])
        for i in range(len(organized_df)):
            if organized_df.loc[i,'year']==year:
                step['args'][1][i] = True # Toggle i'th year to "visible"
        for subtype in subtype_list:
            step['args'][1][(len(organized_df)+subtype_list.index(subtype))] = True
        steps.append(step)

    layout = dict(
            title = 'Global distribution of RSV',
    #         showlegend = True,
            sliders = [dict(
                steps = steps)],
                geo = dict(
                scope='world',
                showland = True,
                landcolor = 'rgb(217, 217, 217)',
                countrywidth=1,
            ),
        )

    fig = dict(data=map_list, layout=layout)
    py.plot(fig)

def main():
    rsv_df = organize_data(DATAFILES)
    YEAR_RANGE = [year for year in range(int(rsv_df.year.min()),int(rsv_df.year.max()))]
    rsv_df = jitter_locs(rsv_df, JITTER_DICT)
    map_rsv(rsv_df, CMAP, YEAR_RANGE, SUBTYPE_LIST)

if __name__ == '__main__':
	start_time = time.time()
	main()
	end_time = time.time()
	print('Program took {0:.3f} minutes to run.'.format((end_time - start_time)/60))
