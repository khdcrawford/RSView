import map_rsv
import pandas as pd
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
import country_converter as coco
import plot_rsv

def country_to_iso3(country):
    if country == 'Global':
        return 'Global'
    else:
        iso3 =  coco.convert(names=country, to='ISO3')
        return iso3

def count_subtypes(rsv_df):
    """
    Restructure the DataFrame so that each row indicates the total number of RSV sequences
    found in each country, each year, for each subtype
    """
    
    #count number of rows(seqs) from each country that are each subtype
    df_group = pd.DataFrame({'count' : rsv_df.groupby(['country', 'subtype']).size()}).reset_index()

    return df_group


def count_subtypes_year(rsv_df):
    """
    Restructure the DataFrame so that each row indicates the total number of RSV sequences
    found in each country, each year, for each subtype
    """
    
    #count number of rows(seqs) from each country that are each subtype
    df_group = pd.DataFrame({'count' : rsv_df.groupby(['country', 'subtype',
                                                           'year']).size()}).reset_index()

    return df_group

def get_ratio_all(df):
    df_ratio = pd.DataFrame(columns=['country', 'ratio', 'count_A', 'count_B', 'count_total'])
    countries = df['country'].unique()
    for i in countries:
        country_df = df[df['country']==i]
        if country_df['subtype'].count() == 1:
            lone_subtype = country_df['subtype'].item()
            if lone_subtype == 'A':
                percent_A = 1
                count_A = country_df['count'].item()
                count_B = 0
                count_total = count_A + count_B

            else:
                percent_A = 0
                count_A = 0
                count_B = country_df['count'].item()
                count_total = count_A + count_B
        else:
            count_A = float(country_df[country_df['subtype'] == 'A']['count'])
            count_B = float(country_df[country_df['subtype'] == 'B']['count'])
            percent_A = count_A/(count_A + count_B)
            count_total = count_A + count_B
        df_ratio = df_ratio.append(pd.Series([i, percent_A, count_A, count_B, count_total], index=df_ratio.columns, name=i))
    return df_ratio


def merge_ratio_health(df_ratio)
	
	df_health_summary = pd.read_csv('./data/health_data_summary.csv')
	df_ratio["iso3"] = [country_to_iso3(x) for x in df_ratio.country.values]
	
	df_merge = pd.merge(df_ratio, df_health_summary, how='inner', on=['iso3'])

	return df_merge



def plot_ratio(df_merge, data_type):

	trace = go.Scatter(
	    x = df_merge['ratio'],
	    y = df_merge[data_type],
	    mode = 'markers',
	    text = df_merge['country_short'] # The hover text goes here... 
	)

	layout= go.Layout(
	    title= 'RSV Subtype Prevalence compared to ' + plot_rsv.DATA_DICT[data_type],
	    hovermode= 'closest',
	    xaxis= dict(
	        title= 'Ratio of Subtype A over Subtype B Sequences Recorded',
	        # ticklen= 5,
	        # zeroline= False,
	        # gridwidth= 2,
	    ),
	    yaxis=dict(
	        title= plot_rsv.DATA_DICT[data_type],
	        # ticklen= 5,
	        # gridwidth= 2,
	    ),
	    showlegend= False
	)

	fig= go.Figure(data=data, layout=layout)
	plot(fig)





def main(data_type):
    """
    Run organize_data, count_types, map_rsv
    """
	rsv_df = map_rsv.organize_data(map_rsv.DATAFILES, map_rsv.GENOTYPE_DICT)
	rsv_df = rsv_df.dropna(subset = ['subtype']).copy()

	rsv_df_count = count_subtypes(rsv_df)

	df_ratio_all = get_ratio_all(rsv_df_count)












