""" Plot the correlation between RSV subtype prevalence and health metrics """

import pandas as pd
import country_converter as coco
from plotly.offline import plot
import plotly.graph_objs as go

import map_rsv
import plot_rsv

from parsearguments import correlation_parser

# import rsview.parsearguments

# from parsearguments import correlationParser

def country_to_iso3(country):
    """ Take user country input and convert to ISO3 to access dataframe """
    if country == 'Global':
        return 'Global'
    iso3 = coco.convert(names=country, to='ISO3')
    return iso3

def count_subtypes(rsv_df):
    """
    Restructure the DataFrame so that each row indicates the total number of RSV sequences
    found in each country, for each subtype (added across all years)
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
    df_group = pd.DataFrame({'count' : rsv_df.groupby(
        ['country', 'subtype', 'year']).size()}).reset_index()

    return df_group

def get_ratio_all(dataframe):
    """
    From counted subtypes, calculate the ratio of subtype A over subtype B in each country
    """
    df_ratio = pd.DataFrame(columns=['country', 'ratio', 'count_A', 'count_B', 'count_total'])
    countries = dataframe['country'].unique()
    for i in countries:
        country_df = dataframe[dataframe['country'] == i]
        if country_df['subtype'].count() == 1:
            lone_subtype = country_df['subtype'].item()
            if lone_subtype == 'A':
                percent_a = 1
                count_a = country_df['count'].item()
                count_b = 0
                count_total = count_a + count_b

            else:
                percent_a = 0
                count_a = 0
                count_b = country_df['count'].item()
                count_total = count_a + count_b
        else:
            count_a = float(country_df[country_df['subtype'] == 'A']['count'])
            count_b = float(country_df[country_df['subtype'] == 'B']['count'])
            percent_a = count_a/(count_a + count_b)
            count_total = count_a + count_b
        df_ratio = df_ratio.append(
            pd.Series([i, percent_a, count_a, count_b, count_total],
                      index=df_ratio.columns, name=i))
    return df_ratio

def get_ratio_year(dataframe):
    """
    From counted subtypes, calculate the ratio of subtype A
    over subtype B in each country, each year
    """
    df_ratio = pd.DataFrame(
        columns=['country', 'year', 'ratio', 'count_A', 'count_B', 'count_total'])
    countries = dataframe['country'].unique()
    for i in countries:
        country_df = dataframe[dataframe['country'] == i]
        years = country_df['year'].unique()
        for j in years:
            years_df = country_df[country_df['year'] == j]
            if years_df['subtype'].count() == 1:
                lone_subtype = years_df['subtype'].item()
                if lone_subtype == 'A':
                    percent_a = 1
                    count_a = years_df['count'].item()
                    count_b = 0
                    count_total = count_a + count_b
                else:
                    percent_a = 0
                    count_a = 0
                    count_b = years_df['count'].item()
                    count_total = count_a + count_b
            else:
                count_a = float(years_df[years_df['subtype'] == 'A']['count'])
                count_b = float(years_df[years_df['subtype'] == 'B']['count'])
                percent_a = count_a/(count_a + count_b)
                count_total = count_a + count_b
            df_ratio = df_ratio.append(
                pd.Series([i, j, percent_a, count_a, count_b, count_total],
                          index=df_ratio.columns, name=i))
    return df_ratio


def merge_ratio_health(df_ratio, datadir):
    """
    Merge the RSV sequence dataframe with calculated ratios with the health dataframe
    """
    df_health_summary = plot_rsv.make_df_health_summary(datadir)
    df_ratio["iso3"] = [country_to_iso3(x) for x in df_ratio.country.values]

    df_merge = pd.merge(df_ratio, df_health_summary, how='inner', on=['iso3'])

    return df_merge

def merge_ratio_health_year(df_ratio, datadir):
    """
    Merge the RSV sequence dataframe with yearly calculated ratios with the health dataframe
    """
    df_health_all = plot_rsv.make_df_health_all(datadir)
    df_ratio["iso3"] = [country_to_iso3(x) for x in df_ratio.country.values]

    df_merge = pd.merge(df_ratio, df_health_all, how='inner', on=['iso3', 'year'])

    return df_merge



def plot_ratio(df_merge, data_type):
    """
    Plot the calculated subtype ratios
    """
    trace1 = go.Scatter(
        x=df_merge['ratio'],
        y=df_merge[data_type],
        mode='markers',
        marker=dict(size=14,
                    line=dict(width=1),
                    color='rgba(204,204,204,1)'
                   ),
        text=df_merge['country_short'] # The hover text goes here
    )

    layout = go.Layout(
        title='RSV Subtype Prevalence compared to ' + plot_rsv.DATA_DICT[data_type],
        hovermode='closest',
        xaxis=dict(
            title='Ratio of Subtype A over Subtype B Sequences Recorded',
            # ticklen= 5,
            # zeroline= False,
            # gridwidth= 2,
        ),
        yaxis=dict(
            title=plot_rsv.DATA_DICT[data_type],
            # ticklen= 5,
            # gridwidth= 2,
        ),
        showlegend=False
    )

    data = [trace1]

    fig = go.Figure(data=data, layout=layout)
    plot(fig)



def plot_ratio_year(df_merge, data_type):
    """
    Plot the calculated subtype ratios for each year
    """

    df_merge['hover'] = df_merge['country_short'] + ', ' + df_merge['year'].astype(int).map(str)

    trace1 = go.Scatter(
        x=df_merge['ratio'],
        y=df_merge[data_type],
        mode='markers',
        marker=dict(size=14,
                    line=dict(width=1),
                    color=df_merge['year'],
                    colorscale='RdBu',
                    showscale=True,
                    # cmin=1990,
                    # cmax=2018
                   ),
        text=df_merge['hover'] # The hover text goes here
    )

    layout = go.Layout(
        title='RSV Subtype Prevalence compared to ' + plot_rsv.DATA_DICT[data_type],
        hovermode='closest',
        xaxis=dict(
            title='Ratio of Subtype A over Subtype B Sequences Recorded',
            # ticklen= 5,
            # zeroline= False,
            # gridwidth= 2,
        ),
        yaxis=dict(
            title=plot_rsv.DATA_DICT[data_type],
            # ticklen= 5,
            # gridwidth= 2,
        ),
        showlegend=False
    )

    data = [trace1]

    fig = go.Figure(data=data, layout=layout)
    plot(fig)



def main(level, data_type, datadir):
    """
    Organize genotype data and plot correlation between subtypes and health metrics
    """
    rsv_df = map_rsv.organize_data(datadir, map_rsv.GENOTYPE_DICT)
    rsv_df = rsv_df.dropna(subset=['subtype']).copy()

    rsv_df_year = rsv_df.dropna(subset=['year']).copy()
    rsv_df_year['year'] = rsv_df_year['year'].astype(int)

    if level == 'all':
        rsv_df_count = count_subtypes(rsv_df)
        df_ratio_all = get_ratio_all(rsv_df_count)
        df_merged = merge_ratio_health(df_ratio_all, datadir)
        plot_ratio(df_merged, data_type)

    else:
        rsv_df_count = count_subtypes_year(rsv_df)
        df_ratio_year = get_ratio_year(rsv_df_count)
        df_merged = merge_ratio_health_year(df_ratio_year, datadir)
        plot_ratio_year(df_merged, data_type)



if __name__ == "__main__":

    ARGS = correlation_parser(plot_rsv.dict_to_help(plot_rsv.DATA_DICT)).parse_args()

    main(ARGS.level, ARGS.data_type, ARGS.datadir)
