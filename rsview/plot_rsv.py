""" Can generate plots of RSV health data on a summary or individual country level """

import pandas as pd
import country_converter as coco
import plotly
import plotly.graph_objs as go
from plotly.offline import plot

from parsearguments import plotParser


plotly.tools.set_credentials_file(username='jillianeb', api_key='eOsTljd6vVMiyuy4Msy0')
plotly.tools.set_config_file(world_readable=False,
                             sharing='private')


def dict_to_help(data_types):
    """ Turns DATA_DICT into readable format for help text """
    output = ""
    for item in data_types:
        output = output + " " + item + ": "
        output = output + data_types[item] + "\n"
    return output

DF_HEALTH_ALL = pd.read_csv('./data/health_data_all.csv')
DF_HEALTH_SUMMARY = pd.read_csv('./data/health_data_summary.csv')

COLOR_DICT = ['rgba(204,204,204,1)'] * len(DF_HEALTH_SUMMARY)

DATA_DICT = {
    ' nnd ':'Total Neonatal Deaths',
    ' pnd ':'Total Post-Neonatal Deaths',
    ' neo9 ':'Neonatal deaths due to Acute Respiratory Infection',
    ' post9 ':'Post-neonatal deaths due to Acute Respiratory Infection',
    ' ufive9 ':'Underfive deaths due to Acute Respiratory Infection',
    ' rneo9 ':'Neonatal death rate from Acute Respiratory Infection (per 1000 live births)',
    ' rpost9 ':'Post-neonatal death rate from Acute Respiratory Infection (per 1000 live births)',
    ' rufive9 ':'Underfive death rate from Acute Respiratory Infection (per 1000 live births)',
    'fneo9':'Percent Neonatal deaths due to Acute Respiratory Infection',
    'fpost9':'Percent Post-neonatal deaths due to Acute Respiratory Infection',
    'fufive9':'Percent Underfive deaths due to Acute Respiratory Infection'}

def input_to_country(country_input):
    """ Takes user input for country and converts it to short country name """
    if country_input in ('Global', 'global'):
        return 'Global'
    country = coco.convert(names=country_input, to='name_short')
    return country


def is_country_present(dataframe, country):
    """ Checks if user input country is a valid country in the dataset """
    iso3 = coco.convert(names=country, to='ISO3')
    # process some error if coco can't convert to a country

    in_dataset = iso3 in dataframe['iso3']
    return in_dataset


def plot_summary(data_type, highlight_country=None):
    """ Plots summary health data. If a highlight_country is specified, it will be highlighted """

    color_highlight = COLOR_DICT.copy()

    df_sorted = DF_HEALTH_SUMMARY.sort_values(data_type)
    df_sorted = df_sorted.reset_index(drop=True)

    if highlight_country is not None:
        country_short = input_to_country(highlight_country)
        highlight_index = df_sorted.index[df_sorted['country_short'] == country_short].values[0]
        #highlight_index = 186
        color_highlight[highlight_index] = 'rgba(222,45,38,0.8)'

    trace1 = go.Bar(
        x=df_sorted['country_short'],
        y=df_sorted[data_type],
        name='neo9',
        marker=dict(
            color=color_highlight
            ),
    )

    data = [trace1]
    layout = go.Layout(
        title=DATA_DICT[data_type] + '<br> Averaged from 2000-2016',
        yaxis=dict(
            #title='Percent',
            titlefont=dict(
                size=16,
                color='rgb(107, 107, 107)'
            ),
            tickfont=dict(
                size=14,
                color='rgb(107, 107, 107)'
            )
        ),
    )

    fig = go.Figure(data=data, layout=layout)

    plot(fig)

    #return py.iplot(fig, filename='stacked-bar', world_readable=True)
    #return plot(fig, filename='stacked-bar', world_readable=True)



def plot_country(data_type, country='Global'):
    """ Plots health data for a specified country over time """

    country_short = input_to_country(country)

    df_country1 = DF_HEALTH_ALL[(DF_HEALTH_ALL['country_short'] == country_short)]

    trace1 = go.Bar(
        x=df_country1['year'],
        y=df_country1[data_type],
        name=DATA_DICT[data_type] + ' in ' + country,
        marker=dict(
            color=COLOR_DICT
            ),
    )

    data = [trace1]
    layout = go.Layout(
        title=DATA_DICT[data_type] + ' in ' + country,
        yaxis=dict(
            #title='Averaged from 2000-2016',
        ),
    )


    fig = go.Figure(data=data, layout=layout)

    plot(fig)

    #return py.iplot(fig, filename='stacked-bar', world_readable=True)


def main(level, data_type, country='Global', highlight_country=None):
    """ Processes user inputs to generate the specified graphs """

    if level != 'country':
        return plot_summary(data_type, highlight_country)
    return plot_country(data_type, country)

if __name__ == "__main__":

    ARGS = plotParser(dict_to_help(DATA_DICT)).parse_args()

    main(
        ARGS.level, data_type=ARGS.data_type, country=ARGS.country,
        highlight_country=ARGS.highlight_country)
