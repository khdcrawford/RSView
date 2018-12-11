""" This module processes the raw health file and generates full and summary CSVs
ready to be loaded into a dataframe for plotting """

import pandas as pd
import numpy as np
import country_converter as coco

from parsearguments import health_parser



def iso3_to_country(iso3):
    """ Take user input and convert it to the short version of the country name """

    if iso3 == 'Global':
        return 'Global'
    country = coco.convert(names=iso3, to='name_short')
    return country

def main(datadir):
    """ Process raw health CSV and generate full and summary CSVs """

    # need to change path to this csv
    df_orig = pd.read_csv(str(datadir) + '/health_data_RAW.csv')

    # change index to iso3 (3 letter country codes)
    df_orig.columns = df_orig.iloc[3]

    # remove spaces from columns names
    df_orig.columns = df_orig.columns.str.strip().str.replace(' ', '')

    # for NaN values in iso3, replace with 'Global'
    df_orig['iso3'] = df_orig['iso3'].replace(np.nan, 'Global', regex=True)

    # copy original dataframe, excluding blank lines
    df_with_na = df_orig[4:3319].copy()
    df_with_na.index = df_with_na['iso3']

    # remove entries with NaN or '-' values
    df_with_na = df_with_na.replace('-', np.nan)
    df_no_na = df_with_na.dropna(axis=0, how='any').fillna(0).copy()

    # remove % sign for fneo9, fpost9 and fufive9
    # '18' = 18%
    column_percents = ['fneo9', 'fpost9', 'fufive9']
    for i in column_percents:
        df_no_na[i] = df_no_na[i].str.rstrip('%').astype('float')


    column_numbers = [
        'nnd', 'pnd', 'neo9', 'post9', 'ufive9', 'rneo9', 'rpost9', 'rufive9']
    for i in column_numbers:
        df_no_na[i] = df_no_na[i].str.strip().str.replace(',', '')

    df_no_na = df_no_na.replace('-', np.nan)
    df_clean = df_no_na.dropna(axis=0, how='any').fillna(0).copy()

    for i in column_numbers:
        df_clean[i] = df_clean[i].astype('float')

    df_clean.rename(columns={'Country/area name': 'country'}, inplace=True)

    # create summary dataframe to groupby iso3 and calculate the mean
    # note: maybe change to median??
    df_summary = df_clean.groupby(df_clean.index).mean().reset_index()

    df_summary["country_short"] = [iso3_to_country(x) for x in df_summary.iso3.values]
    df_clean["country_short"] = [iso3_to_country(x) for x in df_clean.iso3.values]

    outfile_all = str(datadir) + '/health_data_all.csv'
    outfile_summary = str(datadir) + '/health_data_summary.csv'

    # export full and summary dataframes to csv
    df_clean.to_csv(outfile_all)
    df_summary.to_csv(outfile_summary)


if __name__ == "__main__":

    ARGS = health_parser().parse_args()

    main(ARGS.datadir)
