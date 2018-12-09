import pandas as pd
import numpy as np
import country_converter as coco


outfile_all = './data/health_data_all.csv'
outfile_summary = './data/health_data_summary.csv'

def iso3_to_country(iso3):
    if iso3 == 'Global':
        return 'Global'
    else:
        country =  coco.convert(names=iso3, to='name_short')
        return country

def main():

	# need to change path to this csv
	df_orig = pd.read_csv('./data/health_data_RAW.csv')

	# change index to iso3 (3 letter country codes)
	df_orig.columns = df_orig.iloc[3]

	# for NaN values in iso3, replace with 'Global'
	df_orig['iso3'] = df_orig['iso3'].replace(np.nan, 'Global', regex=True)

	# copy original dataframe, excluding blank lines
	df_withNA = df_orig[4:3319].copy()
	df_withNA.index = df_withNA['iso3']

	# remove entries with NaN or '-' values
	df_withNA=df_withNA.replace('-', np.nan)
	df = df_withNA.dropna(axis=0, how='any').fillna(0).copy()

	# remove % sign for fneo9, fpost9 and fufive9
	# '18' = 18%
	column_percents = ['fneo9', 'fpost9', 'fufive9']
	for i in column_percents:
	    df[i] = df[i].str.rstrip('%').astype('float')


	column_numbers = [' nnd ', ' pnd ', ' neo9 ', ' post9 ', ' ufive9 ', ' rneo9 ', ' rpost9 ', ' rufive9 ']
	for i in column_numbers:
	    df[i] = df[i].str.strip().str.replace(',','')


	df=df.replace('-', np.nan)
	df_clean = df.dropna(axis=0, how='any').fillna(0).copy()

	for i in column_numbers:
	    df_clean[i] = df_clean[i].astype('float')

	df_clean.rename(columns={'Country/area name': 'country'}, inplace=True)

	# create summary dataframe to groupby iso3 and calculate the mean 
	# note: maybe change to median??
	df_summary = df_clean.groupby(df_clean.index).mean().reset_index()

	df_summary["country_short"] = [iso3_to_country(x) for x in df_summary.iso3.values]
	df_clean["country_short"] = [iso3_to_country(x) for x in df_clean.iso3.values]

	# export full and summary dataframes to csv
	df_clean.to_csv(outfile_all)
	df_summary.to_csv(outfile_summary)






	