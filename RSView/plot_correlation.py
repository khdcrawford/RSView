import map_rsv
import pandas as pd
from plotly import tools
import plotly.plotly as py
import plotly.graph_objs as go
import country_converter as coco
import plot_rsv
from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import seaborn as sns


from parsearguments import correlationParser

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

def get_ratio_year(df):
	df_ratio = pd.DataFrame(columns=['country', 'year', 'ratio', 'count_A', 'count_B', 'count_total'])
	countries = df['country'].unique()
	for i in countries:
		country_df = df[df['country']==i]
		years = country_df['year'].unique()
		for j in years:
			years_df = country_df[country_df['year']==j]
			if years_df['subtype'].count() == 1:
				lone_subtype = years_df['subtype'].item()
				if lone_subtype == 'A':
					percent_A = 1
					count_A = years_df['count'].item()
					count_B = 0
					count_total = count_A + count_B
				else:
					percent_A = 0
					count_A = 0
					count_B = years_df['count'].item()
					count_total = count_A + count_B
			else:
				count_A = float(years_df[years_df['subtype'] == 'A']['count'])
				count_B = float(years_df[years_df['subtype'] == 'B']['count'])
				percent_A = count_A/(count_A + count_B)
				count_total = count_A + count_B
			df_ratio = df_ratio.append(pd.Series([i, j, percent_A, count_A, count_B, count_total], index=df_ratio.columns, name=i))
	return df_ratio


def merge_ratio_health(df_ratio):
	
	df_health_summary = pd.read_csv('./data/health_data_summary.csv')
	df_ratio["iso3"] = [country_to_iso3(x) for x in df_ratio.country.values]
	
	df_merge = pd.merge(df_ratio, df_health_summary, how='inner', on=['iso3'])

	return df_merge

def merge_ratio_health_year(df_ratio):
	
	df_health_all = pd.read_csv('./data/health_data_all.csv')
	df_ratio["iso3"] = [country_to_iso3(x) for x in df_ratio.country.values]
	
	df_merge = pd.merge(df_ratio, df_health_all, how='inner', on=['iso3','year'])

	return df_merge



def plot_ratio(df_merge, data_type):

	trace1 = go.Scatter(
		x = df_merge['ratio'],
		y = df_merge[data_type],
		mode = 'markers',
		marker= dict(size= 14,
					line= dict(width=1),
					color='rgba(204,204,204,1)'
				   ),
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
			title = plot_rsv.DATA_DICT[data_type],
			# ticklen= 5,
			# gridwidth= 2,
		),
		showlegend = False
	)

	data = [trace1]

	fig= go.Figure(data=data, layout=layout)
	plot(fig)



def plot_ratio_year(df_merge, data_type):

	df_merge['hover'] = df_merge['country_short'] + ', ' + df_merge['year'].astype(int).map(str)

	trace1 = go.Scatter(
		x = df_merge['ratio'],
		y = df_merge[data_type],
		mode = 'markers',
		marker= dict(size= 14,
					line= dict(width=1),
					color=df_merge['year'],
					colorscale='RdBu',
					showscale=True,
					# cmin=1990,
					# cmax=2018
				   ),
		text = df_merge['hover'] # The hover text goes here... 
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
			title = plot_rsv.DATA_DICT[data_type],
			# ticklen= 5,
			# gridwidth= 2,
		),
		showlegend = False
	)

	data = [trace1]

	fig= go.Figure(data=data, layout=layout)
	plot(fig)



def main(level, data_type):
	"""
	Run organize_data, count_types, map_rsv
	"""
	rsv_df = map_rsv.organize_data(map_rsv.DATAFILES, map_rsv.GENOTYPE_DICT)
	rsv_df = rsv_df.dropna(subset = ['subtype']).copy()

	rsv_df_year = rsv_df.dropna(subset = ['year']).copy()
	rsv_df_year['year'] = rsv_df_year['year'].astype(int)

	if level == 'all':
		rsv_df_count = count_subtypes(rsv_df)
		df_ratio_all = get_ratio_all(rsv_df_count)
		df_merged = merge_ratio_health(df_ratio_all)
		plot_ratio(df_merged, data_type)

	else:
		rsv_df_count = count_subtypes_year(rsv_df)
		df_ratio_year = get_ratio_year(rsv_df_count)
		df_merged = merge_ratio_health_year(df_ratio_year)
		plot_ratio_year(df_merged, data_type)

	

if __name__ == "__main__":

	args = correlationParser(plot_rsv.dictToHelp(plot_rsv.DATA_DICT)).parse_args()

	main(args.level, args.data_type)











