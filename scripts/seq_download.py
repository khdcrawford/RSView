from Bio import Entrez
import pandas as pd

Entrez.email = "dusenk@uw.edu"

database = "nuccore"
maxseqs = 1000
query = "human respiratory syncytial virus G"
filetype = "gb"
outmode = "xml"


def getIDs(db, retmax, term):
	search_handle = Entrez.esearch(db=db, retmax=retmax,
		term=term)
	search_record = Entrez.read(search_handle)
	search_handle.close()
	search_IDs = search_record['IdList']
	print('Search returned {0} hits.'.format(len(search_IDs)))
	return(search_IDs)

def getrecords(db, ids, start, rettype, retmode):
	handle = Entrez.efetch(db=db, id=ids, retstart=start, rettype=rettype,
		retmode=retmode)
	parsed = Entrez.parse(handle)
	records = list(parsed)
	handle.close()
	return(records)

def makedf(records):
	date_dict = {}
	except_count = 0
	for i in range(len(records)):
		if i % 100 == 0:
			print('{0} records processed'.format(i))
		sub_dict = {}
		record = records[i]
		features = record['GBSeq_feature-table']
		strain_quals = features[0]['GBFeature_quals']
		for qual in strain_quals:
			qual_dict = dict(qual)
			sub_dict[qual_dict['GBQualifier_name']] = \
					qual_dict['GBQualifier_value']
		try:
			date_dict[sub_dict['collection_date']] = sub_dict
		except:
			except_count += 1
			print(sub_dict)
			continue

	date_df = pd.DataFrame.from_dict(date_dict, orient='index')
	print(date_df.drop(['db_xref', 'mol_type', 'collection_date', 'lat_lon', 'note', 'host', 'isolate', 'isolation_source'], axis=1))
	print(except_count)
	print(len(date_df))

	#print(strain_df.drop(['db_xref', 'mol_type', 'strain', 'lat_lon', 'note', 'host'], axis=1))
	#return strain_dict

def main():
	IDs = getIDs(database, maxseqs, query)
	records = getrecords(database, IDs, 0, filetype, outmode)
	makedf(records)

if __name__ == '__main__':
	main()
