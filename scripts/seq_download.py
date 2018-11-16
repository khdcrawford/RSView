from Bio import Entrez
import pandas as pd

Entrez.email = "dusenk@uw.edu"

database = "nuccore"
maxseqs = 5000
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
	return search_IDs

def gethandle(db, ids, start, rettype, retmode):
	handle = Entrez.efetch(db=db, id=ids, retstart=start, rettype=rettype,
		retmode=retmode)
	return handle


def makedf(handle):
	records = Entrez.parse(handle)
	seqinfo = []
	for record in records:
		sub_dict = {}
		features = record['GBSeq_feature-table']
		strain_quals = features[0]['GBFeature_quals']
		for qual in strain_quals:
			qual_dict = dict(qual)
			sub_dict[qual_dict['GBQualifier_name']] = \
					qual_dict['GBQualifier_value']
		seqinfo.append(sub_dict)
	
	handle.close()

	seqinfo_df = pd.DataFrame(seqinfo)

	return seqinfo_df


def main():
	IDs = getIDs(database, maxseqs, query)
	handle = gethandle(database, IDs, 0, filetype, outmode)
	metadata_df = makedf(handle)
	metadata_df.to_csv('../data/RSVG_gb_metadata.csv')
	print(metadata_df.drop(['db_xref', 'host', 'mol_type'], axis=1))

if __name__ == '__main__':
	main()
