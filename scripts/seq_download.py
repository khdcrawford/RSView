from Bio import Entrez

Entrez.email = "dusenk@uw.edu"

database = "nuccore"
maxseqs = 20
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
	handle = Entrez.efetch(db=db, id=ids, retstart = start, rettype=rettype, 
			retmode=retmode)
	parsed = Entrez.parse(handle)
	records = list(parsed)
	handle.close()
	return(records)

def makedf(records):
	for i in range(len(records)):
		record = records[i]
		features = record['GBSeq_feature-table']
		for feature in features:
			featquals = feature['GBFeature_quals']
			for qual in featquals:
				qual_dict = dict(qual)
				#print(qual_dict['GBQualifier_name'] + ": " + qual_dict['GBQualifier_value'])

def main():
	IDs = getIDs(database, maxseqs, query)
	records = getrecords(database, IDs, 0, filetype, outmode)
	makedf(records)

if __name__ == '__main__':
	main()
