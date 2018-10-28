from Bio import Entrez

Entrez.email = "dusenk@uw.edu"

search_handle = Entrez.esearch(db="nuccore", retmax=100, 
		term="human respiratory syncytial virus G")
search_record = Entrez.read(search_handle)
search_handle.close()
search_IDs = search_record['IdList']

print('Search returned {0} hits.'.format(len(search_IDs)))

handle = Entrez.efetch(db="nuccore", id=search_IDs,
		rettype="gb", retmode="xml")
parsed = Entrez.parse(handle)
records = list(parsed)
handle.close()

for i in range(len(records)):
	record = records[i]
	features = record['GBSeq_feature-table']
	for feature in features:
		featquals = feature['GBFeature_quals']
		for qual in featquals:
			qual_dict = dict(qual)
			#print(qual_dict['GBQualifier_name'] + ": " + qual_dict['GBQualifier_value'])
