from Bio import Entrez
import pandas as pd
import time
import re

Entrez.email = "dusenk@uw.edu"

outfile = './data/RSVG_gb_metadata_5000-10000.csv'

database = "nuccore"
maxseqs = 10000
query = "human respiratory syncytial virus G"
filetype = "gb"
outmode = "xml"
batchsize = 100
begin = 5000

#Make generator and use regex to allow for whitespace between letters and numbers

#Write argparser

GT_LIST = ['GA1', 'GA2', 'GA3', 'GA4', 'GA5', 'GA6', 'GA7', 'NA1', 'NA2', 
		   'SAA1', 'ON1', 'GB1', 'GB2', 'GB3', 'GB4', 'SAB1', 'SAB2', 'SAB3',
		   'SAB4', 'URU1', 'URU2', 'BA1', 'BA2', 'BA3', 'BA4', 'BA5', 'BA6',
		   'BA7', 'BA8', 'BA9', 'BA10', 'BA11', 'BA12', 'THB']


def getIDs(db, retmax, term):
	search_handle = Entrez.esearch(db=db, retmax=retmax,
		term=term)
	search_record = Entrez.read(search_handle)
	search_handle.close()
	search_IDs = search_record['IdList']
	print('Search returned {0} hits.'.format(len(search_IDs)))
	return search_IDs

def gethandle(db, ids, firstseq, numseqs, rettype, retmode):
	handle = Entrez.efetch(db=db, id=ids, retstart=firstseq, retmax=numseqs, 
			 rettype=rettype, retmode=retmode)
	return handle

def find_subtype(meta_dict):
	"""Find subtype
	"""

	subtype = 'NaN'
	if ' A' in meta_dict['organism']:
		subtype = 'A'
	elif ' B' in meta_dict['organism']:
		subtype = 'B'

	elif subtype == 'NaN':
		for val in meta_dict.values():
			if re.search(r'RSVA\b', val) or re.search(r'type: A\b', val) or\
					re.search(r'group: A\b', val) or re.search(r'\bA/', val):
				subtype = 'A'
			elif re.search(r'RSVB\b', val) or re.search(r'type: B\b', val) or\
					re.search(r'group: B\b', val) or re.search(r'\bB/', val):
				subtype = 'B'	

	return subtype


def find_genotype(meta_dict, genotype_list):
	"""Script for extracting genotype data from genbank metadata.

	Args: 
		meta_dict = dictionary of metadata
		genotype_dict = list of possible genotypes
	"""
	genotype = 'NaN'
	for value in meta_dict.values():
		if 'genotype:' in value:
			for gt in genotype_list:
				if gt in value:
					genotype = gt
	return genotype


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

		sub_dict['subtype'] = find_subtype(sub_dict)

		sub_dict['genotype'] = find_genotype(sub_dict, GT_LIST)

		seqinfo.append(sub_dict)
	
	handle.close()

	seqinfo_df = pd.DataFrame(seqinfo)

	return seqinfo_df


def main():

	IDs = getIDs(database, maxseqs, query)
	start = begin
	metadata_frames = []
	while start <= (maxseqs - batchsize):
		handle = gethandle(database, IDs, start, batchsize, filetype, outmode)
		metadata_df = makedf(handle)
		metadata_frames.append(metadata_df)
		start = start + batchsize
		if start % 500 == 0:
			print('Processed {0} seqs'.format(start))
	
	full_df = pd.concat(metadata_frames, ignore_index=True, sort=False)
	#print(full_df.drop('db_xref', axis=1))
	full_df.to_csv(outfile)

if __name__ == '__main__':
	start_time = time.time()
	main()
	end_time = time.time()
	print('Program took {0:.3f} minutes to run.'.format((end_time - start_time)/60))

