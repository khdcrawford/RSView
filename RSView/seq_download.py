from Bio import Entrez
import pandas as pd
import time
import re
import argparse

#Write argparser

Entrez.email = "dusenk@uw.edu"

outfile = './data/RSVG_gb_metadata_15000+.csv'

database = "nuccore"
maxseqs = 20000
query = "human respiratory syncytial virus G"
filetype = "gb"
outmode = "xml"
begin = 15000

batchsize = min(maxseqs - begin, 100)

GTA_LIST = [r'\bGA\s?[0-9]*\b', r'\bNA\s?[0-9]*\b', r'\bSAA\s?[0-9]*\b', 
		   r'\bON\s?[0-9]*\b']

GTB_LIST = [r'\bGB\s?[0-9]*\b', r'\bSAB\s?[0-9]*\b', r'\bURU\s?[0-9]*\b', 
			r'\bBA\s?[0-9]*\b', r'\bBA\s?IV\b', r'\bTHB\b']


def getIDs(db, retmax, term):
	"""Retrieve genbank sequence IDs matching query term. 

	Args:
		`db` (str)
			string specifying genbank database to search
		`retmax` (int)
			maximum number of sequence IDs to return
		`term` (str)
			query term

	Returns:
		`search_IDs` (list)
			list of sequence IDs returned from the search
	"""
	search_handle = Entrez.esearch(db=db, retmax=retmax,
		term=term)
	search_record = Entrez.read(search_handle)
	search_handle.close()
	search_IDs = search_record['IdList']
	return search_IDs

def gethandle(db, ids, firstseq, dload_size, rettype, retmode):
	"""Download Entrez 'handle' for downloading seqs of interest

	Args:
		`db` (str)
			Genbank database containing sequences of interest
		`ids` (list)
			sequence IDs returned from `getIDs` or other Entrez search
		`firstseq` (int)
			index of first sequence to download
		`dload_size` (int)

		See `Entrez.efetch` help for help.

	Return:
		`handle` 
			Entrez object containing sequence information
	"""

	handle = Entrez.efetch(db=db, id=ids, retstart=firstseq, retmax=dload_size,
			 rettype=rettype, retmode=retmode)
	return handle

def find_subtype(meta_dict):
	"""Find subtype from dictionary of sequence metadata.

	Args:
		meta_dict (dict): dictionary of metadata downloaded from genbank

	Returns:
		subtype (str): RSV subtype as one letter string, 'A' or 'B'.
	"""

	subtype = 'NaN'
	if ' A' in meta_dict['organism']:
		subtype = 'A'
	elif ' B' in meta_dict['organism']:
		subtype = 'B'

	elif subtype == 'NaN':
		for val in meta_dict.values():
			if re.search(r'RSV\s?A\b', val) or \
					re.search(r'type: A\b', val) or	\
					re.search(r'group: A\b', val) or re.search(r'\bA\b', val):
				subtype = 'A'
			elif re.search(r'RSV\s?B\b', val) or \
					re.search(r'type: B\b', val) or \
					re.search(r'group: B\b', val) or re.search(r'\bB\b', val):
				subtype = 'B'

	return subtype


def find_genotype(meta_dict, genotype_listA, genotype_listB):
	"""Script for extracting genotype data from genbank metadata.
	
	If the genotype is found, but the subtype is still 'NaN', populate
	subtype data based on genotype. 


	Args: 
		`meta_dict` (dict) 
			dictionary of metadata
		`genotype_listA` (list)
			list of possible genotypes for subtype A
		`genotype_listB` (list)
			list of possible genotypes for subtype B

	Returns:
		`typed_dict` (dict)
			dictionary of metadata with genotype (and missing subtype) data
			filled in.
	"""
	typed_dict = meta_dict
	genotype = 'NaN'
	for value in meta_dict.values():
		if 'genotype:' in value:
			for gt in genotype_listA:
				if re.search(gt, value):
					genotype = re.findall(gt, value)[0]
					if meta_dict['subtype'] == 'NaN':
						typed_dict['subtype'] = 'A'

			for gt in genotype_listB:
				if re.search(gt, value):
					genotype = re.findall(gt, value)[0]
					if meta_dict['subtype'] == 'NaN':
						typed_dict['subtype'] = 'B'

	typed_dict['genotype'] = genotype
	
	return typed_dict


def makedf(handle):
	""" 
	Convert Genbank sequence data into dataframe containing necessary
	metadata.

	Args:
		`handle`
			Entrez object containing information downloaded from GenBank

	Returns:
		`seqinfo_df` (DataFrame)
			pandas DataFrame containing downloaded sequence and metadata
	"""
	records = Entrez.parse(handle)
	seqinfo = []
	for record in records:
		sub_dict = {}
		features = record['GBSeq_feature-table']

		#Retrieve metadata
		strain_quals = features[0]['GBFeature_quals']
		for qual in strain_quals:
			qual_dict = dict(qual)
			if 'GBQualifier_value' in qual_dict.keys():
				sub_dict[qual_dict['GBQualifier_name']] = \
						qual_dict['GBQualifier_value']

		#Retrieve G protein sequence
		for feat_dict in features[1:]:
			if 'GBFeature_quals' in feat_dict.keys():
				for feat_qual in feat_dict['GBFeature_quals']:
					if 'GBQualifier_value' in feat_qual.keys():
						if re.search(r'\bG\b', feat_qual['GBQualifier_value'])\
								or re.search(r'\battachment.*protein\b', 
								             feat_qual['GBQualifier_value']):
							G_quals = feat_dict['GBFeature_quals']
							for q in G_quals:
								if q['GBQualifier_name'] == 'translation':
									sub_dict['G_seq'] = q['GBQualifier_value']

		sub_dict['subtype'] = find_subtype(sub_dict)

		sub_dict = find_genotype(sub_dict, GTA_LIST, GTB_LIST)

		seqinfo.append(sub_dict)
	
	handle.close()

	seqinfo_df = pd.DataFrame(seqinfo)

	return seqinfo_df


def main():
	"""Download sequence data and return dataframe"""

	IDs = getIDs(database, maxseqs, query)
	numseqs = len(IDs)
	print('Search returned {0} hits.'.format(numseqs))

	start = begin
	print('Downloading metadata for seqs: {0} to {1}'.format(start, numseqs))

	metadata_frames = []
	while start <= (numseqs - batchsize):
		handle = gethandle(database, IDs, start, batchsize, filetype, outmode)
		metadata_df = makedf(handle)
		metadata_frames.append(metadata_df)
		start = start + batchsize
		if (start-begin) % 500 == 0:
			print('Processed {0} seqs'.format(start-begin))
	
	if start != numseqs: #Process final seqs
		handle = gethandle(database, IDs, start, numseqs, filetype, outmode)
		metadata_df = makedf(handle)
		metadata_frames.append(metadata_df)
	
	full_df = pd.concat(metadata_frames, ignore_index=True, sort=False)
	assert len(full_df.index) == (numseqs-begin), 'Exported unexpected ' \
			'number of seqs. Expected: {0} Retrieved: {1}'.format(
			(numseqs-begin), len(full_df.index))
	#print(full_df.drop(['db_xref', 'country', 'host'], axis=1))
	full_df.to_csv(outfile)

if __name__ == '__main__':
	start_time = time.time()
	main()
	end_time = time.time()
	print('Program took {0:.3f} minutes to run.'.format((end_time - start_time)/60))

