#!/usr/bin/env python3

import argparse


def load_annotations(file_list):
	presence_matrix = {}
	location_matrix = {}
	sample_names = set()
	for f in file_list:
		sample_name = f.split('/')[-1].replace('.gff', '')
		sample_names.add(sample_name)
		with open(f, 'r') as handle:
			data = handle.read().split('\n')
			for line in data:
				if line == '##FASTA':
					break
				if not line or line[0] == '#':
					continue
				entries = line.split('\t')
				contig = entries[0]
				feature_type = entries[2]
				feature_start = entries[3]
				feature_end = entries[4]
				strand = entries[6]
				metadata = entries[8].split(';')
				gene_id = 'NA'
				gene_name = 'NA'
				gene_product = 'NA'
				ec_id = 'NA'
				for field in metadata:
					name, value = field.split('=')
					if name == 'ID':
						gene_id = value
					elif name == 'gene':
						gene_name = value
					elif name == 'product':
						gene_product = value
					elif name == 'eC_number':
						ec_id = value
				if gene_name == 'NA':
					continue
				if gene_name not in presence_matrix:
					presence_matrix[gene_name] = [{sample_name}, feature_type, gene_id, gene_product, ec_id]
					location_matrix[gene_name] = {sample_name: [contig, feature_start, feature_end, strand]}
				else:
					presence_matrix[gene_name][0].add(sample_name)
					location_matrix[gene_name][sample_name] = [contig, feature_start, feature_end, strand]
	return sample_names, presence_matrix, location_matrix


def output_results(S, P, L, out_dir):
	samples = sorted(list(S))
	with open(out_dir + '/gene_presence_absence_matrix.csv', 'w') as pa_out:
		pa_out.write('GeneName,FeatureType,GeneID,GeneProduct,EC_ID,{}\n'.format(
			','.join(samples)
		))
		for gene_name, values in P.items():
			pa_out.write('{},{},{},{},{}'.format(
				gene_name,
				values[1],
				values[2],
				values[3],
				values[4]
			))
			for sample_name in samples:
				if sample_name in values[0]:
					pa_out.write(',1')
				else:
					pa_out.write(',0')
			pa_out.write('\n')
	with open(out_dir + '/gene_location_matrix.csv', 'w') as l_out:
		l_out.write('GeneName,Sample,Contig,FeatureStart,FeatureEnd,Strand\n')
		for gene_name, subdict in L.items():
			for sample_name in samples:
				if sample_name in subdict:
					l_out.write('{},{},{},{},{},{}\n'.format(
						gene_name,
						sample_name,
						subdict[sample_name][0],
						subdict[sample_name][1],
						subdict[sample_name][2],
						subdict[sample_name][3]
					))


parser = argparse.ArgumentParser('compare_gene_annotations.py')
parser.add_argument('files', nargs='+', help='A list of all files to parse (e.g. my_data/*.gff)')
parser.add_argument('output', type=str, help='Output folder for presence/absence matrix and location data file')


if __name__ == '__main__':
	args = parser.parse_args()
	S, P, L = load_annotations(args.files)
	output_results(S, P, L, args.output)
