#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0
# Kyle M. Scott - kms309@miami.edu

import argparse
import pandas as pd

def parser():
	parser = argparse.ArgumentParser(description='Parse RFMix output')
	parser.add_argument('-f', '--mspFiles', required=True, help='File listing RFMix MSP output files, one per line')
	parser.add_argument('-r', '--regions', required=True, help='File listing regions to parse, tab-separated with columns CHR, SPOS, EPOS')
	parser.add_argument('-s', '--samples', required=False, help='File listing samples to include, one per line')
	parser.add_argument('-o', '--out', required=True, help='Output file')
	args = parser.parse_args()
	return args

def check_first_line(files: list):
	first_line = None
	for file in files:
		with open(file, 'r') as f:
			if first_line is None:
				first_line = f.readline().strip()
			elif first_line != f.readline().strip():
				raise ValueError(f"Subpopulation codes do not match in all files. Ensure RFMix results come from the same run.")
	return

def get_population_mapping(files: list) -> pd.DataFrame:
	file_path = files[0]
	with open(file_path, 'r') as file:
		populations = file.readline().replace('#Subpopulation order/codes: ', '').split('\t')
	populations = [pop.strip().split('=') for pop in populations]
	pop_df = pd.DataFrame(populations, columns=['POP', 'NUM']).astype(str)
	return pop_df

def read_msp(files: list) -> pd.DataFrame:
	dfs = []
	for file in files:
		df = pd.read_csv(file, sep='\t', skiprows=1, dtype=str)
		column_types = {'spos': int, 'epos': int, 'n snps': int, 'sgpos': float, 'egpos': float}
		for col, col_type in column_types.items():
			df[col] = df[col].astype(col_type)
		dfs.append(df)
	columns = [df.columns.tolist() for df in dfs]
	if not all([col == columns[0] for col in columns]):
		raise ValueError(f"Columns do not match in all files. Ensure RFMix results come from the same run.")
	return pd.concat(dfs)

def subset_samples(sample_list: list, msp: pd.DataFrame) -> pd.DataFrame:
	columns = []
	for sample in sample_list:
		columns.extend([f'{sample}.0', f'{sample}.1'])
	columns = ['#chm', 'spos', 'epos', 'n snps', 'sgpos', 'egpos'] + columns
	return msp[columns]

def parse_regions(regions: pd.DataFrame, rfmix: pd.DataFrame, pop_map: pd.DataFrame) -> pd.DataFrame:
	out_df = pd.DataFrame(columns=['ID', 'CHR', 'SPOS', 'EPOS', 'Ancestry', 'Proportion'])
	i = 0
	while i < len(regions):
		region_row = regions.iloc[i]
		if region_row['CHR'] not in rfmix['#chm'].values:
			i += 1
			continue
		working_rfmix = rfmix.copy()
		mask_chr = working_rfmix['#chm'] == region_row['CHR']
		mask1 = mask_chr & (working_rfmix['spos'] < region_row['SPOS']) & (working_rfmix['epos'] > region_row['SPOS'])
		mask2 = mask_chr & (working_rfmix['spos'] > region_row['SPOS']) & (working_rfmix['epos'] < region_row['EPOS'])
		mask3 = mask_chr & (working_rfmix['spos'] < region_row['EPOS']) & (working_rfmix['epos'] > region_row['EPOS'])
		mask4 = mask_chr & (working_rfmix['spos'] < region_row['SPOS']) & (working_rfmix['epos'] > region_row['EPOS'])
		mask = mask1 | mask2 | mask3 | mask4
		working_rfmix.loc[mask1, 'spos'] = region_row['SPOS']
		working_rfmix.loc[mask3, 'epos'] = region_row['EPOS']
		working_rfmix.loc[mask4, 'spos'] = region_row['SPOS']
		working_rfmix.loc[mask4, 'epos'] = region_row['EPOS']
		rfmix_region = working_rfmix[mask]

		pop_dict = pop_map.set_index('NUM')['POP'].to_dict()
		rfmix_region.iloc[:, 6:] = rfmix_region.iloc[:, 6:].replace(pop_dict)

		segment_lengths = {col: [] for col in rfmix_region.columns[6:]}
		total_length = rfmix_region['epos'].max() - rfmix_region['spos'].min()

		for col in rfmix_region.columns[6:]:
			start = rfmix_region.iloc[0]['spos']
			prev_value = rfmix_region.iloc[0][col]
			for j in range(1, len(rfmix_region)):
				if rfmix_region.iloc[j][col] != prev_value:
					end = rfmix_region.iloc[j-1]['epos']
					segment_lengths[col].append([prev_value, ((end - start)/total_length).round(2)])
					start = rfmix_region.iloc[j]['spos']
					prev_value = rfmix_region.iloc[j][col]
			end = rfmix_region.iloc[len(rfmix_region)-1]['epos']
			segment_lengths[col].append([prev_value, ((end - start)/total_length).round(2)])

		processed_data = []
		for key, values in segment_lengths.items():
			id_ = key.rsplit('.', 1)[0]
			ancestry = ':'.join([str(v[0]) for v in values])
			proportion = ':'.join([str(v[1]) for v in values])
			processed_data.append([id_, ancestry, proportion])

		processed_df = pd.DataFrame(processed_data, columns=['ID', 'Ancestry', 'Proportion'])
		processed_df['CHR'] = region_row['CHR']
		processed_df['SPOS'] = rfmix_region['spos'].min()
		processed_df['EPOS'] = rfmix_region['epos'].max()
		processed_df = processed_df.groupby('ID').agg({
			'CHR': 'first',
			'SPOS': 'first',
			'EPOS': 'first',
			'Ancestry': '|'.join,
			'Proportion': '|'.join
		}).reset_index()
		out_df = pd.concat([out_df, processed_df])
		i += 1
	return out_df

def main():
	args = parser()
	files = [line.strip() for line in open(args.mspFiles, 'r')]
	regions_df = pd.read_csv(args.regions, sep='\t')
	check_first_line(files)
	pop_map = get_population_mapping(files)
	msp = read_msp(files)
	if args.samples:
		sample_list = [line.strip() for line in open(args.samples, 'r')]
		msp = subset_samples(sample_list, msp)
	out_df = parse_regions(regions_df, msp, pop_map)
	out_df.to_csv(args.out, sep='\t', index=False)

if __name__ == '__main__':
	main()
