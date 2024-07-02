#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0
# Kyle M. Scott - kms309@miami.edu

import argparse
import io

import pandas as pd
from pysam import bcftools


def parser():
    parser = argparse.ArgumentParser(description='Parse RFMix output')
    parser.add_argument('-f', '--mspFile', required=True, help='Path to RFMix .msp.tsv output file')
    parser.add_argument('-r', '--region', required=True,
                        help='List with region to parse, comma-separated with CHR,SPOS,EPOS')
    parser.add_argument('-m', '--markers', required=False,
                        help='List with markers to append to output, providing genotype and ancestry, comma-separated in format chr:pos')
    parser.add_argument('-v', '--vcf', required=False,
                        help='Path to Phased VCF file, must be indexed. Must be used in conjunction with --markers')
    parser.add_argument('-o', '--out', required=True, help='Output file')
    args = parser.parse_args()
    return args


def get_population_mapping(msp_file: str) -> pd.DataFrame:
    with open(msp_file, 'r') as file:
        populations = file.readline().replace('#Subpopulation order/codes: ', '').split('\t')
    populations = [pop.strip().split('=') for pop in populations]
    pop_df = pd.DataFrame(populations, columns=['POP', 'NUM']).astype(str)
    return pop_df


def read_msp(msp_file: str) -> pd.DataFrame:
    df = pd.read_csv(msp_file, sep='\t', skiprows=1, dtype=str)
    column_types = {'spos': int, 'epos': int, 'n snps': int, 'sgpos': float, 'egpos': float}
    for col, col_type in column_types.items():
        df[col] = df[col].astype(col_type)
    return df


def parse_region(chm: str, spos: int, epos: int, msp: pd.DataFrame, pop_map: pd.DataFrame) -> pd.DataFrame:
    out_df = pd.DataFrame(columns=['ID', 'CHR', 'SPOS', 'EPOS', 'REGION_ANCESTRY', 'REGION_PROPORTION'])
    if chm not in msp['#chm'].values:
        raise ValueError(f"Chromosome {chm} not found in RFMix output. This could be a 'chr' prefix issue.")
    mask_chr = msp['#chm'] == chm
    mask1 = mask_chr & (msp['spos'] < spos) & (msp['epos'] > epos)
    mask2 = mask_chr & (msp['spos'] > spos) & (msp['epos'] < epos)
    mask3 = mask_chr & (msp['spos'] < spos) & (msp['epos'] > epos)
    mask4 = mask_chr & (msp['spos'] < spos) & (msp['epos'] > epos)
    mask = mask1 | mask2 | mask3 | mask4
    msp.loc[mask1, 'spos'] = spos
    msp.loc[mask3, 'epos'] = epos
    msp.loc[mask4, 'spos'] = spos
    msp.loc[mask4, 'epos'] = epos
    msp_region = msp[mask]

    pop_dict = pop_map.set_index('NUM')['POP'].to_dict()
    msp_region.iloc[:, 6:] = msp_region.iloc[:, 6:].replace(pop_dict)

    segment_lengths = {col: [] for col in msp_region.columns[6:]}
    total_length = msp_region['epos'].max() - msp_region['spos'].min()

    for col in msp_region.columns[6:]:
        start = msp_region.iloc[0]['spos']
        prev_value = msp_region.iloc[0][col]
        for j in range(1, len(msp_region)):
            if msp_region.iloc[j][col] != prev_value:
                end = msp_region.iloc[j - 1]['epos']
                segment_lengths[col].append([prev_value, ((end - start) / total_length).round(2)])
                start = msp_region.iloc[j]['spos']
                prev_value = msp_region.iloc[j][col]
        end = msp_region.iloc[len(msp_region) - 1]['epos']
        segment_lengths[col].append([prev_value, ((end - start) / total_length).round(2)])

    processed_data = []
    for key, values in segment_lengths.items():
        id_ = key.rsplit('.', 1)[0]
        ancestry = ':'.join([str(v[0]) for v in values])
        proportion = ':'.join([str(v[1]) for v in values])
        processed_data.append([id_, ancestry, proportion])

    processed_df = pd.DataFrame(processed_data, columns=['ID', 'REGION_ANCESTRY', 'REGION_PROPORTION'])
    processed_df['CHR'] = chm
    processed_df['SPOS'] = msp_region['spos'].min()
    processed_df['EPOS'] = msp_region['epos'].max()
    processed_df = processed_df.groupby('ID').agg({
        'CHR': 'first',
        'SPOS': 'first',
        'EPOS': 'first',
        'REGION_ANCESTRY': '|'.join,
        'REGION_PROPORTION': '|'.join
    }).reset_index()
    out_df = pd.concat([out_df, processed_df])
    return out_df


def get_genotype(vcf_file: str, markers: list) -> pd.DataFrame:
    header = bcftools.view(vcf_file, '-h')
    catch_str = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	"
    header = 'MARKERS\t' + header.split(catch_str)[1]
    query = bcftools.query(vcf_file, '-r', f'{markers}', '-f', '%CHROM:%POS[\t%GT]\n')
    vcf_df = pd.read_csv(io.StringIO(header + query), sep='\t')
    if len(vcf_df) > 1:
        vcf_df["TEMP"] = 0
        vcf_df = vcf_df.groupby('TEMP').agg(lambda x: ','.join(x)).reset_index()
        vcf_df = vcf_df.drop('TEMP', axis=1)
    return vcf_df


def get_marker_ancestry(markers: list, msp: pd.DataFrame, pop_map: pd.DataFrame) -> pd.DataFrame:
    marker_ancestry_df = pd.DataFrame()
    for marker in markers.split(','):
        chm, pos = marker.split(':')
        mask_chr = msp['#chm'] == chm
        mask_pos = (msp['spos'] <= int(pos)) & (msp['epos'] >= int(pos))
        mask = mask_chr & mask_pos
        ancestry = msp[mask].iloc[:, 6:].replace(pop_map.set_index('NUM')['POP'].to_dict())
        ancestry.insert(0, 'MARKERS', marker)
        ids = ancestry.columns[1:].str.split('.').str[0].unique().tolist()
        res = pd.wide_to_long(ancestry, ids, i="MARKERS", j="suffix_key", suffix='.*').reset_index().drop('suffix_key',
                                                                                                          axis=1)
        res = res.groupby('MARKERS').agg(lambda x: '|'.join(x)).reset_index()
        marker_ancestry_df = pd.concat([marker_ancestry_df, res])
    if len(marker_ancestry_df) > 1:
        marker_ancestry_df["TEMP"] = 0
        marker_ancestry_df = marker_ancestry_df.groupby('TEMP').agg(lambda x: ','.join(x)).reset_index()
        marker_ancestry_df = marker_ancestry_df.drop('TEMP', axis=1)
    return marker_ancestry_df


def transpose_marker_results(vcf_df: pd.DataFrame, marker_df: pd.DataFrame) -> pd.DataFrame:
    vcf_df = vcf_df.set_index('MARKERS').reset_index()
    vcf_df = vcf_df.melt(id_vars='MARKERS', var_name='ID', value_name='MARKER_GT')
    marker_df = marker_df.set_index('MARKERS').reset_index()
    marker_df = marker_df.melt(id_vars='MARKERS', var_name='ID', value_name='MARKER_ANCESTRY')
    return pd.merge(vcf_df, marker_df, on=['MARKERS', 'ID'])


def main():
    args = parser()
    pop_map = get_population_mapping(args.mspFile)
    msp = read_msp(args.mspFile)
    region = args.regions.split(',')
    if len(region) != 3:
        raise ValueError('Region must contain exactly three entries and be in format CHR,SPOS,EPOS')
    chm = region[0]
    spos = int(region[1])
    epos = int(region[2])
    anc_prop_df = parse_region(chm, spos, epos, msp, pop_map)
    if args.markers and args.vcf:
        vcf_df = get_genotype(args.vcf, args.markers)
        marker_df = get_marker_ancestry(args.markers, msp, pop_map)
        final_marker_df = transpose_marker_results(vcf_df, marker_df)
        anc_prop_df = pd.merge(anc_prop_df, final_marker_df, on='ID')
    anc_prop_df.to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
