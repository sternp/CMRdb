import argparse
import subprocess
import numpy as np
import pandas as pd
import os
import fastcluster
import scipy
import scipy.cluster.hierarchy as sch
from scipy.sparse import csr_matrix

'''
cat list_of_genomes.txt | parallel -j 32 "paste <(echo {} ) <(fgrep -c '>' {}) <(fgrep -v '>' {}| tr -d 'A|T|G|C|a|t|g|c' |
grep '.' -o | wc -l)" >>  MAG_stats.tsv
'''

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--list_of_genomes', dest='list_of_genomes' , help='list of genomes to run analysis, ideally this script should be run in that folder')
parser.add_argument('--checkm2',dest='checkm2', help='put in a checkm2 file with headers, but remember to put in the fasta extension in the genomes')
parser.add_argument('--cutoff', type=float, default=0.03, help='ANI cutoff value')
parser.add_argument('--threads', type=int, default=1, help='number of threads to use')
#parser.add_argument('--id',dest='id', help="sample ID")
parser.add_argument('--output',dest='output', help="output directory")
args = parser.parse_args()

# Run this if memory is an issue, a 65,000 MAG genome takes up 10GB
#os.system(f'skani sketch -t {args.threads} -o sketched_genomes -l list_of_genomes.txt --slow -m 2000')
#os.system(f'skani search -d sketched_genomes --ql list_of_genomes.txt -t 32 --robust --min-af 50 -s 92 -o skani_search.tsv')
#df0 = pd.read_csv('skani_search.tsv', sep='\t', usecols = ['Ref_file', 'Query_file', 'ANI'])
#df = df0.pivot_table(index='Ref_file', columns='Query_file', values='ANI', fill_value=0)

# Run this if memory is not an issue
#outdir = args.output + '/' + args.id + '/'
#os.system(f'mkdir {args.output}')
#os.mkdir(args.output)
#os.chdir(args.output)
os.system(f'/work/microbiome/sw/skani/skani triangle -l {args.list_of_genomes} -s 92 -m 2000 --slow --min-af 50 --robust --full-matrix -t {args.threads} | tail -n +2 > derep/skani_matrix.txt')
#os.system(f'/work/microbiome/sw/skani/skani triangle -l {args.list_of_genomes} -s 92 -m 2000 --min-af 75 --robust --full-matrix -t {args.threads} | tail -n +2 > skani_matrix.txt')

# Read in dendrogram file as a pandas DataFrame
df = pd.read_csv('derep/skani_matrix.txt', sep='\t', index_col=0, header=None)

df = df / 100
df = 1 - df
#df[df > args.cutoff] = 1
#df[df < args.cutoff] = 0
print(df)


'''
df[df <= args.cutoff] = 0
#df[df > args.cutoff] = 1
# Convert DataFrame to numpy array
'''
#X = df.to_numpy()
X = scipy.sparse.csr_matrix(df.values)
X_dense = X.toarray()
# Calculate pairwise distances
#D = fastcluster.distancematrix(X)

#Z = fastcluster.linkage(X, method='average')

Z = fastcluster.linkage(X_dense, method='average')

print(Z)

# Use fcluster to cluster data
clusters = sch.fcluster(Z, t=args.cutoff, criterion='distance')

# Save cluster assignments to output file
#with open(args.clusters, 'w') as f:
with open(args.output + '_clusters', 'w') as f:
    f.write('Name\tCluster\n')
    output_lines = [f'{df.index[i]}\tcluster_{clusters[i]}\n' for i in range(len(clusters))]
    f.writelines(output_lines)

df_stats = pd.read_csv('/work/microbiome/db/uhgg_v2/MAG_stats.tsv', sep='\t', names=['Name', 'num_contigs', 'amb_bases'])
df_checkm2 = pd.read_csv(args.checkm2, sep='\t', usecols=['Name', 'Completeness', 'Contamination'])
df_clusters = pd.read_csv(args.output + '_clusters', sep='\t', usecols = ['Name', 'Cluster'])

merged_df = pd.merge(df_clusters,df_stats, on='Name')
merged_df = pd.merge(merged_df, df_checkm2, on='Name')
merged_df['park_score'] = merged_df['Completeness'] - 5 * merged_df['Contamination'] - 5 * merged_df['num_contigs'] / 100 - 5 * merged_df['amb_bases'] / 100000
df_max = merged_df.groupby('Cluster')['park_score'].max()
df_result = df_max.to_frame().reset_index().merge(merged_df, on=['Cluster', 'park_score'])
df_result.drop_duplicates(subset=['Name'], keep='first', inplace=True, ignore_index=True)
# join the max new_column values with the original dataframe to get the corresponding rows
#print(df_max)
#print(df_result)
df_result.to_csv(args.output + '_derep_clusters', index=False, sep='\t')
