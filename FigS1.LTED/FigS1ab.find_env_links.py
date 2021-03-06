import sys, argparse

import numpy as np
from scipy.stats import t

def parse_args():
	""" Return dictionary of command line arguments
	"""
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawTextHelpFormatter,
		usage=argparse.SUPPRESS)

	parser.add_argument('--edge', type=str, dest='edge_path', required=True,
			help="""Path to edge file as input in a special space seperated format, e.g. node_1 (pp) node2 = 1.00""")
	parser.add_argument('--otu', type=str, dest='otu_path', required=True,
			help="""Path to otu abundance file in a tab seperated format. Header and row name present as the first row and first column from the 2nd row""")
	parser.add_argument('--env', type=str, dest='env_path', required=True,
			help="""Path to env variable data file in a tab seperated format. Header and row name present as the first row and first column from the 2nd row""")
	parser.add_argument('--r2', type=float, dest='th_r2', default=0.6,
			help="""Specify a r2 value as correlation threshold (Pearson) for identify significant links between env variables and OTUs. All links below this value will not be considered. Default: 0.6""")
	parser.add_argument('--pval', type=float, dest='th_p', default=0.05,
			help="""Specify a p-value as threshold for testing whether the links between env variables and OTU passed significant level of correlation. All links tested with p-value higher than this value will not be considered. Default: 0.05""")
	parser.add_argument('--out', type=str, dest='out', default="/dev/stdout",
			help="""Path to output which includes a succinct summary of valid links. Default: standard output""")
	parser.add_argument('--verbose', action='store_true', dest='verbose', 
			help="""Toggle to output details for each env dependent link. Default: False""")

	return vars(parser.parse_args())

def read_edges(edge_path):
	edges = []

	with open(edge_path, 'r') as fh:
		for nl, line in enumerate(fh):
			items = line.rstrip().split(" ")
			node1 = items[0]
			node2 = items[2]
			strength = float(items[4])

			if (node1 < node2):
				edges.append([node1, node2, strength])	
			else:
				edges.append([node2, node1, strength])	

	return edges 

def read_otus(otu_path, otus_sub=None):
	header = []
	rownames = []
	otu_abunds = []

	with open(otu_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split("\t")
			if len(header) == 0:
				header = items
				continue
			else:
				if otus_sub is None:
					rownames.append(items[0])
					otu_abunds.append([float(item) for item in items[1:]])
				else:
					if items[0] in otus_sub:
						rownames.append(items[0])
						otu_abunds.append([float(item) for item in items[1:]])

	
	np_order = np.argsort(header)
	header_np = np.array(header)[np_order]
	otu_abunds_np = np.array(otu_abunds)[:,np_order]

	return header_np, np.array(rownames), otu_abunds_np

def read_env(env_path):
	header = []
	rownames = []
	env_data = []

	with open(env_path, 'r') as fh:
		for line in fh:
			items = line.rstrip().split("\t")
			if len(header) == 0:
				header = items[1:]
				continue
			else:
				rownames.append(items[0])
				env_data.append([float(item) for item in items[1:]])
	
	np_order = np.argsort(rownames)
	env_data_np = np.array(env_data).transpose()[:,np_order]

	header_np = np.array(rownames)[np_order]
	rownames_np = np.array(header)

	return header_np, rownames_np, env_data_np

def pair_corr(otu_abunds, env_data):
	otu_means = np.mean(otu_abunds, axis=1)[:,np.newaxis]
	otu_stds = np.std(otu_abunds, axis=1)[:,np.newaxis]

	env_means = np.mean(env_data, axis=1)[:,np.newaxis]
	env_stds = np.std(env_data, axis=1)[:,np.newaxis]

	corr_mat = []
	pval_mat = []
	for i in range(0, env_data.shape[0]):
		# Cov here was calculated by dividing N (No. of samples), unbiased estimator should be N-1; However, pearsonr in scipy and cor.test in R both use N as denominator.
		cov_otu_env = np.sum((otu_abunds - otu_means)*(env_data[i,:] - env_means[i]), axis=1)/(otu_abunds.shape[1])
		cov_otu_env = cov_otu_env[:,np.newaxis]
		corr_otu_env = cov_otu_env/(otu_stds * env_stds[i])
		
		# significance level was evaluated by converting r to t using df=N-2 and thus tested in t-distribution
		# formula: r*sqrt(N-2)/sqrt(1-r^2)
		# source: https://stats.libretexts.org/Bookshelves/Introductory_Statistics/Book%3A_Introductory_Statistics_(OpenStax)/12%3A_Linear_Regression_and_Correlation/12.05%3A_Testing_the_Significance_of_the_Correlation_Coefficient
		t_otu_env = corr_otu_env / np.sqrt((1 - np.square(corr_otu_env))/(otu_abunds.shape[1] - 2))
		t_cdf = np.vectorize(t.cdf)

		# a scalar of 2 was mutilplied for two-tailed test, the resulting p-values are the same as what pearsonr and cor.test give
		pval_otu_env = 2*(1-t_cdf(np.abs(t_otu_env), otu_abunds.shape[1] - 2))

		if len(corr_mat) == 0:
			corr_mat = corr_otu_env
			pval_mat = pval_otu_env
		else:
			corr_mat = np.concatenate((corr_mat, corr_otu_env), axis=1)
			pval_mat = np.concatenate((pval_mat, pval_otu_env), axis=1)

	return corr_mat, pval_mat

def search_env_links(edges, otu_inds, envs, corr_mat, pval_mat, corr_th, pval_th):
	env_links = []

	otu_env_mask = (np.abs(corr_mat) >= corr_th) & (pval_mat <= pval_th)

	for edge in edges:
		i = otu_inds[edge[0]]
		j = otu_inds[edge[1]]

		
		if edge[2] > 0:
			row_mask = otu_env_mask[i,:] & otu_env_mask[j,:] & (otu_env_mask[i,:] * otu_env_mask[j,:] > 0)
		else:
			row_mask = otu_env_mask[i,:] & otu_env_mask[j,:] & (otu_env_mask[i,:] * otu_env_mask[j,:] < 0)

		if np.any(row_mask):
			env_links.append(np.concatenate([np.array(edge), envs[row_mask], corr_mat[i,row_mask], pval_mat[i,row_mask], corr_mat[j,row_mask], pval_mat[j,row_mask]]))
	
	return env_links

def main():
	args = parse_args()
	
	edges = read_edges(args['edge_path'])
	
	node_wt_edge = dict()
	for edge in edges:
		node_wt_edge[edge[0]] = 1
		node_wt_edge[edge[1]] = 1

	samples_otu, otu_ids, otu_abunds = read_otus(args['otu_path'], node_wt_edge)
	samples_env, env_vars, env_data = read_env(args['env_path'])

	otu_inds = dict()
	for i, otu_id in enumerate(otu_ids):
		otu_inds[otu_id] = i

	corr_mat, pval_mat = pair_corr(otu_abunds, env_data)
	env_links = search_env_links(edges, otu_inds, env_vars, corr_mat, pval_mat, args['th_r2'], args['th_p'])

	with open(args['out'], 'w') as fh:
		fh.write("total number of links: {}\n".format(len(edges)))
		fh.write("total number of env dependent links: {}\n".format(len(env_links)))

		if args['verbose'] is True:
			for env_link in env_links:
				fh.write("{}\n".format('\t'.join([str(ele) for ele in env_link])))

if __name__ == "__main__":
	main()
