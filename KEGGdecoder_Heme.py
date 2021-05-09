#!/usr/bin/python

'''
Usage: python KEGG-decoder.py <KOALA INPUT> <FUNCTION LIST FORMAT>
Designed to parse through a blastKoala or ghostKoala output to determine
the completeness of various KEGG pathways
Dependencies:
Pandas - http://pandas.pydata.org/pandas-docs/stable/install.html
Seaborn - http://seaborn.pydata.org/installing.html
matplotlib - http://matplotlib.org/users/installing.html
For extended information about KEGG assignments, genes and pathways,
please see accompanying document "KOALA_definitions.txt"
'''

def C5_PPH(ko_match):
#Check for presence of 9 genes
	total = 0
#glutamyl-tRNA reductase, glutamate-1-semialdehyde 2,1-aminomutase
#porphobilinogen synthase, hydroxymethylbilane synthase
#uroporphyrinogen decarboxylase, ferrochelatase
	single_ko = ['K02492', 'K01845', 'K01698', 'K01749', 'K01599', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#coproporphyrinogen III oxidase
	if ('K00228' in ko_match or 'K02495' in ko_match):
		total += 1
#protoporphyrinogen oxidase
	if ('K00230' in ko_match or 'K00231' in ko_match or 'K08973' in ko_match):
		total += 1
	value = float(total)/float(9)
	return {'C5_PPH': float("%.2f" % (value))}
	
def C5_CPH(ko_match):
#Check for presence of 9 genes
	total = 0
#glutamyl-tRNA reductase, glutamate-1-semialdehyde 2,1-aminomutase
#porphobilinogen synthase, hydroxymethylbilane synthase
#uroporphyrinogen decarboxylase, coproporphyrinogen III oxidase, ferrochelatase
	single_ko = ['K02492', 'K01845', 'K01698', 'K01749', 'K01599', 'K00231', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(9)
	return {'C5_CPH': float("%.2f" % (value))}
	
def C5_SIRO(ko_match):
#Check for presence of 11 genes
	total = 0
#glutamyl-tRNA reductase, glutamate-1-semialdehyde 2,1-aminomutase
#porphobilinogen synthase, hydroxymethylbilane synthase
#siroheme decarboxylase, Fe-coproporphyrin III synthase
	single_ko = ['K02492', 'K01845', 'K01698', 'K01749', 'K22225', 'K22226']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#uroporphyrin-III C-methyltransferase
	if ('K00589' in ko_match or 'K02302' in ko_match or 'K02303' in ko_match or 'K02496' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#precorrin-2 dehydrogenase
	if ('K02302' in ko_match or 'K02304' in ko_match):
		total += 1
#sirohydrochlorin ferrochelatase
	if ('K02302' in ko_match or 'K02304' in ko_match or 'K03794' in ko_match):
		total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(11)
	return {'C5_SIRO': float("%.2f" % (value))}

def C4_PPH(ko_match):
#Check for presence of 8 genes
	total = 0
#5-aminolevulinate synthase
#porphobilinogen synthase, hydroxymethylbilane synthase
#uroporphyrinogen decarboxylase, ferrochelatase
	single_ko = ['K00643', 'K01698', 'K01749', 'K01599', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#coproporphyrinogen III oxidase
	if ('K00228' in ko_match or 'K02495' in ko_match):
		total += 1
#protoporphyrinogen oxidase
	if ('K00230' in ko_match or 'K00231' in ko_match or 'K08973' in ko_match):
		total += 1
	value = float(total)/float(8)
	return {'C4_PPH': float("%.2f" % (value))}
	
def C4_CPH(ko_match):
#Check for presence of 8 genes
	total = 0
#5-aminolevulinate synthase
#porphobilinogen synthase, hydroxymethylbilane synthase
#uroporphyrinogen decarboxylase, coproporphyrinogen III oxidase, ferrochelatase
	single_ko = ['K00643', 'K01698', 'K01749', 'K01599', 'K00231', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(8)
	return {'C4_CPH': float("%.2f" % (value))}
	
def C4_SIRO(ko_match):
#Check for presence of 10 genes
	total = 0
#5-aminolevulinate synthase
#porphobilinogen synthase, hydroxymethylbilane synthase
#siroheme decarboxylase, Fe-coproporphyrin III synthase
	single_ko = ['K00643', 'K01698', 'K01749', 'K22225', 'K22226']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#uroporphyrin-III C-methyltransferase
	if ('K00589' in ko_match or 'K02302' in ko_match or 'K02303' in ko_match or 'K02496' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#precorrin-2 dehydrogenase
	if ('K02302' in ko_match or 'K02304' in ko_match):
		total += 1
#sirohydrochlorin ferrochelatase
	if ('K02302' in ko_match or 'K02304' in ko_match or 'K03794' in ko_match):
		total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(10)
	return {'C4_SIRO': float("%.2f" % (value))}
	
def upper_C5(ko_match):
#Check for presence of 2 genes
	total = 0
#glutamyl-tRNA reductase, glutamate-1-semialdehyde 2,1-aminomutase
	single_ko = ['K02492', 'K01845']
	for i in single_ko:
		if i in ko_match:
			total += 1
	value = float(total)/float(2)
	return {'upper_C5': float("%.2f" % (value))}

def upper_C4(ko_match):
#Check for presence of 1 genes
	total = 0
#5-aminolevulinate synthase
	single_ko = ['K00643']
	for i in single_ko:
		if i in ko_match:
			total += 1
	value = float(total)/float(1)
	return {'upper_C4': float("%.2f" % (value))}

def Common(ko_match):
#Check for presence of 3 genes
	total = 0
#porphobilinogen synthase, hydroxymethylbilane synthase
	single_ko = ['K01698', 'K01749']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrinogen-III synthase
	if ('K01719' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
	value = float(total)/float(3)
	return {'Common': float("%.2f" % (value))}

def lower_PPH(ko_match):
#Check for presence of 4 genes
	total = 0
#uroporphyrinogen decarboxylase, ferrochelatase
	single_ko = ['K01599', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#coproporphyrinogen III oxidase
	if ('K00228' in ko_match or 'K02495' in ko_match):
		total += 1
#protoporphyrinogen oxidase
	if ('K00230' in ko_match or 'K00231' in ko_match or 'K08973' in ko_match):
		total += 1
	value = float(total)/float(4)
	return {'lower_PPH': float("%.2f" % (value))}

def lower_CPH(ko_match):
#Check for presence of 4 genes
	total = 0
#uroporphyrinogen decarboxylase, coproporphyrinogen III oxidase, ferrochelatase
	single_ko = ['K01599', 'K00231', 'K01772']
	for i in single_ko:
		if i in ko_match:
			total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(4)
	return {'lower_CPH': float("%.2f" % (value))}
	
def lower_SIRO(ko_match):
#Check for presence of 6 genes
	total = 0
#siroheme decarboxylase, Fe-coproporphyrin III synthase
	single_ko = ['K22225', 'K22226']
	for i in single_ko:
		if i in ko_match:
			total += 1
#uroporphyrin-III C-methyltransferase
	if ('K00589' in ko_match or 'K02302' in ko_match or 'K02303' in ko_match or 'K02496' in ko_match or 'K13542' in ko_match or 'K13543' in ko_match):
		total += 1
#precorrin-2 dehydrogenase
	if ('K02302' in ko_match or 'K02304' in ko_match):
		total += 1
#sirohydrochlorin ferrochelatase
	if ('K02302' in ko_match or 'K02304' in ko_match or 'K03794' in ko_match):
		total += 1
#heme synthase
	if ('K00435' in ko_match or 'K22227' in ko_match):
		total += 1
	value = float(total)/float(6)
	return {'lower_SIRO': float("%.2f" % (value))}	

def default_viz(genome_df, outfile_name):
	import seaborn as sns
	import matplotlib.pyplot as plt
	sns.set(font_scale=1.2)
	sns.set_style({"savefig.dpi": 200})
	ax = sns.heatmap(genome_df, cmap=plt.cm.YlOrRd, linewidths=2,
		linecolor='k', square=True, xticklabels=True,
		yticklabels=True, cbar_kws={"shrink": 0.1})
	ax.xaxis.tick_top()
	#ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
	plt.xticks(rotation=90)
	plt.yticks(rotation=0)
	# get figure (usually obtained via "fig,ax=plt.subplots()" with matplotlib)
	fig = ax.get_figure()
	# specify dimensions and save
	#xLen = len(genome_df.columns.values.tolist())*20
	#yLen = len(genome_df.index.tolist())*20
	fig.set_size_inches(100, 100)
	fig.savefig(outfile_name, bbox_inches='tight', pad_inches=0.1)

def main():
	import os
	import matplotlib
	matplotlib.use('Agg')
	import argparse
	import pandas as pd
	from scipy.cluster import hierarchy
	from scipy.spatial import distance


	parser = argparse.ArgumentParser(description="Accepts KEGG KOALA\
									text file as input. Produces function\
									list and heat map figure.")
	parser.add_argument('-i', '--input', help="Input KOALA file. See documentation\
						for correct format")
	parser.add_argument('-t', '--tangleopt', help="Number of tree iterations for minimizing tangles in tanglegram", default=1000)
	parser.add_argument('-o', '--output', help="List version of the final heat\
						map figure")
	parser.add_argument('-v', '--vizoption', help="Options: static, interactive, tanglegram")
	parser.add_argument('--newick', help="Required input for tanglegram visualization")
	parser.add_argument("-m", "--myorder", help ="Orders output as specified by	user.", default="None")
	args = parser.parse_args()
	arg_dict = vars(args)

	genome_data = {}

	for line in open(str(arg_dict['input']), "r"):
		line = line.rstrip()
		info = line.split()
		if len(info) > 1:
			if info[0].rsplit("_",1)[0] in genome_data.keys():
				genome_data[info[0].rsplit("_",1)[0]].append(info[1])
			else:
				genome_data[info[0].rsplit("_",1)[0]] = [info[1]]

	function_order = ['C5_PPH', 'C5_CPH', 'C5_SIRO', 'C4_PPH', 'C4_CPH', 'C4_SIRO', 'upper_C5', 'upper_C4', 'Common', 'lower_PPH', 'lower_CPH', 'lower_SIRO']


	filehandle = str(arg_dict['output'])
	out_file = open(filehandle, "w")
	out_file.write('Function'+"\t"+str("\t".join(function_order))+"\n")

	for k in genome_data:
		pathway_data = {}
		pathway_data.update(C5_PPH(genome_data[k]))
		pathway_data.update(C5_CPH(genome_data[k]))
		pathway_data.update(C5_SIRO(genome_data[k]))
		pathway_data.update(C4_PPH(genome_data[k]))
		pathway_data.update(C4_CPH(genome_data[k]))
		pathway_data.update(C4_SIRO(genome_data[k]))
		pathway_data.update(upper_C5(genome_data[k]))
		pathway_data.update(upper_C4(genome_data[k]))
		pathway_data.update(Common(genome_data[k]))
		pathway_data.update(lower_PPH(genome_data[k]))
		pathway_data.update(lower_CPH(genome_data[k]))
		pathway_data.update(lower_SIRO(genome_data[k]))

	#    print k, pathway_data

		out_string = str(k)+"\t"
		out_list = [k]
		for i in function_order:
			out_list.append(pathway_data[i])
		out_string = str(out_list).strip('[]')
		tab_string = ""
		for l in out_string:
			if l == "\'":
				continue
			if l == ",":
				tab_string = tab_string + "\t"
			else:
				tab_string = tab_string + l
		out_file.write(tab_string+"\n")
	out_file.close()


	file_in = open(filehandle, "r")
	genome = pd.read_csv(file_in, index_col=0, sep='\t')
	rearrange = False
	if arg_dict["myorder"] != 'None' and os.path.exists(arg_dict["myorder"]):
		rearrange = True
		leaf_order = []
		for line in open(str(arg_dict["myorder"]), "r"):
			line = line.rstrip("\r\n")
			leaf_order.append(line)
		genome = genome.reindex(leaf_order)

	if arg_dict['vizoption'] == 'static':
		from .KEGG_clustering import hClust_euclidean
		if len(genome.index) >= 2 and not rearrange:
			genome = hClust_euclidean(genome)
		default_viz(genome, os.path.splitext(filehandle)[0] + ".svg")
	if arg_dict['vizoption'] == 'interactive':
		from .Plotly_viz import plotly_viz
		plotly_viz(genome, os.path.splitext(filehandle)[0] + ".html")
	if arg_dict['vizoption'] == 'tanglegram':
		from .MakeTanglegram import make_tanglegram
		if len(genome.index) >= 3:
			make_tanglegram(genome, str(arg_dict['newick']), os.path.splitext(filehandle)[0] + ".tanglegram.svg", int(arg_dict["tangleopt"]))
		else:
			raise ValueError("Tanglegram mode requires three or more genomes")


if __name__ == "__main__":
	main()
