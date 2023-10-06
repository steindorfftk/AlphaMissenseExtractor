from gprofiler import GProfiler

def ENST_translator(gene_list):
	gp = GProfiler(return_dataframe=False)
	ENST_dic = {}
	for gene in gene_list:
		if gene not in ENST_dic.keys():
			ENST_dic[gene] = []
		tags = gp.convert(organism='hsapiens',
			query=gene,
			target_namespace='ENST')
		for value in tags:
			ENST_dic[gene].append(value['converted'][4:])
	return ENST_dic				
			
def UNIPROT_translator(gene_list):
	gp = GProfiler(return_dataframe=False)
	UNIPROT_dic = {}
	for gene in gene_list:
		if gene not in UNIPROT_dic.keys():
			UNIPROT_dic[gene] = []
		tags = gp.convert(organism='hsapiens',
			query=gene,
			target_namespace='UNIPROT_GN_ACC')
		for value in tags:
			UNIPROT_dic[gene].append(value['converted'])
	return UNIPROT_dic


def gene_overview(gene_list):
	pathogenicity_score = {}
	with open('AlphaMissense_Data/AlphaMissense_gene_hg38.tsv','r') as text:
		for line in text:
			for gene, ids in gene_list.items():
				for value in ids:
					if value in line:
						pathogenicity_score[gene] = line.split()[1]
	with open('output/pathogenicity_scores.tsv','w') as texto:
		texto.write('Gene , Score')
		for gene, value in pathogenicity_score.items():
			texto.write(gene + ' , ' + value + '\n')


def missense_data_extract(enst_list):
	missense_data = {}
	
	for gene in enst_list.keys():
		missense_data[gene] = [] 
	
	with open('AlphaMissense_Data/AlphaMissense_hg38.tsv','r') as text:
		for line in text:
			for gene, ids in enst_list.items():
				for tag in ids:
					if tag in line:
						missense_data[gene].append(line)
	
	for gene, data in missense_data.items():
		file_path = 'output/' + gene + '_missense.tsv'					
		with open(file_path,'w') as text:
			text.write('CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class \n')
			for value in data:
				text.write(value)

def aminoacid_subst_data_extract(uniprot_list):
	aminoacid_data = {}
	
	for gene in uniprot_list.keys():
		aminoacid_data[gene] = [] 
	
	with open('AlphaMissense_Data/AlphaMissense_aa_substitutions_hg38.tsv','r') as text:
		for line in text:
			for gene, ids in uniprot_list.items():
				for tag in ids:
					if tag in line:
						aminoacid_data[gene].append(line)
	
	for gene, data in aminoacid_data.items():
		file_path = 'output/' + gene + '_substitutions.tsv'				
		with open(file_path,'w') as text:
			text.write('uniprot_id	protein_variant	am_pathogenicity	am_class \n')
			for value in data:
				text.write(value)


			
#Input
input_gene_list = ['IDUA']

#Identificators
enst_ids = ENST_translator(input_gene_list)
uniprot_ids = UNIPROT_translator(input_gene_list)

#Pathogenicity per gene
gene_overview(enst_ids)


#Nucleotide missense variants data extracting
missense_data_extract(enst_ids)

#Nucleotide missense variants exporting (AlphaMissense_hg38.tsv)
aminoacid_subst_data_extract(uniprot_ids)

