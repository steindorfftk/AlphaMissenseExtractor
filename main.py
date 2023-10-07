from gprofiler import GProfiler

def genome_selector():
	while True:
		reference = input('Choose reference genome: (1) - HG19 ; (2) - HG38 ; (3) - Both > ')
		if reference == '1' or reference == '2' or reference == '3':
			break
		else:
			print("Invalid choice. Please choose a reference genome: ")

	if reference == '1':
		genome = 'HG19'
	elif reference == '2':
		genome = 'HG38'
	elif reference == '3':
		genome = 'Both'
	else:
		print('Warning: no reference genome selected - Using both')
		genome = 'Both'
	return genome

def gene_lister():
	user_gene_list = []
	gene = input('Enter a gene and press enter > ').upper()
	user_gene_list.append(gene)

	phrase = 'Gene list: ' + user_gene_list[0]
	print(phrase)

	while True:
		gene = input('Add another gene or send "E" to continue > ').upper()
		if gene == 'E':
			break
		else:
			user_gene_list.append(gene)
			phrase += ' , '
			phrase += gene
			print(phrase)
	print(phrase+'\n\n\n')
	return user_gene_list

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
	if genome == 'HG19':
		with open('AlphaMissense_Data/AlphaMissense_gene_hg19.tsv','r') as text:
			for line in text:
				for gene, ids in gene_list.items():
					for value in ids:
						if value in line:
							pathogenicity_score[gene] = line.split()[1]
		with open('output/pathogenicity_scores_hg19.tsv','w') as texto:
			texto.write('Gene , Score')
			for gene, value in pathogenicity_score.items():
				texto.write(gene + ' , ' + value + '\n')
	elif genome == 'HG38':
		with open('AlphaMissense_Data/AlphaMissense_gene_hg38.tsv','r') as text:
			for line in text:
				for gene, ids in gene_list.items():
					for value in ids:
						if value in line:
							pathogenicity_score[gene] = line.split()[1]
		with open('output/pathogenicity_scores_hg38.tsv','w') as texto:
			texto.write('Gene , Score')
			for gene, value in pathogenicity_score.items():
				texto.write(gene + ' , ' + value + '\n')
	elif genome == 'Both':
		with open('AlphaMissense_Data/AlphaMissense_gene_hg19.tsv','r') as text:
			for line in text:
				for gene, ids in gene_list.items():
					for value in ids:
						if value in line:
							pathogenicity_score[gene] = line.split()[1]
		with open('output/pathogenicity_scores_hg19.tsv','w') as texto:
			texto.write('Gene , Score')
			for gene, value in pathogenicity_score.items():
				texto.write(gene + ' , ' + value + '\n')
		with open('AlphaMissense_Data/AlphaMissense_gene_hg38.tsv','r') as text:
			for line in text:
				for gene, ids in gene_list.items():
					for value in ids:
						if value in line:
							pathogenicity_score[gene] = line.split()[1]
		with open('output/pathogenicity_scores_hg38.tsv','w') as texto:
			texto.write('Gene , Score')
			for gene, value in pathogenicity_score.items():
				texto.write(gene + ' , ' + value + '\n')
					


def missense_data_extract(enst_list):
	missense_data = {}
	
	for gene in enst_list.keys():
		missense_data[gene] = [] 
	if genome == 'HG19':
		with open('AlphaMissense_Data/AlphaMissense_hg19.tsv','r') as text:
			for line in text:
				for gene, ids in enst_list.items():
					for tag in ids:
						if tag in line:
							missense_data[gene].append(line)
		
		for gene, data in missense_data.items():
			file_path = 'output/' + gene + '_missense_hg19.tsv'					
			with open(file_path,'w') as text:
				text.write('CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class \n')
				for value in data:
					text.write(value)
	elif genome == 'HG38':
		with open('AlphaMissense_Data/AlphaMissense_hg38.tsv','r') as text:
			for line in text:
				for gene, ids in enst_list.items():
					for tag in ids:
						if tag in line:
							missense_data[gene].append(line)
		
		for gene, data in missense_data.items():
			file_path = 'output/' + gene + '_missense_hg38.tsv'					
			with open(file_path,'w') as text:
				text.write('CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class \n')
				for value in data:
					text.write(value)
	elif genome == 'Both':
		with open('AlphaMissense_Data/AlphaMissense_hg19.tsv','r') as text:
			for line in text:
				for gene, ids in enst_list.items():
					for tag in ids:
						if tag in line:
							missense_data[gene].append(line)
		
		for gene, data in missense_data.items():
			file_path = 'output/' + gene + '_missense_hg19.tsv'					
			with open(file_path,'w') as text:
				text.write('CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class \n')
				for value in data:
					text.write(value)
		with open('AlphaMissense_Data/AlphaMissense_hg38.tsv','r') as text:
			for line in text:
				for gene, ids in enst_list.items():
					for tag in ids:
						if tag in line:
							missense_data[gene].append(line)
		
		for gene, data in missense_data.items():
			file_path = 'output/' + gene + '_missense_hg38.tsv'					
			with open(file_path,'w') as text:
				text.write('CHROM	POS	REF	ALT	genome	uniprot_id	transcript_id	protein_variant	am_pathogenicity	am_class \n')
				for value in data:
					text.write(value)
	

def aminoacid_subst_data_extract(uniprot_list):
	aminoacid_data = {}
	
	for gene in uniprot_list.keys():
		aminoacid_data[gene] = [] 
	with open('AlphaMissense_Data/AlphaMissense_aa_substitutions.tsv','r') as text:
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
			
#Choose genome
genome = genome_selector()

#Select gene list
input_gene_list = gene_lister()
	

#Identificators
print("Searching for genes' ENSTs...")
enst_ids = ENST_translator(input_gene_list)
print('Done')
print("Searching for genes' UNIPROT IDs...")
uniprot_ids = UNIPROT_translator(input_gene_list)
print('Done')

#a) Pathogenicity per gene
print("Extracting the pathogenicity mean for list of genes...")
gene_overview(enst_ids)
print('Done')

#b) Nucleotide missense variants data extracting
print("Extracting the predictions for missense variants in the cannonical transcript of each gene...")
missense_data_extract(enst_ids)
print('Done')

#c) Aminoacid substitutions 
print("Extracting the predictions for aminoacid substitutions in the cannonical protein sequence of each gene...")
aminoacid_subst_data_extract(uniprot_ids)
print('Done')














