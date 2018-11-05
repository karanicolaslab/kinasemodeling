import os
import glob
import pandas as pd
import csv
from collections import defaultdict
import pyrosetta
pyrosetta.init()

def emboss_needle_search(target_seq_path, template_seq_path):
	
	for template_seq in template_seq_path:
		target_seq_id = os.path.basename(target_seq_path).split('.')[0]
		template_seq_id = os.path.basename(template_seq).split('.')[0]
		print('needle -sid1 {0} -asequence {1} -sid2 {2} -bsequence {3} -gapopen 10.0 -gapextend 0.5 -aformat3 markx3 -outfile {0}_{2}.needle'.format(target_seq_id, target_seq_path, template_seq_id, template_seq))
		os.system('needle -sid1 {0} -asequence {1} -sid2 {2} -bsequence {3} -gapopen 10.0 -gapextend 0.5 -aformat3 markx3 -outfile {0}_{2}.needle'.format(target_seq_id, target_seq_path, template_seq_id, template_seq))
		
def select_top_hits_from_emboss_and_rocs_pdb(emboss_align_file_path, rocs_align_file_path, target_seq_path):

	emboss_result = pd.DataFrame(columns = ('query', 'template', 'length', 'identity', 'similarity', 'gaps', 'score'))
	emboss_ind_title = 1
	for emboss_alignment_file in emboss_align_file_path:
		with open(emboss_alignment_file, 'r') as readFILE:
			for line in readFILE:
				line = line.strip()
				if '# 1:' in line:
					query = line.split()[-1]
				elif '# 2:' in line:
					template = line.split()[-1]
				elif '# Length:' in line:
					length = line.split()[-1]
				elif '# Identity:' in line:
					identity = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
				elif '# Similarity:' in line:
					similarity = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
				elif '# Gaps:' in line:
					gaps = line[24:].replace('(', '').replace(')', '').replace('%', '')
				elif '# Score:' in line:
					score = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
			# print (query, template, length, identity, similarity, gaps, score)
			emboss_result.loc[emboss_ind_title] = query, template, float(length), float(identity), float(similarity), float(gaps), float(score)
			emboss_ind_title += 1
	emboss_result['score_rank'] = emboss_result['score'].rank(method = 'dense', ascending=0) 
	rank_emboss_hits_dict = {}
	for index, row in emboss_result.iterrows():
		rank_emboss_hits_dict[row['template']] = row['score_rank']

	rocs_result = pd.DataFrame(columns = ('rocs_pdb', 'TanimotoCombo'))
	rocs_ind_title = 1
	with open(rocs_align_file_path) as CSVfileREAD:
		reader = csv.DictReader(CSVfileREAD)
		for row in reader:
			row_name_and_shapequery_conc = '{0}_{1}'.format(row['Name'], row['ShapeQuery'])
			pdb_with_chain = '_'.join(row['ShapeQuery'].split('_')[1:3]).upper()
			rocs_result.loc[rocs_ind_title] = row_name_and_shapequery_conc, float(row['TanimotoCombo'])
			rocs_ind_title += 1
	rocs_result['TanimotoCombo_rank'] = rocs_result['TanimotoCombo'].rank(method = 'dense', ascending=0)
	rank_rocs_hits_dict = {}
	for index, row in rocs_result.iterrows():
		rank_rocs_hits_dict[row['rocs_pdb']] = row['TanimotoCombo_rank']


	final_hit_list_df = pd.DataFrame(columns = ('template', 'score_rank', 'rocs_pdb', 'TanimotoCombo_rank', 'final_rank'))
	ind_title = 1
	for template, score_rank in rank_emboss_hits_dict.items():
		for rocs_pdb, TanimotoCombo_rank in rank_rocs_hits_dict.items():
			if len(rocs_pdb.split('_')) == 8:		
				pdb_with_chain = '_'.join(rocs_pdb.split('_')[4:6])
				target_seq_pdb = os.path.basename(target_seq_path).split('.')[0]
				if pdb_with_chain == template and target_seq_pdb != template:
					final_hit_list_df.loc[ind_title] = template, score_rank, pdb_with_chain, TanimotoCombo_rank, score_rank*TanimotoCombo_rank
					ind_title+=1
			elif len(rocs_pdb.split('_')) == 7:
				pdb_with_chain = '_'.join(rocs_pdb.split('_')[3:5])
				target_seq_pdb = os.path.basename(target_seq_path).split('.')[0]
				if pdb_with_chain == template and target_seq_pdb != template:
					final_hit_list_df.loc[ind_title] = template, score_rank, pdb_with_chain, TanimotoCombo_rank, score_rank*TanimotoCombo_rank
					ind_title+=1
	final_hit_list_df = final_hit_list_df.drop_duplicates(subset='rocs_pdb', keep="first")
	final_hit_list_df = final_hit_list_df.sort_values('final_rank', ascending=True)
	final_hit_list_df.iloc[:10].to_csv('template_hits.csv', sep=',', index=False)

def modeling(template_hits, template_pdb_path, alignment_file_path, target_seq_path):
	
	top_hit_template_file_path = []
	with open(template_hits, newline='') as csvFile:
		reader = csv.DictReader(csvFile)
		for row in reader:
			target_seq = os.path.basename(target_seq_path).split('.')[0]
			alignment_file_name = '{}_{}.needle'.format(target_seq, row['template'])
			alignment_file_path_with_name = os.path.join(alignment_file_path, alignment_file_name)
			top_hit_template_file_path.append(alignment_file_path_with_name)

	aligned_seq = defaultdict(list)
	for path in top_hit_template_file_path:
		target_template_file_name = os.path.splitext(os.path.basename(path))[0]
		target_name_fasta_format = '>{} ..'.format('_'.join(target_template_file_name.split('_')[0:2]))
		template_name_fasta_format = '>{} ..'.format('_'.join(target_template_file_name.split('_')[2:]))
		target_aligned_seq = ''
		template_aligned_seq = ''
		with open (path, 'r') as readFile:
			parse = False
			parse2 = False
			for line in readFile:
				line = line.strip()
				if not parse:
					if line.startswith(target_name_fasta_format):
						parse = True
				elif line.startswith(template_name_fasta_format):
					parse = False
				else:
					target_aligned_seq+=line

				if not parse2:
					if line.startswith(template_name_fasta_format):
						parse2 = True
				elif line.startswith('#'):
					parse2 = False
				else:
					template_aligned_seq += line

		aligned_seq[target_template_file_name].append(target_aligned_seq)
		aligned_seq[target_template_file_name].append(template_aligned_seq)
		
	target_seq_for_modeling = {}
	for name, alignment_file in aligned_seq.items():
		top_hits_alignment = '{}\n{}\n{}\n\n'.format(name, alignment_file[0], alignment_file[1])
		with open('top_hits_alignment.txt', 'a') as writeFile:
			writeFile.write(top_hits_alignment)
		target_seq_based_on_temp_pdb = ''
		for i in range(len(alignment_file[0])):
			if not alignment_file[1][i] == '-':
				target_seq_based_on_temp_pdb += alignment_file[0][i]
		target_seq_for_modeling[name]=target_seq_based_on_temp_pdb

	final_target_template_for_modeling = {}
	for target_template, target_final_seq in target_seq_for_modeling.items():
		template_name = '_'.join(target_template.split('_')[2:])
		temp_list_dir = os.listdir(template_pdb_path)
		for template_hit in temp_list_dir:
			if template_name in template_hit:
				final_target_template_for_modeling[template_hit] = target_final_seq

	for template_pdb, target_seq in final_target_template_for_modeling.items():
		output_model_name = '{}_{}.pdb'.format(target_seq_path.split('.')[0], '_'.join(template_pdb.split('_')[0:2]))
		join_apo_dir_path = os.path.join(template_pdb_path, template_pdb)
		pose = pyrosetta.pose_from_file(join_apo_dir_path)
		assert(pose.size() == len(target_seq))
		scorefxn = pyrosetta.get_fa_scorefxn()
		for i in range(len(target_seq)):
			seqpos = i + 1
			name1 = target_seq[i]
			if (name1 == "-"):
				continue
			pyrosetta.rosetta.protocols.toolbox.pose_manipulation.repack_this_residue(seqpos, pose, scorefxn, True, name1)
		pose.dump_pdb(output_model_name)

def protein_ligand_concatenation (template_hits, target_fasta_file, ligands):
	template_id_and_score = {}
	with open(template_hits) as read_CSV:
		reader = csv.DictReader(read_CSV)
		for row in reader:
			template_id_and_score[row['template']]=float(row['final_rank'])
	top_first_hit_pdb = min(template_id_and_score, key=template_id_and_score.get)
	target_pdb_path = os.path.dirname(template_hits)
	target_pdb_name = os.path.basename(target_fasta_file).split('.')[0]
	top_hit_model = '{0}_{1}.pdb'.format(target_pdb_name, top_first_hit_pdb)
	
	for lig in ligands:
		protein_model = os.path.join(target_pdb_path, top_hit_model)
		basename_ligand = os.path.splitext(os.path.basename(lig))[0]
		complex_protein_ligand = '{0}_{1}.pdb'.format(top_hit_model.split('.')[0], basename_ligand)
		print('cat {0} {1} > {2}'.format(protein_model, lig, complex_protein_ligand))
		os.system('cat {0} {1} > {2}'.format(protein_model, lig, complex_protein_ligand))

emboss_needle_search(glob.glob('*.fasta')[0], glob.glob('/Users/kirubapalani/kinase_project/dunbrack_kinase_type_I_library/modeling_studies/modeling_using_protein_comp_structures/fasta_files_generated_using_pose_sequence_function/*.fasta'))
os.system('mkdir protein_comp_modeling')
os.system('mkdir protein_seq_alignment_files')
os.system('for i in *.needle; do mv $i protein_seq_alignment_files; done')
os.system('mv protein_seq_alignment_files protein_comp_modeling;')

select_top_hits_from_emboss_and_rocs_pdb(glob.glob('protein_comp_modeling/protein_seq_alignment_files/*.needle'), 'ROCS/single_report_file_sorted.csv', glob.glob('*.fasta')[0])
os.system('mv template_hits.csv protein_comp_modeling')

modeling('protein_comp_modeling/template_hits.csv', '/Users/kirubapalani/kinase_project/dunbrack_kinase_type_I_library/modeling_studies/modeling_using_protein_comp_structures/apo_pdbs', 'protein_comp_modeling/protein_seq_alignment_files/', glob.glob('*.fasta')[0])
os.system('mv ????_?_????_?.pdb protein_comp_modeling')
os.system('mv top_hits_alignment.txt protein_comp_modeling')

os.system("mkdir protein_ligand_complex_top_1_comp_model")
protein_ligand_concatenation("protein_comp_modeling/template_hits.csv", glob.glob("*.fasta")[0], glob.glob("sdf2params/*.pdb"))
os.system("mv *hits_0001.pdb protein_ligand_complex_top_1_comp_model")
