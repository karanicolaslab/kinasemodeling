import os, glob, pandas as pd, sys, csv

OMEGA = "/usr/local/bin/omega2"
ROCS = "/usr/local/bin/rocs"

def omega_conf(smi,maxconfs=2000):
	smi_prefix = os.path.splitext(os.path.basename(smi))[0]
	os.system ('{0} -in {1} -out {2}_omega.sdf -prefix {2}_omega -warts true -maxconfs {3}'.format(OMEGA, smi, smi_prefix, maxconfs))

def rocs_alignment(conformer, template_database, rocs_maxconfs_output=100):
	sdf_prefix = os.path.basename(os.path.splitext(conformer)[0]).split('_')[0]
	for template in template_database:
		template_id = "_".join(os.path.basename(template).split("_")[0:3])
		os.system ('{0} -dbase {1} -query {2} -prefix {3}_{4}_rocs -oformat sdf -maxconfs 30 -outputquery false -qconflabel title'.format(ROCS, conformer, template, sdf_prefix, template_id))

def combine_report_files(report_file, crystal_structure):
	dataframeList = []
	for report in report_file:
		target_template_file_name = "_".join(os.path.splitext(os.path.basename(report))[0].split("_")[0:5])
		lig_id_three_letter = target_template_file_name.split("_")[0]
		read_rpt = pd.read_csv(report, sep='\t', dtype=str)
		if 'Unnamed: 16' in read_rpt:
			del read_rpt['Unnamed: 16']
		shapequery_changed_table = read_rpt.replace('untitled-query-1', target_template_file_name)
		shapequery_changed_table.to_csv(report, sep='\t', index=None)
		dataframeList.append(shapequery_changed_table)
	single_rpt_csv_file = pd.concat(dataframeList, axis=0)
	crystal_structure_pdb_id = os.path.splitext(os.path.basename(crystal_structure))[0]
	crystal_structure_pdb_id_with_lig = "{0}_{1}_rocs".format(lig_id_three_letter, crystal_structure_pdb_id)
	crystal_temp_removed_boolean = single_rpt_csv_file['ShapeQuery'] != crystal_structure_pdb_id_with_lig
	crystal_temp_removed_table = single_rpt_csv_file[crystal_temp_removed_boolean]
	tanimoto_sorted_table = crystal_temp_removed_table.sort_values(by=['TanimotoCombo'], ascending=False)
	tanimoto_sorted_table.to_csv('single_report_file_sorted.csv', index=None)
	top_100_hit_table = tanimoto_sorted_table.iloc[:100,:]
	top_100_conf_temp_names = top_100_hit_table['Name'].astype(str) + '_' + top_100_hit_table['ShapeQuery'].astype(str)
	top_100_conf_temp_names.to_csv('top_100.txt', sep='\t', index=None)

def sep_hits_from_rocs_sdf_file(top_100_hits_txt_path):
	with open (top_100_hits_txt_path) as open_file:
		dict_sdf_file_name = {}
		for line in open_file:
			if len(line.split('_')) == 8:
				print(line)
				conf_name = "_".join(line.strip().split('_')[0:3])
				temp_name = "_".join(line.strip().split('_')[3:])
				if conf_name not in dict_sdf_file_name:
					dict_sdf_file_name[conf_name] = [temp_name]
				else:
					dict_sdf_file_name[conf_name].append(temp_name)
			elif len(line.split('_')) == 7:
				conf_name = "_".join(line.strip().split('_')[0:2])
				temp_name = "_".join(line.strip().split('_')[2:])
				if conf_name not in dict_sdf_file_name:
					dict_sdf_file_name[conf_name] = [temp_name]
				else:
					dict_sdf_file_name[conf_name].append(temp_name)
		for key, list_value in dict_sdf_file_name.items():
			for single_value in list_value:
				print('sed -n /{0}$/,/\$\$\$\$/p {1}_hits_1.sdf > {0}_{1}_hits.sdf'.format(key, single_value))
				os.system('sed -n /{0}$/,/\$\$\$\$/p {1}_hits_1.sdf > {0}_{1}_hits.sdf'.format(key, single_value))

def sdftoparams(top_hits_sdf_path):
	for file in top_hits_sdf_path:
		out_put_file_name = os.path.splitext(os.path.basename(file))[0]
		print('/usr/bin/python2.7 /Users/kirubapalani/Rosetta/main/source/scripts/python/public/molfile_to_params.py {0}/{1} -p {0}/{2}'.format(os.getcwd(), file, out_put_file_name))
		# os.system('/usr/bin/python2.7 /Users/kirubapalani/Rosetta/main/source/scripts/python/public/molfile_to_params.py {0} -p {1}'.format(file, out_put_file_name))

# omega_conf(glob.glob("*.smi")[0], maxconfs=1000)
# os.system("mkdir OMEGA")
# os.system("mv *omega* OMEGA/")
# rocs_alignment(glob.glob("OMEGA/*.sdf")[0], glob.glob("/Users/kirubapalani/kinase_project/dunbrack_kinase_type_I_library/active_ligands_seperated_from_pdb_one_chain/*"), rocs_maxconfs_output=30)
# os.system("mkdir ROCS")
# os.system("for i in *rocs*; do mv $i ROCS/; done")
# combine_report_files(glob.glob("ROCS/*.rpt"), glob.glob("????_?_???.pdb")[0])
# os.system("mv single_report_file_sorted.csv top_100.txt ROCS/")
# os.system("mkdir top_100_conf")
# os.chdir("ROCS")
# sep_hits_from_rocs_sdf_file("top_100.txt")
# os.system("mv *hits.sdf ../top_100_conf")
# os.chdir("../")
sdftoparams(glob.glob("top_100_conf/*.sdf"))
# os.system("mkdir sdf2params")
# os.system("mv *_0001.pdb *.params sdf2params")
