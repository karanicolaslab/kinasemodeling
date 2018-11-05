import glob, os

def minimization(protein_ligand_path, lig_params_path):
	dict_prtn_lig_with_params = {}
	for p_l_path in protein_ligand_path:
		params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[4:-1])
		params_file_with_ext = "{0}.params".format(params_file)
		dict_prtn_lig_with_params[params_file_with_ext] = p_l_path
	for params in lig_params_path:
		params_basename = os.path.basename(params)
		if params_basename in dict_prtn_lig_with_params:
			protein_ligand = os.path.basename(dict_prtn_lig_with_params[params_basename])
			protein_ligand_without_ext = os.path.splitext(os.path.basename(dict_prtn_lig_with_params[params_basename]))[0]
			# print('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease -database mnt/shared_applications/Rosetta_JKlab/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o mini_{3}.log'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand, protein_ligand_without_ext))
			os.system('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease -database /mnt/shared_applications/Rosetta_JKlab/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o mini_{3}.log'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand, protein_ligand_without_ext))
			
			# print('/Users/kirubapalani/Rosetta/main/source/bin/minimize_ppi.default.macosclangrelease -database /Users/kirubapalani/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1} > mini_{2}.log'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))
			# os.system('/Users/kirubapalani/Rosetta/main/source/bin/minimize_ppi.default.macosclangrelease -database /Users/kirubapalani/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1} > mini_{2}.log'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))

def minimization_score_only(mini_protein_ligand_path, lig_params_path):
	dict_prtn_lig_with_params = {}
	for mini_p_l_path in mini_protein_ligand_path:
		params_file = "_".join(os.path.splitext(os.path.basename(mini_p_l_path))[0].split("_")[5:-1])
		params_file_with_ext = "{0}.params".format(params_file)
		dict_prtn_lig_with_params[params_file_with_ext] = mini_p_l_path
	for params in lig_params_path:
		params_basename = os.path.basename(params)
		if params_basename in dict_prtn_lig_with_params:
			protein_ligand = os.path.basename(dict_prtn_lig_with_params[params_basename])
			protein_ligand_without_ext = os.path.splitext(os.path.basename(dict_prtn_lig_with_params[params_basename]))[0]
			# print('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o {3}.score'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand, protein_ligand_without_ext))
			os.system('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o {3}.score'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand, protein_ligand_without_ext))

			# print('/Users/kirubapalani/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.macosclangrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1} > {2}.score'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))
			# os.system('/Users/kirubapalani/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.macosclangrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1} > {2}.score'.format(dict_prtn_lig_with_params[params_basename], params, protein_ligand_without_ext))

minimization(glob.glob('protein_ligand_complex_top_1_comp_model/*.pdb'), glob.glob('sdf2params/*.params'))
# os.system("mkdir mini_protein_ligand_complex_top_1_comp_model")
# os.system("for i in mini_*.pdb; do echo $i; mv $i mini_protein_ligand_complex_top_1_comp_model; done")
# os.system("for i in *.log; do echo $i; mv $i mini_protein_ligand_complex_top_1_comp_model;done ")
# minimization_score_only(glob.glob('mini_protein_ligand_complex_top_1_comp_model/*.pdb'), glob.glob('sdf2params/*.params'))
# os.system("for i in mini_*.score; do echo $i; mv $i mini_protein_ligand_complex_top_1_comp_model; done")
