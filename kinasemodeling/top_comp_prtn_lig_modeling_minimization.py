import csv
import glob
import os

def minimization(protein_ligand_path, lig_params_path):
	dict_prtn_lig_with_params = {}
	for p_l_path in protein_ligand_path:
		params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[4:-1])
		params_file_with_ext = "{0}.params".format(params_file)
		dict_prtn_lig_with_params[p_l_path] = params_file_with_ext

	for key, value in dict_prtn_lig_with_params.items():
		protein_name = os.path.basename(key)
		protein_ligand_without_ext = protein_name.split('.')[0]
		for lig_params in lig_params_path:
			lig_params_basename = os.path.basename(lig_params)
			if lig_params_basename == value:
				# print('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease -database /mnt/shared_applications/Rosetta_JKlab/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o mini_{3}.log'.format(key, lig_params, protein_name, protein_ligand_without_ext))
				os.system('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/minimize_ppi.default.linuxgccrelease -database /mnt/shared_applications/Rosetta_JKlab/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o mini_{3}.log'.format(key, lig_params, protein_name, protein_ligand_without_ext))


				# print('/Users/kirubapalani/Rosetta/main/source/bin/minimize_ppi.default.macosclangrelease -database /Users/kirubapalani/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1} > mini_{2}.log'.format(key, lig_params, protein_ligand_without_ext))
				# os.system('/Users/kirubapalani/Rosetta/main/source/bin/minimize_ppi.default.macosclangrelease -database /Users/kirubapalani/Rosetta/main/database -ignore_unrecognized_res -s {0} -extra_res_fa {1} > mini_{2}.log'.format(key, lig_params, protein_ligand_without_ext))

def minimization_score_only(mini_protein_ligand_path, lig_params_path):
	dict_prtn_lig_with_params = {}
	for p_l_path in mini_protein_ligand_path:
		params_file = "_".join(os.path.splitext(os.path.basename(p_l_path))[0].split("_")[5:-1])
		params_file_with_ext = "{0}.params".format(params_file)
		dict_prtn_lig_with_params[p_l_path] = params_file_with_ext

	for key, value in dict_prtn_lig_with_params.items():
		protein_name = os.path.basename(key)
		protein_ligand_without_ext = protein_name.split('.')[0]
		for lig_params in lig_params_path:
			lig_params_basename = os.path.basename(lig_params)
			if lig_params_basename == value:
				# print('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o {3}.score'.format(key, lig_params, protein_name, protein_ligand_without_ext))
				os.system('echo \"/mnt/shared_applications/Rosetta_JKlab/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.linuxgccrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1}\" | qsub -N {2} -d $PWD -e error.log -o {3}.score'.format(key, lig_params, protein_name, protein_ligand_without_ext))

				# print('/Users/kirubapalani/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.macosclangrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1} > {2}.score'.format(key, lig_params, protein_ligand_without_ext))
				# os.system('/Users/kirubapalani/Rosetta/main/source/bin/bou-min-ubo-nrg-jump.default.macosclangrelease -ignore_unrecognized_res -s {0} -extra_res_fa {1} > {2}.score'.format(key, lig_params, protein_ligand_without_ext))


minimization(glob.glob('protein_ligand_complex_top_10_comp_models/*.pdb'), glob.glob('lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.params'))
# os.system("mkdir mini_protein_ligand_complex_of_all_models")
# os.system("for i in mini_*.pdb; do echo $i; mv $i mini_protein_ligand_complex_of_all_models; done")
# os.system("for i in *.log; do echo $i; mv $i mini_protein_ligand_complex_of_all_models;done")
# minimization_score_only(glob.glob('mini_protein_ligand_complex_of_all_models/*.pdb'), glob.glob('lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.params'))
# os.system("for i in mini_*.score; do echo $i; mv $i mini_protein_ligand_complex_of_all_models; done")