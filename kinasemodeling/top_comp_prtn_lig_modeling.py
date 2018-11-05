import csv
import glob
import os
from openeye import oechem

def sep_lig_from_first_mini_complex (mini_models_path):
	for path in mini_models_path:
		lig_name = '_'.join(os.path.splitext(os.path.basename(path))[0].split('_')[5:-1])
		print('grep \"HETATM\" {0} > {1}.pdb'.format(path, lig_name))
		os.system('grep \"HETATM\" {0} > {1}.pdb'.format(path, lig_name))

def conv_mini_lig_pdb_to_sdf(mini_lig_pdbs):
	for pdb in mini_lig_pdbs:
		print( pdb, 'pdb_to_sdf using oechem')
		lig_name = pdb.split('.')[0]
		ifs = oechem.oemolistream()
		ofs = oechem.oemolostream()
		if ifs.open(pdb):
			if ofs.open(lig_name + ".sdf"):
				for mol in ifs.GetOEGraphMols():
					oechem.OEWriteMolecule(ofs, mol)

def sdftoparams_top_10(top_hits_sdf_path):	
	for file in top_hits_sdf_path:
		out_put_file_name = os.path.splitext(os.path.basename(file))[0]
		os.system('/usr/bin/python2.7 /Users/kirubapalani/Rosetta/main/source/scripts/python/public/molfile_to_params.py {0} -p {1}'.format(file, out_put_file_name))

def protein_ligand_concatenation_top_10_comp_models(template_hits, target_fasta_file, ligands):
	top_10_comp_models = []
	with open(template_hits) as read_CSV:
		reader = csv.DictReader(read_CSV)
		for row in reader:
			target_pdb_path = os.path.dirname(template_hits)
			target_pdb_name = os.path.basename(target_fasta_file).split('.')[0]
			model = '{0}_{1}.pdb'.format(target_pdb_name, row['template'])
			top_10_comp_models.append(model)
	
	for model in top_10_comp_models:
		for lig in ligands:
			protein_model = os.path.join(target_pdb_path, model)
			basename_ligand = os.path.splitext(os.path.basename(lig))[0]
			complex_protein_ligand = '{0}_{1}.pdb'.format(model.split('.')[0], basename_ligand)
			print('cat {0} {1} > {2}'.format(protein_model, lig, complex_protein_ligand))
			os.system('cat {0} {1} > {2}'.format(protein_model, lig, complex_protein_ligand))

sep_lig_from_first_mini_complex(glob.glob('top_10_mini_models_from_first_comp_model/*.pdb'))
os.system('mkdir lig_from_top_10_min_models_from_first_comp_model')
os.system("mv *rocs_hits* lig_from_top_10_min_models_from_first_comp_model/")

os.chdir('lig_from_top_10_min_models_from_first_comp_model')
mini_lig_pdb_path = glob.glob('*.pdb')
conv_mini_lig_pdb_to_sdf(mini_lig_pdb_path)
os.chdir("../")

sdftoparams_top_10(glob.glob('lig_from_top_10_min_models_from_first_comp_model/*.sdf'))
os.system('mkdir lig_from_top_10_min_models_from_first_comp_model_sdf_to_params')
os.system("mv *rocs_hits* lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/")

protein_ligand_concatenation_top_10_comp_models("protein_comp_modeling/template_hits.csv", glob.glob("*.fasta")[0], glob.glob("lig_from_top_10_min_models_from_first_comp_model_sdf_to_params/*.pdb"))
os.system('mkdir protein_ligand_complex_top_10_comp_models')
os.system("mv *hits_0001.pdb protein_ligand_complex_top_10_comp_models")


