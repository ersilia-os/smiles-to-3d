import os
import csv
import argparse
import subprocess
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit import RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')  

root = os.path.dirname(os.path.abspath(__file__))


def generate_conformers(in_file, out_file, max_time=60, num_confs=10, quiet=False):
    tmp_dir = tempfile.mkdtemp(prefix='smi3d_')
    input_file = os.path.join(tmp_dir, 'input.sdf')
    output_file = os.path.join(tmp_dir, 'output.sdf')
    output_file2 = os.path.join(tmp_dir, 'output2.sdf')

    with open(in_file, 'r') as f:
        smiles_list = []
        id_list = []
        reader = csv.reader(f)
        headers = next(reader)
        
        smiles_col = None
        id_col = None
        for idx, header in enumerate(headers):
            if 'smiles' in header.lower():
                smiles_col = idx
            elif 'id' in header.lower():
                id_col = idx
        
        if smiles_col is None:
            raise ValueError("No column with 'smiles' found in the CSV file.")
        if id_col is None:
            raise ValueError("No column with 'id' found in the CSV file.")
        
        for row in reader:
            smiles_list.append(row[smiles_col])
            id_list.append(row[id_col])

    idxs = []
    mols = []
    for i, (smi, mol_id) in enumerate(zip(smiles_list, id_list)):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.Compute2DCoords(mol)
        if mol is None:
            continue
        idxs += [i]
        mol.SetProp("_MoleculeID", str(mol_id))
        mol.SetProp("_InputIndex", str(i))
        mol.SetProp("_InputSMILES", str(smi))
        mol.SetProp("_InputInChIKey", str(inchi.InchiToInchiKey(inchi.MolToInchi(Chem.MolFromSmiles(smi)))))
        mol.SetProp("_SMILES", str(Chem.MolToSmiles(mol)))
        mol.SetProp("_InChIKey", str(inchi.InchiToInchiKey(inchi.MolToInchi(mol))))
        mols += [mol]

    writer = Chem.SDWriter(input_file)
    for m in mols:
        writer.write(m)
    writer.close()

    if num_confs == 1:
        cmd = "python {0}/tools/gen_3d_structs.py -i {1} -o {2} -t {3}".format(root, input_file, output_file, max_time)
    else:
        cmd = "python {0}/tools/gen_confs.py -i {1} -o {2} -t {3} -n {4}".format(root, input_file, output_file, max_time, num_confs)

    subprocess.Popen(cmd, shell=True).wait()

    # Now, let's read the output SDF file and add _Name property to each molecule conformer
    sdf_supplier = Chem.SDMolSupplier(output_file, removeHs=False)
    conf_counts = {}
    mols_with_names =[]
    for mol in sdf_supplier:
        if mol is None:
            continue
        mol_id = mol.GetProp("_MoleculeID")
        if mol_id not in conf_counts:
            conf_counts[mol_id] = 0
        conf_label = f"{mol_id}_conf{conf_counts[mol_id]}"
        mol.SetProp("_Name", str(conf_label))
        conf_counts[mol_id] += 1
        mols_with_names.append(mol)
    # Rewrite the modified molecules with _Name property to the output SDF file
    sdf_writer = Chem.SDWriter(output_file2)
    for mol in mols_with_names:
        sdf_writer.write(mol)
    sdf_writer.close()

    shutil.copyfile(output_file2, out_file)
    shutil.rmtree(tmp_dir)


def parse_args():
    parser = argparse.ArgumentParser(description='Generates conformers for a given list of SMILES strings.')

    parser.add_argument('-i',
                        dest='in_file',
                        required=True,
                        metavar='<file>',
                        help='Molecule input file in CSV format. A column named "smiles" is required.')
    parser.add_argument('-o',
                        dest='out_file',
                        required=True,
                        metavar='<file>',
                        help='Output file in SDF format.')
    parser.add_argument('-t',
                        dest='max_time',
                        required=False,
                        metavar='<int>',
                        type=int,
                        default=60,
                        help='Max. allowed molecule processing time (default: 60 sec)')
    parser.add_argument('-q',
                        dest='quiet',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Disable progress output (default: false)')
    parser.add_argument('-n',
                        dest='num_confs',
                        required=False,
                        metavar='<int>',
                        type=int,
                        default=10,
                        help='Number of conformers to generate (default: 10)')
    
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    generate_conformers(args.in_file, args.out_file, args.max_time, args.num_confs, args.quiet)