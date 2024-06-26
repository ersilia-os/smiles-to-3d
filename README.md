# SMILES to 3D conformers

Convert a list of SMILES molecules to 3D conformers within reasonable time using the [CPDKit package](https://cdpkit.org/v1.1.1/index.html)

# Installation

We recommend using a conda environment:

```bash
conda create -n smiles3d python=3.11
conda activate smiles3d
pip install rdkit==2023.9.6
pip install cdpkit==1.1.1
```
or directly install in your preferred environment:
```bash
pip install git+https://github.com/ersilia-os/smiles-to-3d.git
```

# Usage

```bash
cd smiles-to-3d
python smiles3d/smiles3d.py -i example/molecules.csv -o example/molecules.sdf
```

You can pass the following specifications:
* Number of conformers (-n): by default 10
* Time to process molecule (-1): by default 60s max

The input file must have one column named "smiles" and one column named "id"

# License
Ersilia code is licensed under a GPLv3 License. The packages used are licensed according to their original licenses.