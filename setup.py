from setuptools import setup, find_packages

setup(
    name="smiles3d",
    version="0.1.0",
    author="Miquel Duran-Frigola, Ersilia Open Source Initiative",
    author_email="miquel@ersilia.io",
    url="https://github.com/ersilia-os/smiles-to-3d",
    description="CDPKit based module to convert lists of smiles to 3D conformers",
    license="GPLv3",
    python_requires=">=3.8",
    install_requires=["rdkit==2023.9.5", "cdpkit==1.1.1"],
    extras_require={},
    packages=find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ),
    keywords="conformers",
    project_urls={"Source Code": "https://github.com/ersilia-os/smiles-to-3d/"},
)