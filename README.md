# Ligand RMSD Calculator

This repository contains a Python script to calculate the Root Mean Square Deviation (RMSD) between two ligand molecules based on their Maximum Common Substructure (MCS). It utilizes the RDKit library for molecular manipulation and MCS finding.
## Features
* Calculates RMSD between two ligand PDB files.
* Uses Maximum Common Substructure (MCS) for robust alignment.
* Optionally renames atoms in the second ligand to match the first ligand's MCS and saves the modified ligand to a new PDB file.

## Requirements
* Python =< 3.9
* numpy
* rdkit

## Installation
This code can be run in any conda environment where both NumPy and RDKit are installed.

### Create a conda environment:
* This will install Numpy
  conda create -n ligrmsd python=3.9 -y

activate ligrmsd

* Install RDKit
  conda install conda-forge::rdkit -y

