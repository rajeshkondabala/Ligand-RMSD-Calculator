import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFMCS
import argparse # Import the argparse module

def CalcLigRMSD(lig1, lig2, rename_lig2=True, output_filename="tmp.pdb"):
    """
    Calculates the Root Mean Square Deviation (RMSD) between two ligands
    based on their Maximum Common Substructure (MCS).

    Args:
        lig1 (rdkit.Chem.Mol): The first ligand molecule.
        lig2 (rdkit.Chem.Mol): The second ligand molecule.
        rename_lig2 (bool): If True, renames atoms in lig2 to match lig1's MCS
                            and saves lig2 to a PDB file. Defaults to True.
        output_filename (str): The filename for the output PDB if rename_lig2 is True.
                               Defaults to "tmp.pdb".

    Returns:
        float: The minimum RMSD between the two ligands based on MCS matches.
    """
    # Remove hydrogens for RMSD calculation as they don't contribute to the heavy atom skeleton
    lig1 = Chem.RemoveHs(lig1)
    lig2 = Chem.RemoveHs(lig2)

    # Get atomic coordinates for both ligands
    coordinates_lig2 = lig2.GetConformer().GetPositions()
    coordinates_lig1 = lig1.GetConformer().GetPositions()

    # Find the Maximum Common Substructure (MCS) between the two ligands
    # This helps in identifying corresponding atoms for alignment
    res = rdFMCS.FindMCS([lig1, lig2])
    ref_mol = Chem.MolFromSmarts(res.smartsString)

    # Get the atom indices in lig1 that correspond to the MCS
    mas1 = list(lig1.GetSubstructMatch(ref_mol))

    # Get all possible sets of atom indices in lig2 that match the MCS
    # uniquify=False ensures we get all matches, not just unique ones
    mas2_list = lig2.GetSubstructMatches(ref_mol, uniquify=False)

    # Extract coordinates for the MCS atoms in lig1
    coordinates_lig1_mcs = coordinates_lig1[mas1]

    list_rmsd = []
    # Iterate through all possible MCS matches in lig2
    for match2 in mas2_list:
        # Extract coordinates for the current MCS match in lig2
        coordinates_lig2_mcs_tmp = coordinates_lig2[list(match2)]
        # Calculate the difference vector between corresponding MCS atoms
        diff = coordinates_lig2_mcs_tmp - coordinates_lig1_mcs
        # Calculate RMSD for the current match: sqrt(sum(diff^2) / num_atoms)
        list_rmsd.append(np.sqrt((diff * diff).sum() / len(coordinates_lig2_mcs_tmp)))

    # The final RMSD is the minimum among all possible MCS alignments
    lig_rmsd = min(list_rmsd)

    # Optional: Rename atoms in lig2 to match lig1's MCS and save to PDB
    if rename_lig2:
        # Get the MCS match in lig2 that resulted in the minimum RMSD
        mas2 = mas2_list[np.argmin(list_rmsd)]
        # Create a mapping from lig2 MCS atom index to lig1 MCS atom index
        correspondence_key2_item1 = dict(zip(mas2, mas1))

        # Get atom names from lig1 for renaming
        atom_names_lig1 = [atom1.GetPDBResidueInfo().GetName() for atom1 in lig1.GetAtoms()]
        # Get the residue name from lig1 to apply to lig2
        lig1_ResName = lig1.GetAtoms()[0].GetPDBResidueInfo().GetResidueName()

        # Iterate through atoms in lig2 and rename them
        for i, atom2 in enumerate(lig2.GetAtoms()):
            # Set the residue name of lig2 atoms to match lig1's residue name
            atom2.GetPDBResidueInfo().SetResidueName(lig1_ResName)
            # If the current atom in lig2 is part of the MCS match, rename it
            if i in correspondence_key2_item1.keys():
                atom2.GetPDBResidueInfo().SetName(atom_names_lig1[correspondence_key2_item1[i]])

        # Save the modified lig2 molecule to a PDB file
        Chem.MolToPDBFile(lig2, output_filename)

    return lig_rmsd

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Calculate RMSD between two ligands based on MCS.")

    # Add arguments
    parser.add_argument('--ligand_1', type=str, required=True,
                        help='Path to the first ligand PDB file.')
    parser.add_argument('--ligand_2', type=str, required=True,
                        help='Path to the second ligand PDB file.')
    parser.add_argument('--output_file', type=str, default="tmp.pdb",
                        help='Optional: Output filename for the renamed ligand 2 PDB file. Defaults to "tmp.pdb".')
    parser.add_argument('--no_rename', action='store_true',
                        help='Optional: Do not rename atoms in ligand 2 or save an output PDB file.')

    # Parse the arguments
    args = parser.parse_args()

    # Load the ligand molecules from PDB files
    try:
        lig1 = Chem.MolFromPDBFile(args.ligand_1)
        lig2 = Chem.MolFromPDBFile(args.ligand_2)
    except Exception as e:
        print(f"Error loading PDB files: {e}")
        print("Please ensure the file paths are correct and the files are valid PDB formats.")
        exit(1)

    if lig1 is None:
        print(f"Could not read ligand_1 from {args.ligand_1}. Please check the file.")
        exit(1)
    if lig2 is None:
        print(f"Could not read ligand_2 from {args.ligand_2}. Please check the file.")
        exit(1)

    # Determine whether to rename/save based on --no_rename flag
    rename_lig2_flag = not args.no_rename

    # Calculate RMSD
    rmsd = CalcLigRMSD(lig1, lig2, rename_lig2=rename_lig2_flag, output_filename=args.output_file)

    print(f"RMSD between {args.ligand_1} and {args.ligand_2}: {rmsd:.4f} Å")

    if rename_lig2_flag:
        print(f"Renamed ligand 2 saved to: {args.output_file}")
