from rdkit import Chem
import itertools

def add_groups(m, mapping):

    # Loop over atoms until there are no wildcard atoms
    while True:

        # Find wildcard atom if available, otherwise exit
        a = None
        for a_ in m.GetAtoms():
            if a_.GetAtomicNum() == 0:
                a = a_
                break
        if a is None:
            break

        # Get appropriate group to substitute in from mapping
        group = mapping[a.GetAtomMapNum()]

        if group.GetNumAtoms() == 1:

            # Directly substitute atom in, if single atom group
            a.SetAtomicNum(group.GetAtomWithIdx(0).GetAtomicNum())
            a.SetAtomMapNum(0)
        else:

            # Set wildcard atoms to having AtomMapNum 1000 for tracking
            a.SetAtomMapNum(1000)

            for a_ in group.GetAtoms():
                if a_.GetAtomicNum() == 0:
                    a_.SetAtomMapNum(1000)

            # Put group and base molecule together and make it editable
            m = Chem.CombineMols(m, group)
            m = Chem.RWMol(m)

            # Find using tracking number the atoms to merge in new molecule
            a1 = None
            a2 = None
            bonds = []
            for a in m.GetAtoms():
                if a.GetAtomMapNum() == 1000:
                    if a1 is None:
                        a1 = a
                    else:
                        a2 = a

            # Find atoms to bind together based on atoms to merge
            b1 = a1.GetBonds()[0]
            start = (b1.GetBeginAtomIdx() if b1.GetEndAtomIdx() == a1.GetIdx()
                else b1.GetEndAtomIdx())

            b2 = a2.GetBonds()[0]
            end = (b2.GetBeginAtomIdx() if b2.GetEndAtomIdx() == a2.GetIdx()
                else b2.GetEndAtomIdx())

            # Add the connection and remove original wildcard atoms
            m.AddBond(start, end, order=Chem.rdchem.BondType.SINGLE)
            m.RemoveAtom(a1.GetIdx())
            m.RemoveAtom(a2.GetIdx())

    return m

import pandas as pd
df = pd.read_csv(r'C:\Users\cmacd\OneDrive - University of Glasgow\Documents\Python Scripts\SMILES Generator\AASC_SMILES_Curated.csv')
corestructurename = "2Nap"

# Create a list of all permutations of rows in the second column (including self-permutations)
permutations = list(itertools.permutations(df['SMILES'], 2)) + [(x, x) for x in df['SMILES']]

# Create a new dataframe with the permutations and respective Column1 values
new_df = pd.DataFrame(permutations, columns=['X', 'Y'])
new_df['AASC1'] = new_df["X"].apply(lambda x: df.loc[df['SMILES'] == x, 'Letter'].iloc[0])
new_df['AASC2'] = new_df["Y"].apply(lambda x: df.loc[df['SMILES'] == x, 'Letter'].iloc[0])

new_df["mol"] = [add_groups(Chem.MolFromSmiles("O=C(NC([*:1])C(NC([*:2])C(O)=O)=O)COC1=CC2=CC=CC=C2C=C1"),
    {1: Chem.MolFromSmiles(r["X"]),
    2: Chem.MolFromSmiles(r["Y"])}) for i, r in new_df.iterrows()]

output = {}
output["Compound"] = [corestructurename + '-' + str(r["AASC1"]) + str(r["AASC2"]) for i, r in new_df.iterrows()]
output["smi"] = [Chem.MolToSmiles(m) for m in new_df["mol"]]
output_df = pd.DataFrame(output)

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(new_df["mol"], molsPerRow=4, useSVG=True)

# Save the output to a CSV file
output_df.to_csv(r"C:\Users\cmacd\OneDrive - University of Glasgow\Documents\Python Scripts\SMILES Generator\Output.csv", index=False)
print("Output saved to output_smiles.csv")