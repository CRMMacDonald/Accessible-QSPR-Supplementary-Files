from  rdkit import Chem, DataStructs
from rdkit import Chem

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
df = pd.read_csv(r"C:\Users\cmacd\OneDrive - University of Glasgow\Documents\Python Scripts\SMILES Generator\AASC_SMILES_Curated.csv")
corestructurename = "2Nap"

df["mol"] = [add_groups(Chem.MolFromSmiles("O=C(NC([*:1])C(O)=O)COC1=CC2=CC=CC=C2C=C1"),
    {1: Chem.MolFromSmiles(r["SMILES"]),
    2: Chem.MolFromSmiles(r["SMILES"])}) for i, r in df.iterrows()]

output = {}
output["Compound"] = [corestructurename + '-' + str(r["Letter"]) for i, r in df.iterrows()]
output["smi"] = [Chem.MolToSmiles(m) for m in df["mol"]]
output_df = pd.DataFrame(output)

from rdkit.Chem import Draw
img = Draw.MolsToGridImage(df["mol"], molsPerRow=4, useSVG=True)

output_df.to_csv(r"C:\Users\cmacd\OneDrive - University of Glasgow\Documents\Python Scripts\SMILES Generator\Output.csv", index=False)
print("Output saved to output_smiles.csv")

