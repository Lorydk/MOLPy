# changelog:
# 1.0 - Initial Release
# 1.1 - Fixed Synonym print if there are none
# 1.2 - Internal Release: first GUI tests
# 2.0 - GUI Release

# MOL test input: C1=CC2=C(C=C1)C(=O)C(=C(C2=O)C3CCC(CC3)C4=CC=C(C=C4)Cl)O

__version__ = "2.0.0"

import pubchempy as pcp
import PySimpleGUI as sg

from rdkit import Chem
from rdkit.Chem import QED, Draw


def get_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None: Chem.Kekulize(mol)
    return mol

def get_cid(smiles):
    cid = pcp.get_compounds(smiles, namespace='smiles')[0]
    if cid is not None: return cid
    
def get_molproperties(mol):
    molProperties = QED.properties(mol)
    if molProperties is not None: return molProperties

def get_qed(mol):
    molQED = QED.qed(mol)
    if molQED is not None: return molQED

def main():
    # Theme
    sg.theme('Dark2')
    
    # Input Column
    mol_input_column = [
        [
            sg.Text("Mol. SMILES:"),
            sg.In(size=(100, 1), tooltip="Write here molecule in SMILES format", key="-MOLSMILES-"),
            sg.Button("OK"),
        ],
        [
            sg.Text("Mol. ALOGP:"),
            sg.In(size=(25, 1), tooltip="Leave empty if unchanged", key="-MOLALOGP-"),
        ]
    ]

    # Data Column
    mol_data_column = [
        [sg.Text("Mol. Data")],
        [sg.Text(key="-CID-")],
        [sg.Text(key="-MF-"), sg.Text(key="-MW-")],
        [sg.Text(key="-SYNONYM-")],
        [sg.Text(key="-INCHIKEY-")],
        [sg.Text(key="-INCHI-")],
        [sg.Text(key="-IUPAC-")],
        [sg.Text(key="-QED-")]
    ]
    
    # Mol Props Column (table)
    mol_props_column = [
        [sg.Text("Mol. Props")],
        [sg.Table(headings=["Prop. Name", "Prop. Value"], values=[], auto_size_columns=True, justification='left', key="-MOLPROPS-")]
    ]
    
    mol_image_column = [
        [sg.Text("Mol. IMAGE (saved as '<SMILES>.png' in the <root>/images folder)")],
        [sg.Image(source="", key="-MOLIMAGE-")]
    ]

    # ----- Full layout -----
    layout = [
        [sg.Column(mol_input_column)],
        [sg.HSeparator()],
        [sg.Column(mol_data_column), sg.VSeparator(), sg.Column(mol_props_column)],
        [sg.HSeparator()],
        [sg.Column(mol_image_column)]
    ]
    
    # Create the window
    window = sg.Window("MolPY POC", layout, finalize=True)

    while True:
        event, values = window.read()
        if event == "Exit" or event == sg.WIN_CLOSED:
            break
        elif event == "OK":
            if values["-MOLSMILES-"]:
                # Main MOL data
                mol = get_mol(values["-MOLSMILES-"])
                
                # Continue printing data only if mol exists
                if mol:
                    if values["-MOLALOGP-"]:
                        mol = Chem.rdmolops.AdjustQueryProperties(mol, ["ALOGP", values["-MOLALOGP-"]])
                    cid = get_cid(values["-MOLSMILES-"])
                    mf = cid.molecular_formula
                    mw = cid.molecular_weight
                    inchikey = cid.inchikey
                    inchi = cid.inchi
                    iupac = cid.iupac_name
                    synonym = cid.synonyms[0] if cid.synonyms else "No Name Found"
                    molProps = get_molproperties(mol)
                    qed = get_qed(mol)
                    
                    molPropsRows = []
                    for propKey in molProps._fields:
                        molPropsRows.append([propKey, getattr(molProps, propKey)])
                    
                    # Create image of current MOL
                    molImg = Draw.MolToImage(mol)
                    molImg.save("./MOLPy/images/" + values["-MOLSMILES-"] + ".png")
                    
                    window["-CID-"].update("CID: " + str(cid))
                    window["-MF-"].update("MF: " + str(mf))
                    window["-MW-"].update("MW: " + str(mw) + " g/mol")
                    window["-SYNONYM-"].update("NAME: " + str(synonym))
                    window["-INCHIKEY-"].update("InChIKey: " + str(inchikey))
                    window["-INCHI-"].update("InChI: " + str(inchi))
                    window["-IUPAC-"].update("IUPAC: " + str(iupac))
                    window["-MOLPROPS-"].update(values=molPropsRows)
                    window["-QED-"].update("QED: " + str(qed))
                    window["-MOLIMAGE-"].update(source=values["-MOLSMILES-"] + ".png")
                else:
                    window["-CID-"].update("Error, SMILES not recognized")
                    window["-MF-"].update("")
                    window["-MW-"].update("")
                    window["-SYNONYM-"].update("")
                    window["-INCHIKEY-"].update("")
                    window["-INCHI-"].update("")
                    window["-IUPAC-"].update("")
                    window["-MOLPROPS-"].update(values=[])
                    window["-QED-"].update("")

    window.close()
    
if __name__ == "__main__":
    main()
