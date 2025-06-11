from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, Draw
from fastapi import HTTPException
import base64
from io import BytesIO
from fastapi.responses import JSONResponse

'''
Typical toxic subgroups, will flag if found in the molecule
'''
TOXIC_SUBSTRUCTURES = {
    "Nitro group": "[N+](=O)[O-]",
    "Azo group": "N=N",
    "Aromatic amine": "cN",
    "Aniline": "c1ccc(cc1)N",
    "Nitrile": "C#N"
}

def predict_properties(smiles: str) -> dict:
    '''
    Logic: predict properties of molecules i.e. is toxic, is soluble, etc. based on functional groups or string patterns
    '''
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise HTTPException(status_code=400, detail="Invalid SMILES string.")
    
    try: 
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        rot_bonds = Lipinski.NumRotatableBonds(mol)
        aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)

        
        # ESOL-based solubility classification
        if logp > 4:
            solubility = "poor"
        elif logp > 2:
            solubility = "moderate"
        else:
            solubility = "high"
        
        # Lipinski's Rule of 5 violations
        '''
        Lipinski's Rule of 5 is a heuristic for molecular absorption in the mody
        '''
        lipinski_violations = sum([
            mw > 500,
            logp > 5,
            hbd > 5,
            hba > 10
        ])
        lipinski_result = {
            "violations": lipinski_violations,
            "is_druglike": lipinski_violations <= 1
        }

        # Molecular visualization as base64 PNG
        img = Draw.MolToImage(mol, size=(250, 250))
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        mol_image_base64 = base64.b64encode(buffer.getvalue()).decode()
        
        # Toxicophore flags
        detected_toxicophores = []
        for name, smarts in TOXIC_SUBSTRUCTURES.items():
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern):
                detected_toxicophores.append(name)

        return JSONResponse({
            "molecular weight": round(mw, 2),
            "logP": round(logp, 2),
            "H_bond_donors": hbd, 
            "H_bond_acceptors": hba,
            "TPSA": round(tpsa, 2),
            "rotatable_bonds": rot_bonds,
            "aromatic_rings": aromatic_rings,
            "esol_solubility": solubility,
            "lipinski": lipinski_result,
            "likely_toxic": bool(detected_toxicophores),
            "toxicophores_detected": detected_toxicophores,
            "validated_smiles": Chem.MolToSmiles(mol),
            "molecule_image_base64": mol_image_base64
        })
    
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"error computing properties: {str(e)}")
    