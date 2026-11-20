
import requests

from rdkit import Chem
from rdkit.Chem import (
    rdMolDescriptors,
    Descriptors,
    Crippen,
    Lipinski
)


def identificar_molecula(smiles):
    """
    Identifica una molécula y calcula sus principales
    descriptores químicos utilizando RDKit.
    """

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("Debes introducir un código SMILES.")

    mol = Chem.MolFromSmiles(smiles.strip())

    if mol is None:
        raise ValueError("El código SMILES no es válido.")

    # Canonicalización de la estructura
    smiles_canonico = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True
    )

    # Identificación de centros quirales
    centros_quirales = Chem.FindMolChiralCenters(
        mol,
        includeUnassigned=True
    )

    # Identificador molecular
    inchi_key = Chem.MolToInchiKey(mol)

    # Descriptores moleculares
    datos = {
        "smiles_original": smiles,
        "smiles_canonico": smiles_canonico,
        "inchi_key": inchi_key,

        "formula": rdMolDescriptors.CalcMolFormula(mol),

        "peso_molecular": Descriptors.MolWt(mol),

        "masa_exacta": rdMolDescriptors.CalcExactMolWt(mol),

        "atomos": mol.GetNumAtoms(),

        "atomos_totales": mol.GetNumAtoms(),

        "atomos_pesados": mol.GetNumHeavyAtoms(),

        "logp": Crippen.MolLogP(mol),

        "donadores_h": Lipinski.NumHDonors(mol),

        "aceptores_h": Lipinski.NumHAcceptors(mol),

        "tpsa": rdMolDescriptors.CalcTPSA(mol),

        "enlaces_rotables": Lipinski.NumRotatableBonds(mol),

        "centros_quirales": len(centros_quirales),

        "detalle_quiralidad": centros_quirales
    }

    return datos


def consultar_pubchem(smiles):
    """
    Consulta PubChem mediante su API PUG REST.

    Devuelve propiedades químicas disponibles.
    No realiza predicciones biológicas.
    """

    if not isinstance(smiles, str) or not smiles.strip():
        return {
            "estado": "error",
            "mensaje": "SMILES vacío."
        }

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
        "compound/smiles/property/"
        "MolecularFormula,MolecularWeight,"
        "InChIKey,IUPACName,XLogP/JSON"
    )

    try:
        response = requests.post(
            url,
            data={"smiles": smiles.strip()},
            timeout=20
        )

        if response.status_code == 404:
            return {
                "estado": "sin_resultados",
                "mensaje": "No se encontró una coincidencia en PubChem."
            }

        response.raise_for_status()

        resultado = response.json()

        propiedades = resultado[
            "PropertyTable"
        ]["Properties"][0]

        return {
            "estado": "exito",
            "datos": propiedades
        }

    except requests.Timeout:
        return {
            "estado": "error",
            "mensaje": "PubChem tardó demasiado en responder."
        }

    except requests.RequestException as e:
        return {
            "estado": "error",
            "mensaje": f"Error de conexión con PubChem: {e}"
        }

    except (KeyError, IndexError, ValueError):
        return {
            "estado": "error",
            "mensaje": "No se pudieron interpretar los datos recibidos."
        }
