
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def analizar_adme(smiles):
    """
    Calcula descriptores moleculares y evalúa criterios
    orientativos de Lipinski y Veber.

    No predice directamente absorción, distribución,
    metabolismo, excreción ni seguridad clínica.
    """

    # Validar SMILES
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        mol = None

    if mol is None:
        return {
            "estado": "error",
            "mensaje": "El SMILES no es válido."
        }

    # Descriptores moleculares
    peso_molecular = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    # Descriptores específicos de Lipinski
    donadores_h = rdMolDescriptors.CalcNumLipinskiHBD(mol)
    aceptores_h = rdMolDescriptors.CalcNumLipinskiHBA(mol)

    enlaces_rotables = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Regla de Lipinski
    criterios_lipinski = {
        "peso_molecular": peso_molecular <= 500,
        "logp": logp <= 5,
        "donadores_h": donadores_h <= 5,
        "aceptores_h": aceptores_h <= 10
    }

    incumplimientos = sum(
        not cumple for cumple in criterios_lipinski.values()
    )

    # Se informa el cumplimiento estricto de los cuatro criterios
    cumple_lipinski = incumplimientos == 0

    # Regla de Veber
    criterios_veber = {
        "enlaces_rotables": enlaces_rotables <= 10,
        "tpsa": tpsa <= 140
    }

    cumple_veber = all(criterios_veber.values())

    # Interpretación de TPSA
    if tpsa <= 90:
        interpretacion_tpsa = "Área polar relativamente baja."
    elif tpsa <= 140:
        interpretacion_tpsa = "Área polar intermedia."
    else:
        interpretacion_tpsa = "Área polar elevada."

    # Interpretación de LogP
    if logp < 0:
        interpretacion_logp = (
            "Predominio de afinidad por la fase acuosa."
        )
    elif logp <= 5:
        interpretacion_logp = (
            "Lipofilicidad dentro del límite de referencia de Lipinski."
        )
    else:
        interpretacion_logp = (
            "Lipofilicidad elevada según el criterio de Lipinski."
        )

    return {
        "estado": "correcto",
        "smiles": smiles,

        "propiedades": {
            "peso_molecular": round(peso_molecular, 3),
            "logp": round(logp, 3),
            "tpsa": round(tpsa, 3),
            "donadores_h": int(donadores_h),
            "aceptores_h": int(aceptores_h),
            "enlaces_rotables": int(enlaces_rotables)
        },

        "lipinski": {
            "criterios": criterios_lipinski,
            "incumplimientos": int(incumplimientos),
            "cumple": cumple_lipinski
        },

        "veber": {
            "criterios": criterios_veber,
            "cumple": cumple_veber
        },

        "interpretacion": {
            "polaridad": interpretacion_tpsa,
            "lipofilicidad": interpretacion_logp
        },

        "limitaciones": [
            "Los descriptores son estimaciones computacionales.",
            "LogP no equivale directamente a solubilidad.",
            "Lipinski y Veber son reglas orientativas.",
            "Cumplir las reglas no demuestra absorción real.",
            "No se predicen concentraciones en sangre ni dosis.",
            "No se evalúa seguridad clínica.",
            "No se modelan directamente metabolismo ni excreción."
        ]
    }

