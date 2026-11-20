
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors


def analizar_adme(smiles):
    """
    Analiza descriptores moleculares relacionados con ADME.

    Importante:
    Los resultados son indicadores fisicoquímicos.
    No representan mediciones experimentales ni predicciones clínicas.
    """

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return {
            "estado": "error",
            "mensaje": "El SMILES no es válido."
        }

    # Descriptores moleculares
    peso_molecular = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    donadores_h = rdMolDescriptors.CalcNumHBD(mol)
    aceptores_h = rdMolDescriptors.CalcNumHBA(mol)
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

    # Criterios de Veber
    criterios_veber = {
        "enlaces_rotables": enlaces_rotables <= 10,
        "tpsa": tpsa <= 140
    }

    cumple_veber = all(criterios_veber.values())

    # Interpretación descriptiva de propiedades
    if tpsa <= 90:
        interpretacion_tpsa = (
            "Área polar relativamente baja."
        )
    elif tpsa <= 140:
        interpretacion_tpsa = (
            "Área polar intermedia."
        )
    else:
        interpretacion_tpsa = (
            "Área polar elevada."
        )

    if logp < 0:
        interpretacion_logp = (
            "Predominio de características hidrofílicas."
        )
    elif logp <= 5:
        interpretacion_logp = (
            "Lipofilicidad dentro del intervalo de referencia de Lipinski."
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
            "incumplimientos": incumplimientos,
            "cumple": incumplimientos == 0
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
            "Los descriptores no miden absorción real.",
            "Lipinski y Veber son reglas orientativas.",
            "No se calculan concentraciones en sangre.",
            "No se predicen dosis ni seguridad clínica.",
            "Metabolismo y excreción requieren modelos y datos adicionales."
        ]
    }

