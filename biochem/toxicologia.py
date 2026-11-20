
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import (
    FilterCatalog,
    FilterCatalogParams
)


def crear_catalogo(tipo):
    """Crea un catálogo de alertas estructurales de RDKit."""

    parametros = FilterCatalogParams()

    if tipo == "BRENK":
        parametros.AddCatalog(
            FilterCatalogParams.FilterCatalogs.BRENK
        )
    elif tipo == "NIH":
        parametros.AddCatalog(
            FilterCatalogParams.FilterCatalogs.NIH
        )

    return FilterCatalog(parametros)


def analizar_toxicologia(smiles):

    # Validar molécula
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
    propiedades = {
        "peso_molecular": round(Descriptors.MolWt(mol), 3),
        "logp": round(Crippen.MolLogP(mol), 3),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "donadores_h": rdMolDescriptors.CalcNumLipinskiHBD(mol),
        "aceptores_h": rdMolDescriptors.CalcNumLipinskiHBA(mol)
    }

    # Catálogos de alertas estructurales
    catalogos = {
        "Brenk": crear_catalogo("BRENK"),
        "NIH": crear_catalogo("NIH")
    }

    alertas = {}

    for nombre, catalogo in catalogos.items():

        resultados = catalogo.GetMatches(mol)

        alertas[nombre] = [
            {
                "nombre": resultado.GetDescription()
            }
            for resultado in resultados
        ]

    total_alertas = sum(
        len(lista) for lista in alertas.values()
    )

    # Interpretación prudente
    if total_alertas == 0:
        interpretacion = (
            "No se encontraron patrones incluidos "
            "en los catálogos utilizados."
        )
    else:
        interpretacion = (
            "Se encontraron patrones estructurales que "
            "merecen revisión. Esto no demuestra toxicidad."
        )

    return {
        "estado": "correcto",
        "smiles": smiles,
        "propiedades": propiedades,
        "alertas": alertas,
        "total_alertas": total_alertas,
        "interpretacion": interpretacion,
        "limitaciones": [
            "Las alertas son patrones estructurales, no diagnósticos.",
            "La ausencia de alertas no demuestra seguridad.",
            "Los resultados dependen de los catálogos utilizados.",
            "No se predicen efectos en personas ni animales.",
            "No se calcula una dosis ni un nivel de toxicidad."
        ]
    }

