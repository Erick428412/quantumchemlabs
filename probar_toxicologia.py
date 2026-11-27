
from biochem.toxicologia import analizar_toxicologia

# Moléculas de prueba
moleculas = {
    "Agua": "O",
    "Etanol": "CCO",
    "Aspirina": "CC(=O)Oc1ccccc1C(=O)O",
    "SMILES incorrecto": "XYZ123"
}

for nombre, smiles in moleculas.items():
    print("\n" + "=" * 45)
    print("MOLÉCULA:", nombre)
    print("SMILES:", smiles)

    resultado = analizar_toxicologia(smiles)

    if resultado["estado"] == "error":
        print("ERROR:", resultado["mensaje"])
        continue

    print("\nPROPIEDADES:")
    for propiedad, valor in resultado["propiedades"].items():
        print(f"{propiedad}: {valor}")

    print("\nALERTAS:")
    for catalogo, alertas in resultado["alertas"].items():
        print(f"\n{catalogo}:")
        if alertas:
            for alerta in alertas:
                print("-", alerta["nombre"])
        else:
            print("Sin patrones encontrados.")

    print("\nTOTAL DE ALERTAS:", resultado["total_alertas"])
    print("INTERPRETACIÓN:", resultado["interpretacion"])
