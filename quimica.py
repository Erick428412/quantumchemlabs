import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import Crippen

def get_mol_svg(mol, width=300, height=200):
    """Genera el string SVG para renderizar el dibujo de la molécula."""
    mol_to_draw = Chem.Mol(mol)
    AllChem.Compute2DCoords(mol_to_draw)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.backgroundColour = (0.12, 0.14, 0.18, 1.0)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_to_draw)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

def generar_isomeros(smiles_input):
    """Procesa el SMILES base y devuelve la lista de isómeros generados por hibridación."""
    mol_base = Chem.MolFromSmiles(smiles_input)
    if mol_base is None:
        return None
    opts = StereoEnumerationOptions()
    opts.maxIsomers = 8
    return list(EnumerateStereoisomers(mol_base, options=opts))

def calcular_propiedades_fisicas(mol_activa):
    """Calcula el punto de ebullición y fusión usando modelos de contribución."""
    mw = Chem.rdMolDescriptors.CalcExactMolWt(mol_activa)
    logp = Crippen.MolLogP(mol_activa)
    h_donors = Chem.rdMolDescriptors.CalcNumHBD(mol_activa)
    
    est_boiling_k = 120 + (8.5 * mw) + (45 * h_donors) - (15 * logp)
    est_boiling_c = est_boiling_k - 273.15
    est_melting_k = est_boiling_k * 0.58 + (20 * h_donors)
    est_melting_c = est_melting_k - 273.15
    
    return est_boiling_c, est_melting_c, mw, logp, h_donors

def calcular_perfil_viscosidad(mw, h_donors, t_min, t_max):
    """Simula la curva de viscosidad usando la ecuación de Arrhenius-Andrade."""
    E_visc = (mw * 0.15 + h_donors * 12.0) * 1000
    R = 8.314
    A_visc = 0.001 / (mw * 0.05)
    
    temp_c_range = np.linspace(t_min, t_max, 100)
    tk_range = temp_c_range + 273.15
    viscosidad_cps = A_visc * np.exp(E_visc / (R * tk_range))
    
    return temp_c_range, viscosidad_cps, E_visc, R, A_visc

def analizar_polaridad_y_polos(mol_with_hs):
    """Calcula cargas de Gasteiger-Marsili y baricentros de polos eléctricos."""
    mol_cargas = Chem.Mol(mol_with_hs)
    AllChem.Compute2DCoords(mol_cargas)
    AllChem.ComputeGasteigerCharges(mol_cargas)
    
    enlaces_polares = 0
    enlaces_no_polares = 0
    for bond in mol_cargas.GetBonds():
        a1 = bond.GetBeginAtom().GetSymbol()
        a2 = bond.GetEndAtom().GetSymbol()
        if (a1 in ['O', 'N', 'F', 'Cl'] and a2 in ['C', 'H']) or (a2 in ['O', 'N', 'F', 'Cl'] and a1 in ['C', 'H']):
            enlaces_polares += 1
        else:
            enlaces_no_polares += 1
            
    conf = mol_cargas.GetConformer()
    pos_x, pos_y, cargas = [], [], []
    for i, atom in enumerate(mol_cargas.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        chg = float(atom.GetProp('_GasteigerCharge'))
        if np.isnan(chg): chg = 0.0
        pos_x.append(pos.x)
        pos_y.append(pos.y)
        cargas.append(chg)
        
    pos_x, pos_y, cargas = np.array(pos_x), np.array(pos_y), np.array(cargas)
    idx_pos, idx_neg = cargas > 0.01, cargas < -0.01
    
    c_pos_x = float(np.mean(pos_x[idx_pos])) if np.any(idx_pos) else 0.0
    c_pos_y = float(np.mean(pos_y[idx_pos])) if np.any(idx_pos) else 0.0
    c_neg_x = float(np.mean(pos_x[idx_neg])) if np.any(idx_neg) else 0.0
    c_neg_y = float(np.mean(pos_y[idx_neg])) if np.any(idx_neg) else 0.0
    
    distancia_polos = float(np.sqrt((c_pos_x - c_neg_x)**2 + (c_pos_y - c_neg_y)**2))
    es_polar = "SÍ" if (distancia_polos > 0.2 and (np.any(idx_pos) and np.any(idx_neg))) else "NO / MUY BAJA"
    
    return enlaces_polares, enlaces_no_polares, es_polar, pos_x, pos_y, idx_pos, idx_neg, c_pos_x, c_pos_y, c_neg_x, c_neg_y, mol_cargas
