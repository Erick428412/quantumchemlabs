from rdkit.Chem import Crippen
import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from qiskit import QuantumCircuit

import quimica
from biochem.adme import analizar_adme


# ==========================================================
# CONFIGURACIÓN GENERAL
# ==========================================================

st.set_page_config(
    page_title="QuantumChemLabs",
    page_icon="🧪",
    layout="wide"
)

st.title("🧬 QuantumChemLabs")
st.write(
    "Plataforma de análisis químico, computación cuántica "
    "y biología computacional."
)


# ==========================================================
# FUNCIONES AUXILIARES
# ==========================================================

def render_mol_local(mol, width=600, height=350):
    """Renderiza una molécula utilizando RDKit."""

    svg = quimica.get_mol_svg(mol, width, height)

    html_content = f"""
    <div style="
        display: flex;
        justify-content: center;
        align-items: center;
        background-color: #0e1117;
        padding: 10px;
        border-radius: 8px;
    ">
        {svg}
    </div>
    """

    components.html(html_content, height=height + 20)


# ==========================================================
# BARRA LATERAL
# ==========================================================

st.sidebar.title("⚙️ Configuración")

st.sidebar.header("Entrada molecular")

smiles_input = st.sidebar.text_input(
    "Ingresa el código SMILES:",
    value="C(O)=C(C)O"
)


# ==========================================================
# PESTAÑAS PRINCIPALES
# ==========================================================

tab_quimica, tab_cuantica, tab_biologia = st.tabs([
    "🧪 Química",
    "⚛️ Cuántica",
    "🧬 Biología"
])


# ==========================================================
# 1. LABORATORIO DE QUÍMICA
# ==========================================================

with tab_quimica:

    st.header("🧪 Laboratorio de Química")

    try:
        isomeros_lista = quimica.generar_isomeros(smiles_input)

        if isomeros_lista is None or len(isomeros_lista) == 0:

            st.error(
                "Código SMILES no válido. "
                "Ingresa una estructura molecular válida."
            )

        else:

            # Reiniciar selección si cambia el SMILES
            if st.session_state.get("ultimo_smiles") != smiles_input:
                st.session_state.ultimo_smiles = smiles_input
                st.session_state.isomero_seleccionado = 0

            # Evitar índices fuera de rango
            if (
                st.session_state.isomero_seleccionado
                >= len(isomeros_lista)
            ):
                st.session_state.isomero_seleccionado = 0

            # --------------------------------------------------
            # SUBPESTAÑAS DE QUÍMICA
            # --------------------------------------------------

            sub_quimica = st.tabs([
                "⚛️ Descriptores Fisicoquímicos",
                "🌡️ Propiedades Físicas",
                "🧭 Polaridad",
                "💧 Solubilidad",
                "⚗️ Reacciones",
                "⚖️ Estequiometría"
            ])

            # --------------------------------------------------
            # MOLÉCULA E ISÓMEROS
            # --------------------------------------------------

            with sub_quimica[0]:

                st.subheader("Isómeros químicos detectados")

                st.write(
                    f"Se generaron **{len(isomeros_lista)}** "
                    "configuraciones moleculares."
                )

                cols = st.columns(min(4, len(isomeros_lista)))

                for idx, iso in enumerate(isomeros_lista):

                    with cols[idx % 4]:

                        svg_iso = quimica.get_mol_svg(
                            iso,
                            width=220,
                            height=160
                        )

                        components.html(
                            f"<div style='text-align:center;'>{svg_iso}</div>",
                            height=170
                        )

                        if st.button(
                            f"Inspeccionar Isómero {idx + 1}",
                            key=f"btn_iso_{idx}"
                        ):
                            st.session_state.isomero_seleccionado = idx
                            st.rerun()

                mol_activa = isomeros_lista[
                    st.session_state.isomero_seleccionado
                ]

                mol_with_hs = Chem.AddHs(mol_activa)

                st.divider()

                st.subheader(
                    "🌐 Estructura seleccionada "
                    f"(Isómero {st.session_state.isomero_seleccionado + 1})"
                )

                render_mol_local(mol_with_hs)

                col1, col2 = st.columns(2)

                with col1:
                    st.subheader("Datos generales")

                    st.write(
                        "**Fórmula química:**",
                        rdMolDescriptors.CalcMolFormula(mol_with_hs)
                    )

                    st.write(
                        "**Peso molecular:**",
                        f"{rdMolDescriptors.CalcExactMolWt(mol_with_hs):.4f} g/mol"
                    )

                    st.write(
                        "**Átomos totales:**",
                        mol_with_hs.GetNumAtoms()
                    )

                with col2:
                    st.subheader("Estereoquímica")

                    centros = Chem.FindMolChiralCenters(
                        mol_activa,
                        includeUnassigned=True
                    )

                    st.write(
                        "**Centros quirales:**",
                        len(centros)
                    )

            # --------------------------------------------------
            # PROPIEDADES FÍSICAS
            # --------------------------------------------------

            
            with sub_quimica[1]:

                st.subheader("⚛️ Descriptores fisicoquímicos")

                peso_molecular = rdMolDescriptors.CalcExactMolWt(mol_activa)
                logp = Crippen.MolLogP(mol_activa)
                tpsa = rdMolDescriptors.CalcTPSA(mol_activa)
                donadores = rdMolDescriptors.CalcNumHBD(mol_activa)
                aceptores = rdMolDescriptors.CalcNumHBA(mol_activa)
                enlaces_rotables = rdMolDescriptors.CalcNumRotatableBonds(mol_activa)

                col1, col2, col3 = st.columns(3)

                col1.metric("Peso molecular", f"{peso_molecular:.3f} g/mol")
                col2.metric("LogP", f"{logp:.3f}")
                col3.metric("TPSA", f"{tpsa:.3f} Å²")

                col4, col5, col6 = st.columns(3)

                col4.metric("Donadores H", donadores)
                col5.metric("Aceptores H", aceptores)
                col6.metric("Enlaces rotables", enlaces_rotables)

                st.info(
                    "Estos son descriptores moleculares calculados "
                    "computacionalmente. No representan mediciones "
                    "experimentales de puntos de fusión, ebullición "
                    "o viscosidad."
                )


            # --------------------------------------------------
            # POLARIDAD
            # --------------------------------------------------

            with sub_quimica[2]:

                st.subheader("🧭 Análisis vectorial de polaridad")

                (
                    enlaces_p_cnt,
                    enlaces_np_cnt,
                    es_polar,
                    pos_x,
                    pos_y,
                    idx_pos,
                    idx_neg,
                    c_pos_x,
                    c_pos_y,
                    c_neg_x,
                    c_neg_y,
                    mol_cargas
                ) = quimica.analizar_polaridad_y_polos(mol_with_hs)

                col1, col2 = st.columns(2)

                col1.write(
                    f"Enlaces covalentes polares: **{enlaces_p_cnt}**"
                )

                col2.write(
                    f"Enlaces covalentes no polares: **{enlaces_np_cnt}**"
                )

                st.divider()

                st.write(
                    f"¿La molécula es globalmente polar?: **{es_polar}**"
                )

                fig_p, ax_p = plt.subplots(figsize=(7, 3.5))

                fig_p.patch.set_facecolor("#0e1117")
                ax_p.set_facecolor("#1e1e1e")

                ax_p.scatter(
                    pos_x,
                    pos_y,
                    c="#4a5568",
                    s=100,
                    alpha=0.5
                )

                for i, atom in enumerate(mol_cargas.GetAtoms()):

                    ax_p.text(
                        pos_x[i] + 0.12,
                        pos_y[i] + 0.12,
                        atom.GetSymbol(),
                        color="#a0aec0",
                        fontsize=9
                    )

                if np.any(idx_pos):

                    ax_p.scatter(
                        c_pos_x,
                        c_pos_y,
                        c="#ff4b4b",
                        s=350,
                        marker="P",
                        label="Polo δ+"
                    )

                if np.any(idx_neg):

                    ax_p.scatter(
                        c_neg_x,
                        c_neg_y,
                        c="#00b4d8",
                        s=350,
                        marker="o",
                        label="Polo δ-"
                    )

                if es_polar == "SÍ":

                    ax_p.arrow(
                        c_pos_x,
                        c_pos_y,
                        c_neg_x - c_pos_x,
                        c_neg_y - c_pos_y,
                        head_width=0.15,
                        head_length=0.15,
                        fc="#ffb703",
                        ec="#ffb703",
                        length_includes_head=True
                    )

                ax_p.set_title(
                    "Mapa de polos eléctricos",
                    color="white"
                )

                ax_p.tick_params(colors="white")
                ax_p.axis("equal")

                st.pyplot(fig_p)
                plt.close(fig_p)

                st.caption(
                    "Representación aproximada en 2D basada en "
                    "cargas parciales de Gasteiger. No es un "
                    "cálculo del momento dipolar molecular 3D."
                )

            # --------------------------------------------------
            # SOLUBILIDAD
            # --------------------------------------------------

            with sub_quimica[3]:

                st.subheader(
                    "💧 Análisis orientativo de lipofilicidad"
                )

                _, _, _, logp_val, _ = (
                    quimica.calcular_propiedades_fisicas(mol_activa)
                )

                if logp_val < 0:

                    perfil_quimico = "Relativamente hidrofílico"
                    desc_solubilidad = (
                        "El descriptor LogP sugiere mayor afinidad "
                        "relativa por medios acuosos."
                    )

                elif logp_val <= 3:

                    perfil_quimico = "Lipofilicidad intermedia"
                    desc_solubilidad = (
                        "Presenta un valor intermedio de LogP."
                    )

                else:

                    perfil_quimico = "Relativamente lipofílico"
                    desc_solubilidad = (
                        "El descriptor LogP sugiere mayor afinidad "
                        "relativa por medios orgánicos."
                    )

                col1, col2 = st.columns(2)

                with col1:

                    st.metric(
                        "Descriptor LogP",
                        f"{logp_val:.3f}"
                    )

                    st.write(f"**Perfil:** {perfil_quimico}")
                    st.write(desc_solubilidad)

                with col2:

                    fig_s, ax_s = plt.subplots(figsize=(6, 2))

                    fig_s.patch.set_facecolor("#0e1117")
                    ax_s.set_facecolor("#1e1e1e")

                    ax_s.axvspan(-3, 0, alpha=0.3, label="Hidrofílico")
                    ax_s.axvspan(0, 3, alpha=0.15, label="Moderado")
                    ax_s.axvspan(3, 6, alpha=0.3, label="Lipofílico")

                    ax_s.axvline(
                        x=logp_val,
                        color="#ff4b4b",
                        linewidth=3,
                        linestyle="--"
                    )

                    ax_s.set_xlim(-3, 6)
                    ax_s.set_xlabel("Escala LogP", color="white")
                    ax_s.get_yaxis().set_visible(False)
                    ax_s.tick_params(colors="white")

                    st.pyplot(fig_s)
                    plt.close(fig_s)

                st.info(
                    "LogP no permite calcular porcentajes reales "
                    "de solubilidad. Para ello se necesitan datos "
                    "experimentales o modelos específicos."
                )

            # --------------------------------------------------
            # REACCIONES
            # --------------------------------------------------

            with sub_quimica[4]:

                st.subheader("⚗️ Representación de reacciones químicas")

                tipo_reaccion = st.selectbox(
                    "Elige una reacción:",
                    [
                        "Combustión completa (con O₂)",
                        "Hidrogenación catalítica"
                    ]
                )

                if st.button(
                    "Mostrar representación",
                    key="btn_reaccion"
                ):

                    if tipo_reaccion == "Combustión completa (con O₂)":

                        st.success("Representación conceptual")

                        st.write(
                            "Productos generales: CO₂ + H₂O"
                        )

                    else:

                        st.info(
                            "La hidrogenación es una reacción "
                            "de adición de hidrógeno a enlaces "
                            "insaturados. Esta versión no simula "
                            "la transformación molecular."
                        )

            # --------------------------------------------------
            # ESTEQUIOMETRÍA
            # --------------------------------------------------

            with sub_quimica[5]:

                st.subheader("⚖️ Calculadora de estequiometría")

                masa_molecular = rdMolDescriptors.CalcExactMolWt(
                    mol_with_hs
                )

                gramos = st.number_input(
                    "Masa en gramos (g):",
                    min_value=0.01,
                    value=10.0,
                    step=0.1
                )

                moles = gramos / masa_molecular

                st.metric(
                    "Cantidad de sustancia",
                    f"{moles:.4f} mol"
                )

    except Exception as e:

        st.error(
            f"Error en el laboratorio de química: {e}"
        )


# ==========================================================
# 2. LABORATORIO DE COMPUTACIÓN CUÁNTICA
# ==========================================================

with tab_cuantica:

    st.header("⚛️ Computación Cuántica")

    st.write(
        "Constructor visual de circuitos cuánticos utilizando Qiskit."
    )

    num_qubits = st.sidebar.slider(
        "Número de cúbits:",
        min_value=1,
        max_value=5,
        value=3
    )

    st.sidebar.subheader("Añadir compuertas")

    agregar_h = st.sidebar.checkbox(
        "Añadir compuertas Hadamard (H)"
    )

    agregar_x = st.sidebar.checkbox(
        "Añadir compuertas Pauli-X"
    )

    agregar_cnot = st.sidebar.checkbox(
        "Añadir compuerta CNOT",
        disabled=(num_qubits < 2)
    )

    agregar_mediciones = st.sidebar.checkbox(
        "Añadir mediciones"
    )

    qc = QuantumCircuit(num_qubits)

    if agregar_h:

        qubits_h = st.sidebar.multiselect(
            "Cúbits para H:",
            list(range(num_qubits)),
            default=[0]
        )

        for q in qubits_h:
            qc.h(q)

    if agregar_x:

        qubits_x = st.sidebar.multiselect(
            "Cúbits para X:",
            list(range(num_qubits))
        )

        for q in qubits_x:
            qc.x(q)

    if agregar_cnot and num_qubits > 1:

        control = st.sidebar.selectbox(
            "Cúbit de control:",
            list(range(num_qubits))
        )

        objetivos = [
            q for q in range(num_qubits)
            if q != control
        ]

        target = st.sidebar.selectbox(
            "Cúbit objetivo:",
            objetivos
        )

        qc.cx(control, target)

    if agregar_mediciones:
        qc.measure_all()

    st.subheader("Circuito resultante")

    try:

        fig = qc.draw(
            output="mpl",
            style="iqp"
        )

        st.pyplot(fig)
        plt.close(fig)

    except Exception as e:

        st.error(
            f"Error al dibujar el circuito: {e}"
        )

        st.info(
            "Comprueba las dependencias de visualización "
            "de Qiskit."
        )

    st.divider()

    st.info(
        "Esta herramienta construye y visualiza circuitos. "
        "Todavía no ejecuta simulaciones cuánticas."
    )


# ==========================================================
# 3. LABORATORIO DE BIOLOGÍA COMPUTACIONAL
# ==========================================================

with tab_biologia:

    st.header("🧬 BioChem Analyzer")

    st.write(
        "Evaluación computacional y documental "
        "de propiedades biológicas moleculares."
    )

    # ------------------------------------------------------
    # IDENTIFICACIÓN MOLECULAR
    # ------------------------------------------------------

    # Se elimina

    # ------------------------------------------------------
    # ADME
    # ------------------------------------------------------

    st.subheader("🧪 Análisis ADME")

    if st.button("Ejecutar análisis ADME", key="btn_adme"):

        try:

            resultado_adme = analizar_adme(smiles_input)

            if resultado_adme["estado"] == "correcto":

                propiedades = resultado_adme["propiedades"]

                st.markdown("### Propiedades moleculares")

                col1, col2, col3 = st.columns(3)

                col1.metric(
                    "Peso molecular",
                    propiedades["peso_molecular"]
                )

                col2.metric(
                    "LogP",
                    propiedades["logp"]
                )

                col3.metric(
                    "TPSA",
                    propiedades["tpsa"]
                )

                col4, col5, col6 = st.columns(3)

                col4.metric(
                    "Donadores H",
                    propiedades["donadores_h"]
                )

                col5.metric(
                    "Aceptores H",
                    propiedades["aceptores_h"]
                )

                col6.metric(
                    "Enlaces rotables",
                    propiedades["enlaces_rotables"]
                )

                st.divider()

                st.markdown("### Evaluación de reglas")

                lipinski = resultado_adme["lipinski"]
                veber = resultado_adme["veber"]

                st.write(
                    "**Lipinski:**",
                    "Cumple" if lipinski["cumple"]
                    else "No cumple todos los criterios"
                )

                st.write(
                    "Incumplimientos:",
                    lipinski["incumplimientos"]
                )

                st.write(
                    "**Veber:**",
                    "Cumple" if veber["cumple"]
                    else "No cumple todos los criterios"
                )

                st.divider()

                st.markdown("### Interpretación")

                st.write(
                    resultado_adme["interpretacion"]["polaridad"]
                )

                st.write(
                    resultado_adme["interpretacion"]["lipofilicidad"]
                )

                with st.expander("Ver criterios y limitaciones"):

                    st.write("Criterios de Lipinski")
                    st.json(lipinski["criterios"])

                    st.write("Criterios de Veber")
                    st.json(veber["criterios"])

                    for limitacion in resultado_adme["limitaciones"]:
                        st.write("•", limitacion)

                st.caption(
                    "Análisis fisicoquímico educativo. "
                    "No constituye una predicción clínica "
                    "de absorción, distribución, metabolismo, "
                    "excreción o seguridad."
                )

            else:

                st.error(resultado_adme["mensaje"])

        except (ValueError, TypeError, KeyError) as e:

            st.error(
                f"Error en el análisis ADME: {e}"
            )

    st.divider()

    # Próximos módulos
    st.subheader("Próximos módulos")

    st.write("🔬 Toxicología computacional — En desarrollo")
    st.write("🎯 Bioactividad molecular — En desarrollo")
    st.write("📚 Literatura científica — En desarrollo")
    st.write("📄 Generación de reportes — En desarrollo")


# ==========================================================
# PIE DE PÁGINA
# ==========================================================

st.sidebar.divider()

st.sidebar.caption(
    "QuantumChemLabs | Plataforma educativa de química "
    "y ciencias computacionales."
)
#Codigo mejorado y completo
