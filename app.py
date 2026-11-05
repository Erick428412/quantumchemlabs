import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from qiskit import QuantumCircuit

# IMPORTAMOS EL CEREBRO QUÍMICO QUE ACABAMOS DE CREAR
import quimica 

# Configuración de la página
st.set_page_config(page_title="Quantum Chemistry Labs", layout="wide")
st.title("Bienvenido a QuantumChemLabs")

def render_mol_local(mol, width=600, height=350):
    svg = quimica.get_mol_svg(mol, width, height)
    html_content = f"""
    <div style="display: flex; justify-content: center; align-items: center; background-color: #0e1117; padding: 10px; border-radius: 8px;">
        {svg}
    </div>
    """
    components.html(html_content, height=height + 20)

# Barra lateral
st.sidebar.header("Entrada de la Molécula")
smiles_input = st.sidebar.text_input("Ingresa el código SMILES base:", value="C(O)=C(C)O")

try:
    isomeros_lista = quimica.generar_isomeros(smiles_input)
    if isomeros_lista is None:
        st.error("Código SMILES no válido. Ingresa una estructura válida.")
    else:
        if "isomero_seleccionado" not in st.session_state:
            st.session_state.isomero_seleccionado = 0
            
        with st.expander("Isómeros Químicos Detectados", expanded=True):
            st.write(f"Se generaron **{len(isomeros_lista)}** configuraciones moleculares:")
            cols = st.columns(min(4, len(isomeros_lista)))
            for idx, iso in enumerate(isomeros_lista):
                col_idx = idx % 4
                with cols[col_idx]:
                    svg_iso = quimica.get_mol_svg(iso, width=220, height=160)
                    st.components.v1.html(f"<div style='text-align:center;'>{svg_iso}</div>", height=170)
                    if st.button(f"Inspeccionar Isómero {idx + 1}", key=f"btn_iso_{idx}"):
                        st.session_state.isomero_seleccionado = idx
                        st.slots = []
                        st.rerun()

        mol_activa = isomeros_lista[st.session_state.isomero_seleccionado]
        mol_with_hs = Chem.AddHs(mol_activa)
        
        st.subheader(f"🌐 Inspección de la Estructura Seleccionada (Isómero {st.session_state.isomero_seleccionado + 1})")
        render_mol_local(mol_with_hs)
        
        st.markdown("---")
        
        # Creación y mapeo de pestañas
        tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
            "📊 Información e Isomería", 
            "🌡️ Propiedades Físicas y Viscosidad",
            "🧭 Polaridad y Polos",
            "💧 Solubilidad y Solventes",
            "⚗️ Simulador de Reacciones", 
            "⚖️ Estequiometría",
            "🧲 Cuantica"
        ])
        
        # --- PESTAÑA 1: INFORMACIÓN E ISOMERÍA ---
        with tab1:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Datos Generales")
                st.write(f"**Fórmula Química:** {Chem.rdMolDescriptors.CalcMolFormula(mol_with_hs)}")
                st.write(f"**Peso Molecular:** {Chem.rdMolDescriptors.CalcExactMolWt(mol_with_hs):.4f} g/mol")
                st.write(f"**Átomos Totales:** {mol_with_hs.GetNumAtoms()}")
            with col2:
                st.subheader("Estereoquímica")
                centros = Chem.FindMolChiralCenters(mol_activa, includeUnassigned=True)
                st.write(f"**Centros Quirales en este Isómero:** {len(centros)}")

        # --- PESTAÑA 2: PROPIEDADES FÍSICAS Y VISCOSIDAD ---
        with tab2:
            st.subheader("Estimación de Propiedades Físicas")
            est_boiling_c, est_melting_c, mw, logp, h_donors = quimica.calcular_propiedades_fisicas(mol_activa)
            
            col_f1, col_f2 = st.columns(2)
            with col_f1: st.metric(label="Punto de Ebullición Estimado", value=f"{est_boiling_c:.1f} °C")
            with col_f2: st.metric(label="Punto de Fusión Estimado", value=f"{est_melting_c:.1f} °C")
                
            st.markdown("---")
            st.subheader("📈 Curva Dinámica de Viscosidad")
            t_min, t_max = float(max(0.0, est_melting_c + 5)), float(est_boiling_c - 5)
            temp_seleccionada_c = st.slider("🌡️ Temperatura de Trabajo (°C):", min_value=t_min, max_value=t_max, value=(t_min + t_max) / 2.0, step=1.0)
            
            temp_c_range, viscosidad_cps, E_visc, R, A_visc = quimica.calcular_perfil_viscosidad(mw, h_donors, t_min, t_max)
            visc_exacta = A_visc * np.exp(E_visc / (R * (temp_seleccionada_c + 273.15)))
            
            st.metric(label=f"Viscosidad a {temp_seleccionada_c:.1f} °C", value=f"{visc_exacta:.4f} cP")
            
            fig, ax = plt.subplots(figsize=(8, 3.5))
            fig.patch.set_facecolor('#0e1117')
            ax.set_facecolor('#1e1e1e')
            ax.plot(temp_c_range, viscosidad_cps, color='#ff4b4b', linewidth=2.5)
            ax.plot(temp_seleccionada_c, visc_exacta, marker='o', markersize=9, color='#00b4d8')
            ax.axvline(x=temp_seleccionada_c, color='#00b4d8', linestyle=':', alpha=0.6)
            ax.axhline(y=visc_exacta, color='#00b4d8', linestyle=':', alpha=0.6)
            ax.set_xlabel('Temperatura (°C)', color='white')
            ax.set_ylabel('Viscosidad (cP)', color='white')
            ax.tick_params(colors='white')
            ax.grid(True, linestyle='--', alpha=0.2, color='gray')
            st.pyplot(fig)

        # --- PESTAÑA 3: POLARIDAD Y POLOS ---
        with tab3:
            st.subheader("🧭 Análisis Vectorial de Polaridad")
            enlaces_p_cnt, enlaces_np_cnt, es_polar, pos_x, pos_y, idx_pos, idx_neg, c_pos_x, c_pos_y, c_neg_x, c_neg_y, mol_cargas = quimica.analizar_polaridad_y_polos(mol_with_hs)
            
            col_b1, col_b2 = st.columns(2)
            with col_b1: st.write(f"• Enlaces **Covalentes Polares**: `{enlaces_p_cnt}`")
            with col_b2: st.write(f"• Enlaces **Covalentes No Polares**: `{enlaces_np_cnt}`")
                
            st.markdown("---")
            st.write(f"¿La molécula es **Globalmente Polar**?: **{es_polar}**")
            
            fig_p, ax_p = plt.subplots(figsize=(7, 3.5))
            fig_p.patch.set_facecolor('#0e1117')
            ax_p.set_facecolor('#1e1e1e')
            ax_p.scatter(pos_x, pos_y, c='#4a5568', s=100, alpha=0.5)
            for i, atom in enumerate(mol_cargas.GetAtoms()):
                ax_p.text(pos_x[i]+0.12, pos_y[i]+0.12, atom.GetSymbol(), color='#a0aec0', fontsize=9)
            
            if np.any(idx_pos): ax_p.scatter(c_pos_x, c_pos_y, c='#ff4b4b', s=350, marker='P', label='Polo δ+')
            if np.any(idx_neg): ax_p.scatter(c_neg_x, c_neg_y, c='#00b4d8', s=350, marker='o', label='Polo δ-')
            
            if es_polar == "SÍ":
                ax_p.arrow(c_pos_x, c_pos_y, c_neg_x - c_pos_x, c_neg_y - c_pos_y, head_width=0.15, head_length=0.15, fc='#ffb703', ec='#ffb703', length_includes_head=True, linewidth=2.5)
                           
            ax_p.set_title('🗺️ Mapa de Polos Eléctricos', color='white', fontweight='bold')
            ax_p.tick_params(colors='white')
            ax_p.axis('equal')
            st.pyplot(fig_p)

        # --- PESTAÑA 4: SOLUBILIDAD Y SOLVENTES ---
        with tab4:
            st.subheader("💧 Predicción de Solubilidad en Sistemas de Solventes")
            _, _, _, logp_val, _ = quimica.calcular_propiedades_fisicas(mol_activa)
            
            if logp_val < 0:
                perfil_quimico = "Altamente Hidrofílico"
                desc_solubilidad = "Se disuelve de forma excelente en agua debido a interacciones dipolo-dipolo fuertes."
            elif 0 <= logp_val <= 3:
                perfil_quimico = "Anfipático / Polaridad Moderada"
                desc_solubilidad = "Solubilidad equilibrada en sistemas acuosos y orgánicos."
            else:
                perfil_quimico = "Altamente Lipofílico / Hidrofóbico"
                desc_solubilidad = "Insoluble o muy poco soluble en agua. Alta afinidad por solventes no polares."

            col_s1, col_s2 = st.columns(2)
            with col_s1:
                st.metric(label="Descriptor LogP (Crippen)", value=f"{logp_val:.3f}")
                st.write(f"**Perfil Químico:** {perfil_quimico}")
                st.write(f"**Comportamiento:** {desc_solubilidad}")
            
            with col_s2:
                fig_s, ax_s = plt.subplots(figsize=(6, 1.8))
                fig_s.patch.set_facecolor('#0e1117')
                ax_s.set_facecolor('#1e1e1e')
                ax_s.axvspan(-3, 0, color='#00b4d8', alpha=0.3, label='Hidrofílico')
                ax_s.axvspan(0, 3, color='#4cc9f0', alpha=0.15, label='Moderado')
                ax_s.axvspan(3, 6, color='#ffb703', alpha=0.3, label='Lipofílico')
                ax_s.axvline(x=logp_val, color='#ff4b4b', linewidth=3, linestyle='--')
                ax_s.plot(logp_val, 0, marker='v', color='#ff4b4b', markersize=12, label='Tu Molécula')
                ax_s.set_xlim(-3, 6)
                ax_s.set_ylim(-0.5, 0.5)
                ax_s.set_xlabel('Escala LogP', color='white', fontsize=9)
                ax_s.get_yaxis().set_visible(False)
                ax_s.tick_params(colors='white', labelsize=8)
                st.pyplot(fig_s)

            st.markdown("---")
            st.write("#### 🧪 Matriz de Solubilidad Estimada por Solvente")
            solv_agua = max(5, min(100, int((1 - (logp_val / 4.0)) * 100))) if logp_val >= 0 else 100
            solv_etanol = max(10, min(100, int((1 - (abs(logp_val - 0.5) / 5.0)) * 100)))
            solv_acetona = max(10, min(100, int((1 - (abs(logp_val - 1.5) / 5.0)) * 100)))
            solv_cloroformo = max(5, min(100, int((1 - (abs(logp_val - 2.5) / 5.0)) * 100)))
            solv_hexano = max(5, min(100, int((logp_val / 5.0) * 100))) if logp_val >= 0 else 5

            col_sol1, col_sol2, col_sol3, col_sol4, col_sol5 = st.columns(5)
            for col, name, val in zip([col_sol1, col_sol2, col_sol3, col_sol4, col_sol5], ["Agua", "Etanol", "Acetona", "Cloroformo", "Hexano"], [solv_agua, solv_etanol, solv_acetona, solv_cloroformo, solv_hexano]):
                with col:
                    st.write(f"{name}")
                    st.progress(val / 100.0)
                    st.caption(f"{val}%")
        # --- PESTAÑA 5: REACCIONES ---
        with tab5:
            st.subheader("Predecir Reacción Química")
            tipo_reaccion = st.selectbox("Elige la reacción:", ["Combustión Completa (con O₂)", "Hidrogenación Catalítica"])
            if st.button("Ejecutar Reacción"):
                if tipo_reaccion == "Combustión Completa (con O₂)":
                    st.success("🔥 ¡Combustión simulada con éxito!")
                    st.write("Productos: $CO_2$ + $H_2O$")
            else:
                st.info("Buscando enlaces insaturados en este isómero específico...")
        # --- PESTAÑA 6: ESTEQUIOMETRÍA ---
        with tab6:
            st.subheader("Calculadora de Laboratorio")
            gramos = st.number_input("Masa en gramos (g):", min_value=0.01, value=10.0)
            st.metric(label="Moles Resultantes", value=f"{gramos / Chem.rdMolDescriptors.CalcExactMolWt(mol_with_hs):.4f} mol")
        with tab7:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Creador de circuitos cuanticos")
                st.write("Crea tu circuito cuantico")
                num_qubits = st.sidebar.slider("Número de Cúbits", min_value=1, max_value=5, value=3)

                # Inicializar el circuito con los cúbits seleccionados
                qc = QuantumCircuit(num_qubits)

                # 2. Controles interactivos en la barra lateral para añadir compuertas
                st.sidebar.subheader("Añadir Compuertas")

                # Compuertas Hadamard (H)
                if st.sidebar.checkbox("Añadir compuertas Hadamard (H)"):
                    qubits_h = st.sidebar.multiselect("¿En qué cúbits aplicar H?", list(range(num_qubits)), default=[0])
                    for q in qubits_h:
                        qc.h(q)

                # Compuertas Pauli-X (NOT)
                if st.sidebar.checkbox("Añadir compuertas Pauli-X (NOT)"):
                    qubits_x = st.sidebar.multiselect("¿En qué cúbits aplicar X?", list(range(num_qubits)))
                    for q in qubits_x:
                        qc.x(q)

                # Compuertas CNOT (Controlled-NOT)
                if num_qubits > 1 and st.sidebar.checkbox("Añadir compuerta CNOT"):
                    control = st.sidebar.selectbox("Cúbit de Control", list(range(num_qubits)), index=0)
                    target = st.sidebar.selectbox("Cúbit Objetivo", [q for q in range(num_qubits) if q != control], index=0)
                    qc.cx(control, target)

                # Medición final
                if st.sidebar.checkbox("Añadir mediciones al final"):
                    qc.measure_all()

                # 3. Renderizar el circuito en la pantalla principal
                st.subheader("Circuito Resultante")

                try:
                    # Generamos la figura de Matplotlib
                    fig = qc.draw(output='mpl', style='iqp')
                    # La mostramos en Streamlit
                    st.pyplot(fig)
                except Exception as e:
                    st.error(f"Hubo un error al dibujar el circuito: {e}")
                    st.info("Asegúrate de tener instalado el paquete de visualización: `pip install qiskit[visualization]`")  
            with col2:
                st.write("Software creado por Erick Lopez")
except Exception as e:
    st.error(f"Error general en el renderizado de la interfaz: {e}")
