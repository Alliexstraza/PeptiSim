import os
os.environ["PORT"] = os.getenv("PORT", "8501")

import streamlit as st
import io
import pandas as pd

# Importa as funções dos seus módulos da pipeline
from pipeline.propriedades import analisar_peptideo
from pipeline.simulador import simular_crescimento_streamlit
# --- NOSSA NOVA IMPORTAÇÃO ---
from pipeline.membrana import analisar_interacao_pdb

# Configuração da página
st.set_page_config(page_title="PeptiTools", layout="wide")

# Inicialização do estado da sessão para todas as variáveis
if 'page' not in st.session_state:
    st.session_state['page'] = 'home'
if 'figura' not in st.session_state:
    st.session_state['figura'] = None
if 'dados_simulacao' not in st.session_state:
    st.session_state['dados_simulacao'] = None
if 'graficos_membrana' not in st.session_state:
    st.session_state['graficos_membrana'] = None

# ==========================
# MENU LATERAL (SIDEBAR)
# ==========================
st.sidebar.title("Navegação")
# --- Menu com a nova página integrada ---
pagina_selecionada = st.sidebar.radio(
    "Escolha a Análise:",
    ["🏠 Home", "🧬 Análise de Sequência", "🔬 Análise de Membrana", "❓ Ajuda", "ℹ️ Sobre"]
)

# Mapeia a seleção do rádio para o nome interno da página
mapa_paginas = {
    "🏠 Home": "home",
    "🧬 Análise de Sequência": "simulador", # Mantive o nome interno 'simulador'
    "🔬 Análise de Membrana": "membrana",
    "❓ Ajuda": "ajuda",
    "ℹ️ Sobre": "sobre"
}
# Atualiza a página no estado da sessão sem causar um rerun imediato
if st.session_state['page'] != mapa_paginas[pagina_selecionada]:
    st.session_state['page'] = mapa_paginas[pagina_selecionada]
    # Limpa os resultados antigos ao trocar de página de análise
    st.session_state['figura'] = None
    st.session_state['dados_simulacao'] = None
    st.session_state['graficos_membrana'] = None
    st.experimental_rerun()

# ==========================
# FUNÇÃO DE RODAPÉ
# ==========================
def mostrar_rodape():
    st.markdown("---")
    st.markdown("### 🔧 Apoio institucional:")
    # Use st.columns para evitar que o rodapé seja muito largo no layout "wide"
    gap1, col1, col2, col3, gap2 = st.columns([1, 0.5, 0.5, 0.5, 1])
    with col1:
        st.image("logo.png", width=90, caption="Laboratório")
    with col2:
        st.image("ufms.png", width=90, caption="UFMS")
    with col3:
        st.image("PMBqBM.png", width=90, caption="PMBqBM")

# ====================================================================
# LÓGICA DAS PÁGINAS
# ====================================================================

# --- PÁGINA HOME ---
if st.session_state['page'] == 'home':
    if os.path.exists("logo.png"):
        st.image("logo.png", width=120)
    st.title("PeptiTools: Plataforma Integrada de Análise de Peptídeos")
    st.markdown("""
    Bem-vindo ao **PeptiTools**! Uma aplicação para acelerar a descoberta e análise de peptídeos com potencial antimicrobiano.
    
    ### Ferramentas Disponíveis:
    Use o menu na barra lateral para navegar entre as análises:
    
    - **🧬 Análise de Sequência:** Calcule propriedades físico-químicas e simule o impacto da inibição de um alvo específico (DNA-girase) no crescimento da *E. coli*.
    - **🔬 Análise de Membrana:** Faça o upload de uma estrutura 3D (PDB) e visualize a energia de interação do seu peptídeo com um modelo de membrana.
    """)
    if st.button("🔬 Iniciar Análise"):
        # Navega para a primeira página de análise
        st.session_state['page'] = 'simulador'
        st.experimental_rerun()
    mostrar_rodape()

# --- PÁGINA DE ANÁLISE DE SEQUÊNCIA (O SEU SIMULADOR ORIGINAL) ---
elif st.session_state['page'] == 'simulador':
    st.title("🧬 Análise de Sequência e Inibição de Alvo")
    seq = st.text_input("Digite a sequência do peptídeo (ex: KLFKFFKFFK):", key="seq_input")
    if seq:
        try:
            props = analisar_peptideo(seq)
            st.subheader("📊 Propriedades físico-químicas")
            st.write(f"**Carga líquida:** {props['carga']}")
            st.write(f"**Hidrofobicidade média:** {props['hidrofobicidade']}")
            st.write(f"**Estabilidade extracelular:** {props['estabilidade_extracelular']}")
            st.write(f"**Permeabilidade membrana externa:** {props['permeabilidade_membrana_externa']}")
            st.write(f"**Translocação membrana citoplasmática:** {props['translocacao_membrana_citoplasmatica']}")
            st.write(f"**Pontuação combinada (0 a 1):** {props['pontuacao_combinada']}")

            st.subheader("🧪 Parâmetro de Docking (vs. DNA-girase)")
            kd_M = st.number_input("Informe o Kd (em M):", min_value=1e-12, format="%.2e", key="kd_input")
            kd_uM = kd_M * 1e6
            st.write(f"Isso equivale a **{kd_uM:.2f} µM**")
            comparar = st.checkbox("Comparar com ácido nalidíxico (Kd = 1.3 µM)?")

            if st.button("Simular impacto metabólico"):
                if kd_M == 0.0:
                    st.warning("Por favor, insira um valor de Kd maior que zero.")
                else:
                    with st.spinner("⏳ Calculando..."):
                        try:
                            fig, df_resultado = simular_crescimento_streamlit(kd_uM, comparar)
                            st.session_state['figura'] = fig
                            st.session_state['dados_simulacao'] = df_resultado
                            st.session_state['page'] = 'resultado'
                            st.experimental_rerun()
                        except Exception as e:
                            st.error(f"❌ Erro durante a simulação:\n\n{e}")
        except Exception as e:
            st.error(f"❌ Erro na análise do peptídeo:\n\n{e}")
    mostrar_rodape()

# --- PÁGINA DE ANÁLISE DE MEMBRANA (A NOSSA NOVA PÁGINA) ---
elif st.session_state['page'] == 'membrana':
    st.title("🔬 Análise de Interação com a Membrana")
    st.markdown("Faça o upload de uma estrutura 3D de peptídeo (`.pdb`) para calcular e visualizar a energia de interação com um modelo de membrana simplificado.")
    
    pdb_file = st.file_uploader("Escolha o seu arquivo PDB", type=["pdb", "ent"], key="pdb_uploader")

    if pdb_file is not None:
        if st.button("Analisar Interação com Membrana"):
            with st.spinner("Calculando energias... Por favor, aguarde."):
                try:
                    # Chama nossa nova função do módulo membrana.py
                    g1_path, g2_path = analisar_interacao_pdb(pdb_file)
                    st.session_state['graficos_membrana'] = [g1_path, g2_path]
                except Exception as e:
                    st.error(f"❌ Erro durante a análise da membrana:\n\n{e}")

    # Mostra os gráficos se eles existirem no estado da sessão
    if st.session_state['graficos_membrana']:
        st.subheader("Resultados da Análise de Energia")
        st.image(st.session_state['graficos_membrana'][0], caption="Energia por Átomo vs Profundidade na Membrana", use_column_width=True)
        st.image(st.session_state['graficos_membrana'][1], caption="Energia Média por Resíduo", use_column_width=True)
        
        if st.button("Analisar Outro Peptídeo"):
            st.session_state['graficos_membrana'] = None
            st.experimental_rerun()
            
    mostrar_rodape()

# --- PÁGINA DE RESULTADO (da simulação original) ---
elif st.session_state['page'] == 'resultado':
    st.subheader("📈 Simulação do crescimento bacteriano")
    fig = st.session_state['figura']
    if fig:
        st.pyplot(fig)
        buf = io.BytesIO()
        fig.savefig(buf, format="png")
        buf.seek(0)
        st.download_button("📥 Baixar gráfico", data=buf, file_name="grafico_simulacao.png", mime="image/png")

    if st.session_state['dados_simulacao'] is not None:
        csv_buf = io.StringIO()
        st.session_state['dados_simulacao'].to_csv(csv_buf, index=False)
        st.download_button(label="📥 Baixar dados (.csv)", data=csv_buf.getvalue(), file_name="dados_simulacao.csv", mime="text/csv")

    if st.button("⬅️ Voltar ao Simulador"):
        st.session_state['page'] = 'simulador'
        st.session_state['figura'] = None
        st.session_state['dados_simulacao'] = None
        st.experimental_rerun()
    mostrar_rodape()

# --- PÁGINA DE AJUDA ---
   elif st.session_state['page'] == 'ajuda':
    st.title("❓ Como usar")
    st.markdown("""
    ### Etapas:
    1. Digite a sequência do peptídeo (ex: KLFKFFKFFK)
    2. Informe o valor de Kd (constante de dissociação), que pode ser obtido via         [**PRODIGY**](https://bianca.science.uu.nl/prodigy/):
        - Faça docking entre o peptídeo e a DNA-girase
        - Submeta a estrutura ao PRODIGY para estimar o Kd

    3. Marque se deseja comparar com o ácido nalidíxico
    4. Clique em **Simular impacto metabólico**

    ### Sobre as pontuações:
    - **Pontuação combinada:** Integra estabilidade extracelular, permeabilidade e translocação
    - Valores mais altos (próximos de 1) indicam maior potencial antimicrobiano
    """)

    mostrar_rodape()

# --- PÁGINA SOBRE ---
   elif st.session_state['page'] == 'sobre':
    st.title("ℹ️ Sobre o projeto")
    st.markdown("""
    Desenvolvido por Jéssica Carretone  
    Sob orientação do Prof. Dr. Malson Lucena

    Esta aplicação integra bioinformática estrutural e simulações metabólicas
    para prever o potencial antimicrobiano de peptídeos.

    ### Tecnologias utilizadas:
    - Python, Streamlit
    - Modelagem metabólica com **COBRApy** e o modelo **iML1515**
    - Análises físico-químicas e docking molecular

    Código-fonte disponível em breve no GitHub.
    """)
    mostrar_rodape()

