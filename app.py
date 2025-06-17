import os
os.environ["PORT"] = os.getenv("PORT", "8501")

import streamlit as st
import matplotlib.pyplot as plt
import io
import pandas as pd

from pipeline.propriedades import analisar_peptideo
from pipeline.simulador import simular_crescimento_streamlit

st.set_page_config(page_title="Análise de Peptídeos", layout="centered")
st.title("🧬 Análise de Peptídeos Antibacterianos")

if 'page' not in st.session_state:
    st.session_state['page'] = 'input'
if 'figura' not in st.session_state:
    st.session_state['figura'] = None
if 'dados_simulacao' not in st.session_state:
    st.session_state['dados_simulacao'] = None

if st.session_state['page'] == 'input':
    seq = st.text_input("Digite a sequência do peptídeo (ex: KLFKFFKFFK):")

    if seq:
        try:
            props = analisar_peptideo(seq)

            st.subheader("🔬 Propriedades físico-químicas")
            st.write(f"**Carga líquida:** {props['carga']}")
            st.write(f"**Hidrofobicidade média:** {props['hidrofobicidade']}")
            st.write(f"**Estabilidade extracelular:** {props['estabilidade_extracelular']}")
            st.write(f"**Permeabilidade membrana externa:** {props['permeabilidade_membrana_externa']}")
            st.write(f"**Translocação membrana citoplasmática:** {props['translocacao_membrana_citoplasmatica']}")
            st.write(f"**Pontuação combinada (0 a 1):** {props['pontuacao_combinada']}")

            st.subheader("🧪 Parâmetro de Docking")
            kd_M = st.number_input("Informe o Kd (em M):", min_value=1e-12, format="%.2e")
            kd_uM = kd_M * 1e6
            st.write(f"🔁 Isso equivale a **{kd_uM:.2f} µM**")

            comparar = st.checkbox("Comparar com ácido nalidíxico (Kd = 1.3 µM)?")

            if st.button("Simular impacto metabólico"):
                if kd_M == 0.0:
                    st.warning("Por favor, insira um valor de Kd maior que zero.")
                else:
                    with st.spinner("🦫 Capivara calculando o impacto..."):
                        try:
                            st.video(open("capivara.mp4", "rb").read())
                            fig, df_resultado = simular_crescimento_streamlit(kd_uM, comparar)
                            st.session_state['figura'] = fig
                            st.session_state['dados_simulacao'] = df_resultado
                            st.session_state['page'] = 'resultado'
                            st.experimental_rerun()
                        except Exception as e:
                            st.error(f"❌ Erro durante a simulação:\n\n{e}")
        except Exception as e:
            st.error(f"❌ Erro na análise do peptídeo:\n\n{e}")

elif st.session_state['page'] == 'resultado':
    st.subheader("📈 Simulação do crescimento bacteriano")
    fig = st.session_state['figura']
    st.pyplot(fig)

    # Exportar gráfico
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    buf.seek(0)
    st.download_button("📥 Baixar gráfico", data=buf, file_name="grafico_simulacao.png", mime="image/png")

    # Exportar CSV
    if st.session_state['dados_simulacao'] is not None:
        csv_buf = io.StringIO()
        st.session_state['dados_simulacao'].to_csv(csv_buf, index=False)
        st.download_button(
            label="📥 Baixar dados (.csv)",
            data=csv_buf.getvalue(),
            file_name="dados_simulacao.csv",
            mime="text/csv"
        )

    if st.button("⬅️ Voltar"):
        st.session_state['page'] = 'input'
        st.session_state['figura'] = None
        st.session_state['dados_simulacao'] = None
        st.experimental_rerun()
