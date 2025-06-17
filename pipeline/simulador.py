import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle

# Carrega modelo pronto do disco
def carregar_modelo():
    with open("modelo_gambiarra.pkl", "rb") as f:
        return pickle.load(f)

def calcular_novo_upper_bound(original_ub, Ki, conc):
    if Ki <= 0:
        return 0
    return original_ub * (Ki / (Ki + conc))

def simular_crescimento_streamlit(kd_uM, comparar=False):
    model = carregar_modelo()
    girase = model.reactions.get_by_id("DNA_GIRASE")
    original_ub = girase.upper_bound

    concs_uM = np.linspace(0.1, 100, 20)
    concs_M = concs_uM * 1e-6

    resultados = []
    dados_csv = []

    inibidores = [("Peptídeo", kd_uM, "blue")]
    if comparar:
        inibidores.append(("Ácido Nalidíxico", 1.3, "green"))

    for nome, Ki_uM, cor in inibidores:
        taxas = []
        for conc_uM, conc_M in zip(concs_uM, concs_M):
            modelo_temp = model.copy()
            Ki_M = Ki_uM * 1e-6
            novo_ub = calcular_novo_upper_bound(original_ub, Ki_M, conc_M)
            modelo_temp.reactions.get_by_id("DNA_GIRASE").upper_bound = novo_ub
            taxa = modelo_temp.slim_optimize()
            taxa = taxa if taxa is not None else 0
            taxas.append(taxa)
            dados_csv.append({
                "Peptídeo": nome,
                "Concentração (µM)": conc_uM,
                "Crescimento (h⁻¹)": taxa
            })

        # Cálculo IC50 corrigido
        taxa_controle = taxas[0]
        meia_taxa = taxa_controle * 0.5

        ic50_val = None
        ic50_label = "Não detectado"

        for i in range(1, len(taxas)):
            t1, t2 = taxas[i-1], taxas[i]
            if (t1 >= meia_taxa >= t2) or (t1 <= meia_taxa <= t2):
                x1, x2 = concs_uM[i-1], concs_uM[i]
                if t2 != t1:  # evitar divisão por zero
                    ic50_val = x1 + (meia_taxa - t1) * (x2 - x1) / (t2 - t1)
                    ic50_label = f"{ic50_val:.2f} µM"
                break

        for linha in dados_csv:
            if linha["Peptídeo"] == nome:
                linha["IC50 estimada (µM)"] = ic50_label

        resultados.append((nome, taxas, cor, ic50_val, ic50_label))

    df_resultado = pd.DataFrame(dados_csv)

    # Gráfico
    fig, ax = plt.subplots(figsize=(8, 6))
    for nome, taxas, cor, ic50_val, ic50_label in resultados:
        ax.plot(concs_uM, taxas, marker='o', label=f"{nome} (IC₅₀ ≈ {ic50_label})", color=cor)
        if ic50_val:
            ax.axvline(ic50_val, linestyle="--", color=cor, alpha=0.5)

    ax.set_xlabel("[Inibidor] (µM)")
    ax.set_ylabel("Crescimento simulado (h⁻¹)")
    ax.set_title("Curva de crescimento com inibição da DNA-girase")
    ax.legend()
    ax.grid(True)

    return fig, df_resultado
