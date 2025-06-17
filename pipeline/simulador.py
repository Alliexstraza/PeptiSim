import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite

def carregar_modelo():
    url = "http://bigg.ucsd.edu/static/models/iML1515.xml"
    filename = "iML1515.xml"
    try:
        model = read_sbml_model(filename)
    except:
        import requests
        r = requests.get(url)
        with open(filename, "wb") as f:
            f.write(r.content)
        model = read_sbml_model(filename)
    return model

def preparar_modelo(model):
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = -10
    model.reactions.get_by_id("EX_o2_e").lower_bound = -20

    try:
        dna_sc = model.metabolites.get_by_id("dna_supercoiling_c")
    except KeyError:
        dna_sc = Metabolite("dna_supercoiling_c", name="DNA Supercoiling (pseudo)", compartment="c")
        model.add_metabolites([dna_sc])

    if "DNA_GIRASE" not in model.reactions:
        girase = Reaction("DNA_GIRASE")
        girase.name = "Simulated DNA Gyrase Activity"
        girase.lower_bound = 0
        girase.upper_bound = 1000
        girase.add_metabolites({dna_sc: 1})
        model.add_reactions([girase])

    if "D_BIOMASS_GYRASE" not in model.reactions:
        biomassa_d = Reaction("D_BIOMASS_GYRASE")
        biomassa_d.name = "Biomass demand for Gyrase Product"
        biomassa_d.lower_bound = 0
        biomassa_d.upper_bound = 1000
        biomassa_d.add_metabolites({dna_sc: -1})
        model.add_reactions([biomassa_d])

    model.objective = "D_BIOMASS_GYRASE"
    return model

def calcular_novo_upper_bound(original_ub, Ki, conc):
    if Ki <= 0:
        return 0
    return original_ub * (Ki / (Ki + conc))

def simular_crescimento_streamlit(kd_uM, comparar=False):
    model = preparar_modelo(carregar_modelo())
    girase = model.reactions.get_by_id("DNA_GIRASE")
    original_ub = girase.upper_bound

    concs_uM = np.linspace(0.1, 100, 20)
    concs_M = concs_uM * 1e-6

    resultados = []

    inibidores = [("Peptídeo", kd_uM, "blue")]
    if comparar:
        inibidores.append(("Ácido Nalidíxico", 1.3, "green"))

    for nome, Ki_uM, cor in inibidores:
        taxas = []
        for conc in concs_M:
            modelo_temp = model.copy()
            Ki_M = Ki_uM * 1e-6
            novo_ub = calcular_novo_upper_bound(original_ub, Ki_M, conc)
            modelo_temp.reactions.get_by_id("DNA_GIRASE").upper_bound = novo_ub
            taxa = modelo_temp.slim_optimize()
            taxas.append(taxa if taxa is not None else 0)

        taxa_controle = taxas[0]
        meia_taxa = taxa_controle * 0.5
        try:
            ic50_val = np.interp(meia_taxa, taxas[::-1], concs_uM[::-1])
            ic50_label = f"{ic50_val:.2f} µM"
        except:
            ic50_val = None
            ic50_label = "Não detectado"

        resultados.append((nome, taxas, cor, ic50_val, ic50_label))

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
        # Monta DataFrame com os dados
    df_lista = []
    for nome, taxas, _, ic50_val, _ in resultados:
        for conc, taxa in zip(concs_uM, taxas):
            df_lista.append({
                "Peptídeo": nome,
                "Concentração (µM)": conc,
                "Crescimento (h⁻¹)": taxa,
                "IC50 estimada (µM)": round(ic50_val, 2) if ic50_val else "N/A"
            })
    df_resultado = pd.DataFrame(df_lista)

    return fig, df_resultado

    return fig
