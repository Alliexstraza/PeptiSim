import numpy as np
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
import copy

# ================================
# 1. Carregar e preparar o modelo
# ================================

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

    # Metabólitos
    atp = model.metabolites.get_by_id("atp_c")
    adp = model.metabolites.get_by_id("adp_c")

    dna_in = model.metabolites.get("dna_substrato_c") or Metabolite("dna_substrato_c", name="DNA Substrato (pseudo)", compartment="c")
    dna_out = model.metabolites.get("dna_supercoiling_c") or Metabolite("dna_supercoiling_c", name="DNA Supercoiling (pseudo)", compartment="c")
    inibidor = model.metabolites.get("inh_peptideo_c") or Metabolite("inh_peptideo_c", name="Peptídeo Inibidor", compartment="c")
    inh_ext = model.metabolites.get("inh_peptideo_e") or Metabolite("inh_peptideo_e", name="Peptídeo Inibidor (ext)", compartment="e")
    complexo = model.metabolites.get("girase_inativa_c") or Metabolite("girase_inativa_c", name="Complexo Inativo", compartment="c")

    for m in [dna_in, dna_out, inibidor, inh_ext, complexo]:
        if m.id not in model.metabolites:
            model.add_metabolites([m])

    if "DNA_GIRASE" not in model.reactions:
        girase = Reaction("DNA_GIRASE")
        girase.name = "Atividade da DNA-girase"
        girase.lower_bound = 0
        girase.upper_bound = 1000
        girase.add_metabolites({
            atp: -2,
            dna_in: -1,
            adp: 2,
            dna_out: 1
        })
        model.add_reactions([girase])

    if "EX_inh_peptideo_e" not in model.reactions:
        ex_inh = Reaction("EX_inh_peptideo_e")
        ex_inh.name = "Troca do inibidor (entrada)"
        ex_inh.lower_bound = 0
        ex_inh.upper_bound = 1000
        ex_inh.add_metabolites({inh_ext: -1})
        model.add_reactions([ex_inh])

    if "TRANS_inh_peptideo" not in model.reactions:
        trans = Reaction("TRANS_inh_peptideo")
        trans.name = "Transporte do inibidor"
        trans.lower_bound = 0
        trans.upper_bound = 1000
        trans.add_metabolites({inh_ext: -1, inibidor: 1})
        model.add_reactions([trans])

    if "SEQUESTRO_GIRASE" not in model.reactions:
        sequestro = Reaction("SEQUESTRO_GIRASE")
        sequestro.name = "Ligação do inibidor à DNA-girase"
        sequestro.lower_bound = 0
        sequestro.upper_bound = 1000
        sequestro.add_metabolites({dna_in: -1, inibidor: -1, complexo: 1})
        model.add_reactions([sequestro])

    biomassa = model.reactions.get_by_id("BIOMASS_Ec_iML1515_core_53p95M")
    if dna_out not in biomassa.metabolites:
        biomassa.add_metabolites({dna_out: -0.01})

    model.objective = biomassa
    return model


# ======================================
# 2. Simulação da curva de crescimento
# ======================================

def calcular_novo_upper_bound(original_ub, Ki, conc):
    if Ki <= 0:
        return 0
    return original_ub * (Ki / (Ki + conc))


def simular_crescimento_streamlit(kd_peptideo_uM, comparar=False):
    model = preparar_modelo(carregar_modelo())
    original_ub = model.reactions.get_by_id("DNA_GIRASE").upper_bound
    concs_uM = np.logspace(-1, 2, 20)
    concs_M = concs_uM * 1e-6
    resultados = []

    # Inibidores (em µM)
    inibidores = []
    if comparar:
        inibidores.append(("Ácido Nalidíxico", 1.3, "green"))
    inibidores.append(("Peptídeo", kd_peptideo_uM, "blue"))

    for nome, Ki_uM, cor in inibidores:
        print(f"🔎 Simulando: {nome} | Ki = {Ki_uM} µM")
        taxas = []
        for conc in concs_M:
            modelo_temp = copy.deepcopy(model)
            Ki_M = Ki_uM * 1e-6
            novo_ub = calcular_novo_upper_bound(original_ub, Ki_M, conc)
            modelo_temp.reactions.get_by_id("DNA_GIRASE").upper_bound = novo_ub
            taxa = modelo_temp.slim_optimize()
            taxas.append(taxa if taxa is not None else 0)
        print(f"→ Crescimento mínimo ({nome}): {min(taxas)}")
        resultados.append((nome, taxas, cor))

    # Gráfico
    fig, ax = plt.subplots(figsize=(8, 6))
    for nome, taxas, cor in resultados:
        ax.plot(concs_uM, taxas, marker='o', label=nome, color=cor)

    ax.set_xscale("log")
    ax.set_xlabel("[Inibidor] (µM)")
    ax.set_ylabel("Crescimento simulado (h⁻¹)")
    ax.set_title("Curva de crescimento com inibição da DNA-girase")
    ax.legend()
    ax.grid(True)
    return fig
