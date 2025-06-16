import re

kd_scale = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

pKa_sidechains = {'D': 3.9, 'E': 4.1, 'C': 8.3, 'Y': 10.1, 'H': 6.0, 'K': 10.5, 'R': 12.5}

def calcular_carga(seq, pH=7):
    carga = 0
    carga += 1 / (1 + 10**(pH - 9))
    carga -= 1 / (1 + 10**(2 - pH))
    for aa in seq:
        if aa in ['K', 'R']:
            carga += 1 / (1 + 10**(pH - pKa_sidechains[aa]))
        elif aa in ['D', 'E', 'C', 'Y']:
            carga -= 1 / (1 + 10**(pKa_sidechains[aa] - pH))
        elif aa == 'H':
            carga += 1 / (1 + 10**(pH - pKa_sidechains[aa]))
    return round(carga, 2)

def calcular_hidrofobicidade(seq):
    valores = [kd_scale.get(aa, 0) for aa in seq]
    return round(sum(valores) / len(valores), 2) if valores else 0

def prever_clivagem_proteases(seq):
    proteases = {
        'Tripsina': [r'K(?=[^P]|$)', r'R(?=[^P]|$)'],
        'Quimotripsina': [r'F(?=[^P]|$)', r'W(?=[^P]|$)', r'Y(?=[^P]|$)'],
        'Elastase': [r'A(?=[^P]|$)', r'G(?=[^P]|$)', r'V(?=[^P]|$)'],
        'Pepsina (pH < 3)': [r'F', r'L']
    }
    resultados = {}
    for enzima, padroes in proteases.items():
        posicoes = []
        for padrao in padroes:
            for match in re.finditer(padrao, seq):
                posicoes.append(match.start() + 1)
        resultados[enzima] = sorted(posicoes)
    return resultados

def pontuar_estabilidade(clivagens_dict):
    pesos = {'Tripsina': 1.0, 'Quimotripsina': 0.8, 'Elastase': 0.6, 'Pepsina (pH < 3)': 0.5}
    total_pontos = sum(len(clivagens_dict.get(e, [])) * peso for e, peso in pesos.items())
    total_max = sum(5 * peso for peso in pesos.values())
    score = 1 - min(total_pontos / total_max, 1)
    return round(score, 3)

def prever_permeabilidade(seq):
    carga = calcular_carga(seq)
    hidrof = calcular_hidrofobicidade(seq)
    tamanho = len(seq)
    if (carga >= 2) and (-1.0 <= hidrof <= 1.5) and (tamanho <= 30):
        return "Alta probabilidade"
    return "Baixa probabilidade"

def prever_translocacao(seq):
    carga = calcular_carga(seq)
    hidrof = calcular_hidrofobicidade(seq)
    tamanho = len(seq)
    if (carga > 0) and (0 <= hidrof <= 2) and (tamanho <= 30):
        return "Alta probabilidade"
    return "Baixa probabilidade"

def pontuacao_combinada(estabilidade, perm, transl):
    p = 1 if perm == "Alta probabilidade" else 0
    t = 1 if transl == "Alta probabilidade" else 0
    return round((estabilidade + p + t) / 3, 3)

def analisar_peptideo(seq):
    seq = seq.upper().strip()
    carga = calcular_carga(seq)
    hidrof = calcular_hidrofobicidade(seq)
    cliv = prever_clivagem_proteases(seq)
    est = pontuar_estabilidade(cliv)
    perm = prever_permeabilidade(seq)
    transl = prever_translocacao(seq)
    score = pontuacao_combinada(est, perm, transl)
    return {
        "carga": carga,
        "hidrofobicidade": hidrof,
        "estabilidade_extracelular": est,
        "permeabilidade_membrana_externa": perm,
        "translocacao_membrana_citoplasmatica": transl,
        "pontuacao_combinada": score
    }
