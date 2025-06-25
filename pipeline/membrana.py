# pipeline/membrana.py

import numpy as np
import os
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

# Parâmetros físicos que definimos
CARGAS_SIMPLIFICADAS = { 'H': 0.1, 'C': -0.1, 'N': -0.3, 'O': -0.5, 'S': -0.2 }

def obter_carga_simplificada(elemento):
    return CARGAS_SIMPLIFICADAS.get(elemento.upper(), 0.0)

def carga_membrana(z):
    if z < -15 or z > 15: return 0.0
    elif -15 <= z < -7.5 or 7.5 < z <= 15: return -0.9
    else: return 0.0

def analisar_interacao_pdb(pdb_file_object, output_dir="static/graficos"):
    """
    Função principal que lê um objeto de arquivo PDB, calcula as energias
    e gera/salva os gráficos de resultado.
    """
    # Garante que o diretório de output exista
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # O Streamlit fornece um objeto de arquivo, precisamos salvá-lo temporariamente
    # para que o PDBParser possa lê-lo pelo caminho do arquivo.
    pdb_path = os.path.join(output_dir, "temp.pdb")
    with open(pdb_path, "wb") as f:
        f.write(pdb_file_object.getbuffer())

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("peptideo", pdb_path)

    # Dicionários e listas para armazenar os resultados
    coordenadas_z = []
    energias_atomicas = []
    energia_por_residuo = defaultdict(list)
    
    # Loop principal de cálculo
    for residue in structure.get_residues():
        if not is_aa(residue, standard=True): continue
        for atom in residue:
            pos = atom.coord
            z = pos[2]
            r = max(abs(z), 2.0)
            q = obter_carga_simplificada(atom.element)
            q_membrana = carga_membrana(z)

            # Simplificação do potencial de Lennard-Jones (apenas a parte radial)
            lj = (1 / r**12) - (1 / r**6)
            coul = (q * q_membrana) / (80.0 * r)
            energia_total_atomo = lj + coul
            
            coordenadas_z.append(z)
            energias_parciais.append(energia_total_atomo)
            energia_por_residuo[residue.get_id()].append(energia_total_atomo)
    
    energia_media_residuo = {rid[1]: np.mean(valores) for rid, valores in energia_por_residuo.items()}

    # --- Geração do Gráfico 1: Energia por Átomo vs Z ---
    plt.figure(figsize=(8, 5))
    plt.scatter(coordenadas_z, energias_parciais, alpha=0.7, color='darkgreen')
    plt.axhline(0, color='gray', linestyle='--')
    plt.xlabel("Coordenada Z (Å)")
    plt.ylabel("Energia Parcial de Interação")
    plt.title("Energia por Átomo vs Profundidade na Membrana")
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    caminho_grafico1 = os.path.join(output_dir, "grafico_energia_vs_z.png")
    plt.savefig(caminho_grafico1)
    plt.close()

    # --- Geração do Gráfico 2: Energia Média por Resíduo ---
    res_indices = list(energia_media_residuo.keys())
    res_energies = list(energia_media_residuo.values())
    
    plt.figure(figsize=(10, 5))
    plt.plot(res_indices, res_energies, marker='o', linestyle='-', color='mediumblue')
    plt.xlabel("Índice do Resíduo")
    plt.ylabel("Energia Média de Interação")
    plt.title("Energia Média por Resíduo")
    plt.xticks(res_indices)
    plt.grid(True)
    plt.tight_layout()
    caminho_grafico2 = os.path.join(output_dir, "grafico_energia_por_residuo.png")
    plt.savefig(caminho_grafico2)
    plt.close()

    # Deleta o arquivo temporário
    os.remove(pdb_path)
    
    # Retorna os caminhos para as imagens salvas
    return caminho_grafico1, caminho_grafico2
