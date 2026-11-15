#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gera_figuras.py — Item 5 do plano de melhoria: visualizações.

Figuras geradas em pipeline/reports/figuras/:
1. volcano_gse160208.png  — genes córtex frontal (CJD vs CT), FDR<0.05 destacado
2. volcano_gse140069.png  — miRNAs sanguíneos (CJD vs CT)
3. timeline_caso_referencia.png      — linha do tempo clínica do caso simulado
4. heatmap_top_genes.png  — top 25 genes × amostras FC

Reusa as funções dos scripts de análise (fonte única de verdade).
"""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analise_gse160208 import (extrair_tabela as ex160, testar_fc,  # noqa: E402
                               q_no_universo, classificar_amostras)
from analise_gse140069 import analisar as analisar140  # noqa: E402

BASE = Path(__file__).resolve().parents[1]
FIGS = BASE / "reports" / "figuras"
FIGS.mkdir(parents=True, exist_ok=True)

AZUL, VERMELHO, CINZA = "#2563eb", "#dc2626", "#9ca3af"


def volcano(nome, pares_delta_p, qmap, rotulos, titulo):
    fig, ax = plt.subplots(figsize=(7, 5))
    for g, d, p in pares_delta_p:
        q = qmap.get(g, 1.0)
        cor = CINZA
        if q < 0.05:
            cor = VERMELHO if d > 0 else AZUL
        ax.scatter(d, -math.log10(max(p, 1e-300)), s=8, c=cor, alpha=0.75)
    pmap = {g: p for g, _, p in pares_delta_p}
    for rot in rotulos:
        g, d = rot[0], rot[1]
        # rótulo na altura do PONTO (p), não do q — antes ficava deslocado
        ax.annotate(g, (d, -math.log10(max(pmap[g], 1e-300))),
                    fontsize=7, alpha=0.9)
    ax.axhline(-math.log10(0.05), ls="--", lw=0.7, c=CINZA)
    ax.set_xlabel("Δ média (CJD − controle)")
    ax.set_ylabel("-log10 p")
    ax.set_title(titulo)
    fig.tight_layout()
    destino = FIGS / nome
    fig.savefig(destino, dpi=150)
    plt.close(fig)
    print(f"[ok] {destino}")


def main() -> None:
    # --- Volcano GSE160208 -------------------------------------------------
    # Mesmos testes e mesmo universo de FDR do relatório (genes estimáveis);
    # os genes sem variância não têm p-valor e não entram na figura.
    amostras, genes, vals, meta, covs = ex160()
    res160 = testar_fc(amostras, genes, vals)
    _gr = classificar_amostras(amostras)
    fc_cjd, fc_ct = _gr["CJD_FC"], _gr["CT_FC"]
    q160 = q_no_universo(res160)
    pares = [(r["gene"], r["delta"], r["p"]) for r in res160 if r["estimavel"]]
    topo_up = sorted((x for x in pares if x[1] > 0), key=lambda t: t[2])[:4]
    topo_dn = sorted((x for x in pares if x[1] < 0), key=lambda t: t[2])[:4]
    volcano("volcano_gse160208.png", pares, q160,
            [t[:3] for t in topo_up + topo_dn],
            "GSE160208 — Córtex frontal sCJD vs. controles (real)")

    # --- Volcano GSE140069 (p do modelo AJUSTADO idade+sexo+RIN) -----------
    # Resultados de analise_gse140069.analisar(): só testes estimáveis.
    r140 = analisar140()
    pares140 = [(m, r140["betaB"][m], r140["pB"][m]) for m in sorted(r140["pB"])]
    qmap140 = r140["qB"]
    topo140 = sorted(pares140, key=lambda t: t[2])[:6]
    volcano("volcano_gse140069.png", pares140, qmap140,
            [(m, l) for m, l, _ in topo140],
            "GSE140069 — sangue sCJD vs. CT (OLS ajustado idade+sexo+RIN)")

    # --- Timeline caso de referência -------------------------------------------------------
    meses = ["M0\ninespecífico", "M1\ncognitivo", "M2\nataxia",
             "M3\nmioclonias", "M4\navançada", "M5\nterminal"]
    dependencia = [10, 35, 60, 80, 95, 100]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(range(len(meses)), dependencia, "o-", c=AZUL, lw=2)
    ax.fill_between(range(len(meses)), dependencia, color=AZUL, alpha=0.12)
    ax.set_xticks(range(len(meses)))
    ax.set_xticklabels(meses, fontsize=8)
    ax.set_ylabel("Dependência de cuidados (%)")
    ax.set_title("Caso Referência (SIMULADO) — progressão típica sCJD MM1")
    ax.set_ylim(0, 105)
    for x, y in zip(range(6), dependencia):
        ax.annotate(f"{y}%", (x, y + 3), ha="center", fontsize=8)
    fig.tight_layout()
    destino = FIGS / "timeline_caso_referencia.png"
    fig.savefig(destino, dpi=150)
    plt.close(fig)
    print(f"[ok] {destino}")

    # --- Heatmap top 25 genes × amostras FC ----------------------------------
    deltas = []
    for gi, g in enumerate(genes):
        xs = [vals[gi][i] for i in fc_cjd]
        ys = [vals[gi][i] for i in fc_ct]
        deltas.append((g, gi, sum(xs)/len(xs) - sum(ys)/len(ys)))
    deltas.sort(key=lambda t: abs(t[2]), reverse=True)
    selecionados = deltas[:25]
    todas_fc = fc_ct + fc_cjd
    matriz = []
    for _, gi, _ in selecionados:
        linha = [vals[gi][i] for i in todas_fc]
        mu = sum(linha) / len(linha)
        sd = (sum((x - mu) ** 2 for x in linha) / (len(linha) - 1)) ** 0.5 or 1.0
        matriz.append([(x - mu) / sd for x in linha])  # z-score por gene
    nomes = [g for g, _, _ in selecionados]
    fig, ax = plt.subplots(figsize=(11, 7))
    im = ax.imshow(matriz, aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(nomes)))
    ax.set_yticklabels(nomes, fontsize=7)
    ax.set_xticks(range(len(todas_fc)))
    ax.set_xticklabels(["CT"]*len(fc_ct) + ["CJD"]*len(fc_cjd),
                       rotation=90, fontsize=6)
    ax.axvline(len(fc_ct)-0.5, c="black", lw=1)
    ax.set_title("Top 25 genes por |Δ| — amostras FC (controles | CJD)")
    fig.colorbar(im, shrink=0.7, label="expressão (z-score por gene)")
    fig.tight_layout()
    destino = FIGS / "heatmap_top_genes.png"
    fig.savefig(destino, dpi=150)
    plt.close(fig)
    print(f"[ok] {destino}")


if __name__ == "__main__":
    main()
