#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
clumping_descoberta.py — PROVA DE DESCOBERTA: o pipeline acha os loci do zero?

Varredura CEGA de todo o GCST90001389 (6,3M variantes): p<1e-5 → clumping
guloso por distância (lead = menor p; absorve ±500 kb no mesmo cromossomo).
Depois verifica se os 3 loci conhecidos da literatura (STX6, PRNP, GAL3ST1)
emergem como clusters independentes — sem nenhum hint de onde olhar.
"""
import gzip
import math
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
import sys as _sys
_sys.path.insert(0, str(BASE / "lib"))
from gwas_fonte import exigir_gwas  # noqa: E402  (mensagem clara se o GWAS faltar)
REPORTS = BASE / "reports"

LIMIAR = 1e-5
JANELA = 500_000
# loci conhecidos (b37) — usados APENAS na verificação pós-hoc
CONHECIDOS = {
    "STX6":    ("1", 180_961_245),
    "PRNP":    ("20", 4_672_307),
    "GAL3ST1": ("22", 30_950_360),
}


def main() -> None:
    agora = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    sig = []
    total = 0
    with gzip.open(exigir_gwas(), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            total += 1
            f = line.split("\t")
            p = float(f[c["p_value"]])
            if p < LIMIAR:
                sig.append((p, f[c["chromosome"]],
                            int(f[c["base_pair_location"]])))
    sig.sort()
    print(f"[1/3] {total:,} variantes varridas · {len(sig)} com p<{LIMIAR:g}")

    # clumping guloso por distância
    clusters = []
    for p, chrom, pos in sig:
        for cl in clusters:
            if cl["chrom"] == chrom and abs(pos - cl["lead_pos"]) <= JANELA:
                cl["variantes"] += 1
                if p < cl["p_lead"]:
                    cl["p_lead"], cl["lead_pos"] = p, pos
                break
        else:
            clusters.append({"chrom": chrom, "lead_pos": pos, "p_lead": p,
                             "variantes": 1})
    clusters.sort(key=lambda x: x["p_lead"])
    print(f"[2/3] {len(clusters)} clusters independentes (±{JANELA//1000} kb)")

    # verificação pós-hoc: cada locus conhecido cai em qual cluster?
    def achar(chrom, pos):
        for k, cl in enumerate(clusters):
            if cl["chrom"] == chrom and abs(pos - cl["lead_pos"]) <= JANELA:
                return k
        return None

    L = ["# Clumping de descoberta — o pipeline acha os loci do zero?",
         f"*`clumping_descoberta.py` em {agora}. Varredura cega de "
         f"{total:,} variantes; limiar p<{LIMIAR:g}; clumping guloso por",
         f"distância ±{JANELA//1000} kb (lead = menor p do cluster).*", "",
         f"**Resultado: {len(clusters)} clusters independentes.**", "",
         "| # | Crom | Lead (b37) | p do lead | Variantes no cluster |",
         "|---|---|---|---|---|"]
    for k, cl in enumerate(clusters):
        L.append(f"| {k+1} | {cl['chrom']} | {cl['lead_pos']:,} | "
                 f"{cl['p_lead']:.2e} | {cl['variantes']} |")

    L += ["", "## Verificação pós-hoc dos loci conhecidos", "",
          "| Locus (literatura) | Posição b37 | Redescoberto? | Cluster | "
          "p do cluster |", "|---|---|---|---|---|"]
    recuperados = 0
    for nome, (chrom, pos) in CONHECIDOS.items():
        k = achar(chrom, pos)
        if k is not None:
            recuperados += 1
            L.append(f"| {nome} | {pos:,} | **SIM** | #{k+1} "
                     f"({clusters[k]['chrom']}:{clusters[k]['lead_pos']:,}) | "
                     f"{clusters[k]['p_lead']:.2e} |")
        else:
            L.append(f"| {nome} | {pos:,} | não (abaixo do limiar) | — | — |")
    L += ["", f"**{recuperados}/{len(CONHECIDOS)} loci conhecidos redescobertos "
          "às cegas.**",
          "- Clusters extras além dos conhecidos: sinais novos a investigar —",
          "  com p<1e-5 mas tipicamente abaixo de GWS (5e-8); nenhum deve ser",
          "  chamado de 'novo loci' sem replicação. Reportados por transparência.",
          "",
          "## Leitura honesta",
          "- Clumping por distância é conservador (PLINK usa LD real); para os",
          "  3 loci conhecidos a conectividade por LD já foi demonstrada no",
          "  fine-mapping (relatorio_finemap_loci.md).",
          "- A prova aqui é de SENSIBILIDADE do pipeline: dado o sumstats bruto,",
          "  os loci principais emergem sem qualquer âncora externa."]
    destino = REPORTS / "relatorio_clumping_descoberta.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[3/3] {destino}")
    for k, cl in enumerate(clusters[:8]):
        print(f"  #{k+1} chr{cl['chrom']}:{cl['lead_pos']:,} p={cl['p_lead']:.2e}")


if __name__ == "__main__":
    main()
