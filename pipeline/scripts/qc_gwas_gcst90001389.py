#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
qc_gwas_gcst90001389.py — QC e triagem independente do sumstats do GWAS de sCJD.

Fonte: GCST90001389 (GWAS Catalog/EBI), build GRCh37, ~6.3M variantes.
Estudo: consórcio MRC Prion Unit, Lancet Neurol 2020 (PMID 32949544),
4.110 casos sCJD × 13.569 controles. Download verificado (HTTP 200).

QC implementada (stdlib puro, streaming):
1. Integridade: linhas totais/malformadas, contagem por cromossomo.
2. Distribuições: MAF (effect_allele_frequency), p-values.
3. Inflação genômica: lambda_GC = mediana(chi2)/0.4549, chi2=(beta/SE)^2.
4. Hits genômicos: p < 5e-8 (limiar consagrado).
5. Locus STX6 (candidato a replicação vs. Brain 2025): janela regional.
   STX6 GRCh37: chr1, ~159.9–160.2 Mb (verificar na saída).
6. Top 20 variantes por p-value com anotação de região citogênica aproximada.

Saída: pipeline/reports/relatorio_qc_gwas_gcst90001389.md
"""
import gzip
import math
import datetime
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
import sys as _sys
_sys.path.insert(0, str(BASE / "lib"))
from gwas_fonte import exigir_gwas, lambda_gc as _lambda_gc, MEDIANA_CHI2_1GL  # noqa: E402
REPORTS = BASE / "reports"

STX6_CHR = "1"
# CORRIGIDO (NCBI Gene, 1q25.3): STX6 GRCh37 = chr1:~180.9-181.0 Mb.
# Versão anterior usava ~160 Mb (janela errada — registrada em memory/mistakes.md).
STX6_START, STX6_END = 180_850_000, 181_050_000
GW_SIG = 5e-8


def main() -> None:
    n_total = n_bad = 0
    por_crom = Counter()
    chi2_vals: list[float] = []
    maf_min, maf_max = 1.0, 0.0
    n_piso_maf = 0
    hits_sig: list[tuple] = []
    stx6_best: tuple | None = None
    top20: list[tuple] = []  # heap implícito por p

    with gzip.open(exigir_gwas(), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for linha in fh:
            partes = linha.rstrip("\n").split("\t")
            n_total += 1
            if len(partes) < len(header):
                n_bad += 1
                continue
            try:
                crom = partes[col["chromosome"]]
                pos = int(partes[col["base_pair_location"]])
                ea = partes[col["effect_allele"]]
                oa = partes[col["other_allele"]]
                eaf = float(partes[col["effect_allele_frequency"]])
                beta = float(partes[col["beta"]])
                se = float(partes[col["standard_error"]])
                p = float(partes[col["p_value"]])
            except (ValueError, IndexError):
                n_bad += 1
                continue
            por_crom[crom] += 1
            if eaf < maf_min:
                maf_min = eaf
            if eaf > maf_max:
                maf_max = eaf
            if eaf == 0.0 or eaf == 1.0:
                n_piso_maf += 1
            if se > 0:
                z2 = (beta / se) ** 2
                chi2_vals.append(z2)
            if p <= 0 or p > 1:
                n_bad += 1
                continue
            registro = (p, crom, pos, oa, ea, eaf, beta, se)
            if p < GW_SIG:
                hits_sig.append(registro)
            if crom == STX6_CHR and STX6_START <= pos <= STX6_END:
                if stx6_best is None or p < stx6_best[0]:
                    stx6_best = registro
            if len(top20) < 20:
                top20.append(registro)
                top20.sort(reverse=True)
            elif p < top20[-1][0]:
                top20[-1] = registro
                top20.sort(reverse=True)

    chi2_vals.sort()
    n_chi = len(chi2_vals)
    mediana_chi2 = chi2_vals[n_chi // 2] if n_chi else float("nan")
    lambda_gc = _lambda_gc(chi2_vals)      # mediana / 0.454936 (constante única)

    agora = datetime.date.today().isoformat()
    L = [
        "# QC independente — GWAS sCJD GCST90001389 (REAIS, 4.110 casos × 13.569 controles)",
        f"*`qc_gwas_gcst90001389.py` em {agora}. Fonte: GWAS Catalog/EBI, GRCh37, "
        "consórcio MRC Prion Unit (PMID 32949544).*",
        "",
        "## 1. Integridade",
        f"- Variantes processadas: **{n_total:,}** · Linhas malformadas: **{n_bad}** "
        f"({100*n_bad/max(1,n_total):.4f}%)",
        "- Variantes por cromossomo (1-22, X): "
        + ", ".join(f"{c}:{por_crom[c]}" for c in
                    sorted(por_crom, key=lambda x: (len(x), x))),
        "",
        "## 2. Distribuições",
        f"- EAF: min={maf_min:.4f}, max={maf_max:.4f}, variantes no piso (0/1): {n_piso_maf}",
        f"- chi2=(beta/SE)^2: mediana={mediana_chi2:.4f}",
        f"- **lambda_GC = {lambda_gc:.4f}** (mediana / {MEDIANA_CHI2_1GL}, sobre todas as variantes) "
        + ("(sem inflação relevante; ≤1.05 é saudável)" if lambda_gc <= 1.05
           else "(⚠ inflação — investigar estratificação/ancestralidade)"),
        "",
        "## 3. Hits genômicos (p < 5e-8)",
        f"- Total: **{len(hits_sig)}**",
    ]
    if hits_sig:
        L += ["", "| p | chr | pos | OA>EA | EAF | beta | SE |", "|---|---|---|---|---|---|---|"]
        for p, c, pos, oa, ea, eaf, beta, se in sorted(hits_sig)[:30]:
            L.append(f"| {p:.2e} | {c} | {pos:,} | {oa}>{ea} | {eaf:.3f} | {beta:.3f} | {se:.3f} |")

    # --- 3b. Os três loci publicados: melhor variante de cada um -----------
    LOCI_PUBLICADOS = [
        ("PRNP", "20", 4_600_000, 4_700_000,
         "região do gene da proteína priônica (chr20p13)"),
        ("STX6", "1", 180_900_000, 181_000_000,
         "dentro do gene STX6 (Ensembl GRCh37: 180,941,861–180,992,047)"),
        ("GAL3ST1", "22", 30_900_000, 31_000_000,
         "borda 5' do gene GAL3ST1 (Ensembl GRCh37: 30,950,622–30,970,574)"),
    ]
    L += ["", "## 3b. Os três loci publicados, reproduzidos",
          "",
          "O artigo original (Jones et al., Lancet Neurol 2020) reporta três "
          "loci com significância genômica: PRNP, STX6 e GAL3ST1. Melhor "
          "variante deste QC em cada um:", "",
          "| Locus | Melhor variante (GRCh37) | p | hits p<5e-8 no locus | Anotação |",
          "|---|---|---|---|---|"]
    n_rep = 0
    for nome, crom, ini, fim, nota in LOCI_PUBLICADOS:
        no_locus = sorted(h for h in hits_sig
                          if h[1] == crom and ini <= h[2] <= fim)
        if no_locus:
            n_rep += 1
            p, c, pos, *_ = no_locus[0]
            L.append(f"| **{nome}** | chr{c}:{pos:,} | {p:.2e} | "
                     f"{len(no_locus)} | {nota} |")
        else:
            L.append(f"| **{nome}** | — | — | 0 | sem hit genômico-significativo |")
    L += ["", f"**{n_rep} de 3 loci publicados reproduzidos** a partir das "
          "estatísticas-resumo depositadas. É uma conferência dos dados "
          "públicos, não uma replicação em coorte independente."]

    L += ["", "## 4. Locus STX6 (candidato a replicação — cf. Brain 2025)",
          f"- Janela: chr{STX6_CHR}:{STX6_START:,}-{STX6_END:,} (GRCh37)"]
    if stx6_best:
        p, c, pos, oa, ea, eaf, beta, se = stx6_best
        L.append(f"- Melhor variante na janela: chr{c}:{pos:,} {oa}>{ea} "
                 f"p={p:.3e} (beta={beta:.3f}, SE={se:.3f}, EAF={eaf:.3f})")
        L.append(f"- **p {'<' if p < GW_SIG else '≥'} 5e-8** — "
                 + ("sinal genômico-significativo nesta coorte 2020, "
                    "consistente com o artigo original (3 loci: PRNP, STX6, GAL3ST1)"
                    if p < GW_SIG
                    else "sem significância genômica nesta coorte 2020"))
    else:
        L.append("- Nenhuma variante na janela (verificar coordenadas).")

    L += ["", "## 5. Top 20 variantes por p-value", "",
          "| p | chr | pos | OA>EA | EAF | beta | SE |", "|---|---|---|---|---|---|---|"]
    for p, c, pos, oa, ea, eaf, beta, se in sorted(top20):
        L.append(f"| {p:.2e} | {c} | {pos:,} | {oa}>{ea} | {eaf:.3f} | {beta:.3f} | {se:.3f} |")

    L += ["", "## Nota de honestidade científica",
          "- QC de primeira passada: sem verificação de strand, sem imputação-info",
          "  (coluna não existe no arquivo), sem clumping por LD (próxima rodada).",
          "- Sem rsIDs no arquivo — coordenadas GRCh37 são a chave primária.",
          "- lambda_GC de sumstats de caso-controle é aproximado (chi2 de z de beta/SE).",
          "- Este QC NÃO é descoberta nova: é verificação independente documentada."]

    destino = REPORTS / "relatorio_qc_gwas_gcst90001389.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    print(f"[ok] variantes={n_total:,} ruins={n_bad} hits_sig={len(hits_sig)} "
          f"lambda_GC={lambda_gc:.4f} STX6_best_p={stx6_best[0]:.2e} "
          f"chr={stx6_best[1]} pos={stx6_best[2]:,}" if stx6_best else "[ok] STX6: sem variante na janela")


if __name__ == "__main__":
    main()
