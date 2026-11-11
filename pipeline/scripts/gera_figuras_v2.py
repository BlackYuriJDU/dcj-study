#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gera_figuras_v2.py — Manhattan (GWAS) + forest (miRNAs brutos vs ajustados)."""
import gzip, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[1]
FIGS = BASE / "reports" / "figuras"
FIGS.mkdir(parents=True, exist_ok=True)

# ---------- Manhattan (p<1e-4 para legibilidade; leads destacados) ----------
LOCI = {"PRNP": ("20", 4_672_307), "GAL3ST1": ("22", 30_950_360),
        "STX6": ("1", 180_961_245)}
xs, ys, cs, labels = [], [], [], []
offset, chr_max, cores = 0, {}, ["#2563eb", "#7c3aed"]
import sys as _sys
_sys.path.insert(0, str(BASE / "lib"))
from gwas_fonte import exigir_gwas  # noqa: E402  (mensagem clara se o GWAS faltar)
with gzip.open(exigir_gwas(), "rt") as fh:
    hdr = fh.readline().split("\t")
    c = {n: i for i, n in enumerate(h) for h in [hdr]}[  # noqa
        "p_value"] if False else {n: i for i, n in enumerate(hdr)}
    c_chr, c_pos = c["chromosome"], c["base_pair_location"]
    atual, off = None, 0
    for line in fh:
        f = line.split("\t")
        p = float(f[c["p_value"]])
        if p >= 1e-4:
            continue
        chrom, pos = f[c_chr], int(f[c_pos])
        if chrom != atual:
            atual, off = chrom, chr_max.get(chrom, offset)
        x = pos + off
        xs.append(x); ys.append(-math.log10(p))
        cs.append(cores[ord(chrom[-1]) % 2])
        chr_max[chrom] = max(chr_max.get(chrom, 0), x + 1e6)
fig, ax = plt.subplots(figsize=(12, 5))
ax.scatter(xs, ys, c=cs, s=4, alpha=0.6)
ax.axhline(-math.log10(5e-8), ls="--", c="#dc2626", lw=1)
ax.text(1, -math.log10(5e-8) + 0.15, "GWS 5×10⁻⁸", c="#dc2626", fontsize=8)
for nome, (chrom, pos) in LOCI.items():
    x = pos + chr_max.get(chrom, 0) - 1e6
    ax.annotate(nome, (x, 9.2), fontsize=10, weight="bold",
                ha="center", c="#111111")
ax.set_xlabel("Genoma (GRCh37)"); ax.set_ylabel("−log10(p)")
ax.set_title("GWAS sCJD — GCST90001389 (p<10⁻⁴ exibidos; QC independente)")
fig.tight_layout(); fig.savefig(FIGS / "manhattan_gwas.png", dpi=150); plt.close(fig)

# ---------- Forest: miRNAs (efeito bruto vs ajustado) ----------
# Calculado dos dados por analise_gse140069 (correção de 2026-11: a versão
# anterior trazia os números digitados à mão, e eles não correspondiam aos
# dados — ex.: d bruto de miR-16-5p é -1,04, não -0,17).
#   d bruto     = d de Cohen em log2(x+1), sem covariáveis;
#   d ajustado  = coeficiente de grupo do modelo ajustado (idade+sexo+RIN)
#                 dividido pelo desvio-padrão residual desse modelo.
import math as _math
_sys.path.insert(0, str(Path(__file__).resolve().parent))
import analise_gse140069 as _m140  # noqa: E402
_r = _m140.analisar()
_d = _r["dados"]
dados = []
for _mi in ("hsa-miR-16-5p", "hsa-miR-93-5p", "hsa-let-7i-5p",
            "hsa-miR-106b-3p", "hsa-miR-500a-3p"):
    _k = _d["mirnas"].index(_mi)
    _a = _r["desenho"].ajustar([_d["vals_log"][_k][i] for i in _d["ok_ajust"]])
    dados.append((_mi, _r["dB"][_mi], _a.beta[1] / _math.sqrt(_a.s2)))
y = list(range(len(dados)))
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.hlines(y, [-1.5]*len(y), [0.4]*len(y), color="#e5e7eb", zorder=1)
ax.axvline(0, c="#9ca3af", lw=1)
ax.scatter([d[1] for d in dados], [i+0.15 for i in y], c="#dc2626", s=45,
           label="d bruto (não ajustado)", zorder=3)
ax.scatter([d[2] for d in dados], [i-0.15 for i in y], c="#2563eb", s=45,
           label="efeito ajustado / DP residual (idade+sexo+RIN)", zorder=3)
ax.set_yticks(y); ax.set_yticklabels([d[0] for d in dados])
ax.invert_yaxis(); ax.set_xlabel("efeito padronizado (CJD − controle), log2(x+1)")
ax.set_title("miRNAs sanguíneos sCJD — efeito bruto vs ajustado")
ax.legend(fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(FIGS / "forest_mirnas.png", dpi=150); plt.close(fig)
print("[ok] manhattan_gwas.png + forest_mirnas.png")
