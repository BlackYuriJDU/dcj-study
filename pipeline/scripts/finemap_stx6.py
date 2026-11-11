#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finemap_stx6.py — Fine-mapping descritivo do locus STX6 (e GAL3ST1/PRNP).

Sumstats REAL GCST90001389 (Lancet Neurol 2020, GRCh37). Sem genótipos
individuais não há r² de LD — portanto isto é fine-mapping DESCRITIVO:
- top variantes regionais por p;
- verificação do SNP índice do Brain 2025 (rs3747957 = chr1:180,953,853 GRCh37);
- direção de efeito comparada (mesmo sinal beta = consistente);
- nota honesta: credible set formal exige painel de LD (próxima rodada).
"""
import gzip
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
import sys as _sys
_sys.path.insert(0, str(BASE / "lib"))
from gwas_fonte import exigir_gwas  # noqa: E402  (mensagem clara se o GWAS faltar)
REPORTS = BASE / "reports"

REGIOES = {
    "STX6": ("1", 180_900_000, 181_000_000, 180_953_853, "rs3747957 (índice Brain 2025)"),
    "GAL3ST1": ("22", 30_900_000, 31_000_000, None, "índice não consultado"),
    "PRNP": ("20", 4_600_000, 4_700_000, None, "região do gene do príon"),
}


def main() -> None:
    resultados = {k: [] for k in REGIOES}
    with gzip.open(exigir_gwas(), "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            crom = p[col["chromosome"]]
            for nome, (chr_alvo, ini, fim, _, _) in REGIOES.items():
                if crom == chr_alvo:
                    try:
                        pos = int(p[col["base_pair_location"]])
                    except ValueError:
                        continue
                    if ini <= pos <= fim:
                        try:
                            rec = (float(p[col["p_value"]]), crom, pos,
                                   p[col["other_allele"]], p[col["effect_allele"]],
                                   float(p[col["effect_allele_frequency"]]),
                                   float(p[col["beta"]]), float(p[col["standard_error"]]))
                        except (ValueError, IndexError):
                            continue
                        resultados[nome].append(rec)

    agora = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    L = [
        "# Fine-mapping descritivo — loci do GWAS sCJD (sumstats 2020, GRCh37)",
        f"*`finemap_stx6.py` em {agora}. Sem LD individual: ranking regional, não credible set formal.*",
        "",
    ]
    for nome, (chr_alvo, ini, fim, snp_idx, rotulo_idx) in REGIOES.items():
        recs = sorted(resultados[nome])
        L += [f"## Locus {nome} (chr{chr_alvo}:{ini:,}-{fim:,})",
              f"- Variantes na região: **{len(recs)}**"]
        if not recs:
            L += [""]
            continue
        lead = recs[0]
        L += [f"- Lead regional: chr{lead[1]}:{lead[2]:,} {lead[3]}>{lead[4]} p={lead[0]:.2e} "
              f"(beta={lead[6]:+.3f}, SE={lead[7]:.3f}, EAF={lead[5]:.3f})"]
        if snp_idx:
            match = [r for r in recs if r[2] == snp_idx]
            if match:
                r = match[0]
                mesmo_sinal = (r[6] * lead[6]) > 0
                L += [f"- **{rotulo_idx}**: presente em chr{r[1]}:{r[2]:,} {r[3]}>{r[4]} "
                      f"p={r[0]:.2e} (beta={r[6]:+.3f}) — "
                      f"{'mesma direção do lead' if mesmo_sinal else 'direção oposta ao lead'} "
                      f"({len(recs)} variantes regionais; rank do rs3747957 por p: "
                      f"{sorted(r2[0] for r2 in recs).index(r[0])+1}º)"]
            else:
                L += [f"- {rotulo_idx}: AUSENTE nos sumstats (verificar build/merge)"]
        L += ["", f"### Top 10 regionais", "",
              "| p | pos | OA>EA | EAF | beta | SE |", "|---|---|---|---|---|---|"]
        for p_v, c, pos, oa, ea, eaf, beta, se in recs[:10]:
            L.append(f"| {p_v:.2e} | {pos:,} | {oa}>{ea} | {eaf:.3f} | {beta:+.3f} | {se:.3f} |")
        L += [""]

    L += ["## Nota de honestidade científica",
          "- Fine-mapping formal (credible set, colocalização eQTL) exige LD entre variantes;",
          "  sem genótipos individuais, este relatório é RANKING DESCRITIVO.",
          "- Comparação com Brain 2025: mesmo lead/efeito = consistência; diferença de p",
          "  esperada (coortes maiores em 2025).",
          "- PRNP: o sinal regional inclui o gene do príon; interpretação biológica",
          "  (códon 129) pertence à literatura, não a este arquivo."]
    destino = REPORTS / "relatorio_finemap_loci.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    for nome, recs in resultados.items():
        if recs:
            print(f"[ok] {nome}: {len(recs)} variantes, lead p={min(r[0] for r in recs):.2e}")


if __name__ == "__main__":
    main()
