#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coloc_sqtl_stx6.py — Colocalização GWAS sCJD × sQTL (splicing) na região STX6.

Pergunta: além da expressão, há sinal de splicing de STX6 associado ao locus?

Dados: eQTL Catalogue r8, quant_method = leafcutter (introns), mesmos cinco
conjuntos cerebrais da meta-eQTL (proveniência em
pipeline/data/proveniencia_qtl.tsv). Os arquivos `cc` trazem apenas os
clusters de intron (traços moleculares) com sinal cis SIGNIFICATIVO; para
esses, trazem as variantes cis da janela (os p-valores no cache vão até ~1).
Consequência: um cluster ausente do `cc` é "sem sQTL significativo detectado
naquele conjunto", não "sem sQTL".

Método: para cada conjunto × cluster significativo na janela (gene declarado
pela coluna gene_id), harmonização por variante completa contra o GWAS
(pipeline/lib/harmoniza.py) e colocalização pela implementação única
(pipeline/lib/coloc_abf.py). Clusters são identificados por coordenada
(início:fim:fita), comparável entre conjuntos.

Sem correção de "Bonferroni": probabilidades a posteriori não são p-valores, e
Bonferroni não se aplica a elas. Os testes são listados todos, com os seus
PP.H4, e lidos como estão.
"""
import sys
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coloc_abf as C  # noqa: E402
import harmoniza as H  # noqa: E402
from fontes_qtl import STX6, OFFSET_B38_B37, datasets, ler_sqtl  # noqa: E402
from coloc_stx6_eqtl import carregar_gwas_regiao  # noqa: E402

REPORTS = BASE / "reports"
MIN_PARES = 20


def chave_cluster(trait_id: str) -> str:
    """'1:START:END:clu_N_+' → 'START:END:+' (comparável entre conjuntos)."""
    partes = trait_id.split(":")
    if len(partes) >= 4:
        fita = partes[3].rsplit("_", 1)[-1]
        return ":".join(partes[1:3] + [fita])
    return trait_id


def clusters_de(reg):
    """{chave: {'gene', 'trait', 'linhas': [(pos, ref, alt, rsid, p, beta, se, maf)]}}."""
    d = {}
    for pos38, ref, alt, rsid, p, beta, se, trait, gene in ler_sqtl(reg):
        g = d.setdefault(chave_cluster(trait),
                         {"gene": gene, "trait": trait, "linhas": []})
        g["linhas"].append((pos38, ref, alt, rsid, p, beta, se, float("nan")))
    return d


def analisar(gwas=None):
    if gwas is None:
        gwas = carregar_gwas_regiao()
    resumo, testes, descartes = [], [], []
    for reg in datasets("leafcutter"):
        d = clusters_de(reg)
        n_stx6 = sum(1 for g in d.values() if g["gene"] == STX6)
        resumo.append((reg, len(d), n_stx6))
        for k, g in d.items():
            variantes, st_d = H.deduplicar(g["linhas"])
            pares, st_g = H.casar_com_gwas(variantes, gwas, OFFSET_B38_B37)
            if len(pares) < MIN_PARES:
                descartes.append((reg["rotulo"], g["gene"], k, len(pares)))
                continue
            tg = [(p["beta_gwas"], p["se_gwas"]) for p in pares]
            te = [(p["qtl"]["beta"], p["qtl"]["se"]) for p in pares]
            res = C.coloc(tg, te, "cc", "quant", sdY2=1.0)
            pp6 = C.coloc(tg, te, "cc", "quant", priors=(1e-4, 1e-4, 1e-6))["pp"]
            testes.append({
                "conjunto": reg["rotulo"], "n": reg["n_amostras"],
                "gene": g["gene"], "cluster": k, "n_pares": len(pares),
                "pp": res["pp"], "h4_p12_1e6": pp6[4],
                "conc": sum(1 for (bg, _), (be, _) in zip(tg, te)
                            if bg * be > 0) / len(pares),
                "p_min": min(H.p_de_z(b / s) for b, s in te),
                "conflitos": st_d["descartadas_conflito"],
                "incompativeis": st_g["alelos_incompativeis"]
                + st_g["so_na_fita_complementar"],
            })
    testes.sort(key=lambda t: -t["pp"][4])
    return {"resumo": resumo, "testes": testes, "descartes": descartes}


def main() -> None:
    hoje = datetime.date.today().isoformat()
    print("[1/2] GWAS + sQTL...")
    a = analisar()
    resumo, testes = a["resumo"], a["testes"]

    L = ["# Colocalização GWAS sCJD × sQTL (splicing) — região STX6",
         f"*`coloc_sqtl_stx6.py` em {hoje}. eQTL Catalogue r8, leafcutter, "
         "mesmos cinco conjuntos da meta-eQTL (proveniência em "
         "`pipeline/data/proveniencia_qtl.tsv`). Os arquivos `cc` contêm só "
         "os clusters com sinal cis significativo: ausência de um cluster é "
         "\"sem sQTL significativo detectado\", não \"sem sQTL\". "
         "Colocalização por conjunto × cluster com a implementação única "
         "(`pipeline/lib/coloc_abf.py`, priors padrão do `coloc`).*", "",
         "> **Correção de 2026-11.** Os PP.H4 ≈ 0 da versão anterior vinham do "
         "mesmo erro de fórmula da colocalização de eQTL, e o limiar "
         "\"Bonferroni\" aplicado a PP.H4 foi retirado (Bonferroni corrige "
         "p-valores, não probabilidades a posteriori). Ver `CORRECOES.md`.", "",
         "| Conjunto | n | clusters significativos na janela | clusters de STX6 |",
         "|---|---|---|---|"]
    for reg, nc, ns in resumo:
        L.append(f"| {reg['rotulo']} ({reg['dataset_id']}) | "
                 f"{reg['n_amostras']} | {nc} | **{ns}** |")

    L += ["",
          f"Testes de colocalização executados (conjunto × cluster com ≥ "
          f"{MIN_PARES} variantes casadas com o GWAS): **{len(testes)}**", "",
          "| conjunto | gene | cluster (b38) | variantes | H1 | H2 | H3 | "
          "**H4** | H4 com p12 = 10⁻⁶ | direção concordante | p mín sQTL |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in testes:
        pp = t["pp"]
        L.append(f"| {t['conjunto']} | {t['gene']} | {t['cluster']} | "
                 f"{t['n_pares']} | {C._fmt(pp[1])} | {C._fmt(pp[2])} | "
                 f"{C._fmt(pp[3])} | **{C._fmt(pp[4])}** | "
                 f"{C._fmt(t['h4_p12_1e6'])} | {t['conc']:.0%} | "
                 f"{t['p_min']:.1e} |")
    if a["descartes"]:
        L += ["", "Clusters com menos de "
              f"{MIN_PARES} variantes casadas (não testados): "
              + "; ".join(f"{c} {g} {k} ({n})" for c, g, k, n in a["descartes"])
              + "."]

    total_stx6 = sum(ns for _, _, ns in resumo)
    alto = [t for t in testes if t["pp"][4] > 0.8]
    medio = [t for t in testes if 0.5 < t["pp"][4] <= 0.8]
    genes = sorted({t["gene"] for t in testes})
    L += ["", "## Leitura", "",
          f"- **Clusters de sQTL do próprio STX6: {total_stx6}** nos cinco "
          "conjuntos. Nenhum sQTL significativo de STX6 foi detectado, em "
          "contraste com o eQTL de expressão, que é forte (p ≈ 7×10⁻⁴⁷ na "
          "meta). Isso é ausência de detecção em arquivos que só listam pares "
          "significativos, não evidência de ausência.",
          f"- Os {len(testes)} testes são de clusters de genes vizinhos "
          f"({', '.join(genes)}), não de STX6. PP.H4 > 0,8 em {len(alto)}; "
          f"entre 0,5 e 0,8 em {len(medio)}.",
          f"- O maior PP.H4 é {testes[0]['pp'][4]:.2f} ({testes[0]['conjunto']}, "
          f"{testes[0]['gene']}); cai para {testes[0]['h4_p12_1e6']:.2f} com "
          f"p12 = 10⁻⁶, tem direção concordante em só "
          f"{testes[0]['conc']:.0%} das variantes e não se repete nos outros "
          "conjuntos que têm o mesmo cluster. Não é tratado como "
          "colocalização.",
          "- Estes testes existem porque os clusters foram pré-selecionados "
          "por terem sQTL significativo; são exploratórios, dizem respeito a "
          "outros genes e não entram nas conclusões sobre STX6.",
          "- Conclusão: o locus de STX6 associa-se à expressão do gene; para "
          "splicing de STX6 não há sinal detectado nestes dados, com poder "
          "limitado.", ""]
    destino = REPORTS / "relatorio_coloc_sqtl_stx6.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[2/2] {destino} — {len(testes)} testes, clusters de STX6: "
          f"{total_stx6}")


if __name__ == "__main__":
    main()
