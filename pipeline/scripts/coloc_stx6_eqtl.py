#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coloc_stx6_eqtl.py — Colocalização por fatores de Bayes aproximados
(Giambartolomei 2014; implementação única em pipeline/lib/coloc_abf.py,
conforme o pacote R `coloc`) entre o GWAS de sCJD (GCST90001389, n=17.679) e os
eQTLs de STX6 em dois tecidos cerebrais do GTEx v10 (eQTL Catalogue r8), um
tecido de cada vez.

Fonte eQTL: caches versionados em pipeline/data/eqtl_QTD*_stx6.tsv; release,
tamanhos amostrais e checksums em pipeline/data/proveniencia_qtl.tsv.

Coordenadas: as posições do eQTL vêm em GRCh38 e são convertidas para GRCh37
por um offset CONSTANTE (-30.864 bp), conferido contra a API de assembly map do
Ensembl ao longo de todo o span de STX6 (valida_offset_grch38_grch37.py →
relatorio_validacao_offset_genomico.md): o deslocamento é exatamente +30.864 bp
em todos os pontos amostrados. Um offset uniforme só é exato onde não há indel
entre os builds; uma reimplementação deve preferir liftOver por chain file.

Harmonização (pipeline/lib/harmoniza.py): identidade completa da variante
(posição, REF, ALT), linhas repetidas por rsID deduplicadas, os dois alelos
conferidos contra o GWAS e sinal do GWAS invertido quando seu alelo de efeito é
o REF do eQTL. Todo descarte é contado no relatório.

Saída: pipeline/reports/relatorio_coloc_stx6.md + figura regional.
"""
import math
import sys
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "lib"))
import coloc_abf as C  # noqa: E402
import harmoniza as H  # noqa: E402
from fontes_qtl import (STX6, REGIAO_B37, OFFSET_B38_B37, datasets,  # noqa: E402,F401
                        ler_eqtl)
from gwas_fonte import carregar_regiao  # noqa: E402

REPORTS = BASE / "reports"
FIGS = REPORTS / "figuras"
N_GWAS = 17_679                   # 4.110 casos + 13.569 controles (GWAS Catalog)
TECIDOS = {"QTD000176": "DLPFC", "QTD000166": "CEREBELO"}


def carregar_gwas_regiao():
    """GWAS da janela de STX6 (GRCh37), indexado por posição → lista de
    registros (nenhum alelo é sobrescrito)."""
    return H.indexar_gwas(carregar_regiao("1", *REGIAO_B37))


def linhas_harmonizacao(st_dedup, st_gwas):
    """Linhas de relatório com o balanço de variantes (o que entrou e o que
    foi descartado, com o motivo)."""
    return [
        f"- Linhas no cache: {st_dedup['linhas']} · variantes distintas "
        f"(posição+REF+ALT): {st_dedup['variantes']} · linhas repetidas "
        f"(mesma variante, outro rsID): {st_dedup['linhas_repetidas_por_rsid']} · "
        f"indels: {st_dedup['indels']} · posições multialélicas: "
        f"{st_dedup['posicoes_multialelicas']} · descartadas por beta/EP "
        f"conflitantes: {st_dedup['descartadas_conflito']}",
        f"- Casadas com o GWAS pelos dois alelos: **{st_gwas['casadas']}** "
        f"(mesmo alelo de efeito: {st_gwas['casadas_mesmo_alelo_de_efeito']}; "
        f"sinal do GWAS invertido: {st_gwas['casadas_sinal_invertido']}; "
        f"palindrômicas A/T ou C/G entre as casadas: "
        f"{st_gwas['palindromicas_casadas']})",
        f"- Não usadas: sem variante do GWAS na posição: "
        f"{st_gwas['sem_variante_no_gwas']} · alelos incompatíveis: "
        f"{st_gwas['alelos_incompativeis']} · só casariam na fita "
        f"complementar: {st_gwas['so_na_fita_complementar']}",
    ]


def main() -> None:
    hoje = datetime.date.today().isoformat()
    print("[1/3] GWAS região STX6...")
    gwas = carregar_gwas_regiao()
    print(f"      {sum(len(v) for v in gwas.values())} variantes GWAS na janela")
    w_cc, w_q = C.prior_w("cc"), C.prior_w("quant", 1.0)

    L = ["# Colocalização STX6: GWAS sCJD × eQTL cerebral (GTEx v10), por tecido",
         f"*`coloc_stx6_eqtl.py` em {hoje}. Método: colocalização por fatores de "
         "Bayes aproximados (Giambartolomei 2014), implementação única em "
         "`pipeline/lib/coloc_abf.py`, conforme o pacote R `coloc`. Prior de "
         f"efeito: variância {w_cc:g} para o GWAS caso-controle e {w_q:g} para o "
         "eQTL quantitativo (0,15²·sdY², com sdY = 1 porque a expressão é "
         "normalizada antes do mapeamento). Priors de hipótese: p1 = p2 = 10⁻⁴, "
         "p12 = 10⁻⁵ (padrão do `coloc`). GWAS n = 17.679 (GCST90001389); eQTL: "
         "eQTL Catalogue r8, caches versionados (ver "
         "`pipeline/data/proveniencia_qtl.tsv`); região chr1:180,9–181,1 Mb "
         "(GRCh37); eQTL em GRCh38 convertido com offset −30.864.*", "",
         "> **Correção de 2026-11.** A versão anterior deste relatório dava "
         "PP.H4 ≈ 0 por um erro de fórmula (o termo de H4 era multiplicado por "
         "p1·p2·p12 em vez de p12). Ver `CORRECOES.md`.", ""]

    regs = {r["dataset_id"]: r for r in datasets("ge")}
    resumo = []
    for qtd, tecido in TECIDOS.items():
        reg = regs[qtd]
        desc = f"{reg['rotulo']} ({qtd}, n = {reg['n_amostras']})"
        print(f"[2/3] eQTL {tecido} — {desc}...")
        variantes, st_d = H.deduplicar(ler_eqtl(reg))
        pares, st_g = H.casar_com_gwas(variantes, gwas, OFFSET_B38_B37)
        tg = [(p["beta_gwas"], p["se_gwas"]) for p in pares]
        te = [(p["qtl"]["beta"], p["qtl"]["se"]) for p in pares]
        res = C.coloc(tg, te, "cc", "quant", sdY2=1.0)
        pp = res["pp"]
        resumo.append((tecido, reg["n_amostras"], len(pares), pp))

        L += [f"## {tecido} — {desc}", ""]
        L += linhas_harmonizacao(st_d, st_g)
        L += ["", "### Posteriores (priors padrão do `coloc`)"]
        L += C.md_tabela_hipoteses(pp, "o GWAS de sCJD", "o eQTL de STX6")
        L += ["",
              f"**PP.H4 = {pp[4]:.4f}** · H3 + H4 (os dois traços associados) = "
              f"{pp[3] + pp[4]:.4f} · H4/(H3+H4) = {C.h4_condicional(pp):.4f} · "
              f"n = {res['n']} variantes",
              "", "### Sensibilidade a p12 (p1 e p2 fixos)"]
        L += C.md_sensibilidade_p12(tg, te, "cc", "quant")
        L += ["", "### Sensibilidade ao prior de efeito W"]
        L += C.md_sensibilidade_w(tg, te)

        # variantes de menor p no GWAS entre as casadas
        ordem = sorted(range(len(pares)), key=lambda i: pares[i]["p_gwas"])
        L += ["", "| rsID | pos b37 | p GWAS | p eQTL | SNP.PP.H4 |",
              "|---|---|---|---|---|"]
        for i in ordem[:10]:
            p = pares[i]
            L.append(f"| {p['qtl']['rsids'][0]} | {p['pos']:,} | "
                     f"{p['p_gwas']:.2e} | {p['qtl']['p']:.2e} | "
                     f"{res['snp_pp_h4'][i]:.4f} |")

        conc = sum(1 for (bg, _), (be, _) in zip(tg, te) if bg * be > 0)
        frac = conc / len(tg) if tg else 0.0
        il = ordem[0]
        z_max = max(abs(b / s) for b, s in te)
        L += ["", "### Concordância de direção",
              f"- Efeitos na mesma direção (GWAS × eQTL, alelo de efeito = ALT): "
              f"**{conc}/{len(tg)} ({frac:.0%})**",
              f"- Lead do GWAS ({pares[il]['qtl']['rsids'][0]}): "
              f"β_GWAS = {tg[il][0]:+.3f} · β_eQTL = {te[il][0]:+.4f} · "
              f"{'mesma direção' if tg[il][0] * te[il][0] > 0 else 'direções opostas'}",
              "", "### Nota de poder",
              f"- Com n = {reg['n_amostras']} doadores, o maior |z| do eQTL de "
              f"STX6 neste tecido é {z_max:.1f} (contra |z| ≈ 5,8 do lead do "
              "GWAS). O sinal de eQTL isolado é fraco, e por isso a posterior "
              "se divide entre H1 (só o GWAS) e H4. Este tecido sozinho é "
              "inconclusivo; a análise com poder é a meta de cinco conjuntos "
              "(`relatorio_coloc_meta_stx6.md`).", ""]

        if tecido == "DLPFC":
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            xs = [p["pos"] for p in pares]
            yg = [-math.log10(p["p_gwas"]) for p in pares]
            ye = [-math.log10(p["qtl"]["p"]) for p in pares]
            fig, ax = plt.subplots(figsize=(9, 5))
            ax.scatter(xs, yg, c="#c0392b", s=14, label="GWAS sCJD")
            ax.scatter(xs, ye, c="#2980b9", s=14, alpha=0.7,
                       label="eQTL STX6 (DLPFC)")
            lead = max(zip(xs, yg), key=lambda t: t[1])
            ax.axvline(lead[0], ls=":", c="grey", alpha=0.6)
            ax.set_xlabel("Posição GRCh37 (chr1)")
            ax.set_ylabel("-log10(p)")
            ax.set_title("STX6: GWAS sCJD × eQTL cerebral (GTEx DLPFC)")
            ax.legend()
            fig.tight_layout()
            FIGS.mkdir(parents=True, exist_ok=True)
            fig.savefig(FIGS / "coloc_stx6_regional.png", dpi=150)
            plt.close(fig)

    L += ["## Leitura", "",
          "| Tecido | n | variantes | H1 | H3 | H4 | H4/(H3+H4) |",
          "|---|---|---|---|---|---|---|"]
    for tecido, n, nv, pp in resumo:
        L.append(f"| {tecido} | {n} | {nv} | {pp[1]:.4f} | {pp[3]:.4f} | "
                 f"{pp[4]:.4f} | {C.h4_condicional(pp):.4f} |")
    L += ["",
          "- Cada tecido do GTEx, sozinho, tem pouco poder para o eQTL de STX6: "
          "a posterior fica dividida entre H1 (associação só com o GWAS) e H4 "
          "(variante causal compartilhada), com H3 pequeno. Nenhum dos dois "
          "tecidos, isolado, sustenta nem refuta a colocalização.",
          "- Os dois tecidos do GTEx vêm, em grande parte, dos mesmos doadores: "
          "não são réplicas independentes.",
          "- H4 alto significaria suporte a uma variante causal compartilhada "
          "sob os pressupostos do modelo (um único sinal causal por traço na "
          "região). Colocalização não testa mediação: não mostra que a "
          "expressão de STX6 esteja no caminho causal entre a variante e a "
          "doença.",
          "- DLPFC e cerebelo não representam todas as regiões afetadas na DCJ.",
          ""]
    destino = REPORTS / "relatorio_coloc_stx6.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[3/3] {destino}")


if __name__ == "__main__":
    main()
