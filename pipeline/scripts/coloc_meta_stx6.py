#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
coloc_meta_stx6.py — Colocalização com poder ampliado: meta-análise por
variância inversa (efeito fixo) dos eQTLs de STX6 em cinco conjuntos cerebrais
do eQTL Catalogue r8 (CommonMind, ROSMAP, BrainSeq, GTEx DLPFC, GTEx cerebelo;
tamanhos amostrais em pipeline/data/proveniencia_qtl.tsv) contra o GWAS de sCJD
(GCST90001389, n=17.679).

Variante a variante — identidade completa (posição, REF, ALT), efeito = ALT:
    β_meta = Σ(β/EP²)/Σ(1/EP²);  EP_meta = 1/√Σ(1/EP²);  mínimo de 3 conjuntos.
Só alelos idênticos são combinados (pipeline/lib/harmoniza.py). Depois,
colocalização (pipeline/lib/coloc_abf.py, conforme o pacote R `coloc`) do
eQTL-meta contra o GWAS, com tabelas de sensibilidade a p12, ao prior de efeito
W, ao modelo de meta-análise e à composição dos conjuntos.

Saída: pipeline/reports/relatorio_coloc_meta_stx6.md
"""
import sys
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coloc_abf as C  # noqa: E402
import harmoniza as H  # noqa: E402
from fontes_qtl import OFFSET_B38_B37, datasets, ler_eqtl  # noqa: E402
from coloc_stx6_eqtl import carregar_gwas_regiao  # noqa: E402

REPORTS = BASE / "reports"
N_GWAS = 17_679
LEAD_GWAS_B37 = 180_961_245
PRIORS_ALTERNATIVOS = (1e-6, 1e-6, 1e-7)   # antes rotulado "conservador"


def carregar_eqtl():
    """{dataset_id: variantes deduplicadas}, [(registro, stats, sdY_estimado)]."""
    por_dataset, resumo = {}, []
    for reg in datasets("ge"):
        variantes, st = H.deduplicar(ler_eqtl(reg))
        por_dataset[reg["dataset_id"]] = variantes
        sdy = C.estimar_sdY([v["se"] for v in variantes.values()],
                            [v["maf"] for v in variantes.values()],
                            reg["n_amostras"])
        resumo.append((reg, st, sdy))
    return por_dataset, resumo


def casar(meta, gwas):
    """Casa a meta com o GWAS; devolve (pares, stats, traço GWAS, traço eQTL)."""
    pares, st = H.casar_com_gwas(meta, gwas, OFFSET_B38_B37)
    tg = [(p["beta_gwas"], p["se_gwas"]) for p in pares]
    te = [(p["qtl"]["beta"], p["qtl"]["se"]) for p in pares]
    return pares, st, tg, te


def h4_de(por_dataset, gwas, aleatorio=False, restringir_a=None,
          min_datasets=H.MIN_DATASETS):
    """(PP.H4, PP.H3, n) de uma configuração de sensibilidade."""
    meta, _ = H.meta_ivw(por_dataset, min_datasets=min_datasets,
                         aleatorio=aleatorio)
    if restringir_a is not None:
        meta = {k: v for k, v in meta.items() if k in restringir_a}
    pares, _, tg, te = casar(meta, gwas)
    if not pares:
        return float("nan"), float("nan"), 0
    pp = C.coloc(tg, te, "cc", "quant")["pp"]
    return pp[4], pp[3], len(pares)


def analisar(gwas=None):
    """Executa a análise inteira e devolve um dict com tudo o que o relatório
    e a validação cruzada precisam."""
    if gwas is None:
        gwas = carregar_gwas_regiao()
    por_dataset, resumo = carregar_eqtl()
    meta, st_meta = H.meta_ivw(por_dataset)
    pares, st_gwas, tg, te = casar(meta, gwas)
    res = C.coloc(tg, te, "cc", "quant", sdY2=1.0)
    chaves = {(p["pos_qtl"], p["ref"], p["alt"]) for p in pares}

    sens = []
    sens.append(("Efeito fixo (análise principal)",) + (res["pp"][4],
                res["pp"][3], len(pares)))
    sens.append(("Efeitos aleatórios (DerSimonian–Laird)",)
                + h4_de(por_dataset, gwas, aleatorio=True))
    ids = list(por_dataset)
    rot = {r["dataset_id"]: r["rotulo"] for r, _, _ in resumo}
    for fora in ids:
        sub = {k: v for k, v in por_dataset.items() if k != fora}
        sens.append((f"Sem {rot[fora]} (4 conjuntos, mínimo 3)",)
                    + h4_de(sub, gwas))
    gtex = [i for i in ids if rot[i].startswith("GTEx")]
    sub = {k: v for k, v in por_dataset.items() if k not in gtex}
    sens.append(("Sem os dois GTEx (3 conjuntos, mínimo 3)",) + h4_de(sub, gwas))
    for so in ids:
        sens.append((f"Só {rot[so]} (variantes da análise principal)",)
                    + h4_de({so: por_dataset[so]}, gwas, restringir_a=chaves,
                            min_datasets=1))
    return {"gwas": gwas, "por_dataset": por_dataset, "resumo": resumo,
            "meta": meta, "st_meta": st_meta, "pares": pares,
            "st_gwas": st_gwas, "tg": tg, "te": te, "res": res, "sens": sens}


def main() -> None:
    hoje = datetime.date.today().isoformat()
    print("[1/3] GWAS + eQTL + meta-análise...")
    a = analisar()
    resumo, meta, pares = a["resumo"], a["meta"], a["pares"]
    st_meta, st_gwas, tg, te, res = (a["st_meta"], a["st_gwas"], a["tg"],
                                     a["te"], a["res"])
    pp = res["pp"]
    n_total = sum(r["n_amostras"] for r, _, _ in resumo)
    w_cc, w_q = res["w1"], res["w2"]

    L = ["# Colocalização STX6 — meta-eQTL de cinco conjuntos cerebrais × GWAS sCJD",
         f"*`coloc_meta_stx6.py` em {hoje}. Meta-análise por variância inversa "
         f"(efeito fixo) de {len(resumo)} conjuntos do eQTL Catalogue r8 "
         f"({n_total:,}".replace(",", ".") + " amostras no total; "
         "proveniência e checksums em `pipeline/data/proveniencia_qtl.tsv`) "
         "contra o GWAS de sCJD (GCST90001389, n = 17.679). Colocalização: "
         "implementação única em `pipeline/lib/coloc_abf.py`, conforme o "
         "pacote R `coloc`.*", "",
         "> **Correção de 2026-11.** As versões até a v0.4 deste relatório "
         "davam PP.H4 ≈ 0 e H3 ≈ 0,995. O número vinha de um erro de fórmula: "
         "o termo de H4 era multiplicado por p1·p2·p12 em vez de p12 (um fator "
         "de 10⁻⁸ com os priors padrão). Com a combinação correta e os mesmos "
         "dados, PP.H4 ≈ 0,98. Ver `CORRECOES.md` e "
         "`relatorio_validacao_coloc_R.md`.", "",
         "## Conjuntos de eQTL", "",
         "| Conjunto | n | linhas no cache | variantes distintas | linhas "
         "repetidas (rsID) | posições multialélicas | indels | sdY estimado |",
         "|---|---|---|---|---|---|---|---|"]
    for reg, st, sdy in resumo:
        L.append(f"| {reg['rotulo']} ({reg['dataset_id']}) | "
                 f"{reg['n_amostras']} | {st['linhas']} | {st['variantes']} | "
                 f"{st['linhas_repetidas_por_rsid']} | "
                 f"{st['posicoes_multialelicas']} | {st['indels']} | "
                 f"{sdy:.2f} |")
    L += ["",
          "Os tamanhos amostrais vêm da metadata oficial r8 (GTEx: 268 e 263; "
          "total 2.156). Eles não entram no fator de Bayes, que é calculado a "
          "partir de β e erro-padrão. *sdY estimado* = "
          "√mediana(EP²·2n·MAF·(1−MAF)), a mesma estimativa do `coloc`; valores "
          "abaixo de 1 são esperados porque o EP reflete a variância residual "
          "depois das covariáveis. A análise usa sdY = 1 (expressão normalizada "
          "antes do mapeamento); a tabela de sensibilidade a W abaixo cobre o "
          "intervalo correspondente.", "",
          "## Harmonização", "",
          f"- Variantes distintas (posição+REF+ALT) em algum conjunto: "
          f"{st_meta['variantes_em_algum_dataset']} · presentes em menos de "
          f"{H.MIN_DATASETS} conjuntos (fora da meta): "
          f"{st_meta['abaixo_do_minimo']} · **na meta: {st_meta['na_meta']}** "
          f"(em {st_meta['posicoes_na_meta']} posições)",
          f"- Posições em que conjuntos diferentes trazem pares REF/ALT "
          f"diferentes: {st_meta['posicoes_com_alelos_discordantes']}. A versão "
          "anterior indexava só pela posição e somava esses alelos; agora cada "
          "alelo é uma variante própria e só alelos idênticos são combinados.",
          f"- Casadas com o GWAS pelos dois alelos: **{st_gwas['casadas']}** "
          f"(mesmo alelo de efeito: {st_gwas['casadas_mesmo_alelo_de_efeito']}; "
          f"sinal do GWAS invertido: {st_gwas['casadas_sinal_invertido']}; "
          f"palindrômicas entre as casadas: {st_gwas['palindromicas_casadas']})",
          f"- Não usadas: sem variante do GWAS na posição: "
          f"{st_gwas['sem_variante_no_gwas']} · alelos incompatíveis: "
          f"{st_gwas['alelos_incompativeis']} · só casariam na fita "
          f"complementar: {st_gwas['so_na_fita_complementar']}", ""]

    L += ["## Posteriores (priors padrão do `coloc`: p1 = p2 = 10⁻⁴, p12 = 10⁻⁵)",
          "",
          f"Prior de efeito: variância {w_cc:g} para o GWAS caso-controle e "
          f"{w_q:g} para o eQTL quantitativo (sdY = 1).", ""]
    L += C.md_tabela_hipoteses(pp, "o GWAS de sCJD", "o eQTL de STX6")
    razao_post = pp[4] / pp[3]
    razao_pri = C.razao_h4_h3_priori(res["n"], *C.PRIORS_PADRAO)
    L += ["",
          f"- **PP.H4 = {pp[4]:.4f}** · PP.H3 = {pp[3]:.4f} · n = {res['n']} "
          "variantes",
          f"- H3 + H4 (os dois traços têm sinal na região) = {pp[3] + pp[4]:.4f}",
          f"- H4/(H3+H4) = {C.h4_condicional(pp):.4f}", "",
          "## Sensibilidade a p12 (p1 e p2 fixos em 10⁻⁴)", "",
          "p12 é o único prior que separa H4 de H3; é a sensibilidade que "
          "importa.", ""]
    L += C.md_sensibilidade_p12(tg, te, "cc", "quant")
    pa = C.coloc(tg, te, "cc", "quant", priors=PRIORS_ALTERNATIVOS)["pp"]
    L += ["",
          f"O cenário p1 = p2 = 10⁻⁶, p12 = 10⁻⁷, que versões anteriores "
          f"chamavam de \"conservador\", dá H4 = {pa[4]:.4f} e H3 = "
          f"{C._fmt(pa[3])}. Ele **não** é mais rigoroso para H4: a razão "
          "p12/(p1·p2) sobe de 10³ para 10⁵, o que favorece a variante "
          "compartilhada. Fica aqui só para registro e não é usado como "
          "evidência.", "",
          "## Sensibilidade ao prior de efeito W", ""]
    L += C.md_sensibilidade_w(tg, te)

    L += ["", "## Sensibilidade ao modelo de meta-análise e à composição", "",
          "| Variação | H3 | H4 | n |", "|---|---|---|---|"]
    for nome, h4, h3, n in a["sens"]:
        L.append(f"| {nome} | {h3:.4f} | {h4:.4f} | {n} |")
    q05 = sum(1 for p in pares if p["qtl"]["p_Q"] < 0.05)
    i50 = sum(1 for p in pares if p["qtl"]["I2"] > 0.5)
    L += ["",
          f"Heterogeneidade entre conjuntos nas {len(pares)} variantes: Q de "
          f"Cochran com p < 0,05 em {q05}; I² > 50% em {i50}. Os dois tecidos "
          "do GTEx compartilham doadores e não são coortes independentes; "
          "sozinhos, cada um dá H4 bem menor (linhas \"Só GTEx\").", ""]

    # variantes
    snp = res["snp_pp_h4"]
    ordem = sorted(range(len(pares)), key=lambda i: -snp[i])
    cred = C.conjunto_credivel(snp, 0.95)
    conc = sum(1 for (bg, _), (be, _) in zip(tg, te) if bg * be > 0)
    L += ["## Variantes", "",
          f"- Mesma direção GWAS × eQTL-meta: **{conc}/{len(pares)} "
          f"({conc / len(pares):.0%})**",
          f"- Conjunto de 95% de credibilidade para a variante compartilhada "
          f"(dado H4): {len(cred)} variantes. A de maior posterior "
          f"(chr1:{pares[ordem[0]]['pos']:,}, SNP.PP.H4 = {snp[ordem[0]]:.4f}) "
          + ("é" if pares[ordem[0]]["pos"] == LEAD_GWAS_B37 else "**não** é")
          + " o lead do GWAS.", "",
          "| pos b37 | REF>ALT | p GWAS | p eQTL-meta | k | I² | SNP.PP.H4 |",
          "|---|---|---|---|---|---|---|"]
    for i in ordem[:10]:
        p = pares[i]
        q = p["qtl"]
        L.append(f"| {p['pos']:,} | {p['ref']}>{p['alt']} | {p['p_gwas']:.2e} | "
                 f"{H.p_de_z(q['beta'] / q['se']):.2e} | {q['k']} | "
                 f"{q['I2']:.0%} | {snp[i]:.4f} |")
    lead = next((p for p in pares if p["pos"] == LEAD_GWAS_B37), None)
    if lead is not None:
        q = lead["qtl"]
        z = q["beta"] / q["se"]
        L += ["",
              f"Lead do GWAS (chr1:{LEAD_GWAS_B37:,}): p GWAS = "
              f"{lead['p_gwas']:.2e}; eQTL-meta z = {z:.2f}, p = "
              f"{H.p_de_z(z):.2e}, k = {q['k']}, I² = {q['I2']:.1%}."]

    L += ["", "## Leitura", "",
          f"- **O que a análise sustenta.** Sob os pressupostos do modelo (no "
          "máximo um sinal causal por traço na região) e com os priors padrão "
          f"do `coloc`, a posterior de uma variante causal compartilhada entre "
          f"o risco de sCJD e a expressão de STX6 no cérebro é {pp[4]:.2f}.",
          f"- **Sensibilidade a p12.** Com p12 = 10⁻⁶ a posterior cai para "
          f"{C.coloc(tg, te, 'cc', 'quant', priors=(1e-4, 1e-4, 1e-6))['pp'][4]:.2f}"
          ". O resultado depende desse prior e deve ser lido junto com a "
          "tabela.",
          "- **Robustez numérica.** O valor muda pouco com o prior de efeito, "
          "com efeitos aleatórios, com a retirada de qualquer conjunto e com a "
          "retirada dos dois GTEx (tabelas acima).",
          f"- **Quanto vem do prior.** O bloco está em LD quase perfeito. Se as "
          f"{res['n']} variantes fossem equivalentes, a razão H4/H3 a priori "
          f"seria p12/(p1·p2·(n−1)) ≈ {razao_pri:.1f}; a razão a posteriori "
          f"observada é ≈ {razao_post:.0f}. Parte da separação entre H3 e H4 "
          "vem do prior, e os dados a ampliam em cerca de "
          f"{razao_post / razao_pri:.0f} vezes.",
          "- **O que a análise não sustenta.** Colocalização não é mediação: "
          "não mostra que a expressão de STX6 esteja no caminho causal entre a "
          "variante e a doença. O modelo não testa mais de um sinal causal por "
          "traço; o LD entre variantes não é modelado; os tecidos são DLPFC e "
          "cerebelo de adultos, não todas as regiões afetadas na DCJ.",
          "- **Não é resultado novo.** A associação entre o locus de STX6 e a "
          "expressão do gene em sCJD já foi publicada (Jones et al. 2020; "
          "estudo multiômico em *Brain* 2025). Este relatório reproduz esse "
          "resultado com dados públicos.", ""]
    destino = REPORTS / "relatorio_coloc_meta_stx6.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[3/3] {destino}")


if __name__ == "__main__":
    main()
