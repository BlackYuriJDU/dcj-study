# -*- coding: utf-8 -*-
"""
test_textos.py — Manuscrito, READMEs e site citam os números que o código
produz AGORA.

Cada valor é recalculado aqui (colocalização a partir dos caches de eQTL e das
390 variantes versionadas do GWAS; expressão e miRNA a partir dos dados do
GEO) e só então procurado no texto. Um número antigo que continue escrito, ou
um número novo digitado errado, derruba o teste. Guardas de redação no fim:
colocalização não é mediação; nada de "réplica exata"; nada de "97%".
"""
import re
import sys
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
sys.path.insert(0, str(BASE / "tests"))
import coloc_abf as C  # noqa: E402
import harmoniza as H  # noqa: E402
import analise_gse140069 as M  # noqa: E402
import analise_gse160208 as G  # noqa: E402
import analise_integrada_cerebro_sangue as I  # noqa: E402
import coloc_meta_stx6 as META  # noqa: E402
import coloc_sqtl_stx6 as SQ  # noqa: E402
from test_coloc_pipeline import mini_gwas  # noqa: E402

SITE = BASE / "Litho Foundation Website"


def ler(rel):
    return (BASE / rel).read_text(encoding="utf-8")


def cientifico(x):
    """3.44e-47 → '3.44×10⁻⁴⁷' (notação do manuscrito)."""
    mant, exp = f"{x:.2e}".split("e")
    sobrescrito = str(int(exp)).translate(str.maketrans("-0123456789",
                                                        "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{mant}×10{sobrescrito}"


@pytest.fixture(scope="module")
def ms():
    return ler("preprint/manuscrito_preprint.md")


@pytest.fixture(scope="module")
def coloc():
    gwas, _ = mini_gwas()
    return META.analisar(gwas), SQ.analisar(gwas)


# ------------------------------------------------------------ colocalização
def test_manuscrito_cita_a_colocalizacao_recalculada(ms, coloc):
    a, _ = coloc
    pp = a["res"]["pp"]
    assert f"**H₃ = {pp[3]:.3f} and H₄ = {pp[4]:.3f}**" in ms
    assert f"H₄ = {pp[4]:.4f}, H₃ = {pp[3]:.4f}" in ms
    s = {p12: x for p12, x in C.sensibilidade_p12(a["tg"], a["te"], "cc", "quant")}
    assert (f"{s[1e-5][4]:.3f} at p₁₂ = 10⁻⁵, {s[5e-6][4]:.3f} at 5×10⁻⁶ and "
            f"{s[1e-6][4]:.3f} at 10⁻⁶") in ms
    assert f"({s[1e-6][4]:.2f} with" in ms            # 0.83 no resumo/nota
    alt = C.coloc(a["tg"], a["te"], "cc", "quant",
                  priors=META.PRIORS_ALTERNATIVOS)["pp"]
    assert f"it gives {alt[4]:.4f} and is not used as evidence" in ms
    sens = {n: h4 for n, h4, _h3, _n in a["sens"]}
    assert f"is {sens['Efeitos aleatórios (DerSimonian–Laird)']:.3f} with " \
           "DerSimonian–Laird" in ms
    assert f"{sens['Sem os dois GTEx (3 conjuntos, mínimo 3)']:.3f} without both GTEx" in ms
    so = {n: h4 for n, h4 in sens.items() if n.startswith("Só GTEx")}
    dl, cb = (so[n] for n in sorted(so, key=lambda n: "cerebelo" in n))
    assert f"give {dl:.2f} (DLPFC, n = 268) and {cb:.2f}" in ms
    conc = sum(1 for (bg, _), (be, _) in zip(a["tg"], a["te"]) if bg * be > 0)
    assert f"({conc}/{len(a['pares'])})" in ms and "89% effect-direction" in ms
    snp = a["res"]["snp_pp_h4"]
    ordem = sorted(range(len(snp)), key=lambda i: -snp[i])
    topo = a["pares"][ordem[0]]
    lead = next(i for i, p in enumerate(a["pares"])
                if p["pos"] == META.LEAD_GWAS_B37)
    assert f"contains {len(C.conjunto_credivel(snp))} variants" in ms
    assert (f"(chr1:{topo['pos']:,}; {snp[ordem[0]]:.3f}) is not the GWAS lead "
            f"({snp[lead]:.3f})") in ms
    n = a["res"]["n"]
    assert f"≈ {C.razao_h4_h3_priori(n, *C.PRIORS_PADRAO):.1f}; the observed " \
           f"posterior ratio is ≈ {pp[4] / pp[3]:.0f}" in ms
    q05 = sum(1 for p in a["pares"] if p["qtl"]["p_Q"] < 0.05)
    i50 = sum(1 for p in a["pares"] if p["qtl"]["I2"] > 0.5)
    assert f"p < 0.05 at {q05} of the {n}" in ms and f"I² > 50% at {i50})" in ms


def test_manuscrito_p_do_eqtl_e_minimo_regional(ms, coloc):
    a, _ = coloc
    ps = [(H.p_de_z(p["qtl"]["beta"] / p["qtl"]["se"]), p["pos"])
          for p in a["pares"]]
    pmin, pos_min = min(ps)
    lead = next(p for p, pos in ps if pos == META.LEAD_GWAS_B37)
    assert cientifico(lead) == "7.60×10⁻⁴⁷" and cientifico(lead) in ms
    assert f"chr1:{pos_min:,} (p = {cientifico(pmin)})" in ms
    outros = sum(1 for p, _ in ps if p <= 3 * pmin) - 1
    assert outros == 13 and "thirteen\nfurther variants" in ms
    # o valor antigo atribuído ao "mínimo regional" não era o mínimo
    assert "6.60×10⁻⁴⁷" not in ms


def test_manuscrito_tamanhos_amostrais(ms):
    assert "(2,156 samples)" in ms and "2,156 brain samples" in ms
    assert "GTEx v10 DLPFC (n=268) and GTEx v10 cerebellum (n=263)" in ms
    for velho in ("n=285", "n=272"):
        assert velho not in ms
    for m in re.finditer(r"2,182", ms):               # só na nota de correção
        assert "not ~2,182" in ms[m.start() - 10: m.end() + 5]


def test_manuscrito_sqtl(ms, coloc):
    _, sq = coloc
    t = sq["testes"]
    assert len(t) == 6 and sum(ns for _, _, ns in sq["resumo"]) == 0
    assert f"The highest is {t[0]['pp'][4]:.2f}" in ms
    assert f"it falls to {t[0]['h4_p12_1e6']:.2f} with p₁₂ = 10⁻⁶" in ms
    assert f"only {t[0]['conc']:.0%} of variants" in ms
    mesmos = [x for x in t[1:] if x["cluster"].split(":")[:2]
              == t[0]["cluster"].split(":")[:2]]
    assert len(mesmos) == 3
    h = sorted((x["pp"][4] for x in mesmos), reverse=True)
    assert f"({h[0]:.2f}, {h[1]:.2f} and below 0.001)" in ms and h[2] < 0.001
    assert "Bonferroni-aware" not in ms
    assert "not evidence of absence" in ms


# ------------------------------------------------------------ expressão e miRNA
def test_manuscrito_expressao_cerebral(ms):
    res = G.testar_fc()
    up, down = G.significativos(res, "estimaveis")
    up_p, down_p = G.significativos(res, "painel")
    n_est = sum(1 for r in res if r["estimavel"])
    assert (f"{len(up) + len(down)} of the {n_est} estimable genes are FDR < 0.05 "
            f"({len(up)} up, {len(down)} down; {len(up_p) + len(down_p)} when") in ms
    a = len(G.degs_criterio_artigo(res, "magnitude_antes_do_fdr"))
    b = len(G.degs_criterio_artigo(res, "fdr_no_universo"))
    assert f"it is {a}, the published count" in ms and f"and {b} when the FDR" in ms
    assert "(CCL4 is the borderline gene)" in ms
    assert "exact replication of the original authors' criteria" not in ms
    assert "replicate exactly" not in ms and "identical count" not in ms


def test_manuscrito_mirna(ms):
    r = M.analisar()
    lfc = r["lfc"]
    up = sum(1 for m in r["sigA"] if lfc[m] > 0)
    up9 = sum(1 for m in r["sigA_939"] if lfc[m] > 0)
    assert (f"yields {len(r['sigA'])} FDR-significant miRNAs ({up}↑/"
            f"{len(r['sigA']) - up}↓; {len(r['sigA_939'])}, {up9}↑/"
            f"{len(r['sigA_939']) - up9}↓, in the v0.4 universe of 939)") in ms
    assert len(r["sigB"]) == len(r["sigB_939"]) == 1
    assert "leaves **1** (hsa-miR-500a-3p) in either universe" in ms
    assert "74 of the 939 miRNAs are constant" in ms
    assert "27 of 57 cases and 46 of 48 controls" in ms
    txt = ler("pipeline/reports/relatorio_gse140069_idade.md")
    linha = re.search(r"\| grupo sozinho \| (\d+) \| base \| (\d+) \|", txt)
    assert f"the group-only model gives {linha.group(1)} FDR-significant" in ms
    for rot, frase in (("**grupo + idade**", "adding age alone leaves {}"),
                       ("**grupo + RIN**", "RNA integrity alone {}"),
                       ("grupo + sexo", "sex alone {}")):
        m = re.search(rf"^\| {re.escape(rot)} \| \**(\d+)\**", txt, re.M)
        assert frase.format(m.group(1)) in ms, rot
    assert "(v0.4 universe: 88, 4, 26, 85 and 1)" in ms


def test_manuscrito_integracao(ms):
    a = I.analisar()
    t = [r for r in a["linhas"] if r["testavel"]]
    assert f"(N = {a['N']}, of which {a['n_up']} are up-regulated)" in ms
    assert (f"they have {t[0]['K']}, {t[1]['K']} and {t[2]['K']} targets on the "
            f"panel, of which {t[0]['k']}, {t[1]['k']} and {t[2]['k']} are "
            f"up-regulated (hypergeometric p = {t[0]['p']:.2f}, {t[1]['p']:.2f} "
            f"and {t[2]['p']:.2f})") in ms
    assert "reflects peripheral processes" not in ms
    assert "all q = 1.0" not in ms


# ------------------------------------------------------------ guardas de redação
def test_manuscrito_v05_com_nota_de_correcao(ms):
    cab = ms[:400]
    assert "Preprint draft v0.5" in cab and "v0.4 —" not in cab
    assert "Arthur Araújo" in cab
    nota = ms[: ms.index("## Abstract")]
    assert "Correction notice (v0.5, replaces v0.4)" in nota
    for trecho in ("p₁·p₂·p₁₂ instead of p₁₂", "two copies of the same wrong",
                   "CORRECOES.md"):
        assert trecho in nota, trecho
    assert (BASE / "CORRECOES.md").exists()


def test_colocalizacao_nao_e_apresentada_como_mediacao(ms):
    assert ms.count("not evidence of mediation") >= 2
    assert "does not\ntest mediation" in ms or "does not test mediation" in ms
    assert "we make no claim of\nnovelty" in ms or "no claim of novelty" in ms
    for proibido in ("expression-mediated", "mediated by STX6 expression",
                     "the mechanism underlying the association",
                     "cannot split H₃", "H₂ and H₃ are hypotheses of two",
                     "H₃ = 0.994997", "prior W = 0.04²"):
        assert proibido not in ms, proibido
    # toda menção ao resultado antigo vem marcada como correção
    for m in re.finditer(r"0\.995", ms):
        janela = ms[max(0, m.start() - 400): m.end() + 200]
        assert "wrong" in janela or "erroneous" in janela


def test_idade_sem_atribuicao_percentual(ms):
    assert "fully explains" not in ms
    assert "variance-inflation" in ms
    for m in re.finditer(r"97%", ms):
        assert "withdrawn" in ms[m.start(): m.end() + 200], \
            "'97%' fora da nota que retira a afirmação"
    assert "do not measure a proportion of the effect attributable to age" in ms


def test_promessas_de_reprodutibilidade(ms):
    en = ler("README.en.md")
    assert "Everything here is reproducible from standard-library Python" not in en
    assert "pure stdlib" not in en
    assert "requirements.txt" in en and "make dados-gwas" in en
    assert "standard-library statistical core" not in ms
    assert "cross-validated against R anchors" not in ms
    assert "obter_gwas.py" in ms
    assert "Pure standard-library Python" not in ler("CITATION.cff")


def test_readmes_citam_os_numeros_corrigidos(coloc):
    a, _ = coloc
    pp = a["res"]["pp"]
    s = {p12: x for p12, x in C.sensibilidade_p12(a["tg"], a["te"], "cc", "quant")}
    en, pt = ler("README.en.md"), ler("README.md")
    assert f"PP.H4 = {pp[4]:.2f}" in en and f"{s[1e-6][4]:.2f} with p12 = 1e-6" in en
    assert f"PP.H4 = {pp[4]:.2f}".replace("0.", "0,") in pt
    assert "2,156 samples" in en and "2.156 amostras" in pt
    assert "441 of 742 estimable genes" in en and "89 → 1 significant" in en
    assert "exact replication" not in en
    assert "CORRECOES.md" in en and "CORRECOES.md" in pt
    assert "Colocalization is not mediation" in en
    assert "Colocalização não é mediação" in pt


def test_site_cita_os_numeros_corrigidos():
    res = G.testar_fc()
    up, down = G.significativos(res, "estimaveis")
    n = len(up) + len(down)
    dados = (SITE / "dados.html").read_text(encoding="utf-8")
    js = (SITE / "app.js").read_text(encoding="utf-8")
    proc = (SITE / "processo.html").read_text(encoding="utf-8")
    idx = (SITE / "index.html").read_text(encoding="utf-8")
    assert f">{n} DEGs</text>" in dados
    for txt in (dados, js, proc, idx):
        baixo = txt.lower()
        for alegacao in ("réplica exata do artigo", "réplica exata:",
                         "r da réplica exata", "réplica exata com r",
                         "exact replication of the gse160208",
                         "r of the exact replication",
                         "exact replication with r"):
            assert alegacao not in baixo, alegacao
        assert "437 DEGs cerebrais" not in txt and "437 brain DEGs" not in txt
        assert "60 miRNAs significantes" not in txt
        assert "60 significant miRNAs" not in txt
        assert "o sangue reflete a periferia" not in txt
        assert "blood reflects the periphery" not in txt
    r = M.analisar()
    assert f"{len(r['sigA'])} miRNAs sem ajuste → {len(r['sigB'])} com ajuste" in dados
    assert f"<b>{len(r['sigA'])}→{len(r['sigB'])}</b>" in idx
