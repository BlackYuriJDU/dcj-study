# -*- coding: utf-8 -*-
"""
test_analises_expressao.py — miRNA de sangue (GSE140069) e expressão cerebral
(GSE160208), RECALCULADOS a partir dos dados versionados.

Cada contagem citada nos relatórios e no manuscrito é obtida aqui rodando o
código de análise; as referências independentes são numpy, scipy e
statsmodels. Âncoras: seções 2.4 e 2.5 do dossiê de auditoria de 2026-11.
"""
import gzip
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
import analise_gse140069 as M  # noqa: E402
import analise_gse160208 as G  # noqa: E402

REPORTS = BASE / "pipeline" / "reports"


# ====================================================== GSE140069 (miRNA)
@pytest.fixture(scope="module")
def mirna():
    return M.analisar()


def test_amostra_e_covariaveis(mirna):
    d = mirna["dados"]
    assert len(d["nomes"]) == 105
    assert sum(d["grupo"]) == 57 and d["grupo"].count(0) == 48
    assert d["sexo"].count(0) == 60 and d["sexo"].count(1) == 45
    assert None not in d["sexo"] and None not in d["idade"] and None not in d["rin"]
    assert len(d["ok_ajust"]) == 105
    idade = np.array(d["idade"])
    g = np.array(d["grupo"])
    assert idade[g == 1].mean() == pytest.approx(66.4, abs=0.05)
    assert idade[g == 0].mean() == pytest.approx(53.6, abs=0.05)
    assert len(d["mirnas"]) == 939 == len(set(d["mirnas"]))
    assert sum(d["det"]) == 269


def test_matriz_de_desenho_posto_e_condicao(mirna):
    d = mirna["dados"]
    X = np.array([[1.0, d["grupo"][i], d["sexo"][i], d["idade"][i], d["rin"][i]]
                  for i in d["ok_ajust"]])
    assert mirna["desenho"].posto == np.linalg.matrix_rank(X) == 5
    assert np.linalg.cond(X) == pytest.approx(471, abs=1)


def test_74_mirnas_constantes_nao_sao_estimaveis(mirna):
    d = mirna["dados"]
    constantes = {m for m, l in zip(d["mirnas"], d["vals"]) if max(l) == min(l)}
    assert len(constantes) == 74
    assert set(mirna["nao_estimavel"]) == constantes
    assert set(mirna["nao_estimavel"].values()) == {"resposta_constante"}
    assert len(mirna["pB"]) == len(mirna["pA"]) == 865
    assert not (constantes & set(mirna["pB"]))


def test_ols_ajustado_bate_com_statsmodels(mirna):
    """β, p de TODOS os 865 testes estimáveis contra statsmodels.OLS."""
    d = mirna["dados"]
    ok = d["ok_ajust"]
    X = np.array([[1.0, d["grupo"][i], d["sexo"][i], d["idade"][i], d["rin"][i]]
                  for i in ok])
    pior_b, pior_p = 0.0, 0.0
    for m, linha in zip(d["mirnas"], d["vals_log"]):
        if m not in mirna["pB"]:
            continue
        r = sm.OLS(np.array([linha[i] for i in ok]), X).fit()
        pior_b = max(pior_b, abs(r.params[1] - mirna["betaB"][m])
                     / max(1e-12, abs(r.params[1])))
        pior_p = max(pior_p, abs(r.pvalues[1] - mirna["pB"][m])
                     / max(1e-300, r.pvalues[1]))
    assert pior_b < 1e-7
    assert pior_p < 1e-6


def test_welch_bate_com_scipy(mirna):
    d = mirna["dados"]
    ic, it = mirna["idx_cjd"], mirna["idx_ct"]
    for m, linha in list(zip(d["mirnas"], d["vals_log"]))[::7]:
        if m not in mirna["pA"]:
            continue
        ref = stats.ttest_ind([linha[i] for i in ic], [linha[i] for i in it],
                              equal_var=False).pvalue
        assert mirna["pA"][m] == pytest.approx(ref, rel=1e-7, abs=1e-300)


def test_fdr_bate_com_statsmodels(mirna):
    ms = sorted(mirna["pB"])
    q_ref = multipletests([mirna["pB"][m] for m in ms], method="fdr_bh")[1]
    for m, q in zip(ms, q_ref):
        assert mirna["qB"][m] == pytest.approx(q, rel=1e-10)


def test_contagens_no_universo_declarado_e_no_legado(mirna):
    """Universo declarado = 865 estimáveis. Legado = 939 (até a v0.4)."""
    det = mirna["detectavel"]
    # legado: os números publicados até a v0.4 (âncoras da seção 2.5)
    assert len(mirna["sigA_939"]) == 84
    assert len(mirna["sigA_939"] & det) == 69
    assert mirna["sigB_939"] == {"hsa-miR-500a-3p"}
    assert mirna["qB_939"]["hsa-miR-500a-3p"] == pytest.approx(0.00754, abs=5e-6)
    # declarado
    assert len(mirna["sigA"]) == 89
    assert len(mirna["sigA"] & det) == 73
    assert mirna["sigB"] == {"hsa-miR-500a-3p"}
    assert mirna["sigB_det"] == {"hsa-miR-1301-3p", "hsa-miR-28-3p",
                                 "hsa-miR-500a-3p", "hsa-miR-502-3p",
                                 "hsa-miR-93-5p"}
    assert mirna["betaB"]["hsa-miR-500a-3p"] == pytest.approx(-7.9358, abs=5e-5)
    assert mirna["pB"]["hsa-miR-500a-3p"] == pytest.approx(8.03e-6, rel=5e-3)
    assert mirna["qB_det"]["hsa-miR-93-5p"] == pytest.approx(0.048, abs=5e-4)


def test_mir500a_e_do_tipo_detectado_nao_detectado(mirna):
    """M5: o β de −7,9 vem da taxa de detecção (27/57 × 46/48), não de uma
    queda de abundância de ~250 vezes."""
    d = mirna["dados"]
    k = d["mirnas"].index("hsa-miR-500a-3p")
    dc = sum(1 for i in mirna["idx_cjd"] if d["vals"][k][i] > M.PISO)
    dt = sum(1 for i in mirna["idx_ct"] if d["vals"][k][i] > M.PISO)
    assert (dc, dt) == (27, 46)
    txt = (REPORTS / "relatorio_gse140069.md").read_text(encoding="utf-8")
    assert "27/57" in txt and "46/48" in txt
    assert "detectado/não detectado" in txt


def test_categorias_inesperadas_levantam_erro(tmp_path, monkeypatch):
    """M3: grupo ou sexo fora do esperado, ou idade não numérica, interrompem
    a análise em vez de virar 0 em silêncio."""
    with gzip.open(M.MATRIX, "rt", encoding="utf-8") as fh:
        original = fh.read()

    def com(troca_de, troca_para):
        assert troca_de in original
        p = tmp_path / "m.txt.gz"
        with gzip.open(p, "wt", encoding="utf-8") as fh:
            fh.write(original.replace(troca_de, troca_para, 1))
        monkeypatch.setattr(M, "MATRIX", p)

    com('"disease status: sCJD"', '"disease status: vCJD"')
    with pytest.raises(ValueError, match="grupo inesperado"):
        M.carregar_covariatas()
    com('"Sex: Female"', '"Sex: unknown"')
    with pytest.raises(ValueError, match="sexo inesperado"):
        M.carregar_covariatas()
    com('"age at sampling: ', '"age at sampling: ~')
    with pytest.raises(ValueError, match="idade não numérico"):
        M.carregar_covariatas()
    monkeypatch.undo()
    assert len(M.carregar_covariatas()) == 105


def test_fdr_rejeita_identificadores_repetidos():
    with pytest.raises(ValueError):
        M.fdr_bh([("a", 0.01), ("a", 0.02)])


def test_relatorio_mirna_cita_as_contagens_recalculadas(mirna):
    txt = (REPORTS / "relatorio_gse140069.md").read_text(encoding="utf-8")
    det = mirna["detectavel"]
    assert "Universo declarado: 865 testes estimáveis" in txt
    assert "74 miRNAs não são estimáveis" in txt
    assert f"| A: Welch log2(x+1), SEM covariáveis | {len(mirna['sigA'])} " in txt
    assert f"| {len(mirna['sigA_939'])} (" in txt
    assert f"| A′: A ∩ filtro de detecção | {len(mirna['sigA'] & det)} | " \
           f"{len(mirna['sigA_939'] & det)} |" in txt
    assert f"| B: OLS ajustado (idade+sexo+RIN) | **{len(mirna['sigB'])}**" in txt


def test_relatorio_idade_e_descritivo():
    """M4: contagens por modelo, sem 'a idade explica 97%'."""
    txt = (REPORTS / "relatorio_gse140069_idade.md").read_text(encoding="utf-8")
    for n in ("| idade sozinha (sem grupo) | 121 | — | 113 |",
              "| RIN sozinho (sem grupo) | 93 | — | 91 |",
              "| grupo sozinho | 92 | base | 88 |",
              "| **grupo + idade** | **4** | **-88** | 4 |",
              "| **grupo + RIN** | **28** | **-64** | 26 |"):
        assert n in txt, n
    assert "1.37" in txt and "1.32" in txt                      # VIF
    assert "não aditivas" in txt
    assert "A idade é o confundidor responsável" not in txt
    assert "accounts for 97%" not in txt
    assert "Essa leitura foi retirada" in txt


# ================================================= GSE160208 (expressão cerebral)
@pytest.fixture(scope="module")
def cerebro():
    amostras, genes, vals, _meta, _covs = G.extrair_tabela()
    return amostras, genes, vals, G.testar_fc(amostras, genes, vals)


def test_composicao_das_amostras(cerebro):
    amostras, genes, _vals, _res = cerebro
    g = G.classificar_amostras(amostras)
    assert {k: len(v) for k, v in g.items()} == \
        {"CJD_FC": 14, "CT_FC": 10, "CJD_CB": 13, "CT_CB": 10}
    assert len(genes) == 800 == len(set(genes))
    with pytest.raises(ValueError):
        G.classificar_amostras(["amostra_sem_grupo"])


def test_58_genes_sem_variancia_nao_sao_estimaveis(cerebro):
    amostras, genes, vals, res = cerebro
    g = G.classificar_amostras(amostras)
    fc = g["CJD_FC"] + g["CT_FC"]
    constantes = {gn for gn, l in zip(genes, vals)
                  if max(l[i] for i in fc) == min(l[i] for i in fc)}
    assert len(constantes) == 58
    assert {r["gene"] for r in res if not r["estimavel"]} == constantes
    assert all(r["p"] is None for r in res if not r["estimavel"])


def test_welch_e_delta_batem_com_scipy(cerebro):
    amostras, genes, vals, res = cerebro
    g = G.classificar_amostras(amostras)
    for r, l in zip(res, vals):
        if not r["estimavel"]:
            continue
        xs, ys = [l[i] for i in g["CJD_FC"]], [l[i] for i in g["CT_FC"]]
        assert r["p"] == pytest.approx(
            stats.ttest_ind(xs, ys, equal_var=False).pvalue, rel=1e-7)
        assert r["delta"] == pytest.approx(np.mean(xs) - np.mean(ys), abs=1e-12)


def test_contagens_fdr_nos_dois_universos(cerebro):
    res = cerebro[3]
    up, down = G.significativos(res, "painel")
    assert (len(up), len(down)) == (314, 123)                 # 437, até a v0.4
    up, down = G.significativos(res, "estimaveis")
    assert (len(up), len(down)) == (315, 126)                 # 441, declarado
    # referência independente para o universo declarado
    est = [r for r in res if r["estimavel"]]
    q = multipletests([r["p"] for r in est], method="fdr_bh")[1]
    assert int((q < 0.05).sum()) == 441


def test_184_depende_da_ordem_dos_filtros(cerebro):
    """E1: 184 só com magnitude antes do FDR (ou sem exigir q); 183 com o FDR
    no universo inteiro. O gene da fronteira é CCL4."""
    res = cerebro[3]
    a = G.degs_criterio_artigo(res, "fdr_no_universo", "estimaveis")
    b = G.degs_criterio_artigo(res, "fdr_no_universo", "painel")
    c = G.degs_criterio_artigo(res, "magnitude_antes_do_fdr")
    assert (len(a), len(b), len(c)) == (183, 183, 184)
    assert c - a == {"CCL4"} and a <= c
    grandes = [r for r in res if r["estimavel"] and abs(r["delta"]) > 1]
    assert len(grandes) == 185
    assert sum(1 for r in grandes if r["p"] < 0.05) == 184
    with pytest.raises(ValueError):
        G.degs_criterio_artigo(res, "ordem_que_nao_existe")
    with pytest.raises(ValueError):
        G.q_no_universo(res, "universo_que_nao_existe")


def test_student_no_painel_da_410(cerebro):
    """Âncora da seção 2.4 (só C): Student + BH sobre 800 → 410."""
    amostras, _genes, vals, res = cerebro
    g = G.classificar_amostras(amostras)
    ps = []
    for r, l in zip(res, vals):
        if r["estimavel"]:
            ps.append(stats.ttest_ind([l[i] for i in g["CJD_FC"]],
                                      [l[i] for i in g["CT_FC"]]).pvalue)
        else:
            ps.append(1.0)
    assert int((multipletests(ps, method="fdr_bh")[1] < 0.05).sum()) == 410


def test_bh_por_posicao_bate_com_statsmodels():
    rng = np.random.default_rng(3)
    ps = rng.random(200).tolist() + [0.01, 0.01, 1.0]
    assert G.bh(ps) == pytest.approx(
        multipletests(ps, method="fdr_bh")[1].tolist(), rel=1e-12)


def test_relatorios_do_cerebro_citam_o_recalculado(cerebro):
    res = cerebro[3]
    up, down = G.significativos(res, "estimaveis")
    up_p, down_p = G.significativos(res, "painel")
    txt = (REPORTS / "relatorio_gse160208.md").read_text(encoding="utf-8")
    assert f"**{len(up) + len(down)}** ({len(up)} ↑, {len(down)} ↓)" in txt
    assert f"**{len(up_p) + len(down_p)}** ({len(up_p)} ↑, {len(down_p)} ↓)" in txt
    assert "não estimáveis: 58" in txt
    assert "EXATO" not in txt and "número idêntico ao publicado" not in txt
    val = (REPORTS / "validacao_cruzada_gse160208_artigo_original.md").read_text(
        encoding="utf-8")
    assert "447 signific" not in val
    assert "Réplica ponto a ponto" not in val
    assert "o número exato publicado" not in val
    script = (BASE / "pipeline" / "scripts" / "analise_gse160208.py").read_text(
        encoding="utf-8")
    assert "virá em versão futura" not in script                # E3
    assert math.isclose(G.DELTA_MIN, 1.0)
