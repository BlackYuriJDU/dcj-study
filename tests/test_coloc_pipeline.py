# -*- coding: utf-8 -*-
"""
test_coloc_pipeline.py — A colocalização de STX6 de ponta a ponta, recalculada.

Não depende do arquivo do GWAS (188 MB, não versionado): as 390 variantes do
GWAS que casam com a meta-eQTL estão em `pipeline/data/stx6_crosscheck_input.tsv`
e servem de mini-GWAS. Os eQTL vêm dos caches versionados. Cada teste roda o
código do pipeline (harmonização → meta → colocalização) e confere:

- as âncoras da seção 2.2 do dossiê de auditoria (sensibilidades do H4);
- a proveniência (checksums, tamanhos amostrais oficiais);
- que os relatórios versionados citam os números que o código produz agora.
"""
import csv
import re
import sys
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
import coloc_abf as C  # noqa: E402
import fontes_qtl as F  # noqa: E402
import harmoniza as H  # noqa: E402
import coloc_meta_stx6 as M  # noqa: E402
import coloc_sqtl_stx6 as SQ  # noqa: E402

ENTRADA = BASE / "pipeline" / "data" / "stx6_crosscheck_input.tsv"
REPORTS = BASE / "pipeline" / "reports"


def mini_gwas():
    """Índice do GWAS restrito às variantes versionadas (alelo de efeito = ALT)."""
    with ENTRADA.open(encoding="utf-8") as fh:
        linhas = list(csv.DictReader(fh, delimiter="\t"))
    return H.indexar_gwas(
        {"pos": int(r["pos_b37"]), "outro": r["ref"], "efeito": r["alt"],
         "beta": float(r["beta_gwas"]), "se": float(r["se_gwas"]),
         "eaf": float(r["maf"]), "p": float(r["p_gwas"])} for r in linhas), linhas


@pytest.fixture(scope="module")
def analise():
    gwas, linhas = mini_gwas()
    return M.analisar(gwas), linhas


# -------------------------------------------------------------- meta-eQTL
def test_entrada_versionada_e_a_meta_recalculada(analise):
    """A entrada versionada tem de ser exatamente o que caches + código dão."""
    a, linhas = analise
    assert len(a["pares"]) == len(linhas) == 390
    for p, r in zip(a["pares"], linhas):
        assert p["pos"] == int(r["pos_b37"])
        assert (p["ref"], p["alt"]) == (r["ref"], r["alt"])
        assert p["qtl"]["beta"] == pytest.approx(float(r["beta_eqtl_meta"]),
                                                 rel=1e-12)
        assert p["qtl"]["se"] == pytest.approx(float(r["se_eqtl_meta"]),
                                               rel=1e-12)
        assert p["qtl"]["k"] == int(r["n_datasets"])


def test_h4_principal(analise):
    a, _ = analise
    pp = a["res"]["pp"]
    assert pp[4] == pytest.approx(0.9796, abs=5e-5)
    assert pp[3] == pytest.approx(0.0204, abs=5e-5)
    assert a["res"]["w1"] == pytest.approx(0.04)
    assert a["res"]["w2"] == pytest.approx(0.0225)


def test_ancoras_de_sensibilidade_da_secao_2_2(analise):
    """Efeito aleatório, leave-one-out, sem GTEx, um conjunto por vez."""
    a, _ = analise
    s = {nome: (h4, n) for nome, h4, _h3, n in a["sens"]}
    assert s["Efeito fixo (análise principal)"][0] == pytest.approx(0.97960, abs=5e-5)
    assert s["Efeitos aleatórios (DerSimonian–Laird)"][0] == pytest.approx(
        0.97547, abs=5e-5)
    loo = [v for k, v in s.items() if k.startswith("Sem ") and "4 conjuntos" in k]
    assert len(loo) == 5
    assert min(h for h, _ in loo) == pytest.approx(0.97940, abs=5e-5)
    assert max(h for h, _ in loo) == pytest.approx(0.97959, abs=5e-5)
    assert sorted(n for _, n in loo) == [365, 373, 374, 380, 380]
    assert s["Sem os dois GTEx (3 conjuntos, mínimo 3)"] == (
        pytest.approx(0.97947, abs=5e-5), 363)
    so = {k: v for k, v in s.items() if k.startswith("Só ")}
    esperado = {"CommonMind": (0.9782, 389), "ROSMAP": (0.9790, 365),
                "BrainSeq": (0.9792, 374), "GTEx_v10 DLPFC": (0.7015, 371),
                "GTEx_v10 cerebelo": (0.5672, 380)}
    for chave, (h4, n) in esperado.items():
        achado = [v for k, v in so.items() if chave in k]
        assert len(achado) == 1
        assert achado[0][0] == pytest.approx(h4, abs=5e-5)
        assert achado[0][1] == n


def test_heterogeneidade_e_direcao(analise):
    a, _ = analise
    pares = a["pares"]
    # 36 e 51 (o dossiê, com a lógica antiga por posição, contava 37 e 52: a
    # diferença é o sítio 180.975.614, cuja heterogeneidade vinha da mistura
    # de alelos)
    assert sum(1 for p in pares if p["qtl"]["p_Q"] < 0.05) == 36
    assert sum(1 for p in pares if p["qtl"]["I2"] > 0.5) == 51
    conc = sum(1 for (bg, _), (be, _) in zip(a["tg"], a["te"]) if bg * be > 0)
    assert conc == 348
    assert a["st_gwas"]["so_na_fita_complementar"] == 0
    assert a["st_gwas"]["alelos_incompativeis"] == 0


def test_sdy_estimado_e_compativel_com_um(analise):
    """A análise supõe sdY = 1 para o eQTL; a estimativa por conjunto fica
    entre 0,6 e 1,05 (EP reflete variância residual após covariáveis)."""
    a, _ = analise
    for reg, _st, sdy in a["resumo"]:
        assert 0.6 < sdy < 1.05, (reg["rotulo"], sdy)


def test_segundo_oraculo_coloc_r_sobre_a_entrada_atual(analise):
    """`pipeline/data/stx6_crosscheck_colocR_output.tsv` é a saída do
    coloc.abf real (R 4.5.3, coloc 5.2.3, 2026-11-11) sobre a entrada ATUAL.
    A implementação do projeto tem de coincidir variante a variante."""
    import crosscheck_coloc_R as X
    a, _ = analise
    r = X.ler_saida_r(X.SAIDA_R)
    assert set(r) == {p["pos"] for p in a["pares"]}
    dif = max(abs(x - r[p["pos"]])
              for x, p in zip(a["res"]["snp_pp_h4"], a["pares"]))
    assert dif < 1e-12
    cab = X.cabecalho_r(X.SAIDA_R)
    m = re.search(r"PP\.H3=([\d.]+) PP\.H4=([\d.]+)", cab)
    assert m, "cabeçalho do script R ausente"
    assert float(m.group(1)) == pytest.approx(a["res"]["pp"][3], abs=1e-6)
    assert float(m.group(2)) == pytest.approx(a["res"]["pp"][4], abs=1e-6)


# ----------------------------------------------------------------- sQTL
def test_sqtl_recalculado():
    gwas, _ = mini_gwas()
    a = SQ.analisar(gwas)
    assert sum(ns for _, _, ns in a["resumo"]) == 0          # nenhum de STX6
    assert len(a["testes"]) == 6
    h4 = sorted((t["pp"][4] for t in a["testes"]), reverse=True)
    assert h4[0] == pytest.approx(0.598, abs=2e-3)
    assert all(x < 0.8 for x in h4)
    assert sum(1 for x in h4 if x > 0.5) == 1
    txt = (BASE / "pipeline" / "scripts" / "coloc_sqtl_stx6.py").read_text(
        encoding="utf-8")
    assert not re.search(r"0\.5\s*/\s*(max\()?n_t", txt), \
        "limiar 'Bonferroni' sobre PP.H4 voltou"


def test_chave_de_cluster_e_comparavel_entre_conjuntos():
    assert SQ.chave_cluster("1:180936513:180938555:clu_46838_+") == \
        "180936513:180938555:+"
    assert SQ.chave_cluster("1:180936513:180938555:clu_51_-") == \
        "180936513:180938555:-"


# ------------------------------------------------------------ proveniência
def test_proveniencia_confere_com_os_caches():
    meta, regs = F.carregar_proveniencia()
    assert meta["release"] == "r8"
    assert meta["metadata_commit"] == "7a5b65b"
    assert len(regs) == 10
    for r in regs:
        cache = F.DATA / F.nome_cache(r)
        if r["linhas"] == 0:
            assert not cache.exists()
            continue
        assert F.sha256_arquivo(cache) == r["sha256"], r["dataset_id"]
        assert len(cache.read_text().splitlines()) == r["linhas"]


def test_tamanhos_amostrais_sao_os_da_metadata_r8():
    n = {r["dataset_id"]: r["n_amostras"] for r in F.datasets("ge")}
    assert n == {"QTD000075": 586, "QTD000434": 560, "QTD000051": 479,
                 "QTD000176": 268, "QTD000166": 263}
    assert sum(n.values()) == 2156
    ns = {r["dataset_id"]: r["n_amostras"] for r in F.datasets("leafcutter")}
    assert sorted(ns.values()) == sorted(n.values())


def test_scripts_nao_digitam_tamanhos_amostrais():
    """Os n antigos do GTEx (285, 272) e o total (~2.182) não podem voltar."""
    for nome in ("coloc_stx6_eqtl.py", "coloc_meta_stx6.py",
                 "coloc_sqtl_stx6.py", "crosscheck_coloc_R.py"):
        txt = (BASE / "pipeline" / "scripts" / nome).read_text(encoding="utf-8")
        assert not re.search(r"\b(285|272)\b", txt), nome
        assert "2.182" not in txt and "2182" not in txt, nome
        assert "d[pos38] =" not in txt, f"{nome}: indexação só por posição"


# ------------------------------------- relatórios citam o que o código produz
def test_relatorio_meta_cita_os_numeros_recalculados(analise):
    a, _ = analise
    pp = a["res"]["pp"]
    txt = (REPORTS / "relatorio_coloc_meta_stx6.md").read_text(encoding="utf-8")
    assert f"**PP.H4 = {pp[4]:.4f}**" in txt
    assert f"PP.H3 = {pp[3]:.4f}" in txt
    for p12, ppx in C.sensibilidade_p12(a["tg"], a["te"], "cc", "quant"):
        assert f"| {C._fmt(ppx[3])} | {C._fmt(ppx[4])} |" in txt, p12
    for nome, h4, h3, n in a["sens"]:
        assert f"| {nome} | {h3:.4f} | {h4:.4f} | {n} |" in txt, nome
    assert "2.156 amostras" in txt
    assert "348/390" in txt
    # descrições corretas das hipóteses (G5)
    assert "associação só com o eQTL de STX6" in txt
    assert "H2 NÃO é hipótese nula" not in txt
    assert "H2+H3+H4" not in txt
    assert "Colocalização não é mediação" in txt
    assert "conservador\"" in txt and "não** é mais rigoroso" in txt


def test_relatorio_validacao_registra_o_que_aconteceu():
    txt = (REPORTS / "relatorio_validacao_coloc_R.md").read_text(encoding="utf-8")
    assert "não era uma validação" in txt
    assert "ACORDO A 6 DECIMAIS" not in txt
    assert "A máquina (ABF + combinação H0..H4) está correta" not in txt
    assert "0.979605" in txt


def test_relatorio_por_tecido_cita_os_numeros_recalculados(analise):
    """GTEx sozinho: H4 = 0,70 (DLPFC) e 0,57 (cerebelo) — os valores do
    relatório por tecido têm de ser os que o código dá para cada tecido."""
    gwas, _ = mini_gwas()
    txt = (REPORTS / "relatorio_coloc_stx6.md").read_text(encoding="utf-8")
    regs = {r["dataset_id"]: r for r in F.datasets("ge")}
    for qtd in ("QTD000176", "QTD000166"):
        variantes, _ = H.deduplicar(F.ler_eqtl(regs[qtd]))
        pares, _ = H.casar_com_gwas(variantes, gwas, F.OFFSET_B38_B37)
        tg = [(p["beta_gwas"], p["se_gwas"]) for p in pares]
        te = [(p["qtl"]["beta"], p["qtl"]["se"]) for p in pares]
        pp = C.coloc(tg, te, "cc", "quant")["pp"]
        # o mini-GWAS só tem as 390 variantes da meta; o relatório usa o GWAS
        # inteiro (372 e 383 variantes). A diferença de 1 a 3 variantes sem
        # sinal não muda o H4 na 3ª casa.
        m = re.search(rf"## \w+ — {re.escape(regs[qtd]['rotulo'])}.*?"
                      r"\*\*PP\.H4 = ([\d.]+)\*\*", txt, re.S)
        assert m, qtd
        assert float(m.group(1)) == pytest.approx(pp[4], abs=2e-3)
