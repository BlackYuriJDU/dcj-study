# -*- coding: utf-8 -*-
"""
test_harmoniza.py — Harmonização por variante completa e meta-análise.

Cobre o achado G6 da auditoria de 2026-11: a indexação só por posição deixava
a última linha vencer e somava efeitos de alelos diferentes em sítios
multialélicos. Os testes sintéticos falham com aquela lógica; os testes sobre
os caches versionados conferem as contagens reais.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import chi2

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
import harmoniza as H  # noqa: E402
import fontes_qtl as F  # noqa: E402


def linha(pos, ref, alt, beta, se, rsid="rs1", p=0.5, maf=0.2):
    return (pos, ref, alt, rsid, p, beta, se, maf)


# ------------------------------------------------------------ normalização
@pytest.mark.parametrize("ent,esp", [
    ((100, "A", "G"), (100, "A", "G")),
    ((100, "a", "g"), (100, "A", "G")),
    ((100, "CTT", "CT"), (100, "CT", "C")),          # deleção: apara sufixo
    ((100, "C", "CT"), (100, "C", "CT")),            # inserção já mínima
    ((100, "AT", "AC"), (101, "T", "C")),            # SNV escrito como MNV
    ((100, "GCA", "GTA"), (101, "C", "T")),
    ((100, "CAA", "CAAA"), (100, "C", "CA")),
])
def test_normalizar_alelos(ent, esp):
    assert H.normalizar_alelos(*ent) == esp


def test_complementar_e_palindromico():
    assert H.complementar("A") == "T"
    assert H.complementar("AC") == "GT"
    assert H.complementar("N") is None
    assert H.palindromico("A", "T") and H.palindromico("G", "C")
    assert not H.palindromico("A", "G")
    assert not H.palindromico("AT", "A")


# ------------------------------------------------------------ deduplicação
def test_deduplicar_linhas_repetidas_por_rsid():
    v, st = H.deduplicar([linha(10, "A", "G", 0.1, 0.02, "rs1"),
                          linha(10, "A", "G", 0.1, 0.02, "rs2"),
                          linha(11, "C", "T", -0.3, 0.05, "rs3")])
    assert set(v) == {(10, "A", "G"), (11, "C", "T")}
    assert v[(10, "A", "G")]["rsids"] == ["rs1", "rs2"]
    assert st["linhas"] == 3 and st["variantes"] == 2
    assert st["linhas_repetidas_por_rsid"] == 1
    assert st["descartadas_conflito"] == 0


def test_deduplicar_descarta_mesma_variante_com_efeitos_diferentes():
    v, st = H.deduplicar([linha(10, "A", "G", 0.1, 0.02),
                          linha(10, "A", "G", 0.4, 0.02)])
    assert v == {}
    assert st["descartadas_conflito"] == 1


def test_sitio_multialelico_nao_e_sobrescrito():
    """Na lógica antiga (`d[pos] = ...`) só o último alelo sobrevivia."""
    v, st = H.deduplicar([linha(10, "G", "T", 0.1, 0.02),
                          linha(10, "G", "C", -0.5, 0.03)])
    assert set(v) == {(10, "G", "T"), (10, "G", "C")}
    assert v[(10, "G", "T")]["beta"] == 0.1
    assert v[(10, "G", "C")]["beta"] == -0.5
    assert st["posicoes"] == 1 and st["posicoes_multialelicas"] == 1


def test_deduplicar_rejeita_erro_padrao_invalido():
    v, st = H.deduplicar([linha(10, "A", "G", 0.1, 0.0),
                          linha(11, "A", "G", float("nan"), 0.1)])
    assert v == {} and st["descartadas_se_invalido"] == 2


# ------------------------------------------------------------- meta-análise
def _ds(*linhas):
    return H.deduplicar(list(linhas))[0]


def test_meta_so_combina_alelos_identicos():
    """Caso do sítio chr1:180.975.614: G>T em dois conjuntos, G>C em três.
    A lógica antiga somava os cinco. Aqui G>C entra com k = 3 e G>T fica
    fora (k = 2 < 3)."""
    por = {
        "d1": _ds(linha(10, "G", "T", 0.50, 0.1)),
        "d2": _ds(linha(10, "G", "T", 0.52, 0.1)),
        "d3": _ds(linha(10, "G", "C", -0.20, 0.1)),
        "d4": _ds(linha(10, "G", "C", -0.22, 0.1)),
        "d5": _ds(linha(10, "G", "C", -0.18, 0.1)),
    }
    meta, st = H.meta_ivw(por)
    assert set(meta) == {(10, "G", "C")}
    assert meta[(10, "G", "C")]["k"] == 3
    assert meta[(10, "G", "C")]["beta"] == pytest.approx(-0.20, abs=1e-12)
    assert st["posicoes_com_alelos_discordantes"] == 1
    assert st["abaixo_do_minimo"] == 1


def test_meta_ivw_bate_com_numpy():
    rng = np.random.default_rng(7)
    betas = rng.normal(0.1, 0.05, 5)
    ses = rng.uniform(0.02, 0.2, 5)
    por = {f"d{i}": _ds(linha(10, "A", "G", float(b), float(s)))
           for i, (b, s) in enumerate(zip(betas, ses))}
    m = H.meta_ivw(por)[0][(10, "A", "G")]
    w = 1 / ses ** 2
    beta = float((betas * w).sum() / w.sum())
    q = float((w * (betas - beta) ** 2).sum())
    assert m["beta"] == pytest.approx(beta, rel=1e-12)
    assert m["se"] == pytest.approx(float(np.sqrt(1 / w.sum())), rel=1e-12)
    assert m["Q"] == pytest.approx(q, rel=1e-10)
    assert m["p_Q"] == pytest.approx(float(chi2.sf(q, 4)), rel=1e-9)
    assert m["I2"] == pytest.approx(max(0.0, (q - 4) / q), rel=1e-10)


def test_meta_efeitos_aleatorios_dersimonian_laird():
    betas = np.array([0.10, 0.50, -0.20, 0.35])
    ses = np.array([0.05, 0.06, 0.05, 0.08])
    por = {f"d{i}": _ds(linha(10, "A", "G", float(b), float(s)))
           for i, (b, s) in enumerate(zip(betas, ses))}
    m = H.meta_ivw(por, aleatorio=True)[0][(10, "A", "G")]
    w = 1 / ses ** 2
    bf = (betas * w).sum() / w.sum()
    q = (w * (betas - bf) ** 2).sum()
    tau2 = max(0.0, (q - 3) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (ses ** 2 + tau2)
    assert tau2 > 0
    assert m["beta"] == pytest.approx(float((betas * wr).sum() / wr.sum()),
                                      rel=1e-12)
    assert m["se"] == pytest.approx(float(np.sqrt(1 / wr.sum())), rel=1e-12)
    fixo = H.meta_ivw(por)[0][(10, "A", "G")]
    assert m["se"] > fixo["se"]


@pytest.mark.parametrize("gl", [1, 2, 3, 4, 5, 8, 9])
def test_chi2_sf_bate_com_scipy(gl):
    for x in (0.0, 0.2, 1.0, 3.84, 9.5, 30.0):
        assert H.chi2_sf(x, gl) == pytest.approx(float(chi2.sf(x, gl)),
                                                 rel=1e-9, abs=1e-15)


# ------------------------------------------------------- casamento com o GWAS
def _g(pos, outro, efeito, beta, se=0.1, eaf=0.3, p=0.01):
    return {"pos": pos, "outro": outro, "efeito": efeito, "beta": beta,
            "se": se, "eaf": eaf, "p": p}


def test_casar_inverte_sinal_quando_o_efeito_do_gwas_e_o_ref():
    var = _ds(linha(100, "A", "G", 0.3, 0.05), linha(200, "C", "T", 0.2, 0.05))
    gw = H.indexar_gwas([_g(90, "A", "G", 0.11), _g(190, "T", "C", 0.22)])
    pares, st = H.casar_com_gwas(var, gw, offset=-10)
    assert [p["pos"] for p in pares] == [90, 190]
    assert pares[0]["beta_gwas"] == pytest.approx(0.11)       # efeito = ALT
    assert pares[1]["beta_gwas"] == pytest.approx(-0.22)      # efeito = REF
    assert pares[1]["freq_alt_gwas"] == pytest.approx(0.7)
    assert st["casadas"] == 2
    assert st["casadas_mesmo_alelo_de_efeito"] == 1
    assert st["casadas_sinal_invertido"] == 1


def test_casar_confere_os_dois_alelos():
    """Mesma posição, alelos diferentes: a lógica por posição casaria."""
    var = _ds(linha(100, "G", "T", 0.3, 0.05))
    gw = H.indexar_gwas([_g(100, "G", "A", 0.5)])
    pares, st = H.casar_com_gwas(var, gw)
    assert pares == [] and st["alelos_incompativeis"] == 1


def test_casar_escolhe_o_alelo_certo_em_sitio_multialelico_do_gwas():
    var = _ds(linha(100, "G", "C", 0.3, 0.05), linha(100, "G", "T", 0.1, 0.05))
    gw = H.indexar_gwas([_g(100, "G", "T", 0.9), _g(100, "G", "C", -0.4)])
    pares, st = H.casar_com_gwas(var, gw)
    d = {(p["ref"], p["alt"]): p["beta_gwas"] for p in pares}
    assert d == {("G", "C"): -0.4, ("G", "T"): 0.9}


def test_casar_nao_usa_fita_complementar_mas_conta():
    var = _ds(linha(100, "A", "G", 0.3, 0.05))
    gw = H.indexar_gwas([_g(100, "T", "C", 0.5)])
    pares, st = H.casar_com_gwas(var, gw)
    assert pares == [] and st["so_na_fita_complementar"] == 1


def test_casar_conta_ausentes_e_palindromicas():
    var = _ds(linha(100, "A", "T", 0.3, 0.05), linha(300, "A", "G", 0.3, 0.05))
    gw = H.indexar_gwas([_g(100, "A", "T", 0.5)])
    pares, st = H.casar_com_gwas(var, gw)
    assert st["casadas"] == 1 and st["palindromicas_casadas"] == 1
    assert st["sem_variante_no_gwas"] == 1


# -------------------------------------------- caches versionados (dados reais)
@pytest.fixture(scope="module")
def eqtl_real():
    por, stats = {}, {}
    for reg in F.datasets("ge"):
        por[reg["dataset_id"]], stats[reg["dataset_id"]] = H.deduplicar(
            F.ler_eqtl(reg))
    return por, stats


def test_caches_contagens_por_conjunto(eqtl_real):
    """Seção 2.2 do dossiê: posições com mais de um par de alelos DENTRO de um
    mesmo conjunto — CommonMind 29, ROSMAP 18, BrainSeq 38, GTEx 0 e 0;
    indels 181, 137, 221, 44, 45."""
    _, st = eqtl_real
    ordem = ["QTD000075", "QTD000434", "QTD000051", "QTD000176", "QTD000166"]
    assert [st[q]["posicoes_multialelicas"] for q in ordem] == [29, 18, 38, 0, 0]
    assert [st[q]["indels"] for q in ordem] == [181, 137, 221, 44, 45]
    assert [st[q]["linhas"] for q in ordem] == [887, 711, 1167, 582, 606]
    assert all(st[q]["descartadas_conflito"] == 0 for q in ordem)
    # linhas repetidas da mesma variante sempre com beta e EP idênticos
    assert [st[q]["linhas_repetidas_por_rsid"] for q in ordem] == \
        [108, 88, 113, 33, 33]


def test_logica_antiga_misturava_18_posicoes(eqtl_real):
    """Reproduz o diagnóstico (não a análise): indexando só pela posição, com
    a última linha vencendo, 651 posições entram na meta e 18 delas combinam
    conjuntos com REF/ALT diferentes."""
    por_pos = {}
    for reg in F.datasets("ge"):
        d = {}
        for pos, ref, alt, _rs, _p, beta, se, _maf in F.ler_eqtl(reg):
            d[pos] = (ref, alt, beta, se)
        por_pos[reg["dataset_id"]] = d
    todas = set().union(*(set(d) for d in por_pos.values()))
    na_meta, misturadas = 0, 0
    for pos in todas:
        vs = [d[pos] for d in por_pos.values()
              if pos in d and d[pos][3] > 0 and d[pos][2] != 0]
        if len(vs) < 3:
            continue
        na_meta += 1
        misturadas += len({(v[0], v[1]) for v in vs}) > 1
    assert (na_meta, misturadas) == (651, 18)


def test_meta_real_por_variante_completa(eqtl_real):
    por, _ = eqtl_real
    meta, st = H.meta_ivw(por)
    assert st["na_meta"] == 674
    assert st["posicoes_na_meta"] == 650
    # o sítio que a lógica antiga misturava (b37 180.975.614): agora G>C, k = 5
    pos38 = 180_975_614 - F.OFFSET_B38_B37
    chaves = [k for k in meta if k[0] == pos38]
    assert chaves == [(pos38, "G", "C")]
    assert meta[chaves[0]]["k"] == 5
    assert meta[chaves[0]]["beta"] == pytest.approx(-0.1806822128302013,
                                                    rel=1e-12)
    # lead do GWAS (b37 180.961.245): k = 5, I² = 48,4%, z = −14,37
    lead = meta[(180_961_245 - F.OFFSET_B38_B37, "G", "A")]
    assert lead["k"] == 5
    assert lead["I2"] == pytest.approx(0.484, abs=5e-4)
    assert lead["beta"] / lead["se"] == pytest.approx(-14.37, abs=5e-3)
    assert H.p_de_z(lead["beta"] / lead["se"]) == pytest.approx(7.60e-47,
                                                                rel=5e-3)


def test_indel_normalizado_nao_perde_variantes(eqtl_real):
    por, st = eqtl_real
    for q, d in por.items():
        for pos, ref, alt in d:
            assert ref == ref.upper() and alt == alt.upper()
            assert not (len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1])
        assert st[q]["variantes"] == len(d)
    assert not any(math.isnan(v["beta"]) for d in por.values()
                   for v in d.values())
