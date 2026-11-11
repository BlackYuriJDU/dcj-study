# -*- coding: utf-8 -*-
"""
test_integracao.py — Integração sangue × cérebro, recalculada.

Achados I1–I5 da auditoria de 2026-11. O teste hipergeométrico é conferido
contra scipy; as âncoras da seção 2.6 do dossiê (valores "repositório" e
"correto") são reproduzidas mudando os parâmetros da análise, e o relatório
versionado tem de citar o que o código produz na configuração corrigida.
"""
import sys
from pathlib import Path

import pytest
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
import analise_integrada_cerebro_sangue as I  # noqa: E402

REL = BASE / "pipeline" / "reports" / "relatorio_integracao_cerebro_sangue.md"
M16, M93, LET7 = "hsa-miR-16-5p", "hsa-miR-93-5p", "hsa-let-7i-5p"
SEM_ALVO = ("hsa-miR-106b-3p", "hsa-miR-500a-3p")


def por_mirna(a):
    return {r["mirna"]: r for r in a["linhas"]}


@pytest.fixture(scope="module")
def corrigida():
    return I.analisar()


@pytest.mark.parametrize("k,K,n,N", [(6, 13, 314, 800), (3, 6, 315, 742),
                                     (0, 5, 100, 800), (5, 5, 100, 800),
                                     (1, 3, 315, 742), (40, 48, 314, 800),
                                     (2, 700, 300, 742)])
def test_hipergeometrico_bate_com_scipy(k, K, n, N):
    assert I.hipergeometrico(k, K, n, N) == pytest.approx(
        float(hypergeom.sf(k - 1, N, K, n)), rel=1e-9, abs=1e-300)


def test_hipergeometrico_rejeita_k_maior_que_o_universo():
    """O cálculo antigo passava K = 48 alvos do genoma com N = 800 — coerente
    só por acaso. K > N tem de ser erro, não um número."""
    with pytest.raises(ValueError):
        I.hipergeometrico(1, 900, 314, 800)


def test_ancoras_do_calculo_antigo_secao_2_6():
    """Reproduz o que o repositório calculava: K = todos os alvos do genoma,
    sem filtro de espécie do alvo, N = 800, n = 314."""
    a = I.analisar(universo="painel", exigir_especie_do_alvo=False,
                   restringir_ao_painel=False)
    assert (a["N"], a["n_up"]) == (800, 314)
    r = por_mirna(a)
    assert [r[m]["alvos_total"] for m in (M16, M93, LET7)] == [48, 16, 11]
    assert [r[m]["alvos_nao_humanos"] for m in (M16, M93, LET7)] == [1, 1, 0]
    assert [r[m]["k"] for m in (M16, M93, LET7)] == [6, 3, 1]
    assert [round(r[m]["esperado"], 1) for m in (M16, M93, LET7)] == [18.8, 6.3, 4.3]
    assert r[M16]["p"] == pytest.approx(1.00, abs=5e-3)
    assert r[M93]["p"] == pytest.approx(0.98, abs=5e-3)
    assert r[LET7]["p"] == pytest.approx(0.996, abs=5e-4)


def test_ancoras_com_k_restrito_ao_painel_secao_2_6():
    """Valores 'correto' do dossiê (ainda com N = 800 e sem filtro de espécie)."""
    a = I.analisar(universo="painel", exigir_especie_do_alvo=False)
    r = por_mirna(a)
    assert [r[m]["alvos_no_painel"] for m in (M16, M93, LET7)] == [13, 6, 3]
    assert [round(r[m]["esperado"], 1) for m in (M16, M93, LET7)] == [5.1, 2.4, 1.2]
    assert r[M16]["p"] == pytest.approx(0.40, abs=5e-3)
    assert r[M93]["p"] == pytest.approx(0.44, abs=5e-3)
    assert r[LET7]["p"] == pytest.approx(0.78, abs=5e-3)


def test_analise_corrigida(corrigida):
    a = corrigida
    assert (a["N"], a["n_up"], a["n_down"]) == (742, 315, 126)
    r = por_mirna(a)
    # I2: espécie do alvo — os dois alvos não humanos saem
    assert [r[m]["alvos_total"] for m in (M16, M93, LET7)] == [47, 15, 11]
    assert all(r[m]["alvos_nao_humanos"] == 0 for m in r)
    # I1: K = alvos no painel testável
    assert [r[m]["K"] for m in (M16, M93, LET7)] == [13, 6, 3]
    assert [r[m]["k"] for m in (M16, M93, LET7)] == [6, 3, 1]
    for m in (M16, M93, LET7):
        assert r[m]["p"] == pytest.approx(
            float(hypergeom.sf(r[m]["k"] - 1, 742, r[m]["K"], 315)), rel=1e-9)
        assert r[m]["esperado"] == pytest.approx(r[m]["K"] * 315 / 742)
        assert r[m]["p"] > 0.4                       # resultado nulo se mantém
    # I3: não testáveis fora do FDR
    assert a["n_testes"] == 3
    for m in SEM_ALVO:
        assert not r[m]["testavel"] and r[m]["p"] is None and "q" not in r[m]
    qs = multipletests([r[m]["p"] for m in (M16, M93, LET7)], method="fdr_bh")[1]
    for m, q in zip((M16, M93, LET7), qs):
        assert r[m]["q"] == pytest.approx(float(q), rel=1e-12)


def test_k_minimo_para_significancia(corrigida):
    r = por_mirna(corrigida)
    for m in (M16, M93, LET7):
        kmin = r[m]["k_min"]
        if kmin is None:
            assert I.hipergeometrico(r[m]["K"], r[m]["K"], 315, 742) >= 0.05
        else:
            assert I.hipergeometrico(kmin, r[m]["K"], 315, 742) < 0.05
            assert I.hipergeometrico(kmin - 1, r[m]["K"], 315, 742) >= 0.05


def test_leitura_do_csv_com_o_modulo_csv():
    """I4: campos com vírgula entre aspas não podem deslocar colunas; o
    resultado coincide com as contagens conhecidas (48, 16, 11 sem filtro de
    espécie do alvo)."""
    txt = (BASE / "pipeline" / "scripts" /
           "analise_integrada_cerebro_sangue.py").read_text(encoding="utf-8")
    assert "csv.DictReader" in txt
    assert '.split(",")' not in txt
    alvos = I.alvos_mirtarbase(exigir_especie_do_alvo=False)
    assert [len(alvos[m]) for m in (M16, M93, LET7)] == [48, 16, 11]
    assert {e for d in alvos.values() for e in d.values()} == {"hsa", "rno", "mmu"}
    humanos = I.alvos_mirtarbase()
    assert {e for d in humanos.values() for e in d.values()} == {"hsa"}


def test_relatorio_cita_o_recalculado_e_nao_fala_em_mecanismo(corrigida):
    txt = REL.read_text(encoding="utf-8")
    r = por_mirna(corrigida)
    assert "**N = 742**" in txt and "**n = 315**" in txt
    for m in (M16, M93, LET7):
        x = r[m]
        assert (f"| {m} | {x['alvos_total']} | {x['K']} | {x['k']} | "
                f"{x['esperado']:.1f} | {x['p']:.2f} | {x['q']:.2f} |") in txt
    for m in SEM_ALVO:
        assert f"| {m} | 0 | 0 | — | — | não testável |" in txt
    assert "Testes realizados: **3**" in txt
    # I5: sem "controle negativo" sem teste e sem leitura de mecanismo
    assert "Controle negativo" not in txt
    assert "biomarcador mecanístico" not in txt
    assert "desmonta inferência causal" not in txt
    assert "pouco poder" in txt
