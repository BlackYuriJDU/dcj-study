# -*- coding: utf-8 -*-
"""
test_ols.py — Regressão por QR com detecção de testes não estimáveis.

Achado M1 da auditoria de 2026-11: o solver antigo (equações normais +
Gauss–Jordan) devolvia, para resposta constante, um "p" que era ruído de
arredondamento e mudava com o tipo numérico da entrada. Os testes de borda
abaixo falham com aquele solver; os de regressão conferem contra numpy e
statsmodels.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import statsmodels.api as sm

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
from ols import Desenho, MOTIVOS, ajustar  # noqa: E402


def desenho_aleatorio(n=60, semente=1):
    rng = np.random.default_rng(semente)
    grupo = (rng.random(n) < 0.5).astype(float)
    sexo = (rng.random(n) < 0.4).astype(float)
    idade = rng.normal(60, 12, n) + 8 * grupo
    rin = rng.normal(6, 1.2, n) - 0.8 * grupo
    X = np.column_stack([np.ones(n), grupo, sexo, idade, rin])
    return X, rng


# ----------------------------------------------------------- casos regulares
def test_bate_com_numpy_e_statsmodels():
    X, rng = desenho_aleatorio()
    d = Desenho(X.tolist())
    assert d.posto == 5 and d.posto_completo and d.gl == 55
    for _ in range(25):
        y = X @ rng.normal(0, 1, 5) + rng.normal(0, 2.0, len(X))
        a = d.ajustar(y.tolist())
        assert a.estimavel and a.motivo == ""
        beta_np = np.linalg.lstsq(X, y, rcond=None)[0]
        m = sm.OLS(y, X).fit()
        assert a.beta == pytest.approx(beta_np.tolist(), rel=1e-9, abs=1e-11)
        assert a.se == pytest.approx(m.bse.tolist(), rel=1e-9)
        assert a.gl == int(m.df_resid)
        assert a.s2 == pytest.approx(float(m.mse_resid), rel=1e-9)


def test_covariaveis_em_escalas_muito_diferentes():
    """Colunas com escalas distantes (idade em dias × indicador 0/1): QR com
    colunas normalizadas não perde o posto nem a precisão."""
    X, rng = desenho_aleatorio(80, 3)
    X[:, 3] *= 365.25 * 1e3
    y = rng.normal(0, 1, len(X))
    a = ajustar(X.tolist(), y.tolist())
    m = sm.OLS(y, X).fit()
    assert a.estimavel
    assert a.se == pytest.approx(m.bse.tolist(), rel=1e-7)
    assert [b / s for b, s in zip(a.beta, a.se)] == pytest.approx(
        m.tvalues.tolist(), rel=1e-7, abs=1e-9)


def test_resposta_quase_constante_ainda_e_estimavel():
    """Um único valor diferente: há variância, o teste existe e bate com a
    referência (caso típico de miRNA detectado em uma amostra)."""
    X, _ = desenho_aleatorio(40, 5)
    y = np.zeros(len(X))
    y[7] = 17.3
    a = ajustar(X.tolist(), y.tolist())
    m = sm.OLS(y, X).fit()
    assert a.estimavel
    assert a.beta == pytest.approx(m.params.tolist(), rel=1e-8, abs=1e-10)
    assert a.se == pytest.approx(m.bse.tolist(), rel=1e-8)


# --------------------------------------------------------- não estimáveis
@pytest.mark.parametrize("valor", [0.0, 1.4426950408889634e-4, 5.0, -3.2e7])
def test_resposta_constante_nao_e_estimavel(valor):
    X, _ = desenho_aleatorio()
    a = ajustar(X.tolist(), [valor] * len(X))
    assert not a.estimavel
    assert a.motivo == "resposta_constante"
    assert a.beta is None and a.se is None


def test_resposta_constante_independe_do_tipo_numerico():
    """O defeito original: os mesmos números como float do Python e como
    numpy.float64 davam p = 0,0203 e p ~ 1e-6. Agora os dois são 'não
    estimável', sem p."""
    X, _ = desenho_aleatorio(105, 11)
    v = math.log2(0.0001 + 1.0)
    a1 = ajustar(X.tolist(), [v] * 105)
    a2 = ajustar([[np.float64(x) for x in l] for l in X],
                 [np.float64(v)] * 105)
    assert (a1.estimavel, a1.motivo) == (a2.estimavel, a2.motivo) == \
        (False, "resposta_constante")


def test_posto_deficiente_coluna_duplicada():
    X, rng = desenho_aleatorio()
    Xd = np.column_stack([X, X[:, 3] * 2.0])            # colinearidade exata
    d = Desenho(Xd.tolist())
    assert d.posto == 5 and not d.posto_completo
    a = d.ajustar(rng.normal(0, 1, len(X)).tolist())
    assert not a.estimavel and a.motivo == "posto_deficiente"


def test_posto_deficiente_grupo_sem_variacao():
    """Todos do mesmo grupo: a coluna do grupo é múltipla do intercepto."""
    n = 30
    X = [[1.0, 1.0, float(i)] for i in range(n)]
    a = ajustar(X, [float(i % 7) for i in range(n)])
    assert not a.estimavel and a.motivo == "posto_deficiente"


def test_sem_graus_de_liberdade():
    X = [[1.0, 0.0], [1.0, 1.0]]
    a = ajustar(X, [3.0, 5.0])
    assert not a.estimavel and a.motivo == "sem_gl"


def test_ajuste_perfeito_nao_gera_estatistica_infinita():
    X, _ = desenho_aleatorio(30, 2)
    y = (X @ np.array([1.0, 2.0, -1.0, 0.5, 3.0])).tolist()
    a = ajustar(X.tolist(), y)
    assert not a.estimavel and a.motivo == "ajuste_perfeito"


def test_valor_nao_finito():
    X, rng = desenho_aleatorio(20, 4)
    y = rng.normal(0, 1, 20).tolist()
    y[3] = float("nan")
    a = ajustar(X.tolist(), y)
    assert not a.estimavel and a.motivo == "nao_finito"


def test_motivos_documentados():
    for m in ("resposta_constante", "posto_deficiente", "sem_gl",
              "ajuste_perfeito", "nao_finito"):
        assert m in MOTIVOS


def test_tamanhos_incompativeis_sao_erro():
    X, _ = desenho_aleatorio(10, 1)
    with pytest.raises(ValueError):
        Desenho(X.tolist()).ajustar([0.0] * 9)
    with pytest.raises(ValueError):
        Desenho([[1.0, 2.0], [1.0]])
