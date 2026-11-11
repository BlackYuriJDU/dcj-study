# -*- coding: utf-8 -*-
"""
ols.py — Mínimos quadrados por fatoração QR (Householder), com detecção
explícita de testes NÃO ESTIMÁVEIS. Só biblioteca padrão.

Por que existe (auditoria 2026-11, achado M1): o solver anterior usava
equações normais com Gauss–Jordan e não detectava resposta constante. Para os
74 miRNAs sem nenhuma variação, o "p" que saía era ruído de arredondamento
(0,0203 numa execução, 5,5×10⁻⁷ em outra com os mesmos números): bastava trocar
o tipo numérico da entrada para criar dezenas de falsos positivos.

Aqui um teste é declarado não estimável, sem p-valor, quando:
  - a resposta é constante (não há variância a explicar);
  - a matriz de desenho tem posto deficiente (colinearidade exata);
  - não sobram graus de liberdade residuais;
  - o ajuste é perfeito (variância residual nula → estatística t indefinida).

Uso:
    d = Desenho(X)                 # fatora uma vez; X = lista de linhas
    a = d.ajustar(y)               # Ajuste(estimavel, motivo, beta, se, gl, ...)
    if a.estimavel: t = a.beta[1] / a.se[1]
"""
import math
from collections import namedtuple

__all__ = ["Desenho", "Ajuste", "ajustar", "MOTIVOS"]

Ajuste = namedtuple("Ajuste", "estimavel motivo beta se gl s2")

MOTIVOS = {
    "resposta_constante": "resposta constante em todas as amostras",
    "posto_deficiente": "matriz de desenho com posto deficiente",
    "sem_gl": "sem graus de liberdade residuais",
    "ajuste_perfeito": "variância residual nula (ajuste perfeito)",
    "nao_finito": "valor não finito na resposta",
}

_TOL_POSTO = 1e-10      # relativo, sobre colunas normalizadas
_TOL_RESID = 1e-20      # variância residual relativa à variância da resposta


class Desenho:
    """Fatoração QR de uma matriz de desenho n×p (colunas normalizadas para
    que a tolerância de posto não dependa das unidades das covariáveis)."""

    def __init__(self, X):
        self.n = len(X)
        self.p = len(X[0]) if X else 0
        if any(len(l) != self.p for l in X):
            raise ValueError("linhas de X com tamanhos diferentes")
        self.X = [[float(v) for v in l] for l in X]
        n, p = self.n, self.p
        # escala das colunas
        self.escala = []
        for j in range(p):
            nrm = math.sqrt(math.fsum(self.X[i][j] ** 2 for i in range(n)))
            self.escala.append(nrm if nrm > 0 else 1.0)
        A = [[self.X[i][j] / self.escala[j] for j in range(p)] for i in range(n)]
        self._v = []                      # vetores de Householder
        diag = []
        for j in range(min(n, p)):
            x = [A[i][j] for i in range(j, n)]
            alfa = math.sqrt(math.fsum(v * v for v in x))
            if alfa == 0.0:
                self._v.append(None)
                diag.append(0.0)
                continue
            if x[0] > 0:
                alfa = -alfa
            v = x[:]
            v[0] -= alfa
            nv = math.sqrt(math.fsum(t * t for t in v))
            if nv == 0.0:
                self._v.append(None)
                diag.append(abs(alfa))
                continue
            v = [t / nv for t in v]
            for c in range(j, p):
                s = 2.0 * math.fsum(v[i - j] * A[i][c] for i in range(j, n))
                for i in range(j, n):
                    A[i][c] -= s * v[i - j]
            self._v.append(v)
            diag.append(abs(A[j][j]))
        self.R = [[A[i][j] if j >= i else 0.0 for j in range(p)]
                  for i in range(min(n, p))]
        maior = max(diag) if diag else 0.0
        self.posto = sum(1 for d in diag if maior > 0 and d > _TOL_POSTO * maior)
        self.posto_completo = (n >= p and self.posto == p)
        self.gl = n - p
        # diagonal de (X'X)^-1 = escala^-1 · R^-1 R^-T · escala^-1
        self._diag_inv = None
        if self.posto_completo:
            Rinv = self._inversa_triangular(self.R)
            self._diag_inv = [
                math.fsum(Rinv[j][k] ** 2 for k in range(j, p))
                / (self.escala[j] ** 2) for j in range(p)]

    @staticmethod
    def _inversa_triangular(R):
        p = len(R)
        inv = [[0.0] * p for _ in range(p)]
        for j in range(p):
            inv[j][j] = 1.0 / R[j][j]
            for i in range(j - 1, -1, -1):
                s = math.fsum(R[i][k] * inv[k][j] for k in range(i + 1, j + 1))
                inv[i][j] = -s / R[i][i]
        return inv

    def _qty(self, y):
        """Primeiras p componentes de Q'y."""
        z = list(y)
        for j, v in enumerate(self._v):
            if v is None:
                continue
            s = 2.0 * math.fsum(v[i - j] * z[i] for i in range(j, self.n))
            for i in range(j, self.n):
                z[i] -= s * v[i - j]
        return z

    def ajustar(self, y):
        """Ajusta y ~ X. Devolve Ajuste; se não estimável, beta/se = None."""
        if len(y) != self.n:
            raise ValueError(f"y tem {len(y)} valores; X tem {self.n} linhas")
        y = [float(v) for v in y]
        nao = lambda m: Ajuste(False, m, None, None, self.gl, None)  # noqa: E731
        if not all(math.isfinite(v) for v in y):
            return nao("nao_finito")
        if max(y) == min(y):
            return nao("resposta_constante")
        if not self.posto_completo:
            return nao("posto_deficiente")
        if self.gl <= 0:
            return nao("sem_gl")
        p = self.p
        z = self._qty(y)
        b = [0.0] * p                       # resolve R b = z[:p]
        for i in range(p - 1, -1, -1):
            s = math.fsum(self.R[i][k] * b[k] for k in range(i + 1, p))
            b[i] = (z[i] - s) / self.R[i][i]
        beta = [b[j] / self.escala[j] for j in range(p)]
        resid = [y[i] - math.fsum(self.X[i][j] * beta[j] for j in range(p))
                 for i in range(self.n)]
        sq = math.fsum(r * r for r in resid)
        media = math.fsum(y) / self.n
        sq_total = math.fsum((v - media) ** 2 for v in y)
        if sq <= _TOL_RESID * sq_total:
            return nao("ajuste_perfeito")
        s2 = sq / self.gl
        se = [math.sqrt(s2 * d) for d in self._diag_inv]
        return Ajuste(True, "", beta, se, self.gl, s2)


def ajustar(X, y):
    """Atalho: fatora X e ajusta um único y."""
    return Desenho(X).ajustar(y)
