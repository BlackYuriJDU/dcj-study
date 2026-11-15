# -*- coding: utf-8 -*-
"""
coloc_abf.py — Implementação ÚNICA de colocalização por fatores de Bayes
aproximados (Giambartolomei et al. 2014), conforme o pacote R `coloc`
(funções `approx.bf.estimates` e `combine.abf`). Só biblioteca padrão.

Todos os scripts do projeto que calculam colocalização ou fator de Bayes de
Wakefield importam daqui. Não há cópias locais (gate em tests/test_coloc_abf.py).

Modelo
------
No máximo uma variante causal por traço na região. Cinco hipóteses:

    H0  nenhuma associação
    H1  associação só com o traço 1
    H2  associação só com o traço 2
    H3  associação com os dois traços, variantes causais distintas
    H4  associação com os dois traços, uma variante causal compartilhada

Priors por variante: p1 (causal só para o traço 1), p2 (só para o traço 2) e
p12 (causal para OS DOIS). p12 já é a probabilidade conjunta; não se
multiplica por p1·p2.

    lH0 = 0
    lH1 = log p1  + logsum(l1)
    lH2 = log p2  + logsum(l2)
    lH3 = log p1 + log p2 + logdiff(logsum(l1) + logsum(l2), logsum(l1 + l2))
    lH4 = log p12 + logsum(l1 + l2)

onde l1, l2 são os log-fatores de Bayes de Wakefield por variante:

    r = W / (W + V),   z = beta / se,   V = se²
    log ABF = ½ · [ log(1 − r) + r · z² ]

Prior de efeito W (variância), como no `coloc`:
    caso-controle  : desvio-padrão 0,20       → W = 0,04
    quantitativo   : desvio-padrão 0,15 · sdY → W = 0,0225 · sdY²

Tudo em espaço logarítmico: regiões com z grande estouram exp() em float.

Histórico (correção de 2026-11, auditoria G1/G3/G4): a versão anterior
multiplicava o termo de H4 por p1·p2·p12, usava sqrt(r) no lugar de
sqrt(1 − r) no fator de Bayes e W = 0,04² rotulado "coloc padrão". Ver
CORRECOES.md.
"""
import math

__all__ = [
    "HIPOTESES", "PRIORS_PADRAO", "W_CASO_CONTROLE", "SD_PRIOR_QUANT",
    "prior_w", "logsum", "logdiff", "log_abf", "abf", "combinar",
    "snp_pp_h4", "coloc", "h4_condicional", "conjunto_credivel",
    "razao_h4_h3_priori", "sensibilidade_p12", "sensibilidade_w",
    "P12_SENSIBILIDADE", "W_SENSIBILIDADE", "estimar_sdY",
    "md_tabela_hipoteses", "md_sensibilidade_p12", "md_sensibilidade_w",
]

HIPOTESES = {
    "H0": "nenhuma associação na região",
    "H1": "associação só com o traço 1 (GWAS)",
    "H2": "associação só com o traço 2 (QTL)",
    "H3": "os dois traços associados, variantes causais distintas",
    "H4": "os dois traços associados, uma variante causal compartilhada",
}

PRIORS_PADRAO = (1e-4, 1e-4, 1e-5)          # p1, p2, p12 — padrão do coloc
W_CASO_CONTROLE = 0.2 ** 2                  # variância; sd = 0,2 (coloc)
SD_PRIOR_QUANT = 0.15                       # × sdY (coloc)
P12_SENSIBILIDADE = (1e-5, 5e-6, 1e-6)      # com p1 = p2 = 1e-4 fixos
W_SENSIBILIDADE = (0.0016, 0.01, 0.0225, 0.04, 0.1)


# --------------------------------------------------------------- prior de efeito
def prior_w(tipo, sdY=1.0):
    """Variância do prior de efeito por tipo de traço (valores do `coloc`)."""
    if tipo == "cc":
        return W_CASO_CONTROLE
    if tipo == "quant":
        return (SD_PRIOR_QUANT * sdY) ** 2
    raise ValueError(f"tipo de traço desconhecido: {tipo!r} (use 'cc' ou 'quant')")


def estimar_sdY(ses, mafs, n):
    """Estimativa de sdY a partir de erros-padrão, MAF e n (como coloc:::sdY.est):
    Var(beta) ≈ sdY² / (2·n·maf·(1−maf))  ⇒  sdY ≈ sqrt(mediana(se²·2n·maf(1−maf))).
    Serve para conferir a suposição sdY = 1 em fenótipos normalizados."""
    vals = sorted(s * s * 2.0 * n * m * (1.0 - m)
                  for s, m in zip(ses, mafs) if s > 0 and 0.0 < m < 1.0)
    if not vals:
        raise ValueError("sem variantes válidas para estimar sdY")
    k = len(vals)
    med = vals[k // 2] if k % 2 else 0.5 * (vals[k // 2 - 1] + vals[k // 2])
    return math.sqrt(med)


# ------------------------------------------------------------------ log-álgebra
def logsum(xs):
    """log(Σ exp(x)) estável."""
    m = max(xs)
    if m == -math.inf:
        return -math.inf
    return m + math.log(sum(math.exp(x - m) for x in xs))


def logdiff(a, b):
    """log(exp(a) − exp(b)) estável; −inf quando a ≤ b (diferença nula ou
    negativa por arredondamento — ex.: região de uma única variante)."""
    if not a > b:
        return -math.inf
    return a + math.log1p(-math.exp(b - a))


# ---------------------------------------------------------------- fator de Bayes
def log_abf(beta, se, w):
    """log do fator de Bayes aproximado de Wakefield (associação vs nulo)."""
    if not se > 0:
        raise ValueError(f"erro-padrão inválido: {se!r}")
    v = se * se
    r = w / (w + v)
    z = beta / se
    return 0.5 * (math.log1p(-r) + r * z * z)


def abf(beta, se, w):
    """Fator de Bayes de Wakefield em escala natural (para fine-mapping de uma
    região; pode estourar com z muito grande — prefira log_abf)."""
    return math.exp(log_abf(beta, se, w))


# ------------------------------------------------------------------- combinação
def _validar(l1, l2):
    if len(l1) != len(l2):
        raise ValueError(f"tamanhos diferentes: {len(l1)} × {len(l2)}")
    if not l1:
        raise ValueError("nenhuma variante")


def combinar(l1, l2, p1, p2, p12):
    """Posteriores [H0, H1, H2, H3, H4] a partir dos log ABF por variante.
    Equivale a `coloc:::combine.abf`."""
    _validar(l1, l2)
    soma = [a + b for a, b in zip(l1, l2)]
    ls1, ls2, ls12 = logsum(l1), logsum(l2), logsum(soma)
    lh = [0.0,
          math.log(p1) + ls1,
          math.log(p2) + ls2,
          math.log(p1) + math.log(p2) + logdiff(ls1 + ls2, ls12),
          math.log(p12) + ls12]
    den = logsum(lh)
    return [math.exp(x - den) for x in lh]


def snp_pp_h4(l1, l2):
    """Posterior de cada variante ser a causal compartilhada, dado H4
    (coluna SNP.PP.H4 do coloc)."""
    _validar(l1, l2)
    soma = [a + b for a, b in zip(l1, l2)]
    den = logsum(soma)
    return [math.exp(x - den) for x in soma]


def coloc(traco1, traco2, tipo1="cc", tipo2="quant", sdY1=1.0, sdY2=1.0,
          priors=PRIORS_PADRAO, w1=None, w2=None):
    """Colocalização de dois traços na mesma lista de variantes harmonizadas.

    traco1, traco2: listas de (beta, se), mesma ordem e mesmo alelo de efeito.
    tipo: 'cc' (caso-controle) ou 'quant'; define o prior de efeito, salvo
    w1/w2 explícitos (usados só em análise de sensibilidade).
    Devolve dict com pp (H0..H4), snp_pp_h4, n, w1, w2, priors, l1, l2.
    """
    w1 = prior_w(tipo1, sdY1) if w1 is None else w1
    w2 = prior_w(tipo2, sdY2) if w2 is None else w2
    l1 = [log_abf(b, s, w1) for b, s in traco1]
    l2 = [log_abf(b, s, w2) for b, s in traco2]
    p1, p2, p12 = priors
    return {"pp": combinar(l1, l2, p1, p2, p12),
            "snp_pp_h4": snp_pp_h4(l1, l2),
            "n": len(l1), "w1": w1, "w2": w2, "priors": tuple(priors),
            "l1": l1, "l2": l2}


# ---------------------------------------------------------------- leituras
def h4_condicional(pp):
    """H4 / (H3 + H4): probabilidade de variante compartilhada DADO que os dois
    traços têm sinal na região."""
    den = pp[3] + pp[4]
    return pp[4] / den if den > 0 else float("nan")


def conjunto_credivel(snp_pp, cobertura=0.95):
    """Índices do menor conjunto de variantes cuja posterior soma ≥ cobertura."""
    ordem = sorted(range(len(snp_pp)), key=lambda i: -snp_pp[i])
    acum, fora = 0.0, []
    for i in ordem:
        fora.append(i)
        acum += snp_pp[i]
        if acum >= cobertura - 1e-12:
            break
    return fora


def razao_h4_h3_priori(n, p1, p2, p12):
    """Razão H4/H3 que o prior sozinho produz quando as n variantes são
    equivalentes (LD perfeito): n·p12 / (n·(n−1)·p1·p2). Diagnóstico de quanto
    da separação H3×H4 vem do prior e quanto vem dos dados."""
    if n < 2:
        return float("inf")
    return p12 / (p1 * p2 * (n - 1))


# ------------------------------------------------------------- sensibilidade
def sensibilidade_p12(traco1, traco2, tipo1="cc", tipo2="quant", sdY1=1.0,
                      sdY2=1.0, p1=1e-4, p2=1e-4, valores=P12_SENSIBILIDADE):
    """[(p12, [H0..H4])] com p1 e p2 FIXOS — a sensibilidade que importa:
    p12 é o único prior que distingue H4 de H3."""
    return [(p12, coloc(traco1, traco2, tipo1, tipo2, sdY1, sdY2,
                        priors=(p1, p2, p12))["pp"]) for p12 in valores]


def sensibilidade_w(traco1, traco2, priors=PRIORS_PADRAO,
                    valores=W_SENSIBILIDADE):
    """[(W, [H0..H4])] com o mesmo W nos dois traços."""
    return [(w, coloc(traco1, traco2, priors=priors, w1=w, w2=w)["pp"])
            for w in valores]


# ------------------------------------------------------------- saída markdown
def _fmt(x):
    """Probabilidade com 4 casas, ou notação científica quando pequena."""
    if x == 0:
        return "0"
    return f"{x:.4f}" if x >= 1e-4 else f"{x:.2e}"


def md_tabela_hipoteses(pp, rotulo1="GWAS sCJD", rotulo2="QTL"):
    desc = {
        "H0": "nenhuma associação na região",
        "H1": f"associação só com {rotulo1}",
        "H2": f"associação só com {rotulo2}",
        "H3": "os dois associados, variantes causais distintas",
        "H4": "os dois associados, uma variante causal compartilhada",
    }
    L = ["| Hipótese | Descrição | Posterior |", "|---|---|---|"]
    for i, h in enumerate(("H0", "H1", "H2", "H3", "H4")):
        neg = "**" if h == "H4" else ""
        L.append(f"| {neg}{h}{neg} | {desc[h]} | {neg}{_fmt(pp[i])}{neg} |")
    return L


def md_sensibilidade_p12(traco1, traco2, tipo1="cc", tipo2="quant",
                         sdY1=1.0, sdY2=1.0):
    L = ["| p12 (p1 = p2 = 10⁻⁴) | p12/(p1·p2) | H3 | H4 | H4/(H3+H4) |",
         "|---|---|---|---|---|"]
    for p12, pp in sensibilidade_p12(traco1, traco2, tipo1, tipo2, sdY1, sdY2):
        nota = " (padrão do coloc)" if p12 == 1e-5 else ""
        L.append(f"| {p12:g}{nota} | {p12 / 1e-8:g} | {_fmt(pp[3])} | "
                 f"{_fmt(pp[4])} | {_fmt(h4_condicional(pp))} |")
    return L


def md_sensibilidade_w(traco1, traco2, priors=PRIORS_PADRAO):
    L = ["| W (igual nos dois traços) | H3 | H4 |", "|---|---|---|"]
    for w, pp in sensibilidade_w(traco1, traco2, priors):
        L.append(f"| {w:g} | {_fmt(pp[3])} | {_fmt(pp[4])} |")
    return L
