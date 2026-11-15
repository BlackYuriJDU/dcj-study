#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analise_gse140069.py — v4 (correção da auditoria de 2026-11: M1, M2, M3, M5).

Dataset REAL GSE140069 — miRNA de sangue total, sCJD vs. controles
(MRC Prion Unit / Simon Mead, Nat Commun 2020, PMID 32769986).

HISTÓRICO DE VERSÕES:
- v1: Welch em escala LINEAR com log2FC por razão de médias — FRÁGIL (abandonada).
- v2: log2(x+1) antes do Welch — correta na máquina, mas SEM covariáveis.
- v3: v2 + regressão linear ajustada por IDADE+SEXO+RIN + filtro de detecção +
  Cohen's d. Solver por equações normais, sem detecção de resposta constante.
- v4 (esta): mesma análise, com
  * regressão por QR (pipeline/lib/ols.py) e testes NÃO ESTIMÁVEIS marcados
    como tal, sem p-valor (74 miRNAs são constantes em todas as amostras; na
    v3 recebiam um "p" que era ruído de arredondamento);
  * universo do FDR declarado: os testes estimáveis. A contagem sobre os 939
    medidos (universo da v3) é reportada ao lado;
  * validação das categorias (grupo, sexo) e das covariáveis numéricas: valor
    inesperado ou ausente levanta erro em vez de virar 0 em silêncio.

Motivação do ajuste: casos 66.4 vs. controles 53.6 anos; RIN 5.59 vs 6.50.
Saída: pipeline/reports/relatorio_gse140069.md
"""
import gzip
import math
import datetime
from pathlib import Path

import sys

import openpyxl

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "lib"))
from ols import Desenho, MOTIVOS  # noqa: E402
XLSX = BASE / "data" / "GSE140069_dados_processados.xlsx"
MATRIX = BASE / "data" / "GSE140069_series_matrix.txt.gz"
REPORTS = BASE / "reports"

PISO = 0.0001
FRACAO_MIN_DETECCAO = 0.25   # filtro de detecção (auditoria M3): ≥25% das amostras acima do piso
NUCLEO_ARTIGO = ("hsa-miR-16-5p", "hsa-miR-93-5p", "hsa-let-7i-5p", "hsa-miR-106b-3p")


# ---------------------------------------------------------------- estatística
def tcdf_p(t: float, df: float) -> float:
    """p bicaudal da distribuição t (beta incompleta regularizada)."""
    def betacf(a, b, x):
        FPMIN = 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c, d = 1.0, 1.0 - qab * x / qap
        d = 1.0 / (d if abs(d) < FPMIN else d)
        h = d
        for m in range(1, 200):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) < FPMIN else d)
            c = 1.0 + aa / c
            c = c if abs(c) >= FPMIN else FPMIN
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) < FPMIN else d)
            c = 1.0 + aa / c
            c = c if abs(c) >= FPMIN else FPMIN
            dele = d * c
            h *= dele
            if abs(dele - 1.0) < 3e-12:
                break
        return h

    def ibeta(a, b, x):
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0
        front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                         + a * math.log(x) + b * math.log(1 - x))
        if x < (a + 1) / (a + b + 2):
            return front * betacf(a, b, x) / a
        return 1.0 - front * betacf(b, a, 1 - x) / b

    return min(1.0, ibeta(df / 2.0, 0.5, df / (df + t * t)))


def welch(xs, ys):
    n1, n2 = len(xs), len(ys)
    m1, m2 = sum(xs) / n1, sum(ys) / n2
    v1 = sum((x - m1) ** 2 for x in xs) / (n1 - 1)
    v2 = sum((y - m2) ** 2 for y in ys) / (n2 - 1)
    se2 = v1 / n1 + v2 / n2
    if se2 == 0:
        return None            # variância nula nos dois grupos: não estimável
    df = se2 ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    t = (m1 - m2) / math.sqrt(se2)
    return t, tcdf_p(abs(t), df)


def cohen_d(x1, x0):
    """d de Cohen (variância agrupada); nan se não houver variância."""
    m1, m0 = sum(x1) / len(x1), sum(x0) / len(x0)
    v1 = sum((v - m1) ** 2 for v in x1) / (len(x1) - 1)
    v0 = sum((v - m0) ** 2 for v in x0) / (len(x0) - 1)
    sp = math.sqrt(((len(x1) - 1) * v1 + (len(x0) - 1) * v0)
                   / (len(x1) + len(x0) - 2))
    return (m1 - m0) / sp if sp > 0 else float("nan")


def ols_grupo(desenho, xs_log, grupo):
    """OLS log2 ~ intercepto + grupo(1=sCJD) + sexo(1=M) + idade + RIN, sobre
    um `Desenho` já fatorado (coluna 1 = grupo).
    Devolve (beta_grupo, p, cohen_d, motivo): se o teste não é estimável,
    beta e p são None e `motivo` diz por quê."""
    a = desenho.ajustar(xs_log)
    if not a.estimavel:
        return None, None, float("nan"), a.motivo
    bg, sg = a.beta[1], a.se[1]
    x1 = [v for i, v in enumerate(xs_log) if grupo[i] == 1]
    x0 = [v for i, v in enumerate(xs_log) if grupo[i] == 0]
    return bg, min(1.0, tcdf_p(abs(bg / sg), a.gl)), cohen_d(x1, x0), ""


def fdr_bh(pares):
    """Benjamini–Hochberg: [(id, p)] → {id: q}. O universo é exatamente a
    lista recebida; ids repetidos são erro (o dict esconderia um teste)."""
    ids = [m for m, _ in pares]
    if len(set(ids)) != len(ids):
        raise ValueError("fdr_bh: identificadores repetidos no universo")
    m = len(pares)
    ordenado = sorted(pares, key=lambda t: t[1])
    prev, out = 1.0, {}
    for k in range(m - 1, -1, -1):
        prev = min(prev, ordenado[k][1] * m / (k + 1))
        out[ordenado[k][0]] = prev
    return out


# ---------------------------------------------------------------- dados
GRUPOS_VALIDOS = {"sCJD": 1, "Control": 0}
SEXOS_VALIDOS = {"Male": 1, "Female": 0}


def _numero(texto, campo, titulo):
    """float(texto) ou None se o campo está declaradamente ausente; qualquer
    outra coisa é erro (antes, um valor não numérico virava None em silêncio)."""
    t = (texto or "").strip()
    if t == "" or t.lower() in ("na", "n/a", "no data", "not available"):
        return None
    try:
        v = float(t)
    except ValueError:
        raise ValueError(f"{campo} não numérico em {titulo!r}: {texto!r}")
    if not math.isfinite(v):
        raise ValueError(f"{campo} não finito em {titulo!r}: {texto!r}")
    return v


def carregar_covariatas():
    """Título de amostra (ex.: 'Control_23463_smallRNASeq') -> covariáveis.
    Chave = !Sample_title do series matrix, que usa a MESMA nomenclatura das
    colunas do xlsx. Categorias fora do esperado levantam erro."""
    covs = {}
    titulos = []
    linhas_cov = []
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith("!Sample_title\t"):
                titulos = [s.strip('"') for s in linha.rstrip("\n").split("\t")[1:]]
            elif linha.startswith("!Sample_characteristics_ch1\t"):
                linhas_cov.append([s.strip('"') for s
                                   in linha.rstrip("\n").split("\t")[1:]])
    for j, titulo in enumerate(titulos):
        d = {}
        for linha in linhas_cov:
            if j >= len(linha) or ":" not in linha[j]:
                continue
            k, v = linha[j].split(":", 1)
            d[k.strip().lower()] = v.strip()
        status = d.get("disease status", "")
        if status not in GRUPOS_VALIDOS:
            raise ValueError(f"grupo inesperado em {titulo!r}: {status!r} "
                             f"(esperado: {sorted(GRUPOS_VALIDOS)})")
        sexo_txt = d.get("sex", "")
        if sexo_txt == "":
            sexo = None                                  # ausência explícita
        elif sexo_txt in SEXOS_VALIDOS:
            sexo = SEXOS_VALIDOS[sexo_txt]
        else:
            raise ValueError(f"sexo inesperado em {titulo!r}: {sexo_txt!r} "
                             f"(esperado: {sorted(SEXOS_VALIDOS)})")
        covs[titulo] = {
            "grupo": GRUPOS_VALIDOS[status],
            "sexo": sexo,
            "idade": _numero(d.get("age at sampling"), "idade", titulo),
            "rin": _numero(d.get("rna integrity number (rin)"), "RIN", titulo),
        }
    return covs


def carregar():
    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    linhas = list(ws.iter_rows(values_only=True))
    cabecalho = [str(c) if c is not None else "" for c in linhas[0]]
    colunas = [(j, nome) for j, nome in enumerate(cabecalho)
               if j >= 5 and "_smallRNASeq" in nome]
    mirnas, vals = [], []
    n_cols = len(colunas)
    for linha in linhas[1:]:
        if not linha or not linha[4]:
            continue
        try:
            nums = [float(linha[j]) for j, _ in colunas]
        except (TypeError, ValueError, IndexError):
            continue
        if len(nums) != n_cols:
            continue
        mirnas.append(str(linha[4]))
        vals.append(nums)
    grupos = [nome.split("_")[0] for _, nome in colunas]
    nomes = [nome for _, nome in colunas]
    return mirnas, grupos, nomes, vals


def preparar():
    """Carrega, valida e devolve tudo o que a análise precisa.
    Levanta erro se: nomes de miRNA repetidos; grupo da coluna da planilha fora
    de {sCJD, Control} ou discordante do series matrix; amostra sem covariáveis."""
    mirnas, grupos, nomes, vals = carregar()
    if len(set(mirnas)) != len(mirnas):
        raise ValueError("nomes de miRNA repetidos na planilha")
    covmap = carregar_covariatas()
    grupo = []
    for g, nome in zip(grupos, nomes):
        if g not in GRUPOS_VALIDOS:
            raise ValueError(f"grupo inesperado na coluna {nome!r}: {g!r}")
        if nome not in covmap:
            raise ValueError(f"amostra sem covariáveis no series matrix: {nome!r}")
        if covmap[nome]["grupo"] != GRUPOS_VALIDOS[g]:
            raise ValueError(f"grupo da planilha e do series matrix divergem "
                             f"em {nome!r}")
        grupo.append(GRUPOS_VALIDOS[g])
    sexo = [covmap[n]["sexo"] for n in nomes]
    idade = [covmap[n]["idade"] for n in nomes]
    rin = [covmap[n]["rin"] for n in nomes]
    # o modelo ajustado usa só amostras com TODAS as covariáveis presentes
    ok_ajust = [i for i in range(len(nomes)) if sexo[i] is not None
                and idade[i] is not None and rin[i] is not None]
    vals_log = [[math.log2(v + 1.0) for v in linha] for linha in vals]
    det = [sum(1 for v in linha if v > PISO) / len(linha) >= FRACAO_MIN_DETECCAO
           for linha in vals]
    return {"mirnas": mirnas, "nomes": nomes, "vals": vals,
            "vals_log": vals_log, "det": det, "grupo": grupo, "sexo": sexo,
            "idade": idade, "rin": rin, "ok_ajust": ok_ajust}


def analisar(d=None):
    """Roda A (Welch) e B (OLS ajustado) e devolve os resultados por miRNA e
    as contagens em cada universo de FDR. Não escreve nada."""
    d = d or preparar()
    mirnas, vals_log, det, grupo = d["mirnas"], d["vals_log"], d["det"], d["grupo"]
    ok = d["ok_ajust"]
    idx_cjd = [i for i, g in enumerate(grupo) if g == 1]
    idx_ct = [i for i, g in enumerate(grupo) if g == 0]
    desenho = Desenho([[1.0, grupo[i], d["sexo"][i], d["idade"][i], d["rin"][i]]
                       for i in ok])
    g_ok = [grupo[i] for i in ok]

    lfc, pA, pB, betaB, dB, motivo = {}, {}, {}, {}, {}, {}
    for k, m in enumerate(mirnas):
        linha = vals_log[k]
        xs = [linha[i] for i in idx_cjd]
        ys = [linha[i] for i in idx_ct]
        lfc[m] = sum(xs) / len(xs) - sum(ys) / len(ys)
        w = welch(xs, ys)
        if w is not None:
            pA[m] = w[1]
        bg, p, dd, mot = ols_grupo(desenho, [linha[i] for i in ok], g_ok)
        dB[m] = dd
        if p is None:
            motivo[m] = mot
        else:
            pB[m], betaB[m] = p, bg
    detectavel = {m for m, ok_ in zip(mirnas, det) if ok_}

    def sig(q):
        return {m for m, v in q.items() if v < 0.05}

    # Universo declarado: testes ESTIMÁVEIS. Universo legado (v3): os 939
    # medidos, com os não estimáveis contando como testes que nunca rejeitam.
    qA = fdr_bh(sorted(pA.items()))
    qB = fdr_bh(sorted(pB.items()))
    qA_939 = fdr_bh([(m, pA.get(m, 1.0)) for m in mirnas])
    qB_939 = fdr_bh([(m, pB.get(m, 1.0)) for m in mirnas])
    qB_det = fdr_bh(sorted((m, p) for m, p in pB.items() if m in detectavel))
    return {
        "dados": d, "desenho": desenho, "idx_cjd": idx_cjd, "idx_ct": idx_ct,
        "lfc": lfc, "pA": pA, "pB": pB, "betaB": betaB, "dB": dB,
        "nao_estimavel": motivo, "detectavel": detectavel,
        "qA": qA, "qB": qB, "qA_939": qA_939, "qB_939": qB_939, "qB_det": qB_det,
        "sigA": sig(qA), "sigB": sig(qB), "sigA_939": sig(qA_939),
        "sigB_939": sig(qB_939), "sigB_det": sig(qB_det),
    }


def main() -> None:
    r = analisar()
    d = r["dados"]
    mirnas, det, vals = d["mirnas"], d["det"], d["vals"]
    idx_cjd, idx_ct = r["idx_cjd"], r["idx_ct"]
    n_est_A, n_est_B = len(r["pA"]), len(r["pB"])
    n_nao = len(r["nao_estimavel"])
    motivos = sorted(set(r["nao_estimavel"].values()))
    sigA, sigB = r["sigA"], r["sigB"]
    lfc, betaB = r["lfc"], r["betaB"]
    up = lambda conj, ef: sum(1 for m in conj if ef[m] > 0)       # noqa: E731
    dn = lambda conj, ef: sum(1 for m in conj if ef[m] < 0)       # noqa: E731
    n_sexo_f = sum(1 for v in d["sexo"] if v == 0)
    n_sexo_m = sum(1 for v in d["sexo"] if v == 1)

    hoje = datetime.date.today().isoformat()
    L = [
        "# Relatório — GSE140069 (SANGUE) — v4 com ajuste de idade/sexo/RIN",
        f"*`analise_gse140069.py` v4 em {hoje}. Mesma análise da v3, com solver "
        "por QR, testes não estimáveis marcados e universo do FDR declarado "
        "(correção da auditoria de 2026-11; ver `CORRECOES.md`).*",
        "",
        f"- Amostras: {len(idx_cjd)} sCJD vs. {len(idx_ct)} controles · miRNAs medidos: {len(mirnas)}",
        f"- Covariáveis do series matrix: idade (casos ~66 vs. controles ~54 anos), "
        f"sexo ({n_sexo_f} F, {n_sexo_m} M, {sum(1 for v in d['sexo'] if v is None)} ausentes), RIN",
        f"- Amostras com sexo, idade e RIN completos (usadas no modelo ajustado): "
        f"{len(d['ok_ajust'])} ({len(d['nomes']) - len(d['ok_ajust'])} excluídas)",
        f"- Matriz de desenho do modelo ajustado: posto {r['desenho'].posto} de "
        f"{r['desenho'].p} colunas",
        f"- Filtro de detecção (≥{int(FRACAO_MIN_DETECCAO*100)}% das amostras acima do piso): "
        f"**{sum(det)} de {len(mirnas)}** miRNAs",
        "",
        "## Universo dos testes (declarado antes do FDR)",
        "",
        f"- **{n_nao} miRNAs não são estimáveis** "
        f"({'; '.join(MOTIVOS[m] for m in motivos)}). Não recebem p-valor e não "
        "entram no FDR.",
        f"- **Universo declarado: {n_est_B} testes estimáveis** (A e B). É o "
        "universo usado nas contagens principais.",
        f"- Universo legado ({len(mirnas)}, usado até a v3): todos os miRNAs "
        "medidos, com os não estimáveis contando como testes que nunca "
        "rejeitam. É válido, porém mais conservador; fica na última coluna "
        "para comparação com os números já publicados.",
        "- Na v3, os miRNAs constantes recebiam do solver antigo um \"p\" que "
        "era ruído de arredondamento: 0,0203 com números do tipo float do "
        "Python e da ordem de 10⁻⁶ a 10⁻⁷ com os mesmos números em "
        "numpy.float64 (o valor exato varia com o ambiente). O resultado "
        "publicado só não mudava porque a execução caiu do lado inofensivo.",
        "",
        "## A vs. B",
        "",
        f"| Análise | FDR<0,05 no universo declarado ({n_est_B}) | FDR<0,05 no universo legado ({len(mirnas)}) | Interpretação |",
        "|---|---|---|---|",
        f"| A: Welch log2(x+1), SEM covariáveis | {len(sigA)} ({up(sigA, lfc)}↑/{dn(sigA, lfc)}↓) | "
        f"{len(r['sigA_939'])} ({up(r['sigA_939'], lfc)}↑/{dn(r['sigA_939'], lfc)}↓) | triagem não ajustada — inflada por idade/RIN |",
        f"| A′: A ∩ filtro de detecção | {len(sigA & r['detectavel'])} | "
        f"{len(r['sigA_939'] & r['detectavel'])} | triagem em miRNAs bem medidos |",
        f"| B: OLS ajustado (idade+sexo+RIN) | **{len(sigB)}** ({up(sigB, betaB)}↑/{dn(sigB, betaB)}↓) | "
        f"**{len(r['sigB_939'])}** | **o número que vale** |",
        f"| A ∩ B | {len(sigA & sigB)} | {len(r['sigA_939'] & r['sigB_939'])} | núcleo robusto às covariáveis |",
        f"| B no universo filtrado por detecção (n={len(r['qB_det'])}) | **{len(r['sigB_det'])}** | — | sensibilidade com correção menor |",
        "",
        "**Veredicto: a maior parte do sinal bruto não sobrevive ao ajuste por "
        "idade, sexo e RIN.** O que resta é a direção (predomínio de miRNAs ↓) "
        "e o núcleo do artigo em significância NOMINAL (não em FDR) — ver "
        "tabela abaixo.",
        "",
        "## Núcleo da assinatura do artigo original (Nat Commun 2020)",
        "",
        f"| miRNA | p nominal (ajustado) | q A ({n_est_A}) | q B ({n_est_B}) | q B (universo filtrado) |",
        "|---|---|---|---|---|",
    ]

    def fq(q, m):
        return f"{q[m]:.2e}" if m in q else "não estimável"

    for m in NUCLEO_ARTIGO:
        pb = f"{r['pB'][m]:.4f}" if m in r["pB"] else "não estimável"
        L.append(f"| {m} | {pb} | {fq(r['qA'], m)} | {fq(r['qB'], m)} | "
                 f"{fq(r['qB_det'], m)} |")

    q93 = r["qB_det"].get("hsa-miR-93-5p")
    L += ["",
          "**Leitura**: os 4 mantêm direção ↓ e p nominal < 0,05; após FDR, "
          + (f"apenas miR-93-5p sobrevive no universo filtrado (q={q93:.3f}). "
             if q93 is not None and q93 < 0.05 and
             sum(1 for m in NUCLEO_ARTIGO if r["qB_det"].get(m, 1) < 0.05) == 1
             else "ver a última coluna. ")
          + "A assinatura publicada é mais frágil sob ajuste padrão do que a "
          "apresentação original sugere. Diferenças plausíveis: Partek GSA "
          "(correção de variância gene-específica) vs. OLS comum, e universo de "
          f"testes (101 no artigo; {len(r['qB_det'])} ou {n_est_B} aqui)."]

    L += ["", "## Significativos no modelo ajustado (B) — com tamanho de efeito",
          "", "| miRNA | β grupo (log2) | p | q(FDR) | d | detectado em sCJD | detectado em controles |",
          "|---|---|---|---|---|---|---|"]
    for m in sorted(sigB, key=lambda x: r["pB"][x])[:15]:
        k = mirnas.index(m)
        dc = sum(1 for i in idx_cjd if vals[k][i] > PISO)
        dt = sum(1 for i in idx_ct if vals[k][i] > PISO)
        L.append(f"| {m} | {betaB[m]:+.2f} | {r['pB'][m]:.2e} | {r['qB'][m]:.2e} | "
                 f"{r['dB'][m]:+.2f} | {dc}/{len(idx_cjd)} | {dt}/{len(idx_ct)} |")
    if "hsa-miR-500a-3p" in sigB:
        k = mirnas.index("hsa-miR-500a-3p")
        dc = sum(1 for i in idx_cjd if vals[k][i] > PISO)
        dt = sum(1 for i in idx_ct if vals[k][i] > PISO)
        L += ["",
              f"**Como ler o β de {betaB['hsa-miR-500a-3p']:.1f} de "
              "hsa-miR-500a-3p.** Não é uma queda de abundância de ~250 vezes. "
              f"O miRNA está acima do piso em {dc} de {len(idx_cjd)} casos "
              f"({dc / len(idx_cjd):.0%}) e em {dt} de {len(idx_ct)} controles "
              f"({dt / len(idx_ct):.0%}); nas demais amostras o valor é o piso "
              f"({PISO}). Em log2(x+1), a distância entre \"detectado\" e "
              "\"no piso\" é de ~17–20 unidades, de modo que o β mede sobretudo "
              "a diferença na taxa de detecção. O achado é do tipo detectado/"
              "não detectado, e um modelo linear sobre log2 não é a descrição "
              "ideal dele."]

    L += ["", "## Notas",
          "- v1 (linear), v2 (log2 sem covariáveis) e v3 (equações normais) estão no histórico do script; esta v4 é a análise vigente.",
          "- O artigo original usou Partek GSA com idade como covariável sobre 101 miRNAs filtrados;",
          f"  aqui o universo é o dos {n_est_B} testes estimáveis (ou os {len(r['qB_det'])} detectáveis) — universos diferentes, declarados.",
          "- A lista ajustada NÃO é 'assinatura': a assinatura validada do artigo são 3 miRNAs com qPCR.",
          "- Sexo codificado M=1; RIN como qualidade de RNA; modelo linear padrão, sem interações.",
          "- Grupo e sexo são validados contra as categorias esperadas; idade e RIN precisam ser numéricos. Valor inesperado interrompe a análise."]
    destino = REPORTS / "relatorio_gse140069.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    print(f"[ok] universo {n_est_B}: A={len(sigA)} A'={len(sigA & r['detectavel'])} "
          f"B={len(sigB)} | universo {len(mirnas)}: A={len(r['sigA_939'])} "
          f"A'={len(r['sigA_939'] & r['detectavel'])} B={len(r['sigB_939'])} | "
          f"B filtrado={len(r['sigB_det'])} | não estimáveis={n_nao}")


if __name__ == "__main__":
    main()
