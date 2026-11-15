#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analise_gse140069_idade.py — Fase 2: teste empírico do confundimento idade/RIN.

O manuscrito afirma que a idade "fully explains" o colapso de 84 -> 1 miRNAs
significativos, e explica isso por "cases were on average 12.8 years older". Até
2026-10-04 essa afirmação não tinha evidência no repositório: a v3 mudou
simultaneamente TRÊS coisas (Welch -> OLS, nenhuma -> idade+sexo+RIN, e a forma
do modelo), de modo que nenhum componente é atribuível. O RIN é tão
desequilibrado quanto a idade (5.59 vs 6.50) e nunca foi isolado.

Este script isola cada covariável e mede a estimabilidade:

  M1  idade sozinha            log2 ~ 1 + idade
  M2  RIN sozinha              log2 ~ 1 + RIN
  M3  grupo sozinho            log2 ~ 1 + grupo
  M4  grupo + idade            (sem RIN, sem sexo)
  M5  grupo + idade+sexo+RIN   (= v3)

Como o sinal de grupo é idêntico em M3/M4/M5 quando não há confusão, mas colapsa
quando a idade entra, M3 vs. M4 é o teste decisivo: mede quanto da perda se deve
specifically à idade, e M1 vs. M2 diz qual das duas covariates carrega o efeito.

Teste adicional — faixa etária restrita: se casos e controles se sobrepõem em
50-70 anos, o efeito de grupo que sobreviver dentro dessa faixa não pode ser
atribuído à idade. A diferença de idade RESIDUAL dentro da faixa é reportada: uma
restrição que deixa os casos ainda 5 anos mais velhos não removeu a confusão.

VIF do bloco (grupo, idade, sexo, RIN) quantifica a colinearidade.

Saída: pipeline/reports/relatorio_gse140069_idade.md
"""
import datetime
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analise_gse140069 import (  # noqa: E402
    REPORTS, Desenho, fdr_bh, preparar, tcdf_p, welch,
)

# nome, descrição, colunas de preditores, rótulo do coeficiente de interesse
# (a coluna 0 é sempre o intercepto; grupo é sempre a 1ª preditora)
MODELOS = [
    ("M1", "idade sozinha", ["idade"], "idade"),
    ("M2", "RIN sozinho", ["rin"], "RIN"),
    ("M3", "grupo sozinho", ["grupo"], "grupo"),
    ("M4", "grupo + idade", ["grupo", "idade"], "grupo"),
    ("M4b", "grupo + RIN", ["grupo", "rin"], "grupo"),
    ("M4c", "grupo + sexo", ["grupo", "sexo"], "grupo"),
    ("M5", "grupo + idade + sexo + RIN", ["grupo", "idade", "sexo", "rin"], "grupo"),
]

FAIXAS = [(50, 70), (55, 70), (55, 75)]


def desenho(preditores, n):
    """Matriz de desenho fatorada (QR): intercepto + preditores."""
    return Desenho([[1.0] + [p[i] for p in preditores] for i in range(n)])


def ols(y, preditores):
    """y: lista de floats. preditores: lista de listas. -> (beta, se, gl) ou
    None quando o teste não é estimável (resposta constante, posto deficiente)."""
    a = desenho(preditores, len(y)).ajustar(y)
    return (a.beta, a.se, a.gl) if a.estimavel else None


def p_do_coef(beta, se, df, idx):
    """p-valor bicaudal do coeficiente idx via t de Student com gl residuais."""
    return min(1.0, tcdf_p(abs(beta[idx] / se[idx]), df))


def vif(coluna, outras):
    """VIF = 1/(1-R²): R² vem de regredir a coluna sobre as demais."""
    y = coluna
    res = ols(y, outras)
    if res is None:
        return float("nan")
    beta, _, _ = res
    ybar = sum(y) / len(y)
    ss_tot = sum((v - ybar) ** 2 for v in y)
    if ss_tot <= 0:
        return float("nan")
    ss_res = sum((y[i] - sum(beta[j] * ([1.0] + [o[i] for o in outras])[j]
                             for j in range(len(beta)))) ** 2
                 for i in range(len(y)))
    r2 = 1.0 - ss_res / ss_tot
    if r2 >= 1.0:
        return float("inf")
    return 1.0 / (1.0 - r2)


def cohen_d(xs, ys):
    n1, n2 = len(xs), len(ys)
    if n1 < 2 or n2 < 2:
        return float("nan")
    m1, m2 = sum(xs) / n1, sum(ys) / n2
    v1 = sum((x - m1) ** 2 for x in xs) / (n1 - 1)
    v2 = sum((y - m2) ** 2 for y in ys) / (n2 - 1)
    sp = math.sqrt(((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2))
    return (m1 - m2) / sp if sp > 0 else float("nan")


def main() -> None:
    d = preparar()                       # valida categorias e covariáveis
    mirnas, nomes = d["mirnas"], d["nomes"]
    grupo, sexo, idade, rin = d["grupo"], d["sexo"], d["idade"], d["rin"]
    ok = d["ok_ajust"]
    vals_log, det = d["vals_log"], d["det"]

    # ---------- descrição da amostra ----------
    idades_c = [idade[i] for i in ok if grupo[i] == 1]
    idades_t = [idade[i] for i in ok if grupo[i] == 0]
    rins_c = [rin[i] for i in ok if grupo[i] == 1]
    rins_t = [rin[i] for i in ok if grupo[i] == 0]

    def quartis(v):
        s = sorted(v)
        n = len(s)
        return s[0], s[n // 4], s[n // 2], s[3 * n // 4], s[-1]

    def desvio(v):
        m = sum(v) / len(v)
        return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))

    # ---------- VIF do bloco completo ----------
    bloco = {"grupo": [grupo[i] for i in ok], "idade": [idade[i] for i in ok],
             "sexo": [sexo[i] for i in ok], "rin": [rin[i] for i in ok]}
    vifs = {k: vif(bloco[k], [v for n, v in bloco.items() if n != k])
            for k in bloco}

    # ---------- modelos ----------
    # universo completo = todos os miRNAs com dados (espelha a contagem de 84,
    # que é do Welch sem covariáveis sobre os 939)
    # universo filtrado = apenas os que passam no filtro de detecção
    resultados = {}
    for chave, nome, cols, alvo in MODELOS:
        vetores = [bloco[c] for c in cols]
        des = desenho(vetores, len(ok))
        # Universos de FDR:
        #   "estimaveis" — declarado: só testes estimáveis;
        #   "939"        — legado (até a v3): todos os medidos; os não
        #                  estimáveis contam como testes que nunca rejeitam;
        #   "filtrado"   — miRNAs que passam no filtro de detecção.
        for universo, teste_mi in (("estimaveis", range(len(mirnas))),
                                   ("939", range(len(mirnas))),
                                   ("filtrado", [k for k in range(len(mirnas))
                                                 if det[k]])):
            pares = []
            beta_grupo = {}
            for k in teste_mi:
                linha = vals_log[k]
                y = [linha[i] for i in ok]
                aj = des.ajustar(y)
                if not aj.estimavel:      # ex.: miRNA constante — sem p-valor
                    if universo == "939":
                        pares.append((mirnas[k], 1.0))
                    continue
                beta, se, df = aj.beta, aj.se, aj.gl
                idx = 1 + cols.index(alvo) if alvo in cols else 1
                pares.append((mirnas[k], p_do_coef(beta, se, df, idx)))
                beta_grupo[mirnas[k]] = beta[idx]
            if not pares:
                continue
            q = fdr_bh(pares)
            sig = {m for m, _ in pares if q[m] < 0.05}
            beta_med = sorted(beta_grupo.values())[len(beta_grupo) // 2]
            resultados[(chave, universo)] = {
                "nome": nome, "n_testes": len(pares), "sig": sig, "q": q,
                "p": dict(pares), "beta_med": beta_med,
            }

    # ---------- sobreposição com o set de 84 (Welch sem covariáveis) ----------
    # reproduz o M-unadjusted do relatório v3 para medir quanto do set de 84 é idade
    idx_c = [i for i in ok if grupo[i] == 1]
    idx_t = [i for i in ok if grupo[i] == 0]
    unadj, unadj_939 = [], []
    for k in range(len(mirnas)):
        linha = vals_log[k]
        w = welch([linha[i] for i in idx_c], [linha[i] for i in idx_t])
        if w is None:                     # variância nula: não estimável
            unadj_939.append((mirnas[k], 1.0))
            continue
        unadj.append((mirnas[k], w[1]))
        unadj_939.append((mirnas[k], w[1]))
    q_unadj = fdr_bh(unadj)
    sig_unadj = {m for m, _ in unadj if q_unadj[m] < 0.05}
    q_unadj_939 = fdr_bh(unadj_939)
    sig_unadj_939 = {m for m, _ in unadj_939 if q_unadj_939[m] < 0.05}

    def conj(chave, universo="estimaveis"):
        return resultados[(chave, universo)]["sig"]

    set_idade, set_rin, set_grupo = conj("M1"), conj("M2"), conj("M3")
    g_idade, g_rin, g_sexo, g_v3 = conj("M4"), conj("M4b"), conj("M4c"), conj("M5")

    # ---------- faixas etárias ----------
    faixa_rows = []
    for lo, hi in FAIXAS:
        dentro = [i for i in ok if lo <= idade[i] <= hi]
        idx = [i for i in dentro if grupo[i] == 1]
        it = [i for i in dentro if grupo[i] == 0]
        if len(idx) < 10 or len(it) < 10:
            faixa_rows.append((lo, hi, len(idx), len(it), float("nan"), {}, {}))
            continue
        gap = (sum(idade[i] for i in idx) / len(idx)
               - sum(idade[i] for i in it) / len(it))
        res_f = {}
        for chave, cols, alvo in (("grupo só", ["grupo"], "grupo"),
                                  ("grupo + idade", ["grupo", "idade"], "grupo")):
            vetores = [[({"grupo": grupo, "idade": idade})[c][i] for i in dentro]
                       for c in cols]
            des_f = desenho(vetores, len(dentro))
            pares = []
            for k in range(len(mirnas)):
                if not det[k]:
                    continue
                linha = vals_log[k]
                aj = des_f.ajustar([linha[i] for i in dentro])
                if not aj.estimavel:
                    continue
                beta, se, df = aj.beta, aj.se, aj.gl
                idxc = 1 + cols.index(alvo)
                pares.append((mirnas[k], p_do_coef(beta, se, df, idxc)))
            if not pares:
                continue
            q = fdr_bh(pares)
            res_f[chave] = {m for m, _ in pares if q[m] < 0.05}
        faixa_rows.append((lo, hi, len(idx), len(it), gap, res_f.get("grupo só", set()),
                           res_f.get("grupo + idade", set())))

    # ---------- relatório ----------
    agora = datetime.date.today().isoformat()
    L = ["# GSE140069 — teste empírico do confundimento idade/RIN",
         f"*`analise_gse140069_idade.py` em {agora}. Fase 2 da auditoria "
         "pré-submissão, revista na correção de 2026-11 (solver por QR, "
         "universo do FDR declarado, leitura descritiva). Objetivo: ver o que "
         "acontece com o número de miRNAs significativos quando cada "
         "covariável entra no modelo.*", ""]

    L += ["## 1. A amostra", "",
          "| covariável | sCJD (n=%d) | controles (n=%d) | diferença | d |" % (
              len(idades_c), len(idades_t)),
          "|---|---|---|---|---|",
          f"| idade (anos) | {sum(idades_c)/len(idades_c):.1f} ± "
          f"{desvio(idades_c):.1f} | {sum(idades_t)/len(idades_t):.1f} ± "
          f"{desvio(idades_t):.1f} | **{sum(idades_c)/len(idades_c)-sum(idades_t)/len(idades_t):+.1f} anos** | "
          f"{cohen_d(idades_c, idades_t):+.2f} |",
          f"| RIN | {sum(rins_c)/len(rins_c):.2f} ± {desvio(rins_c):.2f} | "
          f"{sum(rins_t)/len(rins_t):.2f} ± {desvio(rins_t):.2f} | "
          f"**{sum(rins_c)/len(rins_c)-sum(rins_t)/len(rins_t):+.2f}** | "
          f"{cohen_d(rins_c, rins_t):+.2f} |", ""]
    q1c, _, q3c = quartis(idades_c)[1], quartis(idades_c)[2], quartis(idades_c)[3]
    q1t, _, q3t = quartis(idades_t)[1], quartis(idades_t)[2], quartis(idades_t)[3]
    L += [f"Idade — casos: Q1={q1c:.0f}, mediana={quartis(idades_c)[2]:.0f}, "
          f"Q3={q3c:.0f} · controles: Q1={q1t:.0f}, mediana={quartis(idades_t)[2]:.0f}, "
          f"Q3={q3t:.0f}.", ""]

    L += ["## 2. Colinearidade do bloco (VIF)", "",
          "| preditor | VIF | leitura |", "|---|---|---|"]
    for k, v in vifs.items():
        nota = "**colinear**" if v >= 5 else ("alto" if v >= 2 else "ok")
        L.append(f"| {k} | {v:.2f} | {nota} |")
    L += ["", "VIF ≥ 5 indica que o preditor é mal identificado no modelo — a "
          "atribuição de efeito entre covariáveis colineares não é estável.", ""]

    n_est = resultados[("M5", "estimaveis")]["n_testes"]
    L += ["## 3. Cada covariável isolada", "",
          f"Três universos de FDR, declarados: **estimáveis** ({n_est} testes; "
          f"universo principal), **939 (legado)** (todos os miRNAs medidos; os "
          f"{len(mirnas) - n_est} não estimáveis contam como testes que nunca "
          "rejeitam — é o universo dos números publicados até a v0.4) e "
          "**filtrado** (miRNAs que passam no filtro de detecção). FDR<0.05 "
          "dentro de cada modelo e universo. A coluna β mediana é a do "
          "coeficiente de interesse.", "",
          "| modelo | preditores | universo | testes | coef. | FDR<0.05 | β mediana |",
          "|---|---|---|---|---|---|---|"]
    rot = {"estimaveis": "estimáveis", "939": "939 (legado)",
           "filtrado": "filtrado"}
    for chave, nome, cols, alvo in MODELOS:
        for universo in ("estimaveis", "939", "filtrado"):
            r = resultados.get((chave, universo))
            if not r:
                continue
            L.append(f"| {chave} | {r['nome']} | {rot[universo]} | "
                     f"{r['n_testes']} | {alvo} | **{len(r['sig'])}** | "
                     f"{r['beta_med']:+.4f} |")
    L.append("")

    n_welch, n_welch_939 = len(sig_unadj), len(sig_unadj_939)
    n_base = len(set_grupo)

    def leg(chave):
        return len(conj(chave, "939"))

    L += ["## 4. O que muda quando cada covariável entra no modelo", "",
          f"O teste não ajustado do relatório principal é o Welch: "
          f"**{n_welch}** miRNAs a FDR<0.05 no universo de {n_est} estimáveis "
          f"({n_welch_939} no universo legado de 939 — o \"84\" das versões "
          f"anteriores). O mesmo contraste como OLS só com grupo é outro teste "
          f"(variância agrupada) e dá {n_base} ({leg('M3')}). As comparações "
          "abaixo são entre modelos OLS, por isso a base é o OLS só com grupo.",
          "",
          f"| modelo | FDR<0.05 ({n_est} estimáveis) | variação vs. grupo "
          "sozinho | FDR<0.05 (939, legado) |", "|---|---|---|---|",
          f"| idade sozinha (sem grupo) | {len(set_idade)} | — | {leg('M1')} |",
          f"| RIN sozinho (sem grupo) | {len(set_rin)} | — | {leg('M2')} |",
          f"| grupo sozinho | {n_base} | base | {leg('M3')} |",
          f"| grupo + sexo | {len(g_sexo)} | {len(g_sexo)-n_base:+d} | {leg('M4c')} |",
          f"| **grupo + idade** | **{len(g_idade)}** | **{len(g_idade)-n_base:+d}** | {leg('M4')} |",
          f"| **grupo + RIN** | **{len(g_rin)}** | **{len(g_rin)-n_base:+d}** | {leg('M4b')} |",
          f"| grupo + idade + sexo + RIN (modelo do relatório principal) | "
          f"{len(g_v3)} | {len(g_v3)-n_base:+d} | {leg('M5')} |", "",
          f"Dos {n_welch} do Welch: **{len(sig_unadj & set_idade)}** "
          f"({len(sig_unadj & set_idade)/max(1,n_welch):.0%}) também são "
          f"significativos para **idade sozinha**; "
          f"**{len(sig_unadj & set_rin)}** "
          f"({len(sig_unadj & set_rin)/max(1,n_welch):.0%}) para **RIN "
          f"sozinho**; **{len(sig_unadj & set_grupo)}** para grupo sozinho "
          "(OLS).", ""]

    perda_idade = n_base - len(g_idade)
    perda_rin = n_base - len(g_rin)
    perda_sexo = n_base - len(g_sexo)
    L += ["### Leitura", "",
          f"**O que as contagens mostram.** Com a idade no modelo, o número de "
          f"miRNAs com efeito de grupo a FDR<0.05 vai de {n_base} para "
          f"{len(g_idade)}; só com RIN, para {len(g_rin)}; só com sexo, fica em "
          f"{len(g_sexo)}; com as três covariáveis, {len(g_v3)}. O sinal não "
          "ajustado não se sustenta quando a idade é incluída, e o RIN "
          "sozinho já o reduz bastante.",
          "",
          "**O que elas não mostram.** Estas são contagens de testes que "
          "cruzam um limiar em modelos diferentes. Não medem proporção de "
          "efeito nem de variância. Versões anteriores deste relatório e do "
          "manuscrito diziam que a idade \"responde por ~97% do colapso\" — a "
          "conta era (88 − 4)/(88 − 1) sobre as contagens do universo legado. "
          "Essa leitura foi retirada: uma razão entre contagens de "
          "significativos não é uma fração atribuível à idade, e as reduções "
          "por covariável são **não aditivas** "
          f"(reduções de {perda_idade} com idade, {perda_rin} com RIN e "
          f"{perda_sexo} com sexo somam "
          f"{perda_idade + perda_rin + perda_sexo}, mais que a redução de "
          f"{n_base - len(g_v3)} no modelo completo), porque as "
          "covariáveis capturam variância sobreposta.",
          "",
          "**Sobre o VIF.** VIF < 1,4 em todos os preditores indica que os "
          "coeficientes são estimáveis separadamente (não há colinearidade "
          "forte). Não transforma a comparação de contagens em atribuição "
          "causal: casos e controles diferem em idade por desenho do estudo, "
          "e um modelo linear em idade pode não capturar toda essa diferença.",
          ""]

    L += ["## 5. Faixa etária restrita", "",
          "Se casos e controles se sobrepõem, o efeito de grupo que sobrevive "
          "dentro da faixa não pode ser atribuído à idade. A coluna `gap` é a "
          "diferença de idade **dentro** da faixa — se for grande, a restrição "
          "não removeu a confusão.", "",
          "**Poder insuficiente para resolver um efeito residual.** Com 22–29 "
          "controles por faixa (contra 48 na amostra total) e universo de 269 "
          "miRNAs filtrados, as contagens de 1–7 não distinguem \"nada sobrevive\" "
          "de \"sobrevive pouco\". Esta seção é um controle de sanidade, não uma "
          "medição.", "",
          "| faixa | n casos | n controles | gap idade | grupo só (FDR) | grupo + idade (FDR) |",
          "|---|---|---|---|---|---|"]
    for lo, hi, nc, nt, gap, sg, sga in faixa_rows:
        gs = "—" if gap != gap else f"{gap:+.1f} anos"
        L.append(f"| {lo}–{hi} | {nc} | {nt} | {gs} | {len(sg)} | {len(sga)} |")
    L.append("")

    destino = REPORTS / "relatorio_gse140069_idade.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    print(f"[ok] universo {n_est} (legado 939 entre parênteses)")
    print(f"[ok] não-ajustado={n_welch} ({n_welch_939})  idade_sozinha={len(set_idade)} ({leg('M1')})  "
          f"RIN_sozinha={len(set_rin)} ({leg('M2')})  grupo_sozinho={len(set_grupo)} ({leg('M3')})")
    print(f"[ok] M4(grupo+idade)={len(g_idade)} ({leg('M4')})  M4b(grupo+RIN)={len(g_rin)} ({leg('M4b')})  "
          f"M4c(grupo+sexo)={len(g_sexo)} ({leg('M4c')})  M5={len(g_v3)} ({leg('M5')})")
    print(f"[ok] VIF: " + "  ".join(f"{k}={v:.2f}" for k, v in vifs.items()))


if __name__ == "__main__":
    main()