#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analise_integrada_cerebro_sangue.py — Sangue × cérebro: alvos validados de
miRNAs do sangue entre os genes ↑ do córtex frontal.

Pergunta: os miRNAs reduzidos no SANGUE de pacientes com sCJD (GSE140069) têm
alvos validados sobre-representados entre os genes mais expressos no CÉREBRO
(GSE160208)?

Método:
- alvos: miRTarBase 10.0, suporte "Functional MTI" (evidência forte), com
  miRNA humano E gene-alvo humano (as duas colunas de espécie);
- genes ↑: Welch + FDR<0,05 no córtex frontal, universo declarado em
  analise_gse160208.py (genes estimáveis do painel);
- teste hipergeométrico unicaudal (sobre-representação), em que
      N = genes testáveis do painel,
      K = alvos do miRNA QUE ESTÃO no painel testável,
      n = genes ↑,
      k = alvos ↑.
  K é restrito ao painel: um alvo que o painel não mede não pode aparecer
  entre os genes ↑ e não pertence ao universo.
- miRNA sem nenhum alvo no painel não é testável: fica fora do teste e do FDR.

Correção de 2026-11 (auditoria, I1–I5): a versão anterior usava K = todos os
alvos do miRNA no genoma com N = 800, ignorava a espécie do gene-alvo, punha
no FDR dois miRNAs sem alvo (p = 1) e lia o CSV com split(",").
"""
import csv
import math
import sys
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
REPORTS = BASE / "reports"
MTI = DATA / "mirtarbase_MTI.csv"
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analise_gse160208 as g  # noqa: E402

MIRNAS_SANGUE = [  # 4 do artigo original (todos ↓) + o que sobrevive ao FDR aqui
    "hsa-miR-16-5p", "hsa-miR-93-5p", "hsa-let-7i-5p", "hsa-miR-106b-3p",
    "hsa-miR-500a-3p",
]
SUPORTE = "Functional MTI"
ESPECIE = "hsa"


def carregar_degs_cerebro(universo: str = g.UNIVERSO):
    """(up, down, testaveis): genes ↑ e ↓ a FDR<0,05 e o conjunto de genes do
    universo do teste ("estimaveis": só os testáveis; "painel": os 800)."""
    res = g.testar_fc()
    up, down = g.significativos(res, universo)
    if universo == "estimaveis":
        testaveis = {r["gene"] for r in res if r["estimavel"]}
    else:
        testaveis = {r["gene"] for r in res}
    return up, down, testaveis


def alvos_mirtarbase(mirnas=MIRNAS_SANGUE, exigir_especie_do_alvo: bool = True):
    """{miRNA: {gene: espécie do alvo}} — Functional MTI, miRNA humano.
    Com exigir_especie_do_alvo (padrão), o gene-alvo também tem de ser humano.
    Leitura com o módulo csv (campos com vírgula entre aspas); linhas
    duplicadas do arquivo colapsam no dicionário."""
    out = {m: {} for m in mirnas}
    with open(MTI, encoding="utf-8-sig", errors="replace", newline="") as fh:
        for linha in csv.DictReader(fh):
            mirna = (linha.get("miRNA") or "").strip()
            if mirna not in out:
                continue
            if (linha.get("Species (miRNA)") or "").strip() != ESPECIE:
                continue
            if (linha.get("Support Type") or "").strip() != SUPORTE:
                continue
            especie_alvo = (linha.get("Species (Target Gene)") or "").strip()
            if exigir_especie_do_alvo and especie_alvo != ESPECIE:
                continue
            alvo = (linha.get("Target Gene") or "").strip()
            if alvo:
                out[mirna][alvo] = especie_alvo
    return out


def hipergeometrico(k, K, n, N):
    """P(X ≥ k), X ~ Hipergeométrica(N, K, n) — log-gama, sem dependências."""
    if not (0 <= K <= N and 0 <= n <= N):
        raise ValueError(f"parâmetros inválidos: K={K}, n={n}, N={N}")

    def lnchoose(a, b):
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)
    p = 0.0
    for i in range(max(k, 0, n + K - N), min(K, n) + 1):
        p += math.exp(lnchoose(K, i) + lnchoose(N - K, n - i) - lnchoose(N, n))
    return min(1.0, p)


def bh(ps):
    return g.bh(ps)


def minimo_para_significancia(K, n, N, alfa=0.05):
    """Menor k com P(X ≥ k) < alfa; None se nem k = K alcança."""
    for k in range(0, K + 1):
        if hipergeometrico(k, K, n, N) < alfa:
            return k
    return None


def analisar(universo: str = g.UNIVERSO, exigir_especie_do_alvo: bool = True,
             restringir_ao_painel: bool = True):
    """Roda o teste para cada miRNA. Os dois últimos parâmetros existem só
    para reproduzir, em teste, o cálculo antigo e as âncoras da auditoria."""
    up, down, testaveis = carregar_degs_cerebro(universo)
    alvos = alvos_mirtarbase(exigir_especie_do_alvo=exigir_especie_do_alvo)
    N, n = len(testaveis), len(up)
    linhas = []
    for mirna in MIRNAS_SANGUE:
        todos = set(alvos[mirna])
        no_painel = todos & testaveis
        K = len(no_painel) if restringir_ao_painel else len(todos)
        k = len(no_painel & up)
        testavel = K > 0
        linhas.append({
            "mirna": mirna, "alvos_total": len(todos),
            "alvos_nao_humanos": sum(1 for e in alvos[mirna].values()
                                     if e != ESPECIE),
            "alvos_no_painel": len(no_painel), "K": K, "k": k,
            "k_down": len(no_painel & down),
            "esperado": K * n / N if testavel else None,
            "p": hipergeometrico(k, K, n, N) if testavel else None,
            "p_down": hipergeometrico(len(no_painel & down), K, len(down), N)
            if testavel else None,
            "k_min": minimo_para_significancia(K, n, N) if testavel else None,
            "testavel": testavel,
        })
    testaveis_l = [r for r in linhas if r["testavel"]]
    for r, q in zip(testaveis_l, bh([r["p"] for r in testaveis_l])):
        r["q"] = q
    return {"linhas": linhas, "N": N, "n_up": n, "n_down": len(down),
            "n_testes": len(testaveis_l), "universo": universo}


def main() -> None:
    hoje = datetime.date.today().isoformat()
    print("[1/2] genes ↑ do cérebro e alvos do miRTarBase...")
    a = analisar()
    L = ["# Sangue × cérebro: alvos de miRNAs do sangue entre os genes ↑ do córtex frontal",
         f"*`analise_integrada_cerebro_sangue.py` em {hoje}. Análise "
         "exploratória; nenhum dos dois artigos originais fez esta "
         "integração.*", "",
         "> **Correção de 2026-11.** A versão anterior comparava todos os "
         "alvos de cada miRNA no genoma (48, 16, 11) com um universo de 800 "
         "genes, o que tornava o \"esperado\" e os p-valores incorretos; "
         "também aceitava alvos de rato e camundongo e incluía no FDR dois "
         "miRNAs sem nenhum alvo. A conclusão (nenhuma sobre-representação) "
         "não muda. Ver `CORRECOES.md`.", "",
         "**Desenho**: miRNAs ↓ no sangue (GSE140069) → alvos validados "
         "(miRTarBase 10.0, \"Functional MTI\", miRNA e gene-alvo humanos) → "
         "sobre-representação entre os genes ↑ no córtex frontal (GSE160208, "
         "FDR<0,05). Teste hipergeométrico unicaudal com "
         f"**N = {a['N']}** genes testáveis do painel NanoString e "
         f"**n = {a['n_up']}** genes ↑ ({a['n_down']} ↓). Para cada miRNA, "
         "K = alvos que pertencem ao painel testável.", "",
         "| miRNA (sangue) | alvos validados (humanos) | alvos no painel (K) | "
         "alvos ↑ (k) | esperado (K·n/N) | p | q (BH) | k mínimo para p<0,05 |",
         "|---|---|---|---|---|---|---|---|"]
    for r in a["linhas"]:
        if r["testavel"]:
            kmin = r["k_min"] if r["k_min"] is not None else "inalcançável"
            L.append(f"| {r['mirna']} | {r['alvos_total']} | {r['K']} | "
                     f"{r['k']} | {r['esperado']:.1f} | {r['p']:.2f} | "
                     f"{r['q']:.2f} | {kmin} |")
        else:
            L.append(f"| {r['mirna']} | {r['alvos_total']} | 0 | — | — | "
                     "não testável | — | — |")
    nao = [r["mirna"] for r in a["linhas"] if not r["testavel"]]
    sem_filtro = alvos_mirtarbase(exigir_especie_do_alvo=False)
    excluidos = {m: sorted(f"{gene} ({esp})" for gene, esp in d.items()
                           if esp != ESPECIE)
                 for m, d in sem_filtro.items()}
    excluidos = {m: v for m, v in excluidos.items() if v}
    L += ["",
          f"- Testes realizados: **{a['n_testes']}** (BH sobre esses "
          f"{a['n_testes']}). Não testáveis, fora do FDR: "
          + (", ".join(nao) if nao else "nenhum")
          + " — sem nenhum alvo \"Functional MTI\" humano no painel.",
          "- Alvos excluídos pelo filtro de espécie do gene-alvo: "
          + ("; ".join(f"{m}: {', '.join(v)}" for m, v in excluidos.items())
             if excluidos else "nenhum") + ".", "",
          "## Direção oposta (descritivo)", "",
          "Mesma conta para os genes ↓ do cérebro. Não é um controle negativo "
          "no sentido estrito (um miRNA reduzido não prevê alvos reduzidos, "
          "mas também não os exclui); serve só de comparação.", "",
          "| miRNA | alvos no painel (K) | alvos ↓ | p (hipergeométrico, ↓) |",
          "|---|---|---|---|"]
    for r in a["linhas"]:
        if r["testavel"]:
            L.append(f"| {r['mirna']} | {r['K']} | {r['k_down']} | "
                     f"{r['p_down']:.2f} |")
    testados = [r for r in a["linhas"] if r["testavel"]]
    kmax = max(r["K"] for r in testados)
    kmin_ = min(r["K"] for r in testados)
    L += ["", "## Leitura", "",
          "- **Resultado: nenhuma sobre-representação detectada.** Nos "
          f"{a['n_testes']} miRNAs testáveis, o número de alvos ↑ fica "
          "próximo do esperado ao acaso (p de "
          f"{min(r['p'] for r in testados):.2f} a "
          f"{max(r['p'] for r in testados):.2f}).",
          f"- **O teste tem pouco poder.** O painel mede {a['N']} genes "
          f"testáveis, dos quais {a['n_up']} ({a['n_up'] / a['N']:.0%}) estão "
          f"↑; cada miRNA tem só {kmin_} a {kmax} alvos no painel. A última "
          "coluna da tabela mostra quantos alvos ↑ seriam necessários para "
          "p < 0,05: com tão poucos alvos e uma taxa de base tão alta, só um "
          "enriquecimento quase total seria detectável.",
          "- **O que não se pode concluir.** Um resultado nulo com esse poder "
          "não apoia nem contraria uma relação entre os miRNAs do sangue e a "
          "expressão cerebral. Não dá base para afirmar mecanismo, nem para "
          "afirmar que os miRNAs refletem processos periféricos.",
          "- **Viés do painel.** O painel NanoString é dirigido a "
          "neuroinflamação: os genes medidos não são uma amostra do genoma, e "
          "o resultado vale só para esse universo.", ""]
    destino = REPORTS / "relatorio_integracao_cerebro_sangue.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[2/2] {destino}")


if __name__ == "__main__":
    main()
