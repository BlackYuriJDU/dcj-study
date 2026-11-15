#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analise_gse160208.py — Primeira análise sobre dados REAIS e anonimizados.

Dataset: GSE160208 (GEO/NCBI) — NanoString nCounter, painel neuroinflamação+,
córtex frontal (FC) e cerebelo (CB) de pacientes com sCJD vs. controles.
Fonte: Litman T. et al., Univ. Copenhagen, PMID 33375642. Licença: público NCBI.

Gera pipeline/reports/relatorio_gse160208.md com:
- composição das amostras; 
- top genes diferencialmente expressos (diferença de médias) em FC;
- verificação específica de PRNP;
- estatística inferencial no córtex frontal: teste t de Welch por gene e FDR
  de Benjamini–Hochberg, com o universo dos testes declarado (genes sem
  variância não são estimáveis e não recebem p-valor);
- contagem de genes diferencialmente expressos sob o critério do artigo
  original, com a ORDEM DOS FILTROS como parâmetro explícito (ver ORDEM_FILTROS).
"""
import gzip
import math
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
REPORTS = BASE / "reports"
MATRIX = DATA / "GSE160208_series_matrix.txt.gz"

# --- Critério de gene diferencialmente expresso (DEG) -----------------------
# O número de DEGs depende de DUAS escolhas que o artigo original não fixa de
# forma inequívoca; por isso são parâmetros, não detalhes escondidos no código:
#
#   universo do FDR:
#     "estimaveis" — só os genes com teste estimável (padrão declarado);
#     "painel"     — os 800 genes do painel; os não estimáveis contam como
#                    testes que nunca rejeitam (universo usado até a v0.4).
#   ordem dos filtros:
#     "fdr_no_universo"        — BH sobre o universo inteiro e, depois, o
#                                filtro de magnitude |Δlog2| > DELTA_MIN;
#     "magnitude_antes_do_fdr" — primeiro |Δlog2| > DELTA_MIN e, depois, BH só
#                                sobre os genes que passaram.
ORDENS_FILTROS = ("fdr_no_universo", "magnitude_antes_do_fdr")
UNIVERSOS = ("estimaveis", "painel")
ORDEM_FILTROS = "fdr_no_universo"
UNIVERSO = "estimaveis"
DELTA_MIN = 1.0          # |Δlog2| > 1 = mais de 2 vezes (critério do artigo)


def tcdf_p(t: float, df: float) -> float:
    """p bicaudal da distribuição t via beta incompleta regularizada."""
    def betacf(a: float, b: float, x: float) -> float:
        MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c, d = 1.0, 1.0 - qab * x / qap
        if abs(d) < FPMIN:
            d = FPMIN
        d = 1.0 / d
        h = d
        for m in range(1, MAXIT + 1):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN:
                d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN:
                c = FPMIN
            d = 1.0 / d
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN:
                d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN:
                c = FPMIN
            d = 1.0 / d
            dele = d * c
            h *= dele
            if abs(dele - 1.0) < EPS:
                break
        return h

    def ibeta(a: float, b: float, x: float) -> float:
        if x <= 0:
            return 0.0
        if x >= 1:
            return 1.0
        lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        front = math.exp(lbeta + a * math.log(x) + b * math.log(1 - x))
        if x < (a + 1) / (a + b + 2):
            return front * betacf(a, b, x) / a
        return 1.0 - front * betacf(b, a, 1 - x) / b

    return ibeta(df / 2.0, 0.5, df / (df + t * t))


def welch(xs: list[float], ys: list[float]):
    """Teste t de Welch bicaudal: retorna (t, p), ou None se não estimável."""
    n1, n2 = len(xs), len(ys)
    m1, m2 = sum(xs) / n1, sum(ys) / n2
    v1 = sum((x - m1) ** 2 for x in xs) / (n1 - 1)
    v2 = sum((y - m2) ** 2 for y in ys) / (n2 - 1)
    se2 = v1 / n1 + v2 / n2
    if se2 == 0:
        return None            # variância nula nos dois grupos: não estimável
    t = (m1 - m2) / math.sqrt(se2)
    df = se2 ** 2 / ((v1 / n1) ** 2 / (n1 - 1) + (v2 / n2) ** 2 / (n2 - 1))
    return t, min(1.0, tcdf_p(abs(t), df))


def extrair_tabela() -> tuple[list[str], list[list[float]], dict]:
    """Lê a série matrix e retorna (genes, valores[genes][amostras], metadados)."""
    meta: dict = {}
    genes: list[str] = []
    linhas_vals: list[list[float]] = []
    amostras: list[str] = []          # títulos na ordem das colunas da tabela
    gsms: list[str] = []
    covs_raw: list[list[str]] = []    # cada linha !Sample_characteristics_ch1
    em_tabela = False
    with gzip.open(MATRIX, "rt", encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith("!Series_title"):
                meta["titulo"] = linha.split("\t")[1].strip().strip('"')
            elif linha.startswith("!Sample_title\t"):
                titulos = [s.strip('"') for s in linha.rstrip("\n").split("\t")[1:]]
            elif linha.startswith("!Sample_geo_accession\t"):
                gsms = [s.strip('"') for s in linha.rstrip("\n").split("\t")[1:]]
            elif linha.startswith("!Sample_characteristics_ch1\t"):
                covs_raw.append([s.strip('"') for s
                                 in linha.rstrip("\n").split("\t")[1:]])
            elif linha.startswith("!series_matrix_table_begin"):
                em_tabela = True
            elif linha.startswith("!series_matrix_table_end"):
                break
            elif em_tabela:
                partes = [c.strip().strip('"') for c in linha.rstrip("\n").split("\t")]
                if partes[0] in ("ID_REF", ""):
                    continue  # cabeçalho interno
                try:
                    nums = [float(x) for x in partes[1:]]
                except ValueError:
                    continue  # linha malformada: pula SEM desalinhar
                genes.append(partes[0])
                linhas_vals.append(nums)
    # Mapear colunas GSM -> título clínico (ordem da tabela = ordem dos GSMs)
    mapa = dict(zip(gsms, titulos))
    amostras = [mapa.get(gsm, gsm) for gsm in gsms]
    # Covariáveis: extrair "chave: valor" por amostra
    covs: list[dict[str, str]] = []
    for j in range(len(gsms)):
        d: dict[str, str] = {}
        for linha_c in covs_raw:
            if j < len(linha_c) and ":" in linha_c[j]:
                k, v = linha_c[j].split(":", 1)
                d[k.strip()] = v.strip()
        covs.append(d)
    return amostras, genes, linhas_vals, meta, covs


def bh(ps: list[float]) -> list[float]:
    """Benjamini–Hochberg por posição: q[i] corresponde a ps[i]. O universo é
    exatamente a lista recebida."""
    m = len(ps)
    ordem = sorted(range(m), key=lambda i: ps[i])
    q = [0.0] * m
    prev = 1.0
    for k in range(m - 1, -1, -1):
        prev = min(prev, ps[ordem[k]] * m / (k + 1))
        q[ordem[k]] = prev
    return q


def classificar_amostras(amostras: list[str]) -> dict[str, list[int]]:
    """Índices por grupo a partir do título (ex.: 'FFCJD-12_FC'). Título que
    não se encaixa em doença × região levanta erro (antes era ignorado)."""
    grupos: dict[str, list[int]] = {"CJD_FC": [], "CT_FC": [],
                                    "CJD_CB": [], "CT_CB": []}
    for i, nome in enumerate(amostras):
        doenca = "CJD" if "CJD" in nome else ("CT" if "CT" in nome else None)
        regiao = "FC" if "_FC" in nome else ("CB" if "_CB" in nome else None)
        if doenca is None or regiao is None:
            raise ValueError(f"amostra não classificável pelo título: {nome!r}")
        grupos[f"{doenca}_{regiao}"].append(i)
    return grupos


def testar_fc(amostras=None, genes=None, vals=None):
    """Teste de Welch por gene no córtex frontal (sCJD × controles).
    Devolve lista de dicts: gene, delta (Δlog2 = média sCJD − média controles),
    p (None se não estimável), estimavel."""
    if amostras is None:
        amostras, genes, vals, _meta, _covs = extrair_tabela()
    if len(set(genes)) != len(genes):
        raise ValueError("genes repetidos na matriz")
    grupos = classificar_amostras(amostras)
    out = []
    for gi, g in enumerate(genes):
        xs = [vals[gi][i] for i in grupos["CJD_FC"]]
        ys = [vals[gi][i] for i in grupos["CT_FC"]]
        w = welch(xs, ys)
        out.append({"gene": g, "delta": sum(xs) / len(xs) - sum(ys) / len(ys),
                    "p": None if w is None else w[1],
                    "estimavel": w is not None})
    return out


def q_no_universo(res, universo: str = UNIVERSO) -> dict[str, float]:
    """{gene: q} com BH sobre o universo declarado. No universo "painel" os
    genes não estimáveis entram com p = 1 (e nunca são significativos)."""
    if universo not in UNIVERSOS:
        raise ValueError(f"universo desconhecido: {universo!r}")
    if universo == "estimaveis":
        sel = [r for r in res if r["estimavel"]]
        ps = [r["p"] for r in sel]
    else:
        sel = list(res)
        ps = [r["p"] if r["estimavel"] else 1.0 for r in sel]
    return {r["gene"]: q for r, q in zip(sel, bh(ps))}


def significativos(res, universo: str = UNIVERSO, q_max: float = 0.05):
    """(up, down): genes estimáveis com q < q_max no universo declarado."""
    q = q_no_universo(res, universo)
    sig = [r for r in res if r["estimavel"] and q[r["gene"]] < q_max]
    return ({r["gene"] for r in sig if r["delta"] > 0},
            {r["gene"] for r in sig if r["delta"] < 0})


def degs_criterio_artigo(res, ordem: str = ORDEM_FILTROS,
                         universo: str = UNIVERSO, p_max: float = 0.05,
                         q_max: float = 0.06, delta_min: float = DELTA_MIN):
    """Conjunto de DEGs sob o critério do artigo original (p < 0,05; q ≤ 0,06;
    mais de 2 vezes), com a ordem dos filtros explícita."""
    if ordem not in ORDENS_FILTROS:
        raise ValueError(f"ordem de filtros desconhecida: {ordem!r}")
    est = [r for r in res if r["estimavel"]]
    if ordem == "fdr_no_universo":
        q = q_no_universo(res, universo)
        return {r["gene"] for r in est if r["p"] < p_max
                and q[r["gene"]] <= q_max and abs(r["delta"]) > delta_min}
    grandes = [r for r in est if abs(r["delta"]) > delta_min]
    qs = bh([r["p"] for r in grandes])
    return {r["gene"] for r, q in zip(grandes, qs)
            if r["p"] < p_max and q <= q_max}


def main() -> None:
    amostras, genes, vals, meta, covs = extrair_tabela()

    # --- Classificação das amostras pelo título: FFCJD-*_FC etc. ----------
    grupos = classificar_amostras(amostras)

    # --- Covariáveis por amostra (sexo, códon 129, subtipo) --------------
    def cov(i: int, chave: str) -> str:
        return covs[i].get(chave, "NA") if i < len(covs) else "NA"

    resumo_covs: dict[str, dict[str, int]] = {}
    for k_chave in ("gender", "codon 129", "cjd subtype"):
        contagem: dict[str, int] = {}
        for i in range(len(amostras)):
            contagem[cov(i, k_chave)] = contagem.get(cov(i, k_chave), 0) + 1
        resumo_covs[k_chave] = contagem

    # Subtipos presentes apenas em CJD (controles são NA)
    subtipos_cjd_fc: dict[str, list[int]] = {}
    for i in grupos["CJD_FC"]:
        st = cov(i, "cjd subtype")
        subtipos_cjd_fc.setdefault(st, []).append(i)

    # --- Diferença de médias por gene no córtex frontal -------------------
    def media(gi: int, idxs: list[int]) -> float:
        col = vals[gi]
        sel = [col[i] for i in idxs]
        return sum(sel) / len(sel) if sel else float("nan")

    fc_diffs = []
    for gi, g in enumerate(genes):
        m_cjd, m_ct = media(gi, grupos["CJD_FC"]), media(gi, grupos["CT_FC"])
        fc_diffs.append((g, m_cjd - m_ct, m_cjd, m_ct))
    fc_diffs.sort(key=lambda t: t[1])

    agora = datetime.date.today().isoformat()
    L = [
        "# Relatório — Análise de dados REAIS: GSE160208",
        f"*Gerado por `analise_gse160208.py` em {agora}.*",
        "",
        f"- Dataset: {meta.get('titulo', 'GSE160208')}",
        "- Fonte: GEO/NCBI GSE160208 · PMID 33375642 · Univ. Copenhagen (dados públicos anonimizados)",
        f"- Amostras totais processadas: **{len(amostras)}** · Genes no painel: **{len(genes)}**",
        "",
        "## Composição das amostras",
    ]
    for k in sorted(grupos):
        L.append(f"- {k}: {len(grupos[k])}")

    L += ["", "## Covariáveis disponíveis (metadados reais)", ""]
    for k_chave, contagem in resumo_covs.items():
        partes = ", ".join(f"{v}: {n}" for v, n in sorted(contagem.items()))
        L.append(f"- {k_chave}: {partes}")

    # --- Estratificação por subtipo (MM1 = subtipo do caso de referência) ---------
    mm1_fc = subtipos_cjd_fc.get("MM1", [])
    L += ["", "## Estratificação por subtipo — córtex frontal CJD",
          f"- Grupos CJD-FC por subtipo: "
          + ", ".join(f"{st}: {len(ix)}" for st, ix in sorted(subtipos_cjd_fc.items()))]
    if mm1_fc:
        # Top 5 genes MM1 vs controles (Δ de médias) — foco no subtipo do caso de referência
        diffs_mm1 = []
        for gi, g in enumerate(genes):
            m_mm1 = media(gi, mm1_fc)
            diffs_mm1.append((g, m_mm1 - media(gi, grupos["CT_FC"])))
        diffs_mm1.sort(key=lambda t: t[1])
        L += ["", f"### Subgrupo MM1 (n={len(mm1_fc)}) vs. controles — top 5 up/down",
              "", "| Gene | Δ(MM1−CT) |", "|---|---|"]
        for g, d in reversed(diffs_mm1[-5:]):
            L.append(f"| {g} | +{d:.2f} |")
        for g, d in diffs_mm1[:5]:
            L.append(f"| {g} | {d:.2f} |")

    L += ["", "## Top 10 genes MAIS expressos em CJD (córtex frontal)",
          "", "| Gene | Média CJD | Média CT | Δ |", "|---|---|---|---|"]
    for g, d, mc, mt in reversed(fc_diffs[-10:]):
        L.append(f"| {g} | {mc:.1f} | {mt:.1f} | +{d:.1f} |")

    L += ["", "## Top 10 genes MENOS expressos em CJD (córtex frontal)",
          "", "| Gene | Média CJD | Média CT | Δ |", "|---|---|---|---|"]
    for g, d, mc, mt in fc_diffs[:10]:
        L.append(f"| {g} | {mc:.1f} | {mt:.1f} | {d:.1f} |")

    prnp = [(g, d) for g, d, *_ in fc_diffs if "PRNP" in g.upper()]
    L += ["", "## Verificação específica"]
    L.append(f"- PRNP presente no painel: {'sim → Δ(CJD−CT) = %+.2f' % prnp[0][1] if prnp else 'não (painel é de neuroinflamação)'}")

    # --- Estatística inferencial: Welch t-test + BH-FDR -------------------
    res = testar_fc(amostras, genes, vals)
    n_est = sum(1 for r in res if r["estimavel"])
    nao_est = [r["gene"] for r in res if not r["estimavel"]]
    q_decl = q_no_universo(res, "estimaveis")
    up, down = significativos(res, "estimaveis")
    up_p, down_p = significativos(res, "painel")
    sig = sorted(((r["gene"], r["p"], q_decl[r["gene"]]) for r in res
                  if r["estimavel"] and q_decl[r["gene"]] < 0.05),
                 key=lambda t: t[1])

    L += ["", "## Estatística inferencial (Welch + BH-FDR, córtex frontal)",
          "",
          f"- Genes no painel: {len(res)} · **não estimáveis: {len(nao_est)}** "
          "(mesmo valor em todas as amostras do córtex frontal — o piso do "
          "ensaio; sem variância não há teste nem p-valor).",
          f"- **Universo declarado do FDR: {n_est} genes estimáveis.** "
          f"Significativos com FDR<0.05: **{len(up) + len(down)}** "
          f"({len(up)} ↑, {len(down)} ↓).",
          f"- Universo do painel inteiro ({len(res)} genes, com os não "
          "estimáveis contando como testes que nunca rejeitam; usado até a "
          f"v0.4): **{len(up_p) + len(down_p)}** ({len(up_p)} ↑, "
          f"{len(down_p)} ↓). A diferença entre as duas contagens é só o "
          "tamanho do universo."]
    if sig:
        L += ["", "| Gene | p | q(FDR) |", "|---|---|---|"]
        L += [f"| {g} | {p:.2e} | {q:.2e} |" for g, p, q in sig[:15]]

    # --- Critério do artigo original: a contagem depende da ordem dos filtros
    L += ["", "## Critério do artigo original: 184 genes, e do que esse número depende",
          "",
          "O artigo (Areškevičiūtė et al. 2020) reporta **184** genes "
          "diferencialmente expressos no córtex frontal, com \"p<0.05, "
          "q=0.06, >2-fold change\" (suplemento S1b), calculados no Qlucore. "
          "Aplicando p < 0,05, q ≤ 0,06 e |Δlog2| > 1 aos dados depositados, a "
          "contagem depende da ordem dos filtros e do universo do FDR:", "",
          "| Ordem dos filtros | Universo do FDR | DEGs |", "|---|---|---|"]
    rot_o = {"fdr_no_universo": "FDR no universo inteiro, depois magnitude",
             "magnitude_antes_do_fdr": "magnitude primeiro, depois FDR só nos que passaram"}
    n_grandes = sum(1 for r in res if r["estimavel"] and abs(r["delta"]) > DELTA_MIN)
    for ordem in ORDENS_FILTROS:
        for universo in UNIVERSOS:
            if ordem == "magnitude_antes_do_fdr" and universo == "painel":
                continue                      # o universo passa a ser o subconjunto
            n = len(degs_criterio_artigo(res, ordem, universo))
            uni = (f"{n_est} estimáveis" if universo == "estimaveis"
                   else f"{len(res)} do painel")
            if ordem == "magnitude_antes_do_fdr":
                uni = f"{n_grandes} genes com \\|Δlog2\\| > 1"
            L.append(f"| {rot_o[ordem]} | {uni} | **{n}** |")
    n_nom = sum(1 for r in res if r["estimavel"] and r["p"] < 0.05
                and abs(r["delta"]) > DELTA_MIN)
    L += [f"| só p < 0,05 e \\|Δlog2\\| > 1 (sem exigir q) | — | **{n_nom}** |",
          "",
          "Leitura: a contagem de 184 é reproduzida quando o filtro de "
          "magnitude vem antes do FDR (ou sem exigir q); com o FDR calculado "
          "no universo inteiro, o resultado é 183. O gene que separa as duas "
          "contagens é " + ", ".join(sorted(
              degs_criterio_artigo(res, "magnitude_antes_do_fdr")
              - degs_criterio_artigo(res, "fdr_no_universo"))) + ", um caso-limite "
          "nos três cortes. Isto reproduz a CONTAGEM do artigo sob uma escolha "
          "declarada de ordem de filtros; não é uma réplica exata do "
          "procedimento dos autores: os q-valores do Qlucore não foram "
          "reproduzidos, e p e log2FC são próximos dos oficiais, não idênticos "
          "por construção. A comparação gene a gene com a lista do suplemento "
          "S1b está em `validacao_cruzada_gse160208_artigo_original.md`; o "
          "suplemento não é versionado neste repositório, e essa comparação "
          "não foi refeita na correção de 2026-11.",
          f"A ordem padrão deste script é `{ORDEM_FILTROS}`, no universo "
          f"`{UNIVERSO}`."]

    L += ["", "## Notas",
          "- Welch t-test bicaudal implementado em stdlib; FDR Benjamini–Hochberg.",
          "- Painel dirigido (800 genes neuroinflamatórios), não transcriptoma total.",
          "- Sem correção para covariáveis (idade, PMI) — os metadados brutos não as trazem.",
          "- **\"Perda neuronal\" é inferência deste projeto** (apoiada nos genes ↓ "
          "neurônio-específicos SLC17A6/NEFL/BDNF/TUBB3/GRIN2B), não conclusão do "
          "artigo original, que foca em regionalidade, microglia e células "
          "dendríticas. Ler como \"padrão consistente com disfunção neuronal\", "
          "não como demonstração histológica.",
          "- Genes não estimáveis (variância nula no córtex frontal): "
          + ", ".join(sorted(nao_est)) + "."]

    destino = REPORTS / "relatorio_gse160208.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    print(f"[ok] amostras={len(amostras)} genes={len(genes)} "
          f"CJD_FC={len(grupos['CJD_FC'])} CT_FC={len(grupos['CT_FC'])}")


if __name__ == "__main__":
    main()
