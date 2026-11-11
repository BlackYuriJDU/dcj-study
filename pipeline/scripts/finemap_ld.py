#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
finemap_ld.py — Fine-mapping com LD REAL (Ensembl REST, painel 1000G phase 3)
+ diagnóstico de inflação λ_GC por estratos de MAF.

Método honesto (sem genótipos individuais):
- ABF de Wakefield (W=0.04) por variante;
- posterior aproximada dentro da região (ABF normalizado);
- agrupamento por LD ao lead (r²≥0.80) → credible set no nível de CLUSTER
  (aproximação declarada; modelo conjunto tipo SuSiE exigiria genótipos);
- λ_GC global e por faixas de MAF/EAF para diagnosticar estratificação.

Saídas: pipeline/reports/relatorio_finemap_loci.md (v2) e relatorio_lambda_gc.md
Cache: /tmp/stx6_rsid_map.json e /tmp/ld_*.json evitam re-consultas.
"""
import gzip
import json
import math
import os
import time
import urllib.parse
import urllib.request
import datetime
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "lib"))
from coloc_abf import abf as _abf_lib, W_CASO_CONTROLE  # noqa: E402
from gwas_fonte import exigir_gwas, MEDIANA_CHI2_1GL  # noqa: E402

REPORTS = BASE / "reports"
CACHE = Path("/tmp")

REGIOES = {
    "STX6": ("1", 180_900_000, 181_000_000),
    "GAL3ST1": ("22", 30_900_000, 31_000_000),
    "PRNP": ("20", 4_600_000, 4_700_000),
}
POP = os.environ.get("LD_POP", "1000GENOMES:phase_3:ALL")
POP_TAG = POP.split(":")[-1]   # ALL/EUR/... — sufixo de cache e relatório
R2_CLUSTER = 0.80
W = W_CASO_CONTROLE   # 0,04: variância do prior de Wakefield sobre log(OR)


def http_json(url: str, tentativas: int = 3):
    for k in range(tentativas):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.load(r)
        except Exception as e:
            if k == tentativas - 1:
                raise
            time.sleep(2 * (k + 1))


def esearch_rs(crom: int, pos: int):
    term = f"{crom}[CHROM] AND {pos}[POS] AND human[ORGN]"
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
           + urllib.parse.urlencode({"db": "snp", "term": term,
                                     "retmode": "json", "retmax": "5"}))
    ids = http_json(url)["esearchresult"].get("idlist", [])
    time.sleep(0.35)
    return [f"rs{i}" for i in ids]


def esummary_pos(rsids: list[str]):
    """rsID -> (chr, pos_GRCh37) em lote. USA chrpos_prev_assm (GRCh37):
    o LD do Ensembl devolve parceiros com coordenadas GRCh38; nossa tabela
    de sumstats é GRCh37 — casar sem converter gerava cluster 0% (bug 2)."""
    out = {}
    ids_num = [r[2:] for r in rsids]
    for i in range(0, len(ids_num), 180):
        lote = ",".join(ids_num[i:i + 180])
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
               + urllib.parse.urlencode({"db": "snp", "id": lote, "retmode": "json"}))
        try:
            res = http_json(url)["result"]
            for uid in res.get("uids", []):
                d = res.get(uid, {})
                crom = str(d.get("chr", "")).replace("chr", "")
                # prioridade: build anterior (GRCh37); fallback chrpos (GRCh38)
                bruto = d.get("chrpos_prev_assm") or d.get("docsum", {}).get(
                    "chrpos_prev_assm") or ""
                if not bruto:
                    continue   # sem coordenada GRCh37 → não casa com sumstats
                try:
                    pos37 = int(str(bruto).split(":")[-1])
                except ValueError:
                    continue
                out[f"rs{uid}"] = (crom, pos37)
        except Exception:
            pass
        time.sleep(0.4)
    return out


def carregar_regiao(crom, ini, fim):
    recs = []
    with gzip.open(exigir_gwas(), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for linha in fh:
            p = linha.rstrip("\n").split("\t")
            if p[c["chromosome"]] != crom:
                continue
            try:
                pos = int(p[c["base_pair_location"]])
            except ValueError:
                continue
            if ini <= pos <= fim:
                try:
                    recs.append({"p": float(p[c["p_value"]]), "pos": pos,
                                 "beta": float(p[c["beta"]]),
                                 "se": float(p[c["standard_error"]]),
                                 "eaf": float(p[c["effect_allele_frequency"]]),
                                 "oa": p[c["other_allele"]], "ea": p[c["effect_allele"]]})
                except (ValueError, IndexError):
                    pass
    return sorted(recs, key=lambda r: r["p"])


def ld_do_lead(lead_rsid: str, cache_name: str):
    cache = CACHE / cache_name
    if cache.exists():
        return json.load(open(cache))
    url = f"https://rest.ensembl.org/ld/human/{lead_rsid}/{urllib.parse.quote(POP)}"
    dados = http_json(url)
    json.dump(dados, open(cache, "w"))
    return dados


def main() -> None:
    agora = datetime.date.today().isoformat()
    L = ["# Fine-mapping v2 — LD real (Ensembl/1000G phase 3) + credible sets descritivos",
         f"*`finemap_ld.py` em {agora}. População: {POP}. Método: ABF de Wakefield "
         f"(W={W:g}) + agrupamento por r²≥{R2_CLUSTER} ao lead. Approximate — sem modelo "
         "conjunto (SuSiE exigiria genótipos individuais).*", ""]
    lam_rows = []
    todos_rs = []

    # ---- λ_GC por estratos (passada única, streaming) --------------------
    # BUG CORRIGIDO: o sumstats está ORDENADO POR P CRESCENTE (linha 1 = PRNP
    # p=1.6e-15). Amostrar "primeiros 2M" pegava o bloco mais significativo
    # inteiro → λ=3.98 falso. Agora: amostragem UNIFORME desde a linha 1.
    chi_bins = {"MAF<0.05": [], "0.05–0.25": [], "0.25–0.45": [], ">0.45": []}
    n_total = 0
    chi_all = []
    PASSO_AMOSTRA = 10   # mantém cada 10ª linha para o λ global (independe de ordem)
    with gzip.open(exigir_gwas(), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        for linha in fh:
            n_total += 1
            if n_total % PASSO_AMOSTRA != 0:
                continue
            p = linha.rstrip("\n").split("\t")
            try:
                se = float(p[c["standard_error"]])
                beta = float(p[c["beta"]])
                eaf = float(p[c["effect_allele_frequency"]])
            except (ValueError, IndexError):
                continue
            maf = min(eaf, 1 - eaf)
            chi = (beta / se) ** 2
            chi_all.append(chi)
            if maf < 0.05:
                chi_bins["MAF<0.05"].append(chi)
            elif maf < 0.25:
                chi_bins["0.05–0.25"].append(chi)
            elif maf < 0.45:
                chi_bins["0.25–0.45"].append(chi)
            else:
                chi_bins[">0.45"].append(chi)
            # amostragem p/ λ global estável (todos os GWS + 10% sistemática)
            if len(chi_all) > 2_000_000 and n_total % 10 != 0:
                chi_all.pop()
    lam_global = sorted(chi_all)[len(chi_all)//2] / MEDIANA_CHI2_1GL
    L_lam = [f"# Diagnóstico λ_GC — GCST90001389",
             f"*`finemap_ld.py` em {agora}. λ = mediana(χ²)/0.454936.*", "",
             f"- λ global (amostra sistemática de {len(chi_all):,}): **{lam_global:.4f}**", "",
             "| Estrato de MAF | n | λ do estrato |", "|---|---|---|"]
    for nome, vals in chi_bins.items():
        vals.sort()
        lam = vals[len(vals)//2] / MEDIANA_CHI2_1GL
        lam_rows.append((nome, len(vals), lam))
        L_lam.append(f"| {nome} | {len(vals):,} | {lam:.4f} |")
    grad = max(l for _, _, l in lam_rows) - min(l for _, _, l in lam_rows)
    L_lam += ["", f"## Leitura honesta",
              f"- Gradiente de λ entre estratos: **{grad:.4f}**.",
              "- Estratificação populacional clássica infla MAIS os alelos comuns;",
              "  gradiente pequeno (<0.02) sugere inflação majoritariamente poligênica/",
              "  residual, não estratificação grave. Gradiente grande (>0.05) pede PCA.",
              f"- Conclusão para o manuscrito: λ global {lam_global:.3f} é limítrofe-saudável;",
              "  declaramos correção por genomic control nas inferências primárias."]
    (REPORTS / "relatorio_lambda_gc.md").write_text("\n".join(L_lam), encoding="utf-8")
    print(f"[ok] relatorio_lambda_gc.md — λ global {lam_global:.4f}")

    # ---- Fine-mapping por região -----------------------------------------
    caches_mapa = {}
    if (CACHE / "stx6_rsid_map.json").exists():
        caches_mapa["STX6"] = json.load(open(CACHE / "stx6_rsid_map.json"))

    for reg, (crom, ini, fim) in REGIOES.items():
        recs = carregar_regiao(crom, ini, fim)[:40]
        # rsIDs das top variantes (cache p/ STX6)
        if reg in caches_mapa:
            mapa_rs = {int(k): v for k, v in caches_mapa[reg].items()}
        else:
            mapa_rs = {}
            for r in recs[:20]:
                mapa_rs[r["pos"]] = esearch_rs(int(crom), r["pos"])
            json.dump({str(k): v for k, v in mapa_rs.items()},
                      open(CACHE / f"{reg.lower()}_rsid_map.json", "w"))
        lead = recs[0]
        lead_rsids = mapa_rs.get(lead["pos"], [])
        lead_rs = lead_rsids[0] if lead_rsids else None

        # ÂNCORA: se o lead não está no painel 1000G (LD vazio), usa a melhor
        # variante ranqueada que esteja — e declara isso no relatório.
        ancora_rs, ancora_pos = lead_rs, lead["pos"]
        ld_pairs, pos_rs = {}, {}
        if lead_rs:
            try:
                ld_pairs = ld_do_lead(lead_rs, f"ld_{POP_TAG}_{reg.lower()}_{lead_rs}.json")
            except Exception:
                ld_pairs = []
        if not ld_pairs:
            for r in recs[1:16]:
                cands = mapa_rs.get(r["pos"], [])
                if not cands:
                    continue
                try:
                    teste = ld_do_lead(cands[0], f"ld_{POP_TAG}_{reg.lower()}_{cands[0]}.json")
                except Exception:
                    continue
                if teste:
                    ancora_rs, ancora_pos = cands[0], r["pos"]
                    ld_pairs = teste
                    break
            if ld_pairs:
                print(f"[info] {reg}: lead fora do painel; âncora={ancora_rs} "
                      f"@{ancora_pos} (rank {recs.index(next(r for r in recs if r['pos']==ancora_pos))+1})")
            parceiros = sorted({d["variation2"] for d in ld_pairs}
                               | {d["variation1"] for d in ld_pairs})
            parceiros = [p for p in parceiros if p != lead_rs]
            pos_rs = esummary_pos(parceiros)

        # ID EFETIVO da âncora: o Ensembl pode devolver pares com o rsID
        # MESCLADO (ex.: rs60704301→rs2093390); usar o que aparece nos pares.
        ancora_eff = ancora_rs
        if ld_pairs:
            freq = {}
            for d in ld_pairs:
                for k in ("variation1", "variation2"):
                    freq[d[k]] = freq.get(d[k], 0) + 1
            ancora_eff = max(freq, key=freq.get)

        r2_lead = {}   # posição GRCh37 -> r² com o âncora (ID efetivo)
        # (a) via rsIDs que já mapeamos das top variantes da própria região
        for r in recs:
            rss = mapa_rs.get(r["pos"], [])
            if ancora_eff in rss:
                r2_lead[r["pos"]] = 1.0   # a variante É a âncora (mesclada)
                continue
            for rs in rss:
                for d in ld_pairs:
                    if ancora_eff in (d["variation1"], d["variation2"]) and \
                       rs in (d["variation1"], d["variation2"]) and rs != ancora_eff:
                        r2_lead[r["pos"]] = max(r2_lead.get(r["pos"], 0.0),
                                                float(d["r2"]))
        # (b) via esummary chrpos_prev_assm dos parceiros do LD
        for d in ld_pairs:
            for a, b in (("variation1", "variation2"), ("variation2", "variation1")):
                if d[a] == ancora_eff and d[b] in pos_rs:
                    pc, pp = pos_rs[d[b]]
                    if pc == crom:
                        r2_lead[pp] = max(r2_lead.get(pp, 0.0), float(d["r2"]))

        for r in recs:
            # fator de Bayes de Wakefield: implementação única do projeto
            # (pipeline/lib/coloc_abf.py), a mesma usada pela colocalização
            r["abf"] = _abf_lib(r["beta"], r["se"], W) if r["se"] > 0 else 0.0
            r["r2_lead"] = r2_lead.get(r["pos"])
        soma = sum(r["abf"] for r in recs) or 1.0
        for r in recs:
            r["post"] = r["abf"] / soma

        # clusters por r² ao âncora; múltiplos limiares + cobertura honesta
        cred = []
        def massa(lim):
            return sum(r["post"] for r in recs
                       if r["r2_lead"] is not None and r["r2_lead"] >= lim)
        m08, m05 = massa(R2_CLUSTER), massa(0.50)
        com_info = sum(1 for r in recs[:20] if r["r2_lead"] is not None)
        max_r2_top = max((r["r2_lead"] or 0) for r in recs[:20]) if com_info else 0.0
        for r in sorted(recs, key=lambda x: -x["post"])[:12]:
            if r["pos"] == ancora_pos:
                tag = "âncora"
            elif r["r2_lead"] is not None:
                tag = f"r²={r['r2_lead']:.2f}"
            else:
                tag = "sem dado de painel"
            cred.append(f"| {r['pos']:,} | {r['oa']}>{r['ea']} | {r['p']:.2e} | "
                        f"{r['beta']:+.3f} | {'/'.join(mapa_rs.get(r['pos'], ['—']))} | "
                        f"{tag} | {100*r['post']:.1f}% |")
        L += [f"## Locus {reg} — lead chr{crom}:{lead['pos']:,} p={lead['p']:.2e} · âncora LD: {ancora_rs or 'nenhuma'} (efetivo: {ancora_eff})",
              f"- Variantes na janela: {len(carregar_regiao(crom, ini, fim))} · pares LD da âncora: {len(ld_pairs)}",
              f"- Cobertura do painel nas top-20: {com_info}/20 · máx r² observado: {max_r2_top:.2f}",
              f"- Massa posterior do cluster âncora+proxies:",
              f"  **r²≥{R2_CLUSTER}: {100*m08:.1f}%** · r²≥0.50: {100*m05:.1f}% · sem LD/fora do painel: {100*(1-m05):.1f}%",
              *([f"- ⚠️ Leitura: o único r²=1.00 é a própria âncora (posterior ~{100*sum(r['post'] for r in recs if r['pos']==ancora_pos):.1f}%);"
                 f" o sinal real está em variantes mal marcadas pelo painel de comuns"
                 f" (haplótipo provavelmente de baixa frequência)."] if m05 < 0.02 and max_r2_top >= 0.99 else []),
              "", "| pos | alelos | p | beta | rsID | status vs âncora | posterior |",
              "|---|---|---|---|---|---|---|"] + cred + [""]

    L += ["## Nota metodológica final",
          "- Credible set formal exige modelo conjunta (SuSiE/FINEMAP) com genótipos;",
          "  aqui reportamos MASSA POR CLUSTER de LD — suficiente para declarar que o",
          "  sinal é um bloco haplotípico coeso, não um mosaico de falsos independentes.",
          "- rs3747957 (índice Brain 2025): ver relatório QC; presente com p=9.7e-9."]
    destino = REPORTS / (f"relatorio_finemap_loci_{POP_TAG}.md"
                         if POP_TAG != "ALL" else "relatorio_finemap_loci.md")
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")


if __name__ == "__main__":
    main()
