# -*- coding: utf-8 -*-
"""
harmoniza.py — Harmonização de variantes por identidade COMPLETA
(posição, REF, ALT), meta-análise por variância inversa e casamento com o GWAS.
Só biblioteca padrão.

Por que existe (auditoria 2026-11, achado G6): a versão anterior indexava as
variantes só pela posição (`d[pos38] = ...`), de modo que (a) a última linha de
cada posição vencia em silêncio e (b) a meta-análise somava efeitos de alelos
diferentes em sítios multialélicos. Aqui:

- a chave é (posição, REF, ALT), com alelos em maiúsculas e indels aparados
  (remoção de sufixo e prefixo comuns; sem genoma de referência não é possível
  o alinhamento à esquerda em repetições, o que fica declarado);
- linhas repetidas da mesma variante (uma por rsID) são deduplicadas quando
  beta e erro-padrão coincidem; se divergirem, a variante é descartada e contada;
- a meta-análise só combina alelos idênticos;
- o casamento com o GWAS confere os DOIS alelos e inverte o sinal do GWAS
  quando o alelo de efeito do GWAS é o REF do QTL;
- todo descarte é contado, com o motivo.

Convenção: o alelo de efeito das variantes de QTL é sempre ALT (eQTL Catalogue).
"""
import math
from collections import Counter

__all__ = ["normalizar_alelos", "deduplicar", "meta_ivw", "indexar_gwas",
           "casar_com_gwas", "p_de_z", "complementar", "palindromico",
           "chi2_sf", "MIN_DATASETS"]

MIN_DATASETS = 3
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def p_de_z(z):
    """p bicaudal da normal padrão."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def complementar(alelo):
    """Complemento reverso; None se houver base fora de ACGT."""
    try:
        return "".join(_COMP[b] for b in reversed(alelo))
    except KeyError:
        return None


def palindromico(ref, alt):
    """SNV cuja fita é ambígua (A/T ou C/G)."""
    return len(ref) == 1 and len(alt) == 1 and _COMP.get(ref) == alt


def normalizar_alelos(pos, ref, alt):
    """Representação mínima de (pos, REF, ALT): maiúsculas; apara sufixo comum e
    depois prefixo comum (ajustando a posição), mantendo ≥ 1 base em cada alelo.
    Ex.: (100, 'CTT', 'CT') → (100, 'CT', 'C'); (100, 'AT', 'AC') → (101, 'T', 'C')."""
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


def deduplicar(linhas):
    """linhas: iterável de (pos, ref, alt, rsid, p, beta, se, maf).

    Devolve (variantes, stats):
      variantes = {(pos, ref, alt): {"beta", "se", "maf", "p", "rsids"}}
      stats     = Counter com linhas, variantes, linhas_repetidas_por_rsid,
                  variantes_com_linha_repetida, descartadas_conflito,
                  descartadas_se_invalido, posicoes, posicoes_multialelicas, indels
    """
    stats = Counter()
    bruto = {}
    conflito = set()
    for pos, ref, alt, rsid, p, beta, se, maf in linhas:
        stats["linhas"] += 1
        if not (se > 0) or not math.isfinite(beta):
            stats["descartadas_se_invalido"] += 1
            continue
        chave = normalizar_alelos(pos, ref, alt)
        ant = bruto.get(chave)
        if ant is None:
            bruto[chave] = {"beta": beta, "se": se, "maf": maf, "p": p,
                            "rsids": [rsid], "linhas": 1}
            continue
        ant["linhas"] += 1
        stats["linhas_repetidas_por_rsid"] += 1
        if ant["beta"] != beta or ant["se"] != se:
            conflito.add(chave)
        elif rsid not in ant["rsids"]:
            ant["rsids"].append(rsid)
    for chave in conflito:
        del bruto[chave]
    stats["descartadas_conflito"] = len(conflito)
    stats["variantes"] = len(bruto)
    stats["variantes_com_linha_repetida"] = sum(
        1 for v in bruto.values() if v["linhas"] > 1)
    por_pos = Counter(k[0] for k in bruto)
    stats["posicoes"] = len(por_pos)
    stats["posicoes_multialelicas"] = sum(1 for n in por_pos.values() if n > 1)
    stats["indels"] = sum(1 for k in bruto if len(k[1]) != 1 or len(k[2]) != 1)
    return bruto, stats


def chi2_sf(x, gl):
    """P(χ²_gl > x) para gl inteiro ≥ 1 (forma fechada; stdlib)."""
    if gl < 1:
        return float("nan")
    if x <= 0:
        return 1.0
    if gl % 2 == 0:
        termo = soma = math.exp(-x / 2.0)
        for k in range(1, gl // 2):
            termo *= (x / 2.0) / k
            soma += termo
        return min(1.0, soma)
    soma = math.erfc(math.sqrt(x / 2.0))
    termo = math.sqrt(2.0 * x / math.pi) * math.exp(-x / 2.0)
    for k in range(1, (gl - 1) // 2 + 1):
        soma += termo
        termo *= x / (2 * k + 1)
    return min(1.0, soma)


def meta_ivw(por_dataset, min_datasets=MIN_DATASETS, aleatorio=False):
    """Meta-análise por variância inversa, variante a variante, combinando
    apenas alelos idênticos.

    por_dataset: {id: {(pos, ref, alt): {"beta", "se", ...}}}
    Devolve (meta, stats):
      meta  = {(pos, ref, alt): {"beta", "se", "k", "datasets", "Q", "p_Q",
                                 "I2", "tau2"}}
      stats = Counter: variantes_em_algum_dataset, abaixo_do_minimo, na_meta,
              posicoes_na_meta, posicoes_com_alelos_discordantes
    aleatorio=True usa DerSimonian–Laird (sensibilidade); o padrão é efeito fixo.
    """
    stats = Counter()
    todas = set()
    for d in por_dataset.values():
        todas |= set(d)
    stats["variantes_em_algum_dataset"] = len(todas)
    # posições em que datasets diferentes trazem pares REF/ALT diferentes: é
    # onde a indexação só por posição misturava alelos.
    por_pos = {}
    for pos, ref, alt in todas:
        por_pos.setdefault(pos, set()).add((ref, alt))
    stats["posicoes_com_alelos_discordantes"] = sum(
        1 for s in por_pos.values() if len(s) > 1)
    meta = {}
    for chave in todas:
        ids, betas, ses = [], [], []
        for nome, d in por_dataset.items():
            v = d.get(chave)
            if v is None:
                continue
            ids.append(nome)
            betas.append(v["beta"])
            ses.append(v["se"])
        k = len(betas)
        if k < min_datasets:
            stats["abaixo_do_minimo"] += 1
            continue
        w = [1.0 / s ** 2 for s in ses]
        sw = sum(w)
        beta_f = sum(b * wi for b, wi in zip(betas, w)) / sw
        q = sum(wi * (b - beta_f) ** 2 for b, wi in zip(betas, w))
        c = sw - sum(wi * wi for wi in w) / sw
        tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
        i2 = max(0.0, (q - (k - 1)) / q) if q > 0 else 0.0
        if aleatorio and tau2 > 0:
            w = [1.0 / (s ** 2 + tau2) for s in ses]
            sw = sum(w)
            beta_m = sum(b * wi for b, wi in zip(betas, w)) / sw
        else:
            beta_m = beta_f
        meta[chave] = {"beta": beta_m, "se": math.sqrt(1.0 / sw), "k": k,
                       "datasets": ids, "Q": q, "p_Q": chi2_sf(q, k - 1),
                       "I2": i2, "tau2": tau2}
    stats["na_meta"] = len(meta)
    stats["posicoes_na_meta"] = len({k[0] for k in meta})
    return meta, stats


def indexar_gwas(registros):
    """registros: iterável de dicts com pos, outro, efeito, beta, se, eaf, p
    (posição no MESMO build das chaves depois de aplicado o offset).
    Devolve {pos: [registro, ...]} — nunca descarta um segundo alelo."""
    idx = {}
    for r in registros:
        r = dict(r)
        r["outro"], r["efeito"] = r["outro"].upper(), r["efeito"].upper()
        idx.setdefault(r["pos"], []).append(r)
    return idx


def casar_com_gwas(variantes, gwas_idx, offset=0):
    """Casa variantes de QTL (efeito = ALT) com o GWAS conferindo os dois alelos.

    variantes: {(pos, ref, alt): {"beta", "se", ...}}; a posição do GWAS é
    pos + offset. Devolve (pares, stats). Cada par traz o beta do GWAS já
    orientado para o ALT do QTL (sinal invertido quando o alelo de efeito do
    GWAS é o REF).

    stats (Counter): casadas, casadas_mesmo_alelo_de_efeito,
    casadas_sinal_invertido, sem_variante_no_gwas, alelos_incompativeis,
    so_na_fita_complementar, palindromicas_casadas.
    Variantes que só casariam na fita complementar NÃO são usadas: são contadas
    para diagnóstico (se o número for > 0, a fita do GWAS precisa de revisão).
    """
    stats = Counter()
    pares = []
    for (pos, ref, alt), v in sorted(variantes.items()):
        candidatos = gwas_idx.get(pos + offset)
        if not candidatos:
            stats["sem_variante_no_gwas"] += 1
            continue
        achou = None
        for g in candidatos:
            if g["efeito"] == alt and g["outro"] == ref:
                achou = (g, 1.0)
                break
            if g["efeito"] == ref and g["outro"] == alt:
                achou = (g, -1.0)
                break
        if achou is None:
            cr, ca = complementar(ref), complementar(alt)
            if any({g["efeito"], g["outro"]} == {cr, ca} for g in candidatos):
                stats["so_na_fita_complementar"] += 1
            else:
                stats["alelos_incompativeis"] += 1
            continue
        g, sinal = achou
        stats["casadas"] += 1
        stats["casadas_mesmo_alelo_de_efeito" if sinal > 0
              else "casadas_sinal_invertido"] += 1
        if palindromico(ref, alt):
            stats["palindromicas_casadas"] += 1
        eaf = g.get("eaf")
        pares.append({
            "pos": pos + offset, "pos_qtl": pos, "ref": ref, "alt": alt,
            "beta_gwas": sinal * g["beta"], "se_gwas": g["se"],
            "p_gwas": g["p"],
            "freq_alt_gwas": (eaf if sinal > 0 else 1.0 - eaf)
            if eaf is not None else None,
            "qtl": v,
        })
    return pares, stats
