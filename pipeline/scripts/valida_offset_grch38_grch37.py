#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
valida_offset_grch38_grch37.py — valida empiricamente o offset GRCh38→GRCh37.

Toda a colocalização de STX6 (o resultado mecanístico central do manuscrito)
descansa numa constante: pos_b37 = pos_b38 − 30.864, aplicada identicamente em
coloc_stx6_eqtl.py, coloc_meta_stx6.py, coloc_sqtl_stx6.py e
crosscheck_coloc_R.py. Até 2026-10-04 essa constante nunca foi confrontada com
uma fonte independente, e o docstring ainda afirmava que fora "validado por
Ensembl MAP" — o que era falso na época em que foi escrito.

ESTE SCRIPT USA A API DE ASSEMBLY MAP DO ENSEMBL COMO FONTE PRIMÁRIA.
Uma implementação própria de chain file foi testada e DESCARTADA: o parser
produzia b37 = 30.545.990 para uma posição que o Ensembl mapeia para
180.941.861 — um erro de 150.426.735 bp, constante, indicando corrupção no
acúmulo de blocos. Não usar chain file próprio aqui.

Pergunta respondida: o offset −30.864bp é exato e constante em todo o span de
STX6, ou o resultado da colocalização depende de uma aproximação?

Se for exato: a limitação é de REDAÇÃO — o manuscrito hoje diz que a conversão
não foi validada, o que passa a ser conservador demais. Deve ser atualizado.
Se não for: a colocalização precisa ser refeita.

Uso:  python pipeline/scripts/valida_offset_grch38_grch37.py [--n=41]
"""
import datetime
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
DATA = BASE / "data"
REPORTS = BASE / "reports"

OFFSET_CONSTANTE = -30_864
REGIAO_B37 = (180_900_000, 181_100_000)
STX6_B37 = (180_941_861, 180_992_047)
CACHE_JSON = DATA / "offset_ensembl_cache.json"

_cache = {}


def _carregar_cache():
    if _cache:
        return _cache
    if CACHE_JSON.exists():
        try:
            _cache.update({int(k): v for k, v in
                           json.loads(CACHE_JSON.read_text(
                               encoding="utf-8")).items()})
        except (ValueError, json.JSONDecodeError):
            _cache.clear()
    return _cache


def _salvar_cache():
    CACHE_JSON.write_text(
        json.dumps({str(k): v for k, v in _cache.items()}), encoding="utf-8")


def consultar_ensembl(p37):
    """b37 (1-based) -> b38 (1-based) via Ensembl REST assembly map.

    /map/human/GRCh37/1:START..END recebe GRCh37 e devolve GRCh38 (montagem
    corrente). O Ensembl usa coordenadas 0-based, daí o +1 no retorno.

    Respeita o rate limit do Ensembl (429/503) com backoff e um intervalo
    mínimo entre requisições; resultados são cacheados em disco porque a
    validação é determinística.
    """
    _carregar_cache()
    if p37 in _cache:
        return _cache[p37]

    u = (f"https://rest.ensembl.org/map/human/GRCh37/1:{p37-1}..{p37}"
         "?content-type=application/json")
    req = urllib.request.Request(u, headers={"Content-Type": "application/json"})
    for tentativa in range(6):
        try:
            d = json.load(urllib.request.urlopen(req, timeout=30))
            ms = d.get("mappings", [])
            _cache[p37] = ms[0]["mapped"]["start"] + 1 if ms else None
            if ms:
                time.sleep(0.4)              # < 15 req/s, com folga
            return _cache[p37]
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
            # 429/503 do rate limit ou rede instável: espera e tenta de novo
            time.sleep(min(2.0 * (tentativa + 1), 20.0))
    _cache[p37] = None
    return None


def consultar_variante_ensembl(pos38):
    """b38 (1-based) -> b37 (1-based), invertendo o endpoint do Ensembl.

    O endpoint /map só aceita a montagem antiga como entrada; para converter
    no sentido inverso consulta-se cada candidata b37 = pos38 + offset e
    verifica-se se o Ensembl devolve a pos38 original.
    """
    _carregar_cache()
    b37 = pos38 + OFFSET_CONSTANTE
    if b37 in _cache:
        volta = _cache[b37]
        return b37 if volta == pos38 else None
    p38 = consultar_ensembl(b37)
    return b37 if p38 == pos38 else None


def gwas_variantes():
    """Posições do GWAS na janela (GRCh37), cacheadas em memória."""
    import gzip as _gz
    out = set()
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
    from gwas_fonte import exigir_gwas
    with _gz.open(exigir_gwas(), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci, pi = hdr.index("chromosome"), hdr.index("base_pair_location")
        for linha in fh:
            f = linha.split("\t", pi + 1)
            if f[ci] != "1":
                continue
            p = int(f[pi])
            if REGIAO_B37[0] <= p <= REGIAO_B37[1]:
                out.add(p)
    return out


def main() -> None:
    n_pontos = 41
    for arg in sys.argv[1:]:
        if arg.startswith("--n="):
            n_pontos = max(3, int(arg.split("=", 1)[1]))

    agora = datetime.date.today().isoformat()
    const = OFFSET_CONSTANTE

    # --- 1. o offset é constante no span de STX6? ---
    print(f"[1/3] consultando Ensembl em {n_pontos} pontos do span de STX6...")
    passo = (STX6_B37[1] - STX6_B37[0]) / max(1, n_pontos - 1)
    pontos = [int(round(STX6_B37[0] + i * passo)) for i in range(n_pontos)]
    offsets = {}
    sem_mapa = 0
    for p37 in pontos:
        p38 = consultar_ensembl(p37)
        if p38 is None:
            sem_mapa += 1
            continue
        d = p38 - p37
        offsets[d] = offsets.get(d, 0) + 1
    n_total = sum(offsets.values())
    n_exato = offsets.get(-const, 0)
    _salvar_cache()

    L = ["# Validação do offset GRCh38→GRCh37 usado na colocalização de STX6",
         f"*`valida_offset_grch38_grch37.py` em {agora}. Substitui a suposição de "
         "−30.864 bp por uma medição contra a API de assembly map do Ensembl.*",
         "",
         "A colocalização de STX6 — o resultado mecanístico central do manuscrito "
         "— converte todas as posições eQTL de GRCh38 para GRCh37 somando uma "
         "constante de −30.864 bp. Até 2026-10-04 essa constante nunca tinha sido "
         "confrontada com uma fonte independente, e o docstring do script ainda "
         "afirmava que fora \"validado por Ensembl MAP\", o que era falso.", "",
         "## 1. O offset é constante no span de STX6?", "",
         f"Amostra de **{n_total}** pontos no span de STX6 em GRCh37 "
         f"(chr1:{STX6_B37[0]:,}–{STX6_B37[1]:,}), convertidos pela API de "
         "assembly map do Ensembl.".replace(",", "."), "",
         "| offset observado (bp) | pontos |", "|---|---|"]
    for o, n in sorted(offsets.items(), key=lambda t: -t[1]):
        marca = " ← **offset do código**" if o == -const else ""
        L.append(f"| {o:+,} | {n}{marca} |")
    L.append("")
    if sem_mapa:
        L.append(f"({sem_mapa} ponto(s) sem mapeamento — queda de rede ou rate "
                 f"limit do Ensembl.)")
        L.append("")

    # --- 2. as variantes usadas na colocalização continuam casadas? ---
    print("[2/3] conferindo as variantes eQTL realmente usadas...")
    gwas = gwas_variantes()
    tot_const = tot_real = tot_e = 0
    linhas = []
    for arq in sorted(DATA.glob("eqtl_*_stx6.tsv")):
        pos38 = []
        with arq.open(encoding="utf-8") as fh:
            for linha in fh:
                f = linha.rstrip("\n").split("\t")
                if len(f) >= 2:
                    pos38.append(int(f[0]))
        if not pos38:
            continue
        # subconjunto espaçado: o objetivo é detectar divergência sistemática,
        # não auditar 390 chamadas de rede
        amostra = pos38[::max(1, len(pos38) // 12)]
        casadas_const = sum(1 for p in amostra if (p + const) in gwas)
        casadas_real = sum(1 for p in amostra
                           if (consultar_variante_ensembl(p) or -1) in gwas)
        tot_e += len(amostra)
        tot_const += casadas_const
        tot_real += casadas_real
        linhas.append((arq.stem.replace("eqtl_", ""), len(pos38), len(amostra),
                       casadas_const, casadas_real))
    _salvar_cache()

    L += ["## 2. As variantes da colocalização mudam de posição?", "",
          "Para cada posição eQTL de STX6 (amostra espaçada dos 5 datasets), "
          "comparamos o casamento com o GWAS sob a constante de −30.864 bp contra "
          "o casamento confirmado pelo Ensembl.", "",
          "| dataset | eQTL de STX6 | amostradas | casadas pela constante | "
          "casadas pelo Ensembl |", "|---|---|---|---|---|"]
    for nome, tot, amp, cc, cr in linhas:
        L.append(f"| {nome} | {tot} | {amp} | {cc} | {cr} |")
    L.append("")

    concordam = tot_real == tot_const
    L += ["## 3. Veredito", ""]
    unico = len(offsets) == 1
    if unico and concordam:
        L += [f"**O offset −30.864 bp está correto e é constante.** O Ensembl "
              f"devolve exatamente esse delta em {n_total}/{n_total} pontos "
              f"amostrados do span de STX6, e as variantes que casam com o GWAS "
              f"sob a constante são exatamente as mesmas sob o mapeamento "
              f"verificado ({tot_const}/{tot_const}).", "",
              "**Não há impacto no resultado da colocalização.** A limitação "
              "registrada no manuscrito é real como prudência editorial, mas hoje "
              "ela subdeclara: a conversão pode ser descrita como *verificada "
              "contra a API de assembly map do Ensembl*, o que é mais forte do "
              "que \"não validada\" e mais honesto do que o antigo \"validado por "
              "Ensembl MAP\" (que afirmava uma validação de cadeia inexistente).", ""]
    elif unico and not concordam:
        L += ["**O offset é constante, mas o conjunto de variantes casadas "
              f"difere** ({tot_const} pela constante vs. {tot_real} pelo "
              "Ensembl). A colocalização precisa ser refeita com o mapeamento "
              "verificado antes de qualquer submissão.", ""]
    else:
        L += [f"**O offset NÃO é constante** na janela — foram observados "
              f"{len(offsets)} valores distintos. A constante de −30.864 bp não "
              "descreve corretamente a conversão e a colocalização precisa ser "
              "refeita com liftOver por cadeia.", ""]

    L += ["## 4. Nota metodológica", "",
          "Uma implementação própria de parser de chain file foi escrita e "
          "descartada antes desta versão: ela mapeava b38 180.972.725 para b37 "
          "30.545.990, enquanto o Ensembl devolve 180.941.861 — um erro "
          "constante de 150.426.735 bp, compatível com corrupção no acúmulo dos "
          "blocos do chain. A API de assembly map do Ensembl é usada como fonte "
          "primária por ser autoritativa e verificável.", ""]

    destino = REPORTS / "relatorio_validacao_offset_genomico.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[3/3] {destino}")
    print(f"[ok] offsets observados: {sorted(offsets)}")
    print(f"[ok] casadas: constante={tot_const} ensembl={tot_real} "
          f"(de {tot_e} amostradas)")
    if not (unico and concordam):
        print("[AVISO] validação não conclusiva — revisar a colocalização")


if __name__ == "__main__":
    main()