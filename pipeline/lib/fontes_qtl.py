# -*- coding: utf-8 -*-
"""
fontes_qtl.py — Acesso único aos caches de eQTL/sQTL de STX6 (eQTL Catalogue)
e à sua proveniência versionada (`pipeline/data/proveniencia_qtl.tsv`).

Os tamanhos amostrais NÃO são constantes digitadas nos scripts: vêm do arquivo
de proveniência, que registra release, metadata oficial (commit), janela, URL
de origem e SHA-256 de cada cache (auditoria 2026-11, achados G8 e R4).

Constantes da região (compartilhadas por todos os scripts de colocalização):
    REGIAO_B37      janela de análise em GRCh37
    OFFSET_B38_B37  b37 = b38 + offset (constante verificada para STX6; ver
                    relatorio_validacao_offset_genomico.md)
    JANELA_B38      janela consultada no tabix remoto (GRCh38) quando o cache
                    foi criado; cobre a REGIAO_B37 com margem
"""
import csv
import hashlib
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
PROVENIENCIA = DATA / "proveniencia_qtl.tsv"

STX6 = "ENSG00000135823"
REGIAO_B37 = (180_900_000, 181_100_000)
OFFSET_B38_B37 = -30_864
JANELA_B38 = (180_930_000, 181_132_000)
URL_BASE = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats"

# a janela consultada tem de conter a região de análise convertida para b38
assert JANELA_B38[0] <= REGIAO_B37[0] - OFFSET_B38_B37
assert JANELA_B38[1] >= REGIAO_B37[1] - OFFSET_B38_B37

__all__ = ["STX6", "REGIAO_B37", "OFFSET_B38_B37", "JANELA_B38", "PROVENIENCIA",
           "carregar_proveniencia", "datasets", "ler_eqtl", "ler_sqtl",
           "sha256_arquivo", "url_origem", "nome_cache"]


def sha256_arquivo(caminho):
    h = hashlib.sha256()
    with open(caminho, "rb") as fh:
        for bloco in iter(lambda: fh.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def carregar_proveniencia():
    """(metadados, registros): metadados = {chave: valor} das linhas '# chave: valor';
    registros = lista de dicts, uma por dataset, com n_amostras já inteiro."""
    if not PROVENIENCIA.exists():
        raise SystemExit(f"[ERRO] proveniência ausente: {PROVENIENCIA}")
    meta, corpo = {}, []
    with PROVENIENCIA.open(encoding="utf-8") as fh:
        for linha in fh:
            if linha.startswith("#"):
                if ":" in linha:
                    k, v = linha[1:].split(":", 1)
                    meta[k.strip()] = v.strip()
                continue
            corpo.append(linha)
    registros = list(csv.DictReader(corpo, delimiter="\t"))
    for r in registros:
        r["n_amostras"] = int(r["n_amostras"])
        r["linhas"] = int(r["linhas"])
    return meta, registros


def datasets(quant_method):
    """Registros de proveniência de um método ('ge' = expressão; 'leafcutter' =
    splicing), na ordem do arquivo."""
    _, regs = carregar_proveniencia()
    sel = [r for r in regs if r["quant_method"] == quant_method]
    if not sel:
        raise ValueError(f"nenhum dataset com quant_method={quant_method!r}")
    return sel


def nome_cache(reg):
    prefixo = "eqtl" if reg["quant_method"] == "ge" else "sqtl"
    return f"{prefixo}_{reg['dataset_id']}_stx6.tsv"


def url_origem(reg):
    sufixo = "all" if reg["quant_method"] == "ge" else "cc"
    return (f"{URL_BASE}/{reg['study_id']}/{reg['dataset_id']}/"
            f"{reg['dataset_id']}.{sufixo}.tsv.gz")


def _consultar_remoto(reg):
    """Linhas do cache a partir do tabix remoto (mesma janela e mesmo formato
    do cache original, para que o SHA-256 da proveniência continue conferindo).
    Exige pysam e rede; o índice .tbi é baixado para o diretório corrente."""
    import pysam
    tb = pysam.TabixFile(url_origem(reg))
    linhas = []
    for linha in tb.fetch("1", *JANELA_B38):
        f = linha.split("\t")
        if reg["quant_method"] == "ge":
            if f[0] != STX6:
                continue
            linhas.append("\t".join(map(str, (
                int(f[2]), f[3], f[4], f[18], float(f[8]), float(f[9]),
                float(f[10]), float(f[7])))))
        else:
            if f[1] != "1":
                continue
            linhas.append("\t".join([f[2], f[3], f[4], f[18], f[8], f[9],
                                     f[10], f[0], f[16]]))
    return "\n".join(linhas)


def _texto_cache(reg):
    cache = DATA / nome_cache(reg)
    if cache.exists():
        return cache.read_text()
    if reg["linhas"] == 0:
        return ""                       # dataset sem pares na janela (registrado)
    texto = _consultar_remoto(reg)
    if texto:
        cache.write_text(texto)
    return texto


def ler_eqtl(reg):
    """eQTL de STX6: [(pos_b38, ref, alt, rsid, p, beta, se, maf)].
    Efeito = alelo ALT (convenção do eQTL Catalogue)."""
    out = []
    for linha in _texto_cache(reg).splitlines():
        f = linha.split("\t")
        out.append((int(f[0]), f[1], f[2], f[3], float(f[4]), float(f[5]),
                    float(f[6]), float(f[7])))
    return out


def ler_sqtl(reg):
    """sQTL (leafcutter, arquivo cc = só pares cis SIGNIFICATIVOS):
    [(pos_b38, ref, alt, rsid, p, beta, se, molecular_trait_id, gene_id)]."""
    out = []
    for linha in _texto_cache(reg).splitlines():
        f = linha.split("\t")
        out.append((int(f[0]), f[1], f[2], f[3], float(f[4]), float(f[5]),
                    float(f[6]), f[7], f[8]))
    return out
