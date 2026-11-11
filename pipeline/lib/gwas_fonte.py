# -*- coding: utf-8 -*-
"""
gwas_fonte.py — Fonte única do caminho, da URL e do MD5 do GWAS de sCJD
(GCST90001389, Jones et al. 2020, Lancet Neurol; GRCh37).

O arquivo não é versionado (~188 MB). Scripts que o leem devem chamar
`exigir_gwas()` antes de abrir: se faltar, a mensagem diz como obtê-lo.
"""
import hashlib
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
GWAS_GZ = DATA / "GCST90001389_buildGRCh37.tsv.gz"
_DIR_FTP = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            "GCST90001001-GCST90002000/GCST90001389/")
GWAS_URL = _DIR_FTP + "GCST90001389_buildGRCh37.tsv.gz"
GWAS_MD5_URL = _DIR_FTP + "md5sum.txt"
# Publicado em md5sum.txt no diretório do estudo (conferido em 2026-11-10).
GWAS_MD5 = "071790d80ccae0b41adabc7e0eefbf53"


def md5_arquivo(caminho, bloco=1 << 20):
    """MD5 em fluxo (não carrega o arquivo em memória)."""
    h = hashlib.md5()
    with open(caminho, "rb") as fh:
        while True:
            b = fh.read(bloco)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def exigir_gwas():
    """Devolve o caminho do GWAS ou encerra com instrução clara de obtenção."""
    if not GWAS_GZ.exists():
        raise SystemExit(
            f"\n[ERRO] Arquivo do GWAS ausente: {GWAS_GZ}\n"
            "Ele não é versionado (~188 MB). Para obtê-lo com conferência de "
            "MD5, rode:\n\n"
            "    python3 pipeline/scripts/obter_gwas.py\n\n"
            "(ou `make dados-gwas`). Fonte: GWAS Catalog, estudo GCST90001389.\n")
    return GWAS_GZ


def carregar_regiao(crom, ini, fim):
    """Variantes do GWAS em chr `crom`, [ini, fim] (GRCh37), lidas em fluxo.
    Devolve lista de dicts: pos, outro, efeito, beta, se, eaf, p.
    Nenhuma variante é sobrescrita: sítios com mais de um registro aparecem
    todos (a indexação por alelos fica com pipeline/lib/harmoniza.py)."""
    import gzip
    crom = str(crom)
    prefixo = crom + "\t"
    out = []
    with gzip.open(exigir_gwas(), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        c = {n: i for i, n in enumerate(hdr)}
        if hdr[0] != "chromosome":
            raise SystemExit(f"[ERRO] cabeçalho inesperado no GWAS: {hdr[:3]}")
        for linha in fh:
            if not linha.startswith(prefixo):
                continue
            f = linha.rstrip("\n").split("\t")
            pos = int(f[c["base_pair_location"]])
            if not (ini <= pos <= fim):
                continue
            out.append({"pos": pos,
                        "outro": f[c["other_allele"]],
                        "efeito": f[c["effect_allele"]],
                        "beta": float(f[c["beta"]]),
                        "se": float(f[c["standard_error"]]),
                        "eaf": float(f[c["effect_allele_frequency"]]),
                        "p": float(f[c["p_value"]])})
    return out
