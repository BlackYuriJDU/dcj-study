#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
obter_gwas.py — Obtém as estatísticas-resumo do GWAS de sCJD (GCST90001389,
Jones et al. 2020) do FTP do GWAS Catalog e confere o MD5.

O arquivo (~188 MB) não é versionado. Todos os scripts que dependem dele
chamam `pipeline/lib/gwas_fonte.exigir_gwas()`, que aponta para este script
quando o arquivo falta.

Uso:
    python3 pipeline/scripts/obter_gwas.py            # baixa se faltar; confere MD5
    python3 pipeline/scripts/obter_gwas.py --conferir # só confere o MD5 local

O MD5 de referência está fixado em `pipeline/lib/gwas_fonte.py` e é o mesmo
publicado em `md5sum.txt` no diretório do estudo no GWAS Catalog. O script
também baixa esse `md5sum.txt` e falha se o valor remoto divergir do fixado
(arquivo trocado na origem = investigar antes de usar).
"""
import os
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from gwas_fonte import (GWAS_GZ, GWAS_URL, GWAS_MD5, GWAS_MD5_URL,  # noqa: E402
                        md5_arquivo)


def md5_remoto():
    """MD5 publicado pelo GWAS Catalog para o arquivo (None se indisponível)."""
    try:
        with urllib.request.urlopen(GWAS_MD5_URL, timeout=60) as r:
            for linha in r.read().decode().splitlines():
                campos = linha.split()
                if len(campos) == 2 and campos[1] == GWAS_GZ.name:
                    return campos[0]
    except OSError as e:
        print(f"[aviso] md5sum.txt remoto indisponível ({e}); "
              "uso apenas o MD5 fixado no repositório.")
    return None


def baixar():
    parcial = GWAS_GZ.with_suffix(GWAS_GZ.suffix + ".part")
    GWAS_GZ.parent.mkdir(parents=True, exist_ok=True)
    print(f"[1/3] baixando {GWAS_URL}")
    with urllib.request.urlopen(GWAS_URL, timeout=120) as r, \
            parcial.open("wb") as fh:
        total = int(r.headers.get("Content-Length") or 0)
        lido = 0
        proximo = 0
        while True:
            bloco = r.read(1 << 20)
            if not bloco:
                break
            fh.write(bloco)
            lido += len(bloco)
            if lido >= proximo:
                pct = f" ({lido / total:.0%})" if total else ""
                print(f"      {lido / 1e6:7.1f} MB{pct}", flush=True)
                proximo += 20 << 20
    return parcial


def main() -> int:
    so_conferir = "--conferir" in sys.argv
    remoto = md5_remoto()
    if remoto and remoto != GWAS_MD5:
        print(f"[ERRO] o GWAS Catalog publica MD5 {remoto}, mas o repositório "
              f"fixa {GWAS_MD5}. O arquivo mudou na origem; não prossiga sem "
              "investigar.")
        return 2
    if GWAS_GZ.exists():
        local = md5_arquivo(GWAS_GZ)
        if local == GWAS_MD5:
            print(f"[ok] {GWAS_GZ.name} presente; MD5 {local} confere"
                  + (" (também com o md5sum.txt remoto)." if remoto else "."))
            return 0
        print(f"[ERRO] {GWAS_GZ} existe mas o MD5 é {local} "
              f"(esperado {GWAS_MD5}). Apague-o e rode de novo.")
        return 2
    if so_conferir:
        print(f"[ERRO] {GWAS_GZ} ausente.")
        return 2
    parcial = baixar()
    print("[2/3] conferindo MD5...")
    local = md5_arquivo(parcial)
    if local != GWAS_MD5:
        print(f"[ERRO] MD5 do download é {local} (esperado {GWAS_MD5}); "
              f"arquivo parcial mantido em {parcial} para inspeção.")
        return 2
    os.replace(parcial, GWAS_GZ)
    print(f"[3/3] {GWAS_GZ} ({GWAS_GZ.stat().st_size / 1e6:.1f} MB); "
          f"MD5 {local} confere.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
