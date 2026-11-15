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


TENTATIVAS = 6


def baixar():
    """Baixa para <arquivo>.part, em fluxo. Se a conexão cair (comum num
    arquivo de 188 MB), retoma de onde parou com um cabeçalho Range, até
    TENTATIVAS vezes, com espera crescente. O MD5 final decide se serviu."""
    import time
    parcial = GWAS_GZ.with_suffix(GWAS_GZ.suffix + ".part")
    GWAS_GZ.parent.mkdir(parents=True, exist_ok=True)
    print(f"[1/3] baixando {GWAS_URL}")
    total = None
    for tentativa in range(1, TENTATIVAS + 1):
        ja = parcial.stat().st_size if parcial.exists() else 0
        req = urllib.request.Request(GWAS_URL)
        if ja:
            req.add_header("Range", f"bytes={ja}-")
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                if ja and r.status != 206:      # servidor ignorou o Range
                    ja = 0
                if total is None:
                    resto = int(r.headers.get("Content-Length") or 0)
                    total = ja + resto if resto else 0
                lido, proximo = ja, ja
                with parcial.open("ab" if ja else "wb") as fh:
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
            if not total or parcial.stat().st_size >= total:
                return parcial
            motivo = "conexão encerrada antes do fim"
        except OSError as e:                    # URLError, reset, timeout
            motivo = str(e)
        if tentativa == TENTATIVAS:
            raise SystemExit(f"[ERRO] download falhou após {TENTATIVAS} "
                             f"tentativas ({motivo}). O parcial ficou em "
                             f"{parcial}; rode de novo para retomar.")
        espera = min(2 ** tentativa, 30)
        print(f"      interrompido ({motivo}); retomando em {espera}s "
              f"(tentativa {tentativa + 1}/{TENTATIVAS})", flush=True)
        time.sleep(espera)
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
        parcial.unlink()        # um parcial corrompido não pode ser retomado
        print(f"[ERRO] MD5 do download é {local} (esperado {GWAS_MD5}); "
              "o arquivo baixado foi descartado. Rode de novo.")
        return 2
    os.replace(parcial, GWAS_GZ)
    print(f"[3/3] {GWAS_GZ} ({GWAS_GZ.stat().st_size / 1e6:.1f} MB); "
          f"MD5 {local} confere.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
