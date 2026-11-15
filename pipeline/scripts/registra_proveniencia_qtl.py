#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
registra_proveniencia_qtl.py — Confere (e, com --escrever, regrava) a
proveniência dos caches de eQTL/sQTL em `pipeline/data/proveniencia_qtl.tsv`.

    python3 pipeline/scripts/registra_proveniencia_qtl.py              # confere SHA-256 e nº de linhas
    python3 pipeline/scripts/registra_proveniencia_qtl.py --remoto     # + compara cada cache com o arquivo
                                                                       #   atual do eQTL Catalogue (rede, pysam)
    python3 pipeline/scripts/registra_proveniencia_qtl.py --escrever   # regrava SHA-256/linhas no TSV

Os tamanhos amostrais do TSV vêm de `data_tables/dataset_metadata_r8.tsv` do
repositório eQTL-Catalogue/eQTL-Catalogue-resources (commit registrado no
cabeçalho do TSV). Este script não os altera.
"""
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import fontes_qtl as F  # noqa: E402


def main() -> int:
    escrever = "--escrever" in sys.argv
    remoto = "--remoto" in sys.argv
    meta, regs = F.carregar_proveniencia()
    print(f"release: {meta.get('release')} · metadata: "
          f"{meta.get('metadata_arquivo')} @ {meta.get('metadata_commit')}")
    falhas = 0
    if remoto:                                   # os .tbi vão para um temporário
        os.chdir(tempfile.mkdtemp(prefix="tbi_"))
    novos = {}
    for r in regs:
        cache = F.DATA / F.nome_cache(r)
        if not cache.exists():
            estado = "sem cache (0 linhas registradas)" if r["linhas"] == 0 \
                else "CACHE AUSENTE"
            falhas += r["linhas"] != 0
            sha, n = "-", 0
        else:
            sha = F.sha256_arquivo(cache)
            n = len(cache.read_text().splitlines())
            ok = (sha == r["sha256"] and n == r["linhas"])
            estado = "ok" if ok else "DIVERGE do registrado"
            falhas += (not ok) and (not escrever)
        novos[r["dataset_id"]] = (sha, n)
        linha = f"{r['dataset_id']} {r['rotulo']:<20} n={r['n_amostras']:<4} " \
                f"linhas={n:<5} sha256={sha[:12]}… {estado}"
        if remoto:
            try:
                texto = F._consultar_remoto(r)
                local = cache.read_text() if cache.exists() else ""
                igual = texto == local
                linha += " · remoto: " + ("idêntico" if igual else "DIFERENTE")
                falhas += not igual
            except Exception as e:                       # rede/pysam
                linha += f" · remoto: indisponível ({type(e).__name__})"
        print(linha)
    if escrever:
        linhas = F.PROVENIENCIA.read_text(encoding="utf-8").splitlines()
        cab = None
        for i, l in enumerate(linhas):
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            if cab is None:
                cab = f
                continue
            d = dict(zip(cab, f))
            sha, n = novos[d["dataset_id"]]
            d["sha256"], d["linhas"] = sha, str(n)
            linhas[i] = "\t".join(d[c] for c in cab)
        F.PROVENIENCIA.write_text("\n".join(linhas) + "\n", encoding="utf-8")
        print(f"[ok] regravado: {F.PROVENIENCIA}")
    return 1 if falhas else 0


if __name__ == "__main__":
    sys.exit(main())
