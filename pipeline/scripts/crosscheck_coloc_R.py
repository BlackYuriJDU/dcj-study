#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
crosscheck_coloc_R.py — Exporta a região STX6 harmonizada para conferência
contra o pacote R `coloc` e grava o relatório de validação.

1. Reusa a análise de coloc_meta_stx6.py (mesma harmonização, mesma meta) e
   grava `pipeline/data/stx6_crosscheck_input.tsv`: uma linha por variante,
   com as colunas que `coloc.abf` consome (ver coloc_crosscheck_R.R).
2. Compara a implementação única (pipeline/lib/coloc_abf.py) com o ORÁCULO
   versionado em tests/fixtures/: a entrada de 390 variantes da v0.4 e o
   `SNP.PP.H4` que o `coloc.abf` 5.2.3 real devolveu para ela em 2026-08-29.
   Essa parte dispensa R e dispensa o GWAS.
3. Grava `pipeline/reports/relatorio_validacao_coloc_R.md`.

Com `--so-oraculo`, pula o passo 1 (não precisa do arquivo do GWAS) e usa a
entrada já gravada em pipeline/data/.

Se R e o pacote `coloc` estiverem instalados, `coloc_crosscheck_R.R` roda o
`coloc.abf` sobre a entrada atual e grava
`pipeline/data/stx6_crosscheck_colocR_output.tsv`; este script então inclui
essa comparação no relatório.
"""
import csv
import re
import sys
import datetime
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
RAIZ = BASE.parent
sys.path.insert(0, str(BASE / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import coloc_abf as C  # noqa: E402
import harmoniza as H  # noqa: E402

ENTRADA = BASE / "data" / "stx6_crosscheck_input.tsv"
SAIDA_R = BASE / "data" / "stx6_crosscheck_colocR_output.tsv"
FIX_ENT = RAIZ / "tests" / "fixtures" / "stx6_crosscheck_input_v04.tsv"
FIX_R = RAIZ / "tests" / "fixtures" / "stx6_crosscheck_colocR_output_v04.tsv"
REPORTS = BASE / "reports"
REPORTS_SAIDA = REPORTS          # onde o relatório é gravado (os testes trocam)
COLUNAS = ["pos_b37", "ref", "alt", "beta_gwas", "se_gwas", "p_gwas",
           "beta_eqtl_meta", "se_eqtl_meta", "p_eqtl_meta", "n_datasets", "maf"]


def exportar_entrada():
    """Gera a entrada atual a partir do pipeline (precisa do GWAS)."""
    from coloc_meta_stx6 import analisar
    a = analisar()
    with ENTRADA.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUNAS) + "\n")
        for p in a["pares"]:
            q = p["qtl"]
            f = p["freq_alt_gwas"]
            fh.write("\t".join(str(x) for x in (
                p["pos"], p["ref"], p["alt"], p["beta_gwas"], p["se_gwas"],
                p["p_gwas"], q["beta"], q["se"], H.p_de_z(q["beta"] / q["se"]),
                q["k"], min(f, 1.0 - f))) + "\n")
    print(f"[ok] {ENTRADA}: {len(a['pares'])} variantes harmonizadas "
          "(efeito = ALT nos dois traços)")


def ler_entrada(caminho):
    with open(caminho, encoding="utf-8") as fh:
        linhas = list(csv.DictReader(fh, delimiter="\t"))
    pos = [int(r["pos_b37"]) for r in linhas]
    tg = [(float(r["beta_gwas"]), float(r["se_gwas"])) for r in linhas]
    te = [(float(r["beta_eqtl_meta"]), float(r["se_eqtl_meta"])) for r in linhas]
    return pos, tg, te


def ler_saida_r(caminho):
    """{pos: SNP.PP.H4 do coloc.abf}. Linhas '#' (versões do R) são ignoradas."""
    with open(caminho, encoding="utf-8") as fh:
        corpo = [l for l in fh if not l.startswith("#")]
    return {int(r["pos_b37"]): float(r["PP4_colocR"])
            for r in csv.DictReader(corpo, delimiter="\t")}


def cabecalho_r(caminho):
    """Linha de comentário gravada pelo script R (versões e resumo), ou ''."""
    with open(caminho, encoding="utf-8") as fh:
        primeira = fh.readline().strip()
    return primeira[1:].strip() if primeira.startswith("#") else ""


def _massa_finemap(nome):
    """Linha do STX6 no relatório de fine-mapping (se existir)."""
    p = REPORTS / nome
    if not p.exists():
        return None
    m = re.search(r"## Locus STX6.*?\*\*r²≥0\.8: ([\d.,]+%)\*\*",
                  p.read_text(encoding="utf-8"), re.S)
    return m.group(1) if m else None


def main() -> None:
    if "--so-oraculo" not in sys.argv:
        exportar_entrada()
    hoje = datetime.date.today().isoformat()

    # (1) oráculo congelado: entrada v0.4 + SNP.PP.H4 do coloc.abf real
    pos_f, tg_f, te_f = ler_entrada(FIX_ENT)
    r_f = ler_saida_r(FIX_R)
    res_f = C.coloc(tg_f, te_f, "cc", "quant", sdY2=1.0)
    dif_f = max(abs(a - r_f[p]) for a, p in zip(res_f["snp_pp_h4"], pos_f))
    pf = res_f["pp"]

    # (2) entrada atual (harmonização por variante completa)
    pos_a, tg_a, te_a = ler_entrada(ENTRADA)
    res_a = C.coloc(tg_a, te_a, "cc", "quant", sdY2=1.0)
    pa = res_a["pp"]
    mudou = [p for p, a, b in zip(pos_a, te_a, te_f) if a != b] \
        if pos_a == pos_f else None

    L = ["# Validação da colocalização contra o pacote R `coloc`", "",
         f"*`crosscheck_coloc_R.py` em {hoje}. Compara a implementação única "
         "do projeto (`pipeline/lib/coloc_abf.py`) com a saída do `coloc.abf` "
         "real (R 4.5, coloc 5.2.3, executado em 2026-08-29) guardada em "
         "`tests/fixtures/`.*", "",
         "## O que aconteceu (registro da correção de 2026-11)", "",
         "A versão deste relatório publicada até a v0.4 concluía que \"a "
         "máquina está correta — acordo a 6 decimais\" e que STX6 não "
         "colocalizava (H3 = 0,994997). As duas conclusões estavam erradas:", "",
         "1. **A parte (A) não era uma validação.** Ela reimplementava em R a "
         "mesma fórmula do Python, com o mesmo erro: o termo de H4 "
         "multiplicado por p1·p2·p12 em vez de p12. Duas cópias da mesma "
         "equação errada concordam entre si; a concordância não dizia nada "
         "sobre a equação.",
         "2. **A parte (B) já mostrava o resultado certo.** O `coloc.abf` do "
         "pacote oficial devolvia PP.H4 = 0,9796 para os mesmos 390 pares. A "
         "diferença foi explicada como \"dependência de prior\". Não era: com "
         "o mesmo prior de efeito do código antigo (W = 0,0016) e apenas a "
         "combinação corrigida, PP.H4 = 0,9788.",
         "3. **Havia mais dois erros na mesma função**: o fator de Bayes usava "
         "√r onde Wakefield usa √(1−r), e W = 0,04² era rotulado \"coloc "
         "padrão\" (o padrão do `coloc` é variância 0,04 para caso-controle e "
         "0,0225·sdY² para traço quantitativo). Mudam pouco o H4 desta região, "
         "mas eram erros de fórmula.", "",
         "A fórmula antiga foi removida de todo o código. O que valida a "
         "implementação atual são (i) este oráculo, (ii) casos sintéticos "
         "conferidos contra a enumeração explícita das configurações causais "
         "e (iii) as âncoras numéricas da auditoria — todos em "
         "`tests/test_coloc_abf.py`.", "",
         "| Configuração, entrada de 390 variantes da v0.4 | H2 | H3 | H4 |",
         "|---|---|---|---|",
         "| Publicada até a v0.4 (três erros; valores históricos) | 0,005002 | "
         "0,994997 | 4,6×10⁻⁷ |"]
    w = 0.0016
    pb = C.coloc(tg_f, te_f, w1=w, w2=w)["pp"]
    L += [f"| Fórmula correta com o prior de efeito antigo (W = 0,0016) | "
          f"{pb[2]:.6f} | {pb[3]:.6f} | {pb[4]:.6f} |",
          f"| **Conforme o `coloc` (W = 0,04 / 0,0225)** | {pf[2]:.2e} | "
          f"{pf[3]:.6f} | **{pf[4]:.6f}** |",
          "| `coloc.abf` 5.2.3 real (registrado em 2026-08-29) | 7×10⁻⁶ | "
          "0,0204 | 0,9796 |", "",
          "## Oráculo: `coloc.abf` real × implementação do projeto", "",
          f"- Entrada: `tests/fixtures/stx6_crosscheck_input_v04.tsv` "
          f"({len(pos_f)} variantes).",
          "- Referência: `tests/fixtures/stx6_crosscheck_colocR_output_v04.tsv`, "
          "coluna `PP4_colocR` (posterior por variante, `SNP.PP.H4`).",
          f"- **Diferença máxima, variante a variante: {dif_f:.2e}.**",
          f"- Resumo do projeto: PP.H4 = {pf[4]:.4f}, PP.H3 = {pf[3]:.4f}, "
          f"PP.H2 = {pf[2]:.1e} — igual ao que o pacote registrou "
          "(0,9796 / 0,0204 / 7×10⁻⁶).",
          "- O `coloc.abf` foi chamado com `type = \"cc\"` para o GWAS e "
          "`type = \"quant\"`, `sdY = 1` para o eQTL-meta, com β e variância de "
          "β fornecidos. Nessa forma, N, MAF e a fração de casos não entram "
          "no fator de Bayes.", "",
          "## Entrada atual (harmonização por variante completa)", "",
          f"- `pipeline/data/stx6_crosscheck_input.tsv`: {len(pos_a)} variantes; "
          f"PP.H4 = {pa[4]:.6f}, PP.H3 = {pa[3]:.6f}, PP.H2 = {pa[2]:.2e}."]
    if mudou is not None:
        L.append(
            f"- Mesmas {len(pos_a)} posições da entrada da v0.4. O efeito do "
            f"eQTL-meta mudou em {len(mudou)} variante(s)"
            + (f" (chr1:{', '.join(f'{p:,}' for p in mudou)})" if mudou else "")
            + ": sítio multialélico em que a versão anterior somava efeitos "
            "de alelos diferentes.")
    else:
        L.append("- O conjunto de posições difere da entrada da v0.4; o "
                 "oráculo acima continua valendo para a entrada congelada.")
    if SAIDA_R.exists():
        r_a = ler_saida_r(SAIDA_R)
        if set(r_a) == set(pos_a):
            dif_a = max(abs(a - r_a[p])
                        for a, p in zip(res_a["snp_pp_h4"], pos_a))
            L.append(f"- `coloc.abf` real executado sobre a entrada atual "
                     f"(`pipeline/data/stx6_crosscheck_colocR_output.tsv`; "
                     f"{cabecalho_r(SAIDA_R)}): **diferença máxima por "
                     f"variante {dif_a:.2e}**.")
        else:
            L.append("- `stx6_crosscheck_colocR_output.tsv` existe mas é de "
                     "outra entrada (posições diferentes); rode de novo "
                     "`Rscript pipeline/scripts/coloc_crosscheck_R.R`.")
    else:
        L.append("- O `coloc.abf` do R **não** foi reexecutado sobre a entrada "
                 "atual (R não estava disponível na máquina da correção). Para "
                 "fazê-lo: `Rscript pipeline/scripts/coloc_crosscheck_R.R` e "
                 "depois `python3 pipeline/scripts/crosscheck_coloc_R.py`.")

    m_all = _massa_finemap("relatorio_finemap_loci.md")
    m_eur = _massa_finemap("relatorio_finemap_loci_EUR.md")
    L += ["", "## Fine-mapping: mesmo fator de Bayes", "",
          "`finemap_ld.py` usa a mesma função de fator de Bayes de Wakefield "
          "da colocalização (`coloc_abf.abf`, W = 0,04). A forma já estava "
          "correta nesse script; a troca pela função única não altera os "
          "resultados."]
    if m_all or m_eur:
        L.append(f"Massa posterior do cluster de STX6 (r² ≥ 0,8 com o lead): "
                 f"painel ALL {m_all or 'n/d'}; painel EUR {m_eur or 'n/d'} "
                 "(ver `relatorio_finemap_loci.md` e "
                 "`relatorio_finemap_loci_EUR.md`).")
    L += ["", "## Reprodução", "",
          "- `python3 pipeline/scripts/crosscheck_coloc_R.py --so-oraculo` "
          "(sem GWAS e sem R)",
          "- `python3 -m pytest tests/test_coloc_abf.py -q`",
          "- Opcional, com R e `coloc` instalados: "
          "`Rscript pipeline/scripts/coloc_crosscheck_R.R`", ""]
    destino = REPORTS_SAIDA / "relatorio_validacao_coloc_R.md"
    destino.write_text("\n".join(L), encoding="utf-8")
    print(f"[ok] {destino}")
    print(f"     oráculo: diferença máxima {dif_f:.2e}; H4 (v0.4) = {pf[4]:.6f}; "
          f"H4 (entrada atual) = {pa[4]:.6f}")


if __name__ == "__main__":
    main()
