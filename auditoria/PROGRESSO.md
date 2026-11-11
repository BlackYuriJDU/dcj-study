# Progresso da correção — auditoria 2026-11

Memória entre contextos. Ao retomar: ler este arquivo e `git log` da branch.

- Dossiê: `auditoria/DOSSIE_AUDITORIA_DCJ_STUDY.md`
- Branch de trabalho: `claude/new-session-ukxwtl` (designada pelo ambiente da
  sessão na nuvem; equivale à `correcao/auditoria-2026-11` pedida no prompt —
  ver "Decisões"). Base: `main` em `e0765256110b521220902873d607b4dbe08fd558`.
- Ambiente: `python3 -m venv .venv && .venv/bin/pip install ...` (Python 3.13).
- **Sem `git push`, PR, tag, release, Zenodo ou deploy** (decisões do autor).

## Decisões (cada uma com o custo se estiver errada)

1. **Branch.** O ambiente da sessão exige a branch `claude/new-session-ukxwtl`;
   o prompt pede `correcao/auditoria-2026-11`. Trabalho na primeira e crio a
   segunda como rótulo local no mesmo commit ao final. Custo se errado: um
   `git branch -m`.
2. **Script do Apêndice A não é versionado.** Ele contém, de propósito, uma
   cópia da fórmula antiga (`combina_repo`); o objetivo 1 pede que nenhuma cópia
   permaneça. Continua disponível no texto do dossiê. Custo: nenhum.
3. **Rede.** O FTP do EBI estava bloqueado; criei o ambiente de nuvem
   `dcj-study` liberando `ftp.ebi.ac.uk` (GWAS e eQTL Catalogue) e, depois,
   `rest.ensembl.org` e `eutils.ncbi.nlm.nih.gov` (fine-mapping e offset).

## Onda 0 — Reconhecimento ✅

- `git rev-parse HEAD` → `e0765256110b521220902873d607b4dbe08fd558` (= `main`).
- `python -m pytest tests/ -q -rs` → `37 passed`, 0 pulados.
- Script do Apêndice A sobre `pipeline/data/stx6_crosscheck_input.tsv` e
  `..._colocR_output.tsv`: saída idêntica à esperada, linha a linha
  ([A] H4=4.61894e-07 · [B] 0.978809 · [C] 0.978676 · [D] 0.979605 ·
  p12=5e-6 → 0.960025 · p12=1e-6 → 0.827677 · diferença máxima contra o
  coloc.abf real = 1.11e-16).
- R não está instalado na máquina (`which R` vazio). O oráculo versionado
  dispensa R.
