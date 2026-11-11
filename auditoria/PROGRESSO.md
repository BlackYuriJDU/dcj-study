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

## Onda 1 — Colocalização (G1, G2, G3, G4, G5, G7, T1) ✅

Feito:
- `pipeline/lib/coloc_abf.py`: implementação única, em espaço logarítmico,
  conforme `coloc::combine.abf` / `approx.bf.estimates`; prior de efeito por
  tipo de traço (0,04 caso-controle; 0,0225·sdY² quantitativo); tabelas de
  sensibilidade a p12 e a W; diagnóstico da razão H4/H3 a priori.
- Cópias removidas: `coloc_stx6_eqtl.py`, `coloc_meta_stx6.py`,
  `coloc_sqtl_stx6.py`, `coloc_crosscheck_R.R` (parte A eliminada) e
  `finemap_ld.py` (passa a usar o mesmo fator de Bayes).
- Testes novos: `tests/test_coloc_abf.py` (sintéticos contra enumeração
  explícita, oráculo do coloc.abf real em `tests/fixtures/`, âncoras 2.1).
  O teste antigo `test_coloc_meta_tem_h2_e_validacao` (exigia "0.9950") foi
  substituído.
- `relatorio_validacao_coloc_R.md` reescrito (gerado por
  `crosscheck_coloc_R.py`), dizendo o que aconteceu.
- R instalado num ambiente temporário (micromamba/conda-forge, fora do
  repositório): o `coloc.abf` real (R 4.5.3, coloc 5.2.3) foi executado sobre a
  entrada nova → H4 = 0,979605; diferença máxima por variante 1,6e-15.

Comprovação:
- `python -m pytest tests/test_coloc_abf.py -q` → 36 passed.
- Mutação (não versionada): com H4 × p1·p2·p12 → 19 testes falham; com √r →
  21 falham; com W = 0,0016 → 14 falham.
- `python pipeline/scripts/crosscheck_coloc_R.py --so-oraculo`

Decisões:
4. **Scripts reescritos uma vez só.** A troca da fórmula (onda 1) e a
   harmonização por variante completa (G6, onda 2) mexem nas mesmas funções;
   por isso `harmoniza.py`, `fontes_qtl.py`, `gwas_fonte.py` e a proveniência
   entram já no commit da onda 1. Custo: a fronteira entre os commits das ondas
   1 e 2 não é a do prompt.
5. **`pipeline/data/stx6_crosscheck_colocR_output.tsv`** deixou de ser a saída
   de 2026-08-29 (movida, intacta, para `tests/fixtures/…_v04.tsv`) e passou a
   ser a saída do coloc.abf sobre a entrada nova. É arquivo derivado.
6. **sdY = 1** para o eQTL (expressão normalizada); a estimativa por conjunto
   (0,67–1,00) está no relatório da meta e a sensibilidade a W cobre o efeito.

Divergências em relação ao dossiê (detalhes em CORRECOES.md):
- Linhas repetidas por rsID: 108/88/113/33/33 (dossiê: 86/69/91/24/24) —
  contagem de linhas excedentes × outra definição; posições multialélicas e
  indels coincidem (29/18/38/0/0; 181/137/221/44/45).
- Heterogeneidade nas 390: Q p<0,05 em 36 e I²>50% em 51 (dossiê 37 e 52): a
  diferença é o sítio 180.975.614, cuja heterogeneidade vinha da mistura de
  alelos.
- Sítio 180.975.614: com a chave completa, G>C existe nos cinco conjuntos
  (k = 5); o dossiê descreve "G>T em dois e G>C em três", que é o que a lógica
  "última linha vence" enxergava.
