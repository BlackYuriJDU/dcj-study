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

## Onda 2 — Dados genéticos (R1, G6, G8, R4, G10, seção 2.7) ✅

Feito:
- R1: `pipeline/scripts/obter_gwas.py` (URL do GWAS Catalog + MD5
  `071790d80ccae0b41adabc7e0eefbf53`, conferido também contra o `md5sum.txt`
  remoto) e `pipeline/lib/gwas_fonte.py` (`exigir_gwas()` com mensagem clara).
  Os sete scripts que abriam o arquivo direto passam por `exigir_gwas()`.
- G6: `pipeline/lib/harmoniza.py` — chave (posição, REF, ALT), indels aparados,
  deduplicação por rsID, meta só de alelos idênticos, casamento com o GWAS
  pelos dois alelos com inversão de sinal, descartes contados.
- G8/R4: `pipeline/data/proveniencia_qtl.tsv` (release r8, metadata commit
  7a5b65b, n oficiais 268/263, SHA-256 de cada cache). Os scripts leem os n
  dali. Verificação remota em 2026-11-10: os 9 caches são idênticos, byte a
  byte, ao que o eQTL Catalogue serve hoje para a mesma janela
  (`registra_proveniencia_qtl.py --remoto`).
- G10: sQTL recalculado; "Bonferroni" sobre PP.H4 retirado. Maior PP.H4 entre
  os 6 testes em genes vizinhos: 0,60 (antes ≈ 0 por causa do erro G1).
- Seção 2.7 reexecutada com o GWAS real: 6.314.492 variantes, 0 malformadas,
  41 hits p<5e-8, lead STX6 chr1:180.961.245 p=7,51e-9, massa de fine-mapping
  do STX6 90,5% (painéis ALL e EUR, via Ensembl/NCBI).
- D3 (adiantado): constante única 0,454936 para o λGC → 1,0587 nos dois
  relatórios (antes 1,0588 no QC, que usava 0,4549).

Comprovação:
- `python pipeline/scripts/obter_gwas.py --conferir`
- `python pipeline/scripts/registra_proveniencia_qtl.py`
- `python -m pytest tests/test_harmoniza.py tests/test_coloc_pipeline.py -q`
- Relatórios regenerados: coloc por tecido, meta, sQTL, QC, clumping,
  fine-mapping (ALL/EUR), λGC. Clumping e fine-mapping saíram idênticos aos
  versionados, salvo a linha de data.

Sensibilidades da seção 2.2 (todas reproduzidas): DL 0,9755; leave-one-out
0,9794–0,9796 (n 365–380); sem GTEx 0,9795 (n 363); só CommonMind/ROSMAP/
BrainSeq 0,9782/0,9790/0,9792; só GTEx DLPFC 0,7015; só GTEx cerebelo 0,5672.

Achado novo (não estava no dossiê): vários relatórios versionados tinham
trechos escritos à mão DEPOIS de gerados (QC seção 3b, nota "Reconciliação" e
"perda neuronal" no GSE160208). Regenerar os apagava. Os trechos foram
passados para os scripts, com os números calculados. Um deles estava errado:
"30 dos 41 hits no bloco de PRNP" — são 16 (PRNP) + 23 (STX6) + 2 (GAL3ST1).

## Onda 3 — miRNA e expressão (M1, M2, M3, M5, E1, E2, E3) ✅

Feito:
- M1: `pipeline/lib/ols.py` (QR de Householder; resposta constante, posto
  deficiente, gl ≤ 0 e ajuste perfeito → "não estimável", sem p).
  Reproduzido o defeito antigo: os 74 miRNAs constantes recebiam p = 0,0203
  com float do Python e ~9e-6 com numpy.float64 (o dossiê cita 5,5e-7; o valor
  depende do ambiente, o que é o próprio defeito).
- M2: universo do FDR declarado = testes estimáveis (865 miRNAs; 742 genes).
  O universo legado (939; 800) é reportado ao lado.
- M3: grupo/sexo validados, idade/RIN numéricos, ids únicos no FDR; erro no
  inesperado.
- M5: miR-500a-3p detectado em 27/57 casos × 46/48 controles — efeito do tipo
  detectado/não detectado; dito no relatório.
- E1/E2/E3: ordem dos filtros e universo como parâmetros em
  `analise_gse160208.py`; "réplica exata" retirado; 437 × 447 resolvido.

Números (universo declarado | legado):
- miRNA A (Welch): 89 | 84 · A′: 73 | 69 · B (ajustado): 1 | 1 · B filtrado: 5
- idade: idade sozinha 121 | 113 · RIN 93 | 91 · grupo 92 | 88 ·
  grupo+idade 4 | 4 · grupo+RIN 28 | 26 · grupo+sexo 94 | 85 · completo 1 | 1
- cérebro FDR<0,05: 441 (315↑/126↓) | 437 (314↑/123↓); 58 genes sem variância
- critério do artigo: 183 (FDR no universo) ou 184 (magnitude antes do FDR);
  o gene da fronteira é CCL4.

Divergência: `validacao_cruzada_gse160208` citava "447" para o universo de 742
genes; o recalculado é 441. Não reproduzi 447 em nenhuma configuração.

Pendente (motivo concreto): identidade gene a gene com a lista S1b e os
q-valores do Qlucore — o suplemento do artigo não está no repositório.

Comprovação: `python -m pytest tests/test_ols.py tests/test_analises_expressao.py -q`

## Onda 4 — Integração sangue × cérebro (I1 a I5) ✅

Feito: K restrito ao painel testável; espécie do alvo humana; não testáveis
fora do BH (3 testes em vez de 5); leitura com `csv`; "controle negativo" e
leitura mecanística retirados; poder declarado (k mínimo para p<0,05).

Números: antes (esperado; p) 18,8/6,3/4,3; 1,00/0,98/0,996 → âncora "correto"
do dossiê (N=800, n=314) 5,1/2,4/1,2; 0,40/0,44/0,78 → final (N=742, n=315,
alvos humanos) 5,5/2,5/1,3; 0,50/0,51/0,81. Os dois alvos não humanos
(Ccnd1 de rato, Sp7 de camundongo) não estavam no painel: K não muda (13/6/3).

Comprovação: `python -m pytest tests/test_integracao.py -q`

Decisão 7: ondas 3 e 4 feitas em sequência, sem subagentes — a onda 4 importa
funções que a onda 3 reescreve; rodar em paralelo criaria conflito no mesmo
arquivo. Subagente fica para a verificação independente (onda 7).

## Onda 5 — Reprodutibilidade (R2, R3, R5, R6, T2) ✅ (R6: ver resultado abaixo)

Feito:
- R2: `requirements.txt` com versões fixadas; CI instala a partir dele
  (Python 3.13).
- R3: Makefile com `.venv`, `dados-gwas`, `conferir-dados`,
  `relatorios-offline`, `relatorios-gwas`, `relatorios-rede`, `relatorios`,
  `pdf`; `clean` não toca em dados versionados.
- R5: dez `.tbi` removidos da raiz; `*.tbi` ignorado.
- T2: testes de infraestrutura (`tests/test_reprodutibilidade.py`), incluindo
  um que roda cada script offline e exige saída igual ao relatório versionado.
- `obter_gwas.py` retoma downloads interrompidos (Range).
- `finemap_ld.py`: cache configurável por `DCJ_CACHE` (padrão `/tmp`).

R6 — regeneração a partir de um clone limpo (outro diretório, `.venv` novo):
- 1ª tentativa: `make install` e `make test` ok; `make dados-gwas` falhou por
  reset de conexão (o script não retomava); os relatórios offline regeneraram,
  e a comparação mostrou dois relatórios versionados defasados em relação aos
  scripts. Os dois problemas foram corrigidos (retomada; teste de defasagem).
- 2ª tentativa: ver "Resultado da 2ª regeneração" mais abaixo.

Achados novos tratados nesta onda: gráfico de floresta dos miRNAs com valores
digitados à mão e errados (passou a ser calculado); rótulos dos volcanos na
altura de q em vez de p.

Decisão 8: uma tentativa de limpar caches em `/tmp` foi barrada pelo
ambiente (diretório compartilhado). Em vez de apagar, o cache do
fine-mapping virou parâmetro (`DCJ_CACHE`).

## Onda 6 — Textos (G9, M4, D1, D2, D3, T3) ✅

Feito:
- Manuscrito v0.5 (2026-11-11) com nota de correção no início; reescritos:
  resumo, §2 (integridade), §3.4, §3.5, §4.1 (tabela de estratos), §4.5,
  §4.7 inteiro, §4.8, discussão (2º e 4º pontos e parágrafo metodológico),
  limitações, §8 e rodapé. PDF regenerado (o gerador ganhou suporte a bloco de
  citação e a itens de lista).
- `README.md`, `README.en.md`, `CITATION.cff` (resumo), site (pt e en: 441
  DEGs, contagem de 184 sob ordem de filtros declarada, 89→1 miRNAs,
  integração sem leitura de mecanismo).
- `CORRECOES.md` na raiz.
- T3: `test_preprint_v04_e_assinado`, `test_confundimento_idade_atribuido` e
  `test_relatorio_gse140069_e_v3` substituídos; `tests/test_textos.py`
  recalcula cada número antes de procurá-lo no manuscrito, nos READMEs e no
  site.

Achados novos tratados: o manuscrito dizia que o menor p do eQTL-meta era
6,60e-47 em 180.949.780 (é 3,44e-47 em 180.951.692); tabela de λ por estrato
com n em dobro; ENSG00000230124 era citado como RP5-1180C10.2 (o Ensembl hoje
devolve ACBD6; o texto passou a citar o identificador).

Não tratado (vai para o relatório final): referência 2 do manuscrito com
título/autores divergentes do registro do próprio repositório; mapa do
`README.en.md` lista arquivos não versionados; versão/DOI em `CITATION.cff` e
a frase "release v1.0.0-preprint-clean" continuam os da v0.4 (publicar é
decisão do autor).

Comprovação: `python -m pytest tests/ -q -rs` → 196 passed, 0 pulados.
