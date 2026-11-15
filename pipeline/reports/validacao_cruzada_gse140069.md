# Validação cruzada — nossa análise independente × artigo original (GSE140069)

*Validação executada em 2026-08-24 · revisor: Jarvis (ox-alpha) · método: leitura integral do texto PMC + Supplementary Data 1 oficial (xlsx dos autores) + re-execução controlada do nosso pipeline.*

---

## (i) O que o artigo diz

**Referência:** Norsworthy PJ, Pal S, Alibhai Z, et al. **"A blood miRNA signature associates with sporadic Creutzfeldt-Jakob disease diagnosis."** *Nature Communications* 11:3960 (2020). DOI [10.1038/s41467-020-17655-x](https://doi.org/10.1038/s41467-020-17655-x) · PMID 32769986.
**URLs verificadas:** [PMC7414116 (texto completo)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7414116/) · [nature.com](https://www.nature.com/articles/s41467-020-17655-x) · [Supplementary Data 1 (tabela DE, xlsx baixado)](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-020-17655-x/MediaObjects/41467_2020_17655_MOESM3_ESM.xlsx) · preprint medRxiv 2020.01.08.19015214.

### Assinatura final proposta
Painel de **3 miRNAs DOWNREGULADOS** no sangue total, validados e replicados por qPCR:

| miRNA | FC (descoberta, seq) | q (BH) | FC (replicação, qPCR) | p (qPCR) |
|---|---|---|---|---|
| hsa-miR-16-5p | −2.76 | 1.48×10⁻⁴ | −1.87 | 0.023 |
| hsa-miR-93-5p | −2.34 | 6.48×10⁻⁴ | −1.74 | 0.023 |
| hsa-let-7i-5p | −2.49 | 1.06×10⁻² | −1.82 | 0.014 |

Contexto da descoberta: **101 miRNAs testados** (filtro de cobertura média ≥5000 reads); **4 significativos** FDR<0.05 — os 3 acima + **hsa-miR-106b-3p** (FC −1.74, q=8.4×10⁻³; não replicável por qPCR por expressão muito baixa). **hsa-let-7d-3p** UP limítrofe (+1.94, q=0.053) e **NÃO replicou** no qPCR (FC +1.14, p=0.114). Outros 30 miRNAs de baixa abundância testados em subconjunto de 36 amostras: nenhum significativo.

### Métodos
- Small RNA-seq (TruSeq), Bowtie2/hg38, miRBase v21, normalização Cufflinks (total hits → FPKM).
- Estatística: **Partek Gene Specific Analysis (GSA)** — regressão/ANOVA com correção de variância gene-específica, **IDADE como covariável**, BH-FDR. **Não é DESeq2 nem limma.**
- Sexo NÃO ajustado; RIN tratado por análise de sensibilidade (excluir RIN<4 não alterou resultados); excluir os 3 casos "prováveis" não alterou; estratificação por códon 129 e tipo PrPSc sem efeito; sem normalização por contagem de células (sCJD não as altera).
- Replicação: coorte independente qPCR, 29 sCJD vs 30 controles (miScript/RNU6-2+snRNAs; Mann-Whitney). Alvos mRNA upregulados: CCND3, CDKN1A, ZFP36, NAP1L1 (p=0.032 cada); RNF44 ns.
- Clínica: sem correlação com idade de início/duração/MRC score-slope; longitudinal (21 pacientes): taxa de queda dos miRNAs ≠ taxa de progressão.

### Performance do classificador (ROC/Z-scores, SPSS — sem ML)
- sCJD vs controles (descoberta): AUC individual 0.736–0.762; **combinado (3 miRNAs) AUC 0.788**.
- Alzheimer vs controles: combinado AUC 0.860.
- **sCJD vs Alzheimer: AUC 0.924 combinado, especificidade 100%** (individuais 0.897–0.934). Sensibilidade no ótimo de Youden está no Suppl Table 6 (não citada no texto principal).

---

## (ii) Tabela de concordância (nossa v2 correta, log2(x+1)+Welch+BH × Supplementary Data 1 deles)

### Núcleo da assinatura deles — TODOS capturados por nós

| miRNA | Eles: FC / q (n=101) | Nós: log2FC / q (n=939) | Direção | Signif.? |
|---|---|---|---|---|
| hsa-miR-16-5p | −2.76 / 1.5e−4 | −1.35 / 4e−4 | ✔ | ✔✔ |
| hsa-miR-93-5p | −2.34 / 6.5e−4 | −1.23 / 5e−4 | ✔ | ✔✔ |
| hsa-miR-106b-3p | −1.74 / 8.4e−3 | −0.96 / 1.0e−3 (**nosso hit nº1**) | ✔ | ✔✔ |
| hsa-let-7i-5p | −2.49 / 1.1e−2 | −1.43 / 1.3e−3 | ✔ | ✔✔ |
| hsa-let-7d-3p (limítrofe UP) | +1.94 / 0.053 | +0.86 / 0.018 | ✔ | ✔(deles ns por pouco) |

Controles negativos também convergem: **miR-25-3p** ns nos dois (q 0.125 vs 0.124!) · **miR-484** ns nos dois (hemólise afastada).

### Nossos top hits fora do núcleo deles

| Nosso hit | Eles | Interpretação |
|---|---|---|
| miR-142-5p ↓ | NA tabela: FC −2.19, q=0.55 (ns) | direção idêntica, só a significância diverge |
| miR-4732-3p ↑ | NA tabela: FC +1.37, q=0.09 (quase!) | quase-significativo para eles — apoio parcial |
| miR-532-5p ↓, miR-92b-3p ↑, miR-320b ↑, miR-671-3p ↑, miR-4732-5p ↑ | NA tabela, todos ns | mesma direção, significância divergente |
| miR-106b-5p ↓, miR-29a-3p ↓, miR-221-3p ↓, miR-17-3p ↓, miR-22-3p ↓, miR-486-3p ↑, miR-423-5p ↑ | **FILTRADOS FORA** (não estão nos 101 testados) | hipóteses novas; o artigo não as contradiz |

### Métricas globais de concordância (n=101 miRNAs comuns)
- Direção do efeito concordante: **81/101 (80%)**
- Correlação de Pearson entre log2FCs: **r = +0.64**
- Significativos dentro dos 101 comuns: eles 4, nós 18 → nosso teste é sistematicamente mais liberal no MESMO conjunto.
- Sinal global: dominância de DOWN nos dois (eles 4/4 hits; nós 52/60 na v1, 74/84 na v2).

---

## (iii) Divergências e causas prováveis

1. **Universo testado (maior causa):** eles filtraram cobertura média ≥5000 reads → 101 miRNAs; nós testamos as 939 linhas do xlsx sem filtro. Múltiplas comparações 939×101 penalizam o nosso BH — e ainda assim achamos mais significativos, porque:
2. **Motor estatístico:** Partek GSA (regressão com idade + correção conservadora de variância gene-específica) vs Welch t-test simples sem covariáveis. O GSA é deliberadamente conservador.
3. **Escala/transformação:** eles FPKM→GSA; nós v2 log2(x+1). O pseudocount sobre valores-piso (0.0001) infla |log2FC| de miRNAs "ligado/desligado" (ex.: nosso miR-29a −7.55 em log2(x+1), mas −2.09 como razão de médias aritméticas).
4. **Definição de significância é igual (FDR<0.05, BH)** — a diferença vem de 1–3, não do critério.
5. **Idade:** covariável neles, ausente na nossa análise (limitação já auto-declarada no relatório).

---

## (iv) Veredicto honesto sobre a nossa análise

### O que está VALIDADO
O núcleo do artigo **reproduz integralmente na nossa análise independente**: os 4 miRNAS significantes da descoberta estão entre os nossos significativos, mesma direção, magnitudes quase idênticas (diferença ≤0.35 em log2); let-7d-3p limítrofe-up nos dois; os dois controles negativos deles (miR-25-3p, miR-484) são não-significativos nos dois; correlação de efeitos r=+0.64. Pipeline diferente (xlsx processado × counts brutos), mesmo sinal — isso é genuína replicação analítica.

### ERROS/EXAGEROS nossos (corrigir)
1. **CRÍTICO — inconsistência interna no `relatorio_gse140069.md`:** os números publicados (60 sig; log2FC da tabela) vieram do script **v1**: Welch em escala LINEAR + "log2FC" calculado como razão de médias aritméticas — duas escalas diferentes misturadas no mesmo relatório. O script atual (v2, com log2(x+1) antes do teste) produz **84 significativos (10↑/74↓)** com rankings diferentes (ex.: miR-320b cai para q=0.12; aparecem miR-500a/miR-29c/miR-144 no topo-down). **O relatório precisa ser regenerado com o v2.** *(RESOLVIDO 2026-08-24 14:29 — regenerado como v3 com OLS idade+sexo+RIN; ver ADDENDUM ao fim.)*
2. Chamar nossos 60/84 de resultado comparável ao "do artigo" induz a erro: eles testaram 101 pós-filtro; nós, 939 sem filtro. Comparação justa exige aplicar filtro de detecção.
3. Nossos top hits fora dos 101 deles (miR-29a, miR-221, miR-106b-5p etc.) são **hipóteses novas não testadas por eles** — não podem ser apresentados como confirmação ou contradição.
4. Sem ajuste por idade (eles ajustaram); pseudocount 1 sobre piso 0.0001 infla efeitos de não-detecção — preferível filtrar ou usar método de contagens.
5. Ponto forte a manter: honestidade das notas de limitações no relatório original estava correta ("sem correção idade/sexo/RIN nesta rodada").

### Recomendações práticas
- Regenerar o relatório com o script v2 e reportar lado a lado: (a) conjunto filtrado ≈ aos 101 deles; (b) conjunto completo como triagem exploratória.
- Adicionar idade como covariável (está no series matrix do GEO) e sensibilidade sem amostras de RIN<4, espelhando o artigo.
- Se formos citar performance, usar os AUCs deles com a ressalva de que são ROC empíricos em coortes pequenas e que a comparação sCJD×AD foi indireta (controles comuns).

---
*Fontes primárias: PMC7414116 (texto completo, acesso aberto CC-BY 4.0); Supplementary Data 1 (MOESM3, xlsx oficial dos autores, 101 miRNAs); Supplementary Table 6 referenciada mas não necessária ao veredicto.*

---

## ADDENDUM (2026-08-29) — estado final pós-correção

Este documento foi escrito às 13:28 de 24/08, ANTES da regeneração do relatório (14:29).
Estado entregue e verificado por re-execução em 29/08:

- `analise_gse140069.py` **v3** (log2(x+1) + OLS com idade+sexo+RIN + filtro de detecção ≥25% + Cohen's d) → `relatorio_gse140069.md` v3.
- Números finais: **84** significativos não-ajustado (FDR<0.05) → **69** com filtro de detecção → **1** após ajuste de covariáveis → **5** no universo filtrado (n=269). Núcleo do artigo: direção ↓ preservada nos 4; apenas **miR-93-5p** sobrevive ao FDR no universo filtrado (q=0.048).
- Re-execução em 29/08 reproduz o relatório byte-a-byte (única diferença: timestamp) e o `volcano_gse140069.png` é md5-idêntico à regeneração (c92c5b029fc2007d33543d85c73ecb32) — figura gerada pela mesma lógica v3.
- Os itens 1–4 da seção (iv) acima foram executados; o item 5 (nota de limitação) permanece correto como escrito.

## Atualização de 2026-11 (script v4, correção da auditoria)

- `analise_gse140069.py` **v4**: mesma análise da v3, com regressão por QR, testes não estimáveis marcados e universo do FDR declarado. 74 dos 939 miRNAs são constantes em todas as amostras e não têm teste; o universo declarado passa a ser o dos 865 estimáveis.
- Números vigentes (universo de 865; entre parênteses, o universo antigo de 939, que é o dos números acima): **89** (84) significativos sem ajuste → **73** (69) com filtro de detecção → **1** (1) após ajuste de covariáveis → **5** no universo filtrado (n=269).
- A afirmação de reprodução "byte a byte" acima vale para a v3 de 29/08; o relatório atual é o da v4 e é conferido contra a saída do script por `tests/test_reprodutibilidade.py`.
- O sobrevivente do modelo ajustado, hsa-miR-500a-3p, é um efeito do tipo detectado/não detectado (acima do piso em 27 de 57 casos e 46 de 48 controles).
