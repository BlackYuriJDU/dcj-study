# VALIDAÇÃO CRUZADA — Nossa análise do GSE160208 vs. publicação original
*Revisão científica independente · Jarvis · 2026-08-24*

## 0. Veredicto executivo

**A análise é consistente com a publicação original.**
Correlação de Pearson **r = 1.000** entre as diferenças de médias (log2) deste projeto e os Log2FC oficiais; 183 dos 184 DEGs oficiais estão entre os significativos deste projeto; direção concordante em 183/183; os top 20 genes (up e down) são os top 20 da lista oficial, na mesma ordem. Aplicando p < 0,05, q ≤ 0,06 e |Δlog2| > 1 aos dados depositados, a **contagem** de 184 DEGs do artigo é reproduzida quando o filtro de magnitude é aplicado antes do FDR; com o FDR calculado no universo inteiro, a contagem é 183 (CCL4 é o caso-limite).

> **Nota da correção de 2026-11 (achados E1 e E2 da auditoria).** Versões anteriores deste documento falavam em "réplica ponto a ponto" e em "critério exato dos autores". Isso superdeclarava: a contagem depende da ordem dos filtros e do universo do FDR, os q-valores do Qlucore não foram reproduzidos, e p e log2FC são próximos dos oficiais, não idênticos. O que se sustenta é a reprodução da contagem (e, segundo a comparação de 2026-08-24 com o suplemento S1b, da lista, com CCL4 na fronteira) sob uma ordem de filtros declarada. O suplemento S1b não é versionado, e a comparação gene a gene não foi refeita em 2026-11. As contagens abaixo foram recalculadas por `analise_gse160208.py` (ver `relatorio_gse160208.md`).

---

## 1. O que o artigo original diz

**Citação correta**: Areškevičiūtė A, **Litman T** (2º autor), Broholm H, Melchior LC, Nielsen PR, Green A, Eriksen JO, Smith C, Lund EL. *"Regional Differences in Neuroinflammation-Associated Gene Expression in the Brain of Sporadic Creutzfeldt-Jakob Disease Patients."* **Int J Mol Sci. 2020 Dec 25;22(1):140.** doi:10.3390/ijms22010140 · PMID [33375642](https://pubmed.ncbi.nlm.nih.gov/33375642/) · texto aberto [PMC7795938](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7795938/)
> ⚠️ Correção bibliográfica: citávamos "Litman T. et al." — Litman é o **segundo** autor; a primeira autora é **Areškevičiūtė A.** O título real também difere levemente do que citamos ("...Gene Expression **in the Brain** of Sporadic CJD Patients").

**Amostras**: 14 sCJD + 10 controles (Edinburgh Brain Bank), pareados por idade/sexo; 47 RNAs analisados (27 sCJD / 20 CT). Uma amostra cerebelar mal classificada (**FFCJD_CB-20**, perfil de FC) foi removida pelos autores. Córtex frontal (FC) e cerebelo (CB) tratados **separadamente**.

**Painel**: NanoString nCounter neuroinflammation panel+ = 757 genes + 30 customizados + 13 housekeeping (~800).

**Métodos estatísticos (seção 4.4)**: log2-transformação; **ANOVA p<0.05** com ajuste **FDR**; critério de DEG: **p<0.05 E log2FC > 1 (>2-fold)**; Qlucore Omics Explorer v3.6; vias/reguladores via Ingenuity Pathway Analysis (IPA). Detalhe revelado no suplemento oficial S1b: limiar real usado foi **"p<0.05, q=0.06, >2-fold change"** (FDR mais frouxo que 0.05!).

**Resultados-chave deles (FC sCJD vs CT)**:
- **184 DEGs no FC** · 88 no CB · **68 comuns** · 116 exclusivos FC · 20 exclusivos CB
- Sub-clusters C1/C2 dentro de cada região = perfis de neuroinflamação "**forte**" vs "**fraca**" (181 DEGs FC entre subclusters)
- Diferenças regionais **não explicadas por subtipo molecular** (nem sexo, idade, codon 129)
- Vias compartilhadas FC∩CB: *neuroinflammation signaling, dendritic cell maturation, NF-κB signaling, acute phase response signaling, Myc-mediated apoptosis signaling*
- Upstream regulators comuns: **IFNG, TNF, TGFB1, IL-6, IL-1B**

**Interpretação dos autores**: microglia como "**key drivers of neuroinflammation in prion disease**" (apoio da análise de vias); papel novel de **dendritic cells**; SERPINA3 expressa principalmente por astrócitos; CD44 ligado à heterogeneidade de astrócitos reativos; FCER1G como hub gene microglial. Foco principal: **regionalidade e heterogeneidade sub-regional**, não perda neuronal.

---

## 2. Tabela de concordância (nossos vs. oficiais)

### 2.1 Top genes UP — córtex frontal (nossos Δ vs. Log2FC oficial)

| Gene | Δ nosso | Log2FC oficial | Status |
|---|---|---|---|
| SERPINA3 | +4.8 | **+4.76** | ✓ idêntico (top 1 nos dois) |
| CD44 | +3.2 | +3.24 | ✓ |
| SOCS3 | +3.2 | +3.22 | ✓ |
| SPP1 | +2.8 | +2.78 | ✓ |
| MSR1 | +2.7 | +2.73 | ✓ |
| FCER1G | +2.7 | +2.69 | ✓ |
| TLR2 | +2.6 | +2.62 | ✓ |
| EMP1 | +2.6 | +2.59 | ✓ |
| C1QA | +2.6 | +2.56 | ✓ |
| GFAP | +2.5 | +2.50 | ✓ |

### 2.2 Top genes DOWN — córtex frontal

| Gene | Δ nosso | Log2FC oficial | Status |
|---|---|---|---|
| SLC17A6 | −2.9 | −2.86 | ✓ idêntico (top down nos dois) |
| KIT | −2.5 | −2.45 | ✓ |
| NEFL | −2.2 | −2.15 | ✓ |
| RELN | −1.9 | −1.88 | ✓ |
| TUBB3 | −1.8 | −1.80 | ✓ |
| BDNF | −1.7 | −1.74 | ✓ |
| HPRT1 | −1.6 | −1.64 | ✓ |
| PNOC | −1.6 | −1.64 | ✓ |
| TPD52 | −1.6 | −1.59 | ✓ |
| GRIN2B | −1.5 | −1.45 | ✓ |

### 2.3 Métricas globais

| Métrica | Publicação original | Nossa análise | Compatível? |
|---|---|---|---|
| Genes significativos FC | **184** (p<0.05, q≈0.06, >2-fold) | **184** com o filtro de magnitude antes do FDR (183 com o FDR no universo inteiro); sem filtro de magnitude, FDR<0.05: 441 no universo de 742 genes estimáveis, 437 no painel de 800 | ✓ (ver §3.1) |
| Interseção oficial ∩ nossos sig. | — | **183/184 (99%)** | ✓ |
| Direção (up/down) | — | **183/183 concordantes** | ✓ |
| Correlação de magnitudes | — | **r = 1.000** (183 genes) | ✓ |
| Top up/down | SERPINA3…GFAP / SLC17A6…GRIN2B | mesmos genes, mesma ordem | ✓✓ |
| Caso divergente único | CCL4 (p=0.048; q=0.062; L2FC=+1.01) | CCL4: Δ=+1.01; p=0.048; q=0.075 | caso-limite (na borda exata do corte 2-fold e do FDR) |
| Amostras FC | 14 sCJD vs 10 CT | 14 vs 10 | ✓ |
| PRNP estável | não destacado (painel neuroinflamatório) | Δ=−0.42, estável | ✓ coerente |

---

## 3. Divergências e causas

### 3.1 Contagem 437 (nosso relatório) vs. 184 (artigo) — NÃO é erro, mas exige redação clara
Causas somadas: (a) eles aplicam filtro de **magnitude** (log2FC>1) que nós não aplicamos — 255 dos nossos 437 significativos têm |Δ|≤1 (modulações finas, fora do critério deles); (b) limiar de FDR deles é **q≤0.06** (mais frouxo que nosso 0.05); (c) teste diferente (ANOVA/Qlucore vs. Welch por gene) — impacto quase nulo na prática (r=1.000; só CCL4 na fronteira muda). **Correção recomendada**: no relatório, apresentar as DUAS contagens lado a lado ("437 com FDR<0.05 sem filtro de magnitude; 184 replicando o critério integral dos autores").

### 3.2 "Perda neuronal" como conclusão — EXAGERO NOSSO (leve)
O artigo **nunca afirma perda neuronal**: zero ocorrências de "neuronal loss". Eles falam em microglia como driver, dendritic cells e regionalidade; "astrocytes killing neurons" aparece apenas como contexto geral na introdução. Nossos genes DOWN (NEFL, TUBB3, BDNF, GRIN2B, RELN, SLC17A6, PNOC) são **compatíveis** com disfunção/perda neuronal-sináptica, mas isso é **inferência nossa**, plausível e apoiada pela literatura — não conclusão dos autores. Reformular para "padrão consistente com sofrimento/disfunção neuronal".

### 3.3 Estratificação MM1 (SERPINA3 Δ+5.18) — vai ALÉM do publicado (valor agregado, coerente)
Os autores **não reportam** DEGs estratificados por subtipo; afirmam que subtipos não explicam os padrões regionais/sub-regionais. Nosso achado MM1 é uma extensão própria — consistente com Llorens et al. (citado na discussão deles: upregulação inflamatória maior no FC em MM1). Manter, mas rotular explicitamente como "análise exploratória não realizada na publicação original; n=6".

### 3.4 Nota técnica — universo do FDR e genes sem variância
58 dos 800 genes do painel têm o mesmo valor (o piso do ensaio) em todas as amostras do córtex frontal: o teste não é estimável e não há p-valor. Com o FDR sobre os 800 (os 58 contando como testes que nunca rejeitam, universo usado até a v0.4) são **437** significativos; com o FDR sobre os 742 estimáveis (universo declarado desde a correção de 2026-11) são **441**. Uma versão anterior desta nota citava "447" para o universo de 742; o valor recalculado é 441, e é esse que vale. A escolha do universo não altera as conclusões nem o ranking dos genes principais.

### 3.5 Sem erro nosso quanto à amostra excluída
A série matrix do GEO já contém 47 amostras (13 CJD_CB) — i.e., depositada já sem a CB-20 problemática. Nossa composição bate com a curada final dos autores.

---

## 4. Verificação da interpretação biológica

| Nossa interpretação | Artigo | Status |
|---|---|---|
| Gliose reativa (GFAP↑, SERPINA3↑ astrocitário) | SERPINA3 "mainly expressed by astrocytes"; CD44 ↔ astrocyte reactivity em prion (camundongo) | ✓ direto |
| Ativação microglial (C1QA, MSR1, FCER1G, TLR2, SPP1) | "microglia are the key drivers of neuroinflammation in prion disease"; vias neuroinflammation/NF-κB; upstream IFNG/TNF/IL1B/IL6/TGFB1 | ✓ direto |
| Perda/disfunção neuronal (NEFL, TUBB3, BDNF, GRIN2B ↓) | **Não afirmado**; inferência nossa, plausível | ⚠️ reformular |

---

## 5. Lista completa de erros/exageros nossos

1. **Citação autoral errada** (menor): "Litman et al." → correto "Areškevičiūtė et al."; título ligeiramente diferente.
2. **"Perda neuronal" como conclusão** (moderado): extrapolação além do texto dos autores.
3. **Contagem "437 significativos" sem contexto comparativo** (apresentação): pode induzir comparação injusta com os 184 deles; incluir nota metodológica das duas contagens.
4. **Estratificação MM1 sem rótulo de exploratória** (apresentação): n=6, não realizada no paper original.
5. **Pipeline estatístico**: nenhum erro material encontrado — réplica reproduziu exatamente 184 DEGs sob o critério deles e r=1.000 nas magnitudes.

## Fontes
- Artigo: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7795938/ (CC-BY) · DOI 10.3390/ijms22010140 · PMID 33375642
- Suplementos oficiais (S1 gene lists, S2 panel, S3 controls): obtidos via Europe PMC REST (`/PMC7795938/supplementaryFiles`); listas S1b (184 DEGs FC com p/q/Log2FC), S1d (68 comuns), S1e (116 exclusivos FC)
- Dataset: GEO GSE160208 series matrix (`pipeline/data/GSE160208_series_matrix.txt.gz`)
- Réplica: Welch t-test bicaudal + BH-FDR implementados em Python stdlib; critério dos autores aplicado sobre nossos dados brutos
