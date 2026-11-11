# Relatório — Análise de dados REAIS: GSE160208
*Gerado por `analise_gse160208.py` em 2026-11-11.*

- Dataset: Gene expression in the brain of sporadic Creutzfeldt-Jakob disease patients (CJD), and normal controls (CT)
- Fonte: GEO/NCBI GSE160208 · PMID 33375642 · Univ. Copenhagen (dados públicos anonimizados)
- Amostras totais processadas: **47** · Genes no painel: **800**

## Composição das amostras
- CJD_CB: 13
- CJD_FC: 14
- CT_CB: 10
- CT_FC: 10

## Covariáveis disponíveis (metadados reais)

- gender: F: 23, M: 24
- codon 129: MM: 22, MV: 14, VV: 11
- cjd subtype: MM1: 12, MM1+2: 2, MV1: 4, MV2: 2, NA: 20, VV2: 7

## Estratificação por subtipo — córtex frontal CJD
- Grupos CJD-FC por subtipo: MM1: 6, MM1+2: 1, MV1: 2, MV2: 1, VV2: 4

### Subgrupo MM1 (n=6) vs. controles — top 5 up/down

| Gene | Δ(MM1−CT) |
|---|---|
| SERPINA3 | +5.18 |
| SOCS3 | +3.75 |
| CD44 | +3.66 |
| FCER1G | +3.48 |
| MSR1 | +3.46 |
| SLC17A6 | -3.08 |
| KIT | -2.72 |
| NEFL | -2.55 |
| BDNF | -2.03 |
| TPD52 | -2.00 |

## Top 10 genes MAIS expressos em CJD (córtex frontal)

| Gene | Média CJD | Média CT | Δ |
|---|---|---|---|
| SERPINA3 | 10.4 | 5.6 | +4.8 |
| CD44 | 9.0 | 5.7 | +3.2 |
| SOCS3 | 8.1 | 4.9 | +3.2 |
| SPP1 | 12.0 | 9.3 | +2.8 |
| MSR1 | 8.0 | 5.3 | +2.7 |
| FCER1G | 8.5 | 5.8 | +2.7 |
| TLR2 | 7.6 | 5.0 | +2.6 |
| EMP1 | 9.0 | 6.5 | +2.6 |
| C1QA | 8.1 | 5.6 | +2.6 |
| GFAP | 16.6 | 14.1 | +2.5 |

## Top 10 genes MENOS expressos em CJD (córtex frontal)

| Gene | Média CJD | Média CT | Δ |
|---|---|---|---|
| SLC17A6 | 5.5 | 8.4 | -2.9 |
| KIT | 5.6 | 8.1 | -2.5 |
| NEFL | 11.1 | 13.3 | -2.2 |
| RELN | 7.5 | 9.4 | -1.9 |
| TUBB3 | 5.5 | 7.3 | -1.8 |
| BDNF | 5.4 | 7.2 | -1.7 |
| HPRT1 | 10.0 | 11.6 | -1.6 |
| PNOC | 5.4 | 7.1 | -1.6 |
| TPD52 | 8.8 | 10.4 | -1.6 |
| GRIN2B | 8.8 | 10.2 | -1.5 |

## Verificação específica
- PRNP presente no painel: sim → Δ(CJD−CT) = -0.42

## Estatística inferencial (Welch + BH-FDR, córtex frontal)

- Genes no painel: 800 · **não estimáveis: 58** (mesmo valor em todas as amostras do córtex frontal — o piso do ensaio; sem variância não há teste nem p-valor).
- **Universo declarado do FDR: 742 genes estimáveis.** Significativos com FDR<0.05: **441** (315 ↑, 126 ↓).
- Universo do painel inteiro (800 genes, com os não estimáveis contando como testes que nunca rejeitam; usado até a v0.4): **437** (314 ↑, 123 ↓). A diferença entre as duas contagens é só o tamanho do universo.

| Gene | p | q(FDR) |
|---|---|---|
| DLG4 | 1.03e-07 | 3.11e-05 |
| EHMT2 | 1.53e-07 | 3.11e-05 |
| OSMR | 1.67e-07 | 3.11e-05 |
| MYD88 | 1.68e-07 | 3.11e-05 |
| DLX1 | 3.46e-07 | 5.13e-05 |
| TRAF3 | 4.40e-07 | 5.39e-05 |
| TLR2 | 5.09e-07 | 5.39e-05 |
| CSF1 | 7.95e-07 | 5.63e-05 |
| C1QA | 8.57e-07 | 5.63e-05 |
| IFITM3 | 9.72e-07 | 5.63e-05 |
| PMP22 | 1.05e-06 | 5.63e-05 |
| SPP1 | 1.05e-06 | 5.63e-05 |
| PLXDC2 | 1.06e-06 | 5.63e-05 |
| KIT | 1.08e-06 | 5.63e-05 |
| SERPINA3 | 1.14e-06 | 5.63e-05 |

## Critério do artigo original: 184 genes, e do que esse número depende

O artigo (Areškevičiūtė et al. 2020) reporta **184** genes diferencialmente expressos no córtex frontal, com "p<0.05, q=0.06, >2-fold change" (suplemento S1b), calculados no Qlucore. Aplicando p < 0,05, q ≤ 0,06 e |Δlog2| > 1 aos dados depositados, a contagem depende da ordem dos filtros e do universo do FDR:

| Ordem dos filtros | Universo do FDR | DEGs |
|---|---|---|
| FDR no universo inteiro, depois magnitude | 742 estimáveis | **183** |
| FDR no universo inteiro, depois magnitude | 800 do painel | **183** |
| magnitude primeiro, depois FDR só nos que passaram | 185 genes com \|Δlog2\| > 1 | **184** |
| só p < 0,05 e \|Δlog2\| > 1 (sem exigir q) | — | **184** |

Leitura: a contagem de 184 é reproduzida quando o filtro de magnitude vem antes do FDR (ou sem exigir q); com o FDR calculado no universo inteiro, o resultado é 183. O gene que separa as duas contagens é CCL4, um caso-limite nos três cortes. Isto reproduz a CONTAGEM do artigo sob uma escolha declarada de ordem de filtros; não é uma réplica exata do procedimento dos autores: os q-valores do Qlucore não foram reproduzidos, e p e log2FC são próximos dos oficiais, não idênticos por construção. A comparação gene a gene com a lista do suplemento S1b está em `validacao_cruzada_gse160208_artigo_original.md`; o suplemento não é versionado neste repositório, e essa comparação não foi refeita na correção de 2026-11.
A ordem padrão deste script é `fdr_no_universo`, no universo `estimaveis`.

## Notas
- Welch t-test bicaudal implementado em stdlib; FDR Benjamini–Hochberg.
- Painel dirigido (800 genes neuroinflamatórios), não transcriptoma total.
- Sem correção para covariáveis (idade, PMI) — os metadados brutos não as trazem.
- **"Perda neuronal" é inferência deste projeto** (apoiada nos genes ↓ neurônio-específicos SLC17A6/NEFL/BDNF/TUBB3/GRIN2B), não conclusão do artigo original, que foca em regionalidade, microglia e células dendríticas. Ler como "padrão consistente com disfunção neuronal", não como demonstração histológica.
- Genes não estimáveis (variância nula no córtex frontal): BIK, BLK, C6, CCL7, CCR2, CCR5, CD200R, CD244, CD300LF, CD3D, CD3E, CD3G, CD69, CD70, CEACAM3, CPA3, CRYBA4, CSF2RB, CTSE, CTSW, CXCL10, CXCL9, DUOXA1, EOMES, FCRLA, FCRLB, GZMA, GZMB, HCAR2, IL1A, IL1RN, IL21R, IL27, IRF4, KIR3DL1, KLRD1, LCN2, LY9, MARCO, MMP12, MS4A1, MS4A2, NCR1, RAD51, S1PR4, SH2D1A, SLAMF9, SPIB, TBX21, TCL1A, TNFRSF17, TNFRSF4, TP73, TRAT1, TTR, UTY, XCL1, ZBP1.