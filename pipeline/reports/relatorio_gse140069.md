# Relatório — GSE140069 (SANGUE) — v4 com ajuste de idade/sexo/RIN
*`analise_gse140069.py` v4 em 2026-11-11. Mesma análise da v3, com solver por QR, testes não estimáveis marcados e universo do FDR declarado (correção da auditoria de 2026-11; ver `CORRECOES.md`).*

- Amostras: 57 sCJD vs. 48 controles · miRNAs medidos: 939
- Covariáveis do series matrix: idade (casos ~66 vs. controles ~54 anos), sexo (60 F, 45 M, 0 ausentes), RIN
- Amostras com sexo, idade e RIN completos (usadas no modelo ajustado): 105 (0 excluídas)
- Matriz de desenho do modelo ajustado: posto 5 de 5 colunas
- Filtro de detecção (≥25% das amostras acima do piso): **269 de 939** miRNAs

## Universo dos testes (declarado antes do FDR)

- **74 miRNAs não são estimáveis** (resposta constante em todas as amostras). Não recebem p-valor e não entram no FDR.
- **Universo declarado: 865 testes estimáveis** (A e B). É o universo usado nas contagens principais.
- Universo legado (939, usado até a v3): todos os miRNAs medidos, com os não estimáveis contando como testes que nunca rejeitam. É válido, porém mais conservador; fica na última coluna para comparação com os números já publicados.
- Na v3, os miRNAs constantes recebiam do solver antigo um "p" que era ruído de arredondamento: 0,0203 com números do tipo float do Python e da ordem de 10⁻⁶ a 10⁻⁷ com os mesmos números em numpy.float64 (o valor exato varia com o ambiente). O resultado publicado só não mudava porque a execução caiu do lado inofensivo.

## A vs. B

| Análise | FDR<0,05 no universo declarado (865) | FDR<0,05 no universo legado (939) | Interpretação |
|---|---|---|---|
| A: Welch log2(x+1), SEM covariáveis | 89 (12↑/77↓) | 84 (10↑/74↓) | triagem não ajustada — inflada por idade/RIN |
| A′: A ∩ filtro de detecção | 73 | 69 | triagem em miRNAs bem medidos |
| B: OLS ajustado (idade+sexo+RIN) | **1** (0↑/1↓) | **1** | **o número que vale** |
| A ∩ B | 1 | 1 | núcleo robusto às covariáveis |
| B no universo filtrado por detecção (n=269) | **5** | — | sensibilidade com correção menor |

**Veredicto: a maior parte do sinal bruto não sobrevive ao ajuste por idade, sexo e RIN.** O que resta é a direção (predomínio de miRNAs ↓) e o núcleo do artigo em significância NOMINAL (não em FDR) — ver tabela abaixo.

## Núcleo da assinatura do artigo original (Nat Commun 2020)

| miRNA | p nominal (ajustado) | q A (865) | q B (865) | q B (universo filtrado) |
|---|---|---|---|---|
| hsa-miR-16-5p | 0.0060 | 3.84e-04 | 3.06e-01 | 1.46e-01 |
| hsa-miR-93-5p | 0.0007 | 4.87e-04 | 1.29e-01 | 4.80e-02 |
| hsa-let-7i-5p | 0.0404 | 1.16e-03 | 4.83e-01 | 2.47e-01 |
| hsa-miR-106b-3p | 0.0079 | 9.34e-04 | 3.26e-01 | 1.46e-01 |

**Leitura**: os 4 mantêm direção ↓ e p nominal < 0,05; após FDR, apenas miR-93-5p sobrevive no universo filtrado (q=0.048). A assinatura publicada é mais frágil sob ajuste padrão do que a apresentação original sugere. Diferenças plausíveis: Partek GSA (correção de variância gene-específica) vs. OLS comum, e universo de testes (101 no artigo; 269 ou 865 aqui).

## Significativos no modelo ajustado (B) — com tamanho de efeito

| miRNA | β grupo (log2) | p | q(FDR) | d | detectado em sCJD | detectado em controles |
|---|---|---|---|---|---|---|
| hsa-miR-500a-3p | -7.94 | 8.03e-06 | 6.95e-03 | -1.20 | 27/57 | 46/48 |

**Como ler o β de -7.9 de hsa-miR-500a-3p.** Não é uma queda de abundância de ~250 vezes. O miRNA está acima do piso em 27 de 57 casos (47%) e em 46 de 48 controles (96%); nas demais amostras o valor é o piso (0.0001). Em log2(x+1), a distância entre "detectado" e "no piso" é de ~17–20 unidades, de modo que o β mede sobretudo a diferença na taxa de detecção. O achado é do tipo detectado/não detectado, e um modelo linear sobre log2 não é a descrição ideal dele.

## Notas
- v1 (linear), v2 (log2 sem covariáveis) e v3 (equações normais) estão no histórico do script; esta v4 é a análise vigente.
- O artigo original usou Partek GSA com idade como covariável sobre 101 miRNAs filtrados;
  aqui o universo é o dos 865 testes estimáveis (ou os 269 detectáveis) — universos diferentes, declarados.
- A lista ajustada NÃO é 'assinatura': a assinatura validada do artigo são 3 miRNAs com qPCR.
- Sexo codificado M=1; RIN como qualidade de RNA; modelo linear padrão, sem interações.
- Grupo e sexo são validados contra as categorias esperadas; idade e RIN precisam ser numéricos. Valor inesperado interrompe a análise.