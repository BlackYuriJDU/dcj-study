# GSE140069 — teste empírico do confundimento idade/RIN
*`analise_gse140069_idade.py` em 2026-11-11. Fase 2 da auditoria pré-submissão, revista na correção de 2026-11 (solver por QR, universo do FDR declarado, leitura descritiva). Objetivo: ver o que acontece com o número de miRNAs significativos quando cada covariável entra no modelo.*

## 1. A amostra

| covariável | sCJD (n=57) | controles (n=48) | diferença | d |
|---|---|---|---|---|
| idade (anos) | 66.4 ± 8.1 | 53.6 ± 15.3 | **+12.8 anos** | +1.07 |
| RIN | 5.59 ± 1.32 | 6.50 ± 1.26 | **-0.92** | -0.71 |

Idade — casos: Q1=63, mediana=66, Q3=72 · controles: Q1=47, mediana=58, Q3=64.

## 2. Colinearidade do bloco (VIF)

| preditor | VIF | leitura |
|---|---|---|
| grupo | 1.37 | ok |
| idade | 1.32 | ok |
| sexo | 1.01 | ok |
| rin | 1.15 | ok |

VIF ≥ 5 indica que o preditor é mal identificado no modelo — a atribuição de efeito entre covariáveis colineares não é estável.

## 3. Cada covariável isolada

Três universos de FDR, declarados: **estimáveis** (865 testes; universo principal), **939 (legado)** (todos os miRNAs medidos; os 74 não estimáveis contam como testes que nunca rejeitam — é o universo dos números publicados até a v0.4) e **filtrado** (miRNAs que passam no filtro de detecção). FDR<0.05 dentro de cada modelo e universo. A coluna β mediana é a do coeficiente de interesse.

| modelo | preditores | universo | testes | coef. | FDR<0.05 | β mediana |
|---|---|---|---|---|---|---|
| M1 | idade sozinha | estimáveis | 865 | idade | **121** | -0.0144 |
| M1 | idade sozinha | 939 (legado) | 939 | idade | **113** | -0.0144 |
| M1 | idade sozinha | filtrado | 269 | idade | **100** | -0.0646 |
| M2 | RIN sozinho | estimáveis | 865 | RIN | **93** | +0.1248 |
| M2 | RIN sozinho | 939 (legado) | 939 | RIN | **91** | +0.1248 |
| M2 | RIN sozinho | filtrado | 269 | RIN | **112** | +0.7153 |
| M3 | grupo sozinho | estimáveis | 865 | grupo | **92** | -0.4034 |
| M3 | grupo sozinho | 939 (legado) | 939 | grupo | **88** | -0.4034 |
| M3 | grupo sozinho | filtrado | 269 | grupo | **93** | -1.6473 |
| M4 | grupo + idade | estimáveis | 865 | grupo | **4** | -0.3943 |
| M4 | grupo + idade | 939 (legado) | 939 | grupo | **4** | -0.3943 |
| M4 | grupo + idade | filtrado | 269 | grupo | **8** | -1.0401 |
| M4b | grupo + RIN | estimáveis | 865 | grupo | **28** | -0.3957 |
| M4b | grupo + RIN | 939 (legado) | 939 | grupo | **26** | -0.3957 |
| M4b | grupo + RIN | filtrado | 269 | grupo | **45** | -1.1257 |
| M4c | grupo + sexo | estimáveis | 865 | grupo | **94** | -0.4012 |
| M4c | grupo + sexo | 939 (legado) | 939 | grupo | **85** | -0.4012 |
| M4c | grupo + sexo | filtrado | 269 | grupo | **90** | -1.6211 |
| M5 | grupo + idade + sexo + RIN | estimáveis | 865 | grupo | **1** | -0.3369 |
| M5 | grupo + idade + sexo + RIN | 939 (legado) | 939 | grupo | **1** | -0.3369 |
| M5 | grupo + idade + sexo + RIN | filtrado | 269 | grupo | **5** | -0.6715 |

## 4. O que muda quando cada covariável entra no modelo

O teste não ajustado do relatório principal é o Welch: **89** miRNAs a FDR<0.05 no universo de 865 estimáveis (84 no universo legado de 939 — o "84" das versões anteriores). O mesmo contraste como OLS só com grupo é outro teste (variância agrupada) e dá 92 (88). As comparações abaixo são entre modelos OLS, por isso a base é o OLS só com grupo.

| modelo | FDR<0.05 (865 estimáveis) | variação vs. grupo sozinho | FDR<0.05 (939, legado) |
|---|---|---|---|
| idade sozinha (sem grupo) | 121 | — | 113 |
| RIN sozinho (sem grupo) | 93 | — | 91 |
| grupo sozinho | 92 | base | 88 |
| grupo + sexo | 94 | +2 | 85 |
| **grupo + idade** | **4** | **-88** | 4 |
| **grupo + RIN** | **28** | **-64** | 26 |
| grupo + idade + sexo + RIN (modelo do relatório principal) | 1 | -91 | 1 |

Dos 89 do Welch: **63** (71%) também são significativos para **idade sozinha**; **43** (48%) para **RIN sozinho**; **86** para grupo sozinho (OLS).

### Leitura

**O que as contagens mostram.** Com a idade no modelo, o número de miRNAs com efeito de grupo a FDR<0.05 vai de 92 para 4; só com RIN, para 28; só com sexo, fica em 94; com as três covariáveis, 1. O sinal não ajustado não se sustenta quando a idade é incluída, e o RIN sozinho já o reduz bastante.

**O que elas não mostram.** Estas são contagens de testes que cruzam um limiar em modelos diferentes. Não medem proporção de efeito nem de variância. Versões anteriores deste relatório e do manuscrito diziam que a idade "responde por ~97% do colapso" — a conta era (88 − 4)/(88 − 1) sobre as contagens do universo legado. Essa leitura foi retirada: uma razão entre contagens de significativos não é uma fração atribuível à idade, e as reduções por covariável são **não aditivas** (−88 com idade, −64 com RIN e +2 com sexo somam 150, mais que os 91 perdidos no modelo completo), porque as covariáveis capturam variância sobreposta.

**Sobre o VIF.** VIF < 1,4 em todos os preditores indica que os coeficientes são estimáveis separadamente (não há colinearidade forte). Não transforma a comparação de contagens em atribuição causal: casos e controles diferem em idade por desenho do estudo, e um modelo linear em idade pode não capturar toda essa diferença.

## 5. Faixa etária restrita

Se casos e controles se sobrepõem, o efeito de grupo que sobrevive dentro da faixa não pode ser atribuído à idade. A coluna `gap` é a diferença de idade **dentro** da faixa — se for grande, a restrição não removeu a confusão.

**Poder insuficiente para resolver um efeito residual.** Com 22–29 controles por faixa (contra 48 na amostra total) e universo de 269 miRNAs filtrados, as contagens de 1–7 não distinguem "nada sobrevive" de "sobrevive pouco". Esta seção é um controle de sanidade, não uma medição.

| faixa | n casos | n controles | gap idade | grupo só (FDR) | grupo + idade (FDR) |
|---|---|---|---|---|---|
| 50–70 | 37 | 29 | +2.7 anos | 7 | 5 |
| 55–70 | 34 | 22 | +0.7 anos | 3 | 4 |
| 55–75 | 44 | 24 | +2.0 anos | 1 | 1 |
