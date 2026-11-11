# Colocalização STX6 — meta-eQTL de cinco conjuntos cerebrais × GWAS sCJD
*`coloc_meta_stx6.py` em 2026-11-11. Meta-análise por variância inversa (efeito fixo) de 5 conjuntos do eQTL Catalogue r8 (2.156 amostras no total; proveniência e checksums em `pipeline/data/proveniencia_qtl.tsv`) contra o GWAS de sCJD (GCST90001389. n = 17.679). Colocalização: implementação única em `pipeline/lib/coloc_abf.py`. conforme o pacote R `coloc`.*

> **Correção de 2026-11.** As versões até a v0.4 deste relatório davam PP.H4 ≈ 0 e H3 ≈ 0,995. O número vinha de um erro de fórmula: o termo de H4 era multiplicado por p1·p2·p12 em vez de p12 (um fator de 10⁻⁸ com os priors padrão). Com a combinação correta e os mesmos dados, PP.H4 ≈ 0,98. Ver `CORRECOES.md` e `relatorio_validacao_coloc_R.md`.

## Conjuntos de eQTL

| Conjunto | n | linhas no cache | variantes distintas | linhas repetidas (rsID) | posições multialélicas | indels | sdY estimado |
|---|---|---|---|---|---|---|---|
| CommonMind DLPFC (QTD000075) | 586 | 887 | 779 | 108 | 29 | 181 | 0.80 |
| ROSMAP DLPFC (QTD000434) | 560 | 711 | 623 | 88 | 18 | 137 | 0.68 |
| BrainSeq DLPFC (QTD000051) | 479 | 1167 | 1054 | 113 | 38 | 221 | 0.67 |
| GTEx_v10 DLPFC (QTD000176) | 268 | 582 | 549 | 33 | 0 | 44 | 0.81 |
| GTEx_v10 cerebelo (QTD000166) | 263 | 606 | 573 | 33 | 0 | 45 | 1.00 |

Os tamanhos amostrais vêm da metadata oficial r8 (GTEx: 268 e 263; total 2.156). Eles não entram no fator de Bayes, que é calculado a partir de β e erro-padrão. *sdY estimado* = √mediana(EP²·2n·MAF·(1−MAF)), a mesma estimativa do `coloc`; valores abaixo de 1 são esperados porque o EP reflete a variância residual depois das covariáveis. A análise usa sdY = 1 (expressão normalizada antes do mapeamento); a tabela de sensibilidade a W abaixo cobre o intervalo correspondente.

## Harmonização

- Variantes distintas (posição+REF+ALT) em algum conjunto: 1111 · presentes em menos de 3 conjuntos (fora da meta): 437 · **na meta: 674** (em 650 posições)
- Posições em que conjuntos diferentes trazem pares REF/ALT diferentes: 40. A versão anterior indexava só pela posição e somava esses alelos; agora cada alelo é uma variante própria e só alelos idênticos são combinados.
- Casadas com o GWAS pelos dois alelos: **390** (mesmo alelo de efeito: 390; sinal do GWAS invertido: 0; palindrômicas entre as casadas: 59)
- Não usadas: sem variante do GWAS na posição: 284 · alelos incompatíveis: 0 · só casariam na fita complementar: 0

## Posteriores (priors padrão do `coloc`: p1 = p2 = 10⁻⁴, p12 = 10⁻⁵)

Prior de efeito: variância 0.04 para o GWAS caso-controle e 0.0225 para o eQTL quantitativo (sdY = 1).

| Hipótese | Descrição | Posterior |
|---|---|---|
| H0 | nenhuma associação na região | 2.81e-46 |
| H1 | associação só com o GWAS de sCJD | 8.52e-43 |
| H2 | associação só com o eQTL de STX6 | 7.04e-06 |
| H3 | os dois associados, variantes causais distintas | 0.0204 |
| **H4** | os dois associados, uma variante causal compartilhada | **0.9796** |

- **PP.H4 = 0.9796** · PP.H3 = 0.0204 · n = 390 variantes
- H3 + H4 (os dois traços têm sinal na região) = 1.0000
- H4/(H3+H4) = 0.9796

## Sensibilidade a p12 (p1 e p2 fixos em 10⁻⁴)

p12 é o único prior que separa H4 de H3; é a sensibilidade que importa.

| p12 (p1 = p2 = 10⁻⁴) | p12/(p1·p2) | H3 | H4 | H4/(H3+H4) |
|---|---|---|---|---|
| 1e-05 (padrão do coloc) | 1000 | 0.0204 | 0.9796 | 0.9796 |
| 5e-06 | 500 | 0.0400 | 0.9600 | 0.9600 |
| 1e-06 | 100 | 0.1723 | 0.8277 | 0.8277 |

O cenário p1 = p2 = 10⁻⁶, p12 = 10⁻⁷, que versões anteriores chamavam de "conservador", dá H4 = 0.9998 e H3 = 0.0002. Ele **não** é mais rigoroso para H4: a razão p12/(p1·p2) sobe de 10³ para 10⁵, o que favorece a variante compartilhada. Fica aqui só para registro e não é usado como evidência.

## Sensibilidade ao prior de efeito W

| W (igual nos dois traços) | H3 | H4 |
|---|---|---|
| 0.0016 | 0.0212 | 0.9787 |
| 0.01 | 0.0205 | 0.9795 |
| 0.0225 | 0.0204 | 0.9796 |
| 0.04 | 0.0204 | 0.9796 |
| 0.1 | 0.0204 | 0.9796 |

## Sensibilidade ao modelo de meta-análise e à composição

| Variação | H3 | H4 | n |
|---|---|---|---|
| Efeito fixo (análise principal) | 0.0204 | 0.9796 | 390 |
| Efeitos aleatórios (DerSimonian–Laird) | 0.0245 | 0.9755 | 390 |
| Sem CommonMind DLPFC (4 conjuntos, mínimo 3) | 0.0206 | 0.9794 | 365 |
| Sem ROSMAP DLPFC (4 conjuntos, mínimo 3) | 0.0205 | 0.9795 | 380 |
| Sem BrainSeq DLPFC (4 conjuntos, mínimo 3) | 0.0205 | 0.9795 | 380 |
| Sem GTEx_v10 DLPFC (4 conjuntos, mínimo 3) | 0.0204 | 0.9796 | 373 |
| Sem GTEx_v10 cerebelo (4 conjuntos, mínimo 3) | 0.0205 | 0.9795 | 374 |
| Sem os dois GTEx (3 conjuntos, mínimo 3) | 0.0205 | 0.9795 | 363 |
| Só CommonMind DLPFC (variantes da análise principal) | 0.0218 | 0.9782 | 389 |
| Só ROSMAP DLPFC (variantes da análise principal) | 0.0210 | 0.9790 | 365 |
| Só BrainSeq DLPFC (variantes da análise principal) | 0.0208 | 0.9792 | 374 |
| Só GTEx_v10 DLPFC (variantes da análise principal) | 0.0463 | 0.7015 | 371 |
| Só GTEx_v10 cerebelo (variantes da análise principal) | 0.0513 | 0.5672 | 380 |

Heterogeneidade entre conjuntos nas 390 variantes: Q de Cochran com p < 0,05 em 36; I² > 50% em 51. Os dois tecidos do GTEx compartilham doadores e não são coortes independentes; sozinhos, cada um dá H4 bem menor (linhas "Só GTEx").

## Variantes

- Mesma direção GWAS × eQTL-meta: **348/390 (89%)**
- Conjunto de 95% de credibilidade para a variante compartilhada (dado H4): 17 variantes. A de maior posterior (chr1:180,951,692, SNP.PP.H4 = 0.1048) **não** é o lead do GWAS.

| pos b37 | REF>ALT | p GWAS | p eQTL-meta | k | I² | SNP.PP.H4 |
|---|---|---|---|---|---|---|
| 180,951,692 | A>G | 9.99e-09 | 3.44e-47 | 5 | 44% | 0.1048 |
| 180,956,015 | A>G | 7.73e-09 | 6.93e-47 | 5 | 47% | 0.0672 |
| 180,956,985 | G>A | 9.60e-09 | 5.79e-47 | 5 | 45% | 0.0656 |
| 180,961,245 | G>A | 7.51e-09 | 7.60e-47 | 5 | 48% | 0.0634 |
| 180,952,516 | C>A | 9.71e-09 | 6.17e-47 | 5 | 45% | 0.0610 |
| 180,949,780 | T>G | 9.09e-09 | 6.60e-47 | 5 | 47% | 0.0600 |
| 180,958,946 | A>G | 9.02e-09 | 6.82e-47 | 5 | 48% | 0.0592 |
| 180,956,905 | A>G | 9.60e-09 | 6.75e-47 | 5 | 47% | 0.0560 |
| 180,953,038 | A>C | 9.73e-09 | 6.86e-47 | 5 | 47% | 0.0545 |
| 180,954,130 | C>T | 9.74e-09 | 6.86e-47 | 5 | 47% | 0.0544 |

Lead do GWAS (chr1:180,961,245): p GWAS = 7.51e-09; eQTL-meta z = -14.37, p = 7.60e-47, k = 5, I² = 48.4%.

## Leitura

- **O que a análise sustenta.** Sob os pressupostos do modelo (no máximo um sinal causal por traço na região) e com os priors padrão do `coloc`, a posterior de uma variante causal compartilhada entre o risco de sCJD e a expressão de STX6 no cérebro é 0.98.
- **Sensibilidade a p12.** Com p12 = 10⁻⁶ a posterior cai para 0.83. O resultado depende desse prior e deve ser lido junto com a tabela.
- **Robustez numérica.** O valor muda pouco com o prior de efeito, com efeitos aleatórios, com a retirada de qualquer conjunto e com a retirada dos dois GTEx (tabelas acima).
- **Quanto vem do prior.** O bloco está em LD quase perfeito. Se as 390 variantes fossem equivalentes, a razão H4/H3 a priori seria p12/(p1·p2·(n−1)) ≈ 2.6; a razão a posteriori observada é ≈ 48. Parte da separação entre H3 e H4 vem do prior, e os dados a ampliam em cerca de 19 vezes.
- **O que a análise não sustenta.** Colocalização não é mediação: não mostra que a expressão de STX6 esteja no caminho causal entre a variante e a doença. O modelo não testa mais de um sinal causal por traço; o LD entre variantes não é modelado; os tecidos são DLPFC e cerebelo de adultos, não todas as regiões afetadas na DCJ.
- **Não é resultado novo.** A associação entre o locus de STX6 e a expressão do gene em sCJD já foi publicada (Jones et al. 2020; estudo multiômico em *Brain* 2025). Este relatório reproduz esse resultado com dados públicos.
