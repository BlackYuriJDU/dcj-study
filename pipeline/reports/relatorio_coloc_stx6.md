# Colocalização STX6: GWAS sCJD × eQTL cerebral (GTEx v10), por tecido
*`coloc_stx6_eqtl.py` em 2026-11-11. Método: colocalização por fatores de Bayes aproximados (Giambartolomei 2014), implementação única em `pipeline/lib/coloc_abf.py`, conforme o pacote R `coloc`. Prior de efeito: variância 0.04 para o GWAS caso-controle e 0.0225 para o eQTL quantitativo (0,15²·sdY², com sdY = 1 porque a expressão é normalizada antes do mapeamento). Priors de hipótese: p1 = p2 = 10⁻⁴, p12 = 10⁻⁵ (padrão do `coloc`). GWAS n = 17.679 (GCST90001389); eQTL: eQTL Catalogue r8, caches versionados (ver `pipeline/data/proveniencia_qtl.tsv`); região chr1:180,9–181,1 Mb (GRCh37); eQTL em GRCh38 convertido com offset −30.864.*

> **Correção de 2026-11.** A versão anterior deste relatório dava PP.H4 ≈ 0 por um erro de fórmula (o termo de H4 era multiplicado por p1·p2·p12 em vez de p12). Ver `CORRECOES.md`.

## DLPFC — GTEx_v10 DLPFC (QTD000176, n = 268)

- Linhas no cache: 582 · variantes distintas (posição+REF+ALT): 549 · linhas repetidas (mesma variante, outro rsID): 33 · indels: 44 · posições multialélicas: 0 · descartadas por beta/EP conflitantes: 0
- Casadas com o GWAS pelos dois alelos: **372** (mesmo alelo de efeito: 372; sinal do GWAS invertido: 0; palindrômicas A/T ou C/G entre as casadas: 59)
- Não usadas: sem variante do GWAS na posição: 177 · alelos incompatíveis: 0 · só casariam na fita complementar: 0

### Posteriores (priors padrão do `coloc`)
| Hipótese | Descrição | Posterior |
|---|---|---|
| H0 | nenhuma associação na região | 8.32e-05 |
| H1 | associação só com o GWAS de sCJD | 0.2521 |
| H2 | associação só com o eQTL de STX6 | 1.55e-05 |
| H3 | os dois associados, variantes causais distintas | 0.0463 |
| **H4** | os dois associados, uma variante causal compartilhada | **0.7015** |

**PP.H4 = 0.7015** · H3 + H4 (os dois traços associados) = 0.7478 · H4/(H3+H4) = 0.9381 · n = 372 variantes

### Sensibilidade a p12 (p1 e p2 fixos)
| p12 (p1 = p2 = 10⁻⁴) | p12/(p1·p2) | H3 | H4 | H4/(H3+H4) |
|---|---|---|---|---|
| 1e-05 (padrão do coloc) | 1000 | 0.0463 | 0.7015 | 0.9381 |
| 5e-06 | 500 | 0.0714 | 0.5403 | 0.8833 |
| 1e-06 | 100 | 0.1257 | 0.1903 | 0.6023 |

### Sensibilidade ao prior de efeito W
| W (igual nos dois traços) | H3 | H4 |
|---|---|---|
| 0.0016 | 0.0361 | 0.2260 |
| 0.01 | 0.0423 | 0.6137 |
| 0.0225 | 0.0463 | 0.7015 |
| 0.04 | 0.0507 | 0.7160 |
| 0.1 | 0.0573 | 0.6900 |

| rsID | pos b37 | p GWAS | p eQTL | SNP.PP.H4 |
|---|---|---|---|---|
| rs11586493 | 180,961,245 | 7.51e-09 | 1.71e-03 | 0.0577 |
| rs7553330 | 180,956,015 | 7.73e-09 | 1.71e-03 | 0.0561 |
| rs12754041 | 180,957,962 | 8.30e-09 | 1.71e-03 | 0.0524 |
| rs4111520 | 180,958,946 | 9.02e-09 | 1.71e-03 | 0.0484 |
| rs6680541 | 180,949,780 | 9.09e-09 | 1.71e-03 | 0.0482 |
| rs1411478 | 180,962,282 | 9.34e-09 | 1.71e-03 | 0.0469 |
| rs6425658 | 180,956,985 | 9.60e-09 | 1.71e-03 | 0.0456 |
| rs6425657 | 180,956,905 | 9.60e-09 | 1.71e-03 | 0.0456 |
| rs12744212 | 180,952,516 | 9.71e-09 | 1.71e-03 | 0.0452 |
| rs7543927 | 180,953,038 | 9.73e-09 | 1.71e-03 | 0.0451 |

### Concordância de direção
- Efeitos na mesma direção (GWAS × eQTL, alelo de efeito = ALT): **303/372 (81%)**
- Lead do GWAS (rs11586493): β_GWAS = -0.149 · β_eQTL = -0.2189 · mesma direção

### Nota de poder
- Com n = 268 doadores, o maior |z| do eQTL de STX6 neste tecido é 3.8 (contra |z| ≈ 5,8 do lead do GWAS). O sinal de eQTL isolado é fraco, e por isso a posterior se divide entre H1 (só o GWAS) e H4. Este tecido sozinho é inconclusivo; a análise com poder é a meta de cinco conjuntos (`relatorio_coloc_meta_stx6.md`).

## CEREBELO — GTEx_v10 cerebelo (QTD000166, n = 263)

- Linhas no cache: 606 · variantes distintas (posição+REF+ALT): 573 · linhas repetidas (mesma variante, outro rsID): 33 · indels: 45 · posições multialélicas: 0 · descartadas por beta/EP conflitantes: 0
- Casadas com o GWAS pelos dois alelos: **383** (mesmo alelo de efeito: 383; sinal do GWAS invertido: 0; palindrômicas A/T ou C/G entre as casadas: 59)
- Não usadas: sem variante do GWAS na posição: 190 · alelos incompatíveis: 0 · só casariam na fita complementar: 0

### Posteriores (priors padrão do `coloc`)
| Hipótese | Descrição | Posterior |
|---|---|---|
| H0 | nenhuma associação na região | 0.0001 |
| H1 | associação só com o GWAS de sCJD | 0.3813 |
| H2 | associação só com o eQTL de STX6 | 1.71e-05 |
| H3 | os dois associados, variantes causais distintas | 0.0514 |
| **H4** | os dois associados, uma variante causal compartilhada | **0.5672** |

**PP.H4 = 0.5672** · H3 + H4 (os dois traços associados) = 0.6186 · H4/(H3+H4) = 0.9169 · n = 383 variantes

### Sensibilidade a p12 (p1 e p2 fixos)
| p12 (p1 = p2 = 10⁻⁴) | p12/(p1·p2) | H3 | H4 | H4/(H3+H4) |
|---|---|---|---|---|
| 1e-05 (padrão do coloc) | 1000 | 0.0514 | 0.5672 | 0.9169 |
| 5e-06 | 500 | 0.0717 | 0.3959 | 0.8466 |
| 1e-06 | 100 | 0.1050 | 0.1159 | 0.5247 |

### Sensibilidade ao prior de efeito W
| W (igual nos dois traços) | H3 | H4 |
|---|---|---|
| 0.0016 | 0.0359 | 0.1581 |
| 0.01 | 0.0440 | 0.4435 |
| 0.0225 | 0.0514 | 0.5670 |
| 0.04 | 0.0602 | 0.6033 |
| 0.1 | 0.0823 | 0.5870 |

| rsID | pos b37 | p GWAS | p eQTL | SNP.PP.H4 |
|---|---|---|---|---|
| rs11586493 | 180,961,245 | 7.51e-09 | 2.37e-03 | 0.0687 |
| rs7553330 | 180,956,015 | 7.73e-09 | 2.90e-03 | 0.0572 |
| rs12754041 | 180,957,962 | 8.30e-09 | 2.90e-03 | 0.0534 |
| rs4111520 | 180,958,946 | 9.02e-09 | 2.37e-03 | 0.0576 |
| rs6680541 | 180,949,780 | 9.09e-09 | 3.52e-03 | 0.0425 |
| rs1411478 | 180,962,282 | 9.34e-09 | 2.37e-03 | 0.0558 |
| rs6425658 | 180,956,985 | 9.60e-09 | 2.90e-03 | 0.0465 |
| rs6425657 | 180,956,905 | 9.60e-09 | 2.90e-03 | 0.0465 |
| rs12744212 | 180,952,516 | 9.71e-09 | 2.90e-03 | 0.0460 |
| rs7543927 | 180,953,038 | 9.73e-09 | 2.90e-03 | 0.0459 |

### Concordância de direção
- Efeitos na mesma direção (GWAS × eQTL, alelo de efeito = ALT): **298/383 (78%)**
- Lead do GWAS (rs11586493): β_GWAS = -0.149 · β_eQTL = -0.2635 · mesma direção

### Nota de poder
- Com n = 263 doadores, o maior |z| do eQTL de STX6 neste tecido é 3.8 (contra |z| ≈ 5,8 do lead do GWAS). O sinal de eQTL isolado é fraco, e por isso a posterior se divide entre H1 (só o GWAS) e H4. Este tecido sozinho é inconclusivo; a análise com poder é a meta de cinco conjuntos (`relatorio_coloc_meta_stx6.md`).

## Leitura

| Tecido | n | variantes | H1 | H3 | H4 | H4/(H3+H4) |
|---|---|---|---|---|---|---|
| DLPFC | 268 | 372 | 0.2521 | 0.0463 | 0.7015 | 0.9381 |
| CEREBELO | 263 | 383 | 0.3813 | 0.0514 | 0.5672 | 0.9169 |

- Cada tecido do GTEx, sozinho, tem pouco poder para o eQTL de STX6: a posterior fica dividida entre H1 (associação só com o GWAS) e H4 (variante causal compartilhada), com H3 pequeno. Nenhum dos dois tecidos, isolado, sustenta nem refuta a colocalização.
- Os dois tecidos do GTEx vêm, em grande parte, dos mesmos doadores: não são réplicas independentes.
- H4 alto significaria suporte a uma variante causal compartilhada sob os pressupostos do modelo (um único sinal causal por traço na região). Colocalização não testa mediação: não mostra que a expressão de STX6 esteja no caminho causal entre a variante e a doença.
- DLPFC e cerebelo não representam todas as regiões afetadas na DCJ.
