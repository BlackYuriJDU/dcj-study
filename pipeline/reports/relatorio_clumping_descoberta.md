# Clumping de descoberta — o pipeline acha os loci do zero?
*`clumping_descoberta.py` em 2026-11-11 03:28. Varredura cega de 6,314,492 variantes; limiar p<1e-05; clumping guloso por
distância ±500 kb (lead = menor p do cluster).*

**Resultado: 35 clusters independentes.**

| # | Crom | Lead (b37) | p do lead | Variantes no cluster |
|---|---|---|---|---|
| 1 | 20 | 4,672,307 | 1.62e-15 | 18 |
| 2 | 22 | 30,950,360 | 6.18e-10 | 3 |
| 3 | 1 | 180,961,245 | 7.51e-09 | 50 |
| 4 | 16 | 15,539,902 | 5.73e-08 | 23 |
| 5 | 22 | 29,983,139 | 2.26e-07 | 8 |
| 6 | 7 | 92,659,158 | 3.00e-07 | 17 |
| 7 | 15 | 66,915,596 | 1.39e-06 | 7 |
| 8 | 7 | 148,700,849 | 1.66e-06 | 25 |
| 9 | 21 | 19,560,949 | 1.68e-06 | 3 |
| 10 | 10 | 57,025,110 | 3.06e-06 | 5 |
| 11 | 22 | 33,952,852 | 3.24e-06 | 1 |
| 12 | 6 | 161,746,022 | 3.46e-06 | 8 |
| 13 | 1 | 57,633,496 | 4.73e-06 | 6 |
| 14 | 3 | 49,468,689 | 5.04e-06 | 2 |
| 15 | 7 | 75,445,774 | 5.21e-06 | 4 |
| 16 | 20 | 1,542,232 | 5.27e-06 | 6 |
| 17 | 6 | 96,440,475 | 5.59e-06 | 2 |
| 18 | 9 | 27,161,517 | 5.99e-06 | 1 |
| 19 | 12 | 23,189,810 | 6.07e-06 | 1 |
| 20 | 1 | 28,632,870 | 6.23e-06 | 1 |
| 21 | 15 | 70,720,034 | 6.38e-06 | 3 |
| 22 | 10 | 12,849,532 | 6.58e-06 | 4 |
| 23 | 18 | 25,959,904 | 6.83e-06 | 2 |
| 24 | 7 | 50,454,652 | 7.29e-06 | 2 |
| 25 | 17 | 32,003,819 | 7.36e-06 | 1 |
| 26 | 14 | 87,346,119 | 7.39e-06 | 1 |
| 27 | 15 | 33,410,453 | 7.95e-06 | 1 |
| 28 | 13 | 35,022,460 | 8.53e-06 | 1 |
| 29 | 1 | 231,685,298 | 8.58e-06 | 1 |
| 30 | 3 | 116,381,899 | 8.64e-06 | 1 |
| 31 | 2 | 79,163,654 | 8.67e-06 | 1 |
| 32 | 1 | 243,083,639 | 8.98e-06 | 1 |
| 33 | 13 | 104,138,193 | 9.15e-06 | 1 |
| 34 | 8 | 96,256,866 | 9.42e-06 | 1 |
| 35 | 4 | 40,451,693 | 9.67e-06 | 1 |

## Verificação pós-hoc dos loci conhecidos

| Locus (literatura) | Posição b37 | Redescoberto? | Cluster | p do cluster |
|---|---|---|---|---|
| STX6 | 180,961,245 | **SIM** | #3 (1:180,961,245) | 7.51e-09 |
| PRNP | 4,672,307 | **SIM** | #1 (20:4,672,307) | 1.62e-15 |
| GAL3ST1 | 30,950,360 | **SIM** | #2 (22:30,950,360) | 6.18e-10 |

**3/3 loci conhecidos redescobertos às cegas.**
- Clusters extras além dos conhecidos: sinais novos a investigar —
  com p<1e-5 mas tipicamente abaixo de GWS (5e-8); nenhum deve ser
  chamado de 'novo loci' sem replicação. Reportados por transparência.

## Leitura honesta
- Clumping por distância é conservador (PLINK usa LD real); para os
  3 loci conhecidos a conectividade por LD já foi demonstrada no
  fine-mapping (relatorio_finemap_loci.md).
- A prova aqui é de SENSIBILIDADE do pipeline: dado o sumstats bruto,
  os loci principais emergem sem qualquer âncora externa.