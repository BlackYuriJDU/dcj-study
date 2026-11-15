# QC independente — GWAS sCJD GCST90001389 (REAIS, 4.110 casos × 13.569 controles)
*`qc_gwas_gcst90001389.py` em 2026-11-11. Fonte: GWAS Catalog/EBI, GRCh37, consórcio MRC Prion Unit (PMID 32949544).*

## 1. Integridade
- Variantes processadas: **6,314,492** · Linhas malformadas: **0** (0.0000%)
- Variantes por cromossomo (1-22, X): 1:485963, 2:540025, 3:460593, 4:474861, 5:421728, 6:444501, 7:366860, 8:359438, 9:268595, 10:330488, 11:315896, 12:309747, 13:244942, 14:199871, 15:175829, 16:175809, 17:150391, 18:182103, 19:116331, 20:136760, 21:84540, 22:69221

## 2. Distribuições
- EAF: min=0.0050, max=0.5000, variantes no piso (0/1): 0
- chi2=(beta/SE)^2: mediana=0.4817
- **lambda_GC = 1.0587** (mediana / 0.454936, sobre todas as variantes) (⚠ inflação — investigar estratificação/ancestralidade)

## 3. Hits genômicos (p < 5e-8)
- Total: **41**

| p | chr | pos | OA>EA | EAF | beta | SE |
|---|---|---|---|---|---|---|
| 1.62e-15 | 20 | 4,672,307 | C>T | 0.300 | -0.219 | 0.027 |
| 1.65e-15 | 20 | 4,675,155 | A>T | 0.299 | -0.218 | 0.027 |
| 2.36e-15 | 20 | 4,672,275 | A>G | 0.299 | -0.218 | 0.027 |
| 2.68e-15 | 20 | 4,680,251 | A>G | 0.328 | -0.211 | 0.027 |
| 3.70e-14 | 20 | 4,672,816 | A>G | 0.331 | -0.201 | 0.027 |
| 7.89e-14 | 20 | 4,671,225 | T>G | 0.347 | -0.198 | 0.027 |
| 1.26e-12 | 20 | 4,684,286 | T>A | 0.258 | -0.208 | 0.029 |
| 4.26e-12 | 20 | 4,677,369 | A>G | 0.244 | -0.202 | 0.029 |
| 2.39e-11 | 20 | 4,675,980 | G>A | 0.237 | -0.197 | 0.030 |
| 3.43e-10 | 20 | 4,675,589 | T>C | 0.204 | -0.194 | 0.031 |
| 3.83e-10 | 20 | 4,670,510 | G>A | 0.248 | -0.185 | 0.030 |
| 6.18e-10 | 22 | 30,950,360 | T>C | 0.325 | -0.169 | 0.027 |
| 7.35e-10 | 20 | 4,671,381 | T>G | 0.287 | -0.173 | 0.028 |
| 7.95e-10 | 20 | 4,668,328 | T>G | 0.252 | -0.182 | 0.030 |
| 8.60e-10 | 22 | 30,953,295 | C>T | 0.314 | -0.168 | 0.027 |
| 2.79e-09 | 20 | 4,669,534 | G>A | 0.280 | -0.169 | 0.028 |
| 2.80e-09 | 20 | 4,667,829 | T>C | 0.369 | -0.157 | 0.026 |
| 7.51e-09 | 1 | 180,961,245 | G>A | 0.420 | -0.149 | 0.026 |
| 7.73e-09 | 1 | 180,956,015 | A>G | 0.420 | -0.149 | 0.026 |
| 8.30e-09 | 1 | 180,957,962 | T>G | 0.419 | -0.149 | 0.026 |
| 9.02e-09 | 1 | 180,958,946 | A>G | 0.419 | -0.149 | 0.026 |
| 9.09e-09 | 1 | 180,949,780 | T>G | 0.420 | -0.149 | 0.026 |
| 9.34e-09 | 1 | 180,962,282 | A>G | 0.420 | -0.148 | 0.026 |
| 9.60e-09 | 1 | 180,956,985 | G>A | 0.419 | -0.148 | 0.026 |
| 9.60e-09 | 1 | 180,956,905 | A>G | 0.419 | -0.148 | 0.026 |
| 9.71e-09 | 1 | 180,952,516 | C>A | 0.420 | -0.149 | 0.026 |
| 9.73e-09 | 1 | 180,953,038 | A>C | 0.420 | -0.148 | 0.026 |
| 9.74e-09 | 1 | 180,953,853 | A>G | 0.419 | -0.148 | 0.026 |
| 9.74e-09 | 1 | 180,954,130 | C>T | 0.419 | -0.148 | 0.026 |
| 9.77e-09 | 1 | 180,954,089 | A>G | 0.419 | -0.148 | 0.026 |

## 3b. Os três loci publicados, reproduzidos

O artigo original (Jones et al., Lancet Neurol 2020) reporta três loci com significância genômica: PRNP, STX6 e GAL3ST1. Melhor variante deste QC em cada um:

| Locus | Melhor variante (GRCh37) | p | hits p<5e-8 no locus | Anotação |
|---|---|---|---|---|
| **PRNP** | chr20:4,672,307 | 1.62e-15 | 16 | região do gene da proteína priônica (chr20p13) |
| **STX6** | chr1:180,961,245 | 7.51e-09 | 23 | dentro do gene STX6 (Ensembl GRCh37: 180,941,861–180,992,047) |
| **GAL3ST1** | chr22:30,950,360 | 6.18e-10 | 2 | borda 5' do gene GAL3ST1 (Ensembl GRCh37: 30,950,622–30,970,574) |

**3 de 3 loci publicados reproduzidos** a partir das estatísticas-resumo depositadas. É uma conferência dos dados públicos, não uma replicação em coorte independente.

## 4. Locus STX6 (candidato a replicação — cf. Brain 2025)
- Janela: chr1:180,850,000-181,050,000 (GRCh37)
- Melhor variante na janela: chr1:180,961,245 G>A p=7.512e-09 (beta=-0.149, SE=0.026, EAF=0.420)
- **p < 5e-8** — sinal genômico-significativo nesta coorte 2020, consistente com o artigo original (3 loci: PRNP, STX6, GAL3ST1)

## 5. Top 20 variantes por p-value

| p | chr | pos | OA>EA | EAF | beta | SE |
|---|---|---|---|---|---|---|
| 1.62e-15 | 20 | 4,672,307 | C>T | 0.300 | -0.219 | 0.027 |
| 1.65e-15 | 20 | 4,675,155 | A>T | 0.299 | -0.218 | 0.027 |
| 2.36e-15 | 20 | 4,672,275 | A>G | 0.299 | -0.218 | 0.027 |
| 2.68e-15 | 20 | 4,680,251 | A>G | 0.328 | -0.211 | 0.027 |
| 3.70e-14 | 20 | 4,672,816 | A>G | 0.331 | -0.201 | 0.027 |
| 7.89e-14 | 20 | 4,671,225 | T>G | 0.347 | -0.198 | 0.027 |
| 1.26e-12 | 20 | 4,684,286 | T>A | 0.258 | -0.208 | 0.029 |
| 4.26e-12 | 20 | 4,677,369 | A>G | 0.244 | -0.202 | 0.029 |
| 2.39e-11 | 20 | 4,675,980 | G>A | 0.237 | -0.197 | 0.030 |
| 3.43e-10 | 20 | 4,675,589 | T>C | 0.204 | -0.194 | 0.031 |
| 3.83e-10 | 20 | 4,670,510 | G>A | 0.248 | -0.185 | 0.030 |
| 6.18e-10 | 22 | 30,950,360 | T>C | 0.325 | -0.169 | 0.027 |
| 7.35e-10 | 20 | 4,671,381 | T>G | 0.287 | -0.173 | 0.028 |
| 7.95e-10 | 20 | 4,668,328 | T>G | 0.252 | -0.182 | 0.030 |
| 8.60e-10 | 22 | 30,953,295 | C>T | 0.314 | -0.168 | 0.027 |
| 2.79e-09 | 20 | 4,669,534 | G>A | 0.280 | -0.169 | 0.028 |
| 2.80e-09 | 20 | 4,667,829 | T>C | 0.369 | -0.157 | 0.026 |
| 7.51e-09 | 1 | 180,961,245 | G>A | 0.420 | -0.149 | 0.026 |
| 7.73e-09 | 1 | 180,956,015 | A>G | 0.420 | -0.149 | 0.026 |
| 8.30e-09 | 1 | 180,957,962 | T>G | 0.419 | -0.149 | 0.026 |

## Nota de honestidade científica
- QC de primeira passada: sem verificação de strand, sem imputação-info
  (coluna não existe no arquivo), sem clumping por LD (próxima rodada).
- Sem rsIDs no arquivo — coordenadas GRCh37 são a chave primária.
- lambda_GC de sumstats de caso-controle é aproximado (chi2 de z de beta/SE).
- Este QC NÃO é descoberta nova: é verificação independente documentada.