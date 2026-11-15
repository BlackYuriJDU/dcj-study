# Colocalização GWAS sCJD × sQTL (splicing) — região STX6
*`coloc_sqtl_stx6.py` em 2026-11-11. eQTL Catalogue r8, leafcutter, mesmos cinco conjuntos da meta-eQTL (proveniência em `pipeline/data/proveniencia_qtl.tsv`). Os arquivos `cc` contêm só os clusters com sinal cis significativo: ausência de um cluster é "sem sQTL significativo detectado", não "sem sQTL". Colocalização por conjunto × cluster com a implementação única (`pipeline/lib/coloc_abf.py`, priors padrão do `coloc`).*

> **Correção de 2026-11.** Os PP.H4 ≈ 0 da versão anterior vinham do mesmo erro de fórmula da colocalização de eQTL, e o limiar "Bonferroni" aplicado a PP.H4 foi retirado (Bonferroni corrige p-valores, não probabilidades a posteriori). Ver `CORRECOES.md`.

| Conjunto | n | clusters significativos na janela | clusters de STX6 |
|---|---|---|---|
| CommonMind DLPFC (QTD000079) | 586 | 1 | **0** |
| ROSMAP DLPFC (QTD000438) | 560 | 1 | **0** |
| BrainSeq DLPFC (QTD000055) | 479 | 2 | **0** |
| GTEx_v10 DLPFC (QTD000180) | 268 | 0 | **0** |
| GTEx_v10 cerebelo (QTD000170) | 263 | 2 | **0** |

Testes de colocalização executados (conjunto × cluster com ≥ 20 variantes casadas com o GWAS): **6**

| conjunto | gene | cluster (b38) | variantes | H1 | H2 | H3 | **H4** | H4 com p12 = 10⁻⁶ | direção concordante | p mín sQTL |
|---|---|---|---|---|---|---|---|---|---|---|
| GTEx_v10 cerebelo | ENSG00000135835 | 180936513:180938555:+ | 383 | 0.0762 | 0.0001 | 0.3258 | **0.5979** | 0.1294 | 21% | 1.6e-07 |
| CommonMind DLPFC | ENSG00000135835 | 180936513:180938555:+ | 389 | 0.4522 | 0.0001 | 0.3555 | **0.1920** | 0.0232 | 35% | 7.3e-06 |
| BrainSeq DLPFC | ENSG00000135835 | 180936513:180938555:+ | 376 | 0.6020 | 0.0001 | 0.3508 | **0.0468** | 0.0049 | 32% | 7.3e-06 |
| GTEx_v10 cerebelo | ENSG00000230124 | 180270423:180273338:- | 383 | 0.9154 | 1.24e-05 | 0.0375 | **0.0467** | 0.0049 | 48% | 2.9e-03 |
| BrainSeq DLPFC | ENSG00000116260 | 180175366:180178794:+ | 376 | 0.9339 | 8.95e-06 | 0.0271 | **0.0387** | 0.0040 | 48% | 6.4e-03 |
| ROSMAP DLPFC | ENSG00000135835 | 180936513:180938555:- | 365 | 0.0012 | 0.0003 | 0.9984 | **6.74e-05** | 6.74e-06 | 62% | 3.8e-09 |

## Leitura

- **Clusters de sQTL do próprio STX6: 0** nos cinco conjuntos. Nenhum sQTL significativo de STX6 foi detectado, em contraste com o eQTL de expressão, que é forte (p ≈ 7×10⁻⁴⁷ na meta). Isso é ausência de detecção em arquivos que só listam pares significativos, não evidência de ausência.
- Os 6 testes são de clusters de genes vizinhos (ENSG00000116260, ENSG00000135835, ENSG00000230124), não de STX6. PP.H4 > 0,8 em 0; entre 0,5 e 0,8 em 1.
- O maior PP.H4 é 0.60 (GTEx_v10 cerebelo, ENSG00000135835); cai para 0.13 com p12 = 10⁻⁶, tem direção concordante em só 21% das variantes e não se repete nos outros conjuntos que têm o mesmo cluster. Não é tratado como colocalização.
- Estes testes existem porque os clusters foram pré-selecionados por terem sQTL significativo; são exploratórios, dizem respeito a outros genes e não entram nas conclusões sobre STX6.
- Conclusão: o locus de STX6 associa-se à expressão do gene; para splicing de STX6 não há sinal detectado nestes dados, com poder limitado.
