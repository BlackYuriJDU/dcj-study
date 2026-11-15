# Sangue × cérebro: alvos de miRNAs do sangue entre os genes ↑ do córtex frontal
*`analise_integrada_cerebro_sangue.py` em 2026-11-11. Análise exploratória; nenhum dos dois artigos originais fez esta integração.*

> **Correção de 2026-11.** A versão anterior comparava todos os alvos de cada miRNA no genoma (48, 16, 11) com um universo de 800 genes, o que tornava o "esperado" e os p-valores incorretos; também aceitava alvos de rato e camundongo e incluía no FDR dois miRNAs sem nenhum alvo. A conclusão (nenhuma sobre-representação) não muda. Ver `CORRECOES.md`.

**Desenho**: miRNAs ↓ no sangue (GSE140069) → alvos validados (miRTarBase 10.0, "Functional MTI", miRNA e gene-alvo humanos) → sobre-representação entre os genes ↑ no córtex frontal (GSE160208, FDR<0,05). Teste hipergeométrico unicaudal com **N = 742** genes testáveis do painel NanoString e **n = 315** genes ↑ (126 ↓). Para cada miRNA, K = alvos que pertencem ao painel testável.

| miRNA (sangue) | alvos validados (humanos) | alvos no painel (K) | alvos ↑ (k) | esperado (K·n/N) | p | q (BH) | k mínimo para p<0,05 |
|---|---|---|---|---|---|---|---|
| hsa-miR-16-5p | 47 | 13 | 6 | 5.5 | 0.50 | 0.76 | 9 |
| hsa-miR-93-5p | 15 | 6 | 3 | 2.5 | 0.51 | 0.76 | 6 |
| hsa-let-7i-5p | 11 | 3 | 1 | 1.3 | 0.81 | 0.81 | inalcançável |
| hsa-miR-106b-3p | 0 | 0 | — | — | não testável | — | — |
| hsa-miR-500a-3p | 0 | 0 | — | — | não testável | — | — |

- Testes realizados: **3** (BH sobre esses 3). Não testáveis, fora do FDR: hsa-miR-106b-3p, hsa-miR-500a-3p — sem nenhum alvo "Functional MTI" humano no painel.
- Alvos excluídos pelo filtro de espécie do gene-alvo: hsa-miR-16-5p: Ccnd1 (rno); hsa-miR-93-5p: Sp7 (mmu).

## Direção oposta (descritivo)

Mesma conta para os genes ↓ do cérebro. Não é um controle negativo no sentido estrito (um miRNA reduzido não prevê alvos reduzidos, mas também não os exclui); serve só de comparação.

| miRNA | alvos no painel (K) | alvos ↓ | p (hipergeométrico, ↓) |
|---|---|---|---|
| hsa-miR-16-5p | 13 | 1 | 0.91 |
| hsa-miR-93-5p | 6 | 1 | 0.67 |
| hsa-let-7i-5p | 3 | 1 | 0.43 |

## Leitura

- **Resultado: nenhuma sobre-representação detectada.** Nos 3 miRNAs testáveis, o número de alvos ↑ fica próximo do esperado ao acaso (p de 0.50 a 0.81).
- **O teste tem pouco poder.** O painel mede 742 genes testáveis, dos quais 315 (42%) estão ↑; cada miRNA tem só 3 a 13 alvos no painel. A última coluna da tabela mostra quantos alvos ↑ seriam necessários para p < 0,05: com tão poucos alvos e uma taxa de base tão alta, só um enriquecimento quase total seria detectável.
- **O que não se pode concluir.** Um resultado nulo com esse poder não apoia nem contraria uma relação entre os miRNAs do sangue e a expressão cerebral. Não dá base para afirmar mecanismo, nem para afirmar que os miRNAs refletem processos periféricos.
- **Viés do painel.** O painel NanoString é dirigido a neuroinflamação: os genes medidos não são uma amostra do genoma, e o resultado vale só para esse universo.
