# Fine-mapping v2 — LD real (Ensembl/1000G phase 3) + credible sets descritivos
*`finemap_ld.py` em 2026-11-11. População: 1000GENOMES:phase_3:ALL. Método: ABF de Wakefield (W=0.04) + agrupamento por r²≥0.8 ao lead. Approximate — sem modelo conjunto (SuSiE exigiria genótipos individuais).*

## Locus STX6 — lead chr1:180,961,245 p=7.51e-09 · âncora LD: rs11586493 (efetivo: rs11586493)
- Variantes na janela: 162 · pares LD da âncora: 184
- Cobertura do painel nas top-20: 20/20 · máx r² observado: 1.00
- Massa posterior do cluster âncora+proxies:
  **r²≥0.8: 90.5%** · r²≥0.50: 90.5% · sem LD/fora do painel: 9.5%

| pos | alelos | p | beta | rsID | status vs âncora | posterior |
|---|---|---|---|---|---|---|
| 180,961,245 | G>A | 7.51e-09 | -0.149 | rs11586493 | âncora | 5.8% |
| 180,956,015 | A>G | 7.73e-09 | -0.149 | rs1404986578/rs1389132786/rs7553330 | r²=0.99 | 5.6% |
| 180,957,962 | T>G | 8.30e-09 | -0.149 | rs1293366509/rs74227337/rs61433244/rs12754041 | r²=0.99 | 5.2% |
| 180,958,946 | A>G | 9.02e-09 | -0.149 | rs4366283/rs4111520 | r²=0.99 | 4.8% |
| 180,949,780 | T>G | 9.09e-09 | -0.149 | rs6680541 | r²=0.97 | 4.8% |
| 180,962,282 | A>G | 9.34e-09 | -0.148 | rs2525668307/rs1654947403/rs59443232/rs1411478 | r²=1.00 | 4.7% |
| 180,956,985 | G>A | 9.60e-09 | -0.148 | rs58087663/rs56685500/rs6425658 | r²=0.99 | 4.6% |
| 180,956,905 | A>G | 9.60e-09 | -0.148 | rs57262540/rs6425657 | r²=0.99 | 4.6% |
| 180,952,516 | C>A | 9.71e-09 | -0.149 | rs58426280/rs57054063/rs56521928/rs12744212 | r²=0.99 | 4.5% |
| 180,953,038 | A>C | 9.73e-09 | -0.148 | rs1014849176/rs57149666/rs7543927 | r²=0.99 | 4.5% |
| 180,953,853 | A>G | 9.74e-09 | -0.148 | rs1654776103/rs3747957 | r²=0.99 | 4.5% |
| 180,954,130 | C>T | 9.74e-09 | -0.148 | rs961691834/rs4652548 | r²=0.99 | 4.5% |

## Locus GAL3ST1 — lead chr22:30,950,360 p=6.18e-10 · âncora LD: rs386462923 (efetivo: rs8142452)
- Variantes na janela: 322 · pares LD da âncora: 61
- Cobertura do painel nas top-20: 4/20 · máx r² observado: 1.00
- Massa posterior do cluster âncora+proxies:
  **r²≥0.8: 0.0%** · r²≥0.50: 0.0% · sem LD/fora do painel: 100.0%
- ⚠️ Leitura: o único r²=1.00 é a própria âncora (posterior ~0.0%); o sinal real está em variantes mal marcadas pelo painel de comuns (haplótipo provavelmente de baixa frequência).

| pos | alelos | p | beta | rsID | status vs âncora | posterior |
|---|---|---|---|---|---|---|
| 30,950,360 | T>C | 6.18e-10 | -0.169 | rs2517624666/rs58431710/rs2267158 | r²=0.45 | 57.8% |
| 30,953,295 | C>T | 8.60e-10 | -0.168 | rs2040917680/rs61593263/rs17858302/rs17856591/rs17845430 | r²=0.38 | 42.2% |
| 30,949,820 | A>G | 3.28e-05 | -0.112 | rs386462923/rs60067782/rs8142452 | âncora | 0.0% |
| 30,916,518 | C>T | 1.86e-04 | +0.112 | rs2040410816/rs886617 | sem dado de painel | 0.0% |
| 30,915,804 | A>G | 1.90e-04 | +0.112 | rs11912737 | sem dado de painel | 0.0% |
| 30,915,186 | A>G | 1.93e-04 | +0.112 | rs2147213388/rs77170318/rs9606749 | sem dado de painel | 0.0% |
| 30,917,077 | T>C | 2.27e-04 | +0.111 | rs571918479/rs1859479 | sem dado de painel | 0.0% |
| 30,916,443 | C>T | 2.32e-04 | +0.111 | rs886616 | sem dado de painel | 0.0% |
| 30,916,380 | C>A | 2.33e-04 | +0.111 | rs1459415747/rs886615 | sem dado de painel | 0.0% |
| 30,916,877 | A>G | 2.37e-04 | +0.110 | rs2412990 | sem dado de painel | 0.0% |
| 30,912,328 | C>G | 2.56e-04 | +0.110 | rs56612612/rs9606747 | sem dado de painel | 0.0% |
| 30,911,996 | C>T | 2.58e-04 | +0.110 | rs9608967 | sem dado de painel | 0.0% |

## Locus PRNP — lead chr20:4,672,307 p=1.62e-15 · âncora LD: rs60704301 (efetivo: rs2093390)
- Variantes na janela: 337 · pares LD da âncora: 106
- Cobertura do painel nas top-20: 18/20 · máx r² observado: 1.00
- Massa posterior do cluster âncora+proxies:
  **r²≥0.8: 58.9%** · r²≥0.50: 100.0% · sem LD/fora do painel: 0.0%

| pos | alelos | p | beta | rsID | status vs âncora | posterior |
|---|---|---|---|---|---|---|
| 4,672,307 | C>T | 1.62e-15 | -0.219 | rs60704301/rs4254562/rs2093390 | âncora | 29.8% |
| 4,675,155 | A>T | 1.65e-15 | -0.218 | rs17249667/rs6037932 | r²=0.87 | 29.1% |
| 4,672,275 | A>G | 2.36e-15 | -0.218 | rs1025221938/rs6052770/rs2093391 | r²=0.73 | 20.6% |
| 4,680,251 | A>G | 2.68e-15 | -0.211 | rs52800775/rs17858648/rs17850971/rs1799990 | r²=0.68 | 18.3% |
| 4,672,816 | A>G | 3.70e-14 | -0.201 | rs57766978/rs6052771 | r²=0.71 | 1.4% |
| 4,671,225 | T>G | 7.89e-14 | -0.198 | rs58875214/rs6052769 | r²=0.55 | 0.7% |
| 4,684,286 | T>A | 1.26e-12 | -0.208 | rs6084836 | r²=0.48 | 0.0% |
| 4,677,369 | A>G | 4.26e-12 | -0.202 | rs17328364/rs6116475 | r²=0.56 | 0.0% |
| 4,675,980 | G>A | 2.39e-11 | -0.197 | rs6052772 | r²=0.70 | 0.0% |
| 4,675,589 | T>C | 3.43e-10 | -0.194 | rs17249737/rs13045348 | r²=0.40 | 0.0% |
| 4,670,510 | G>A | 3.83e-10 | -0.185 | rs4815729 | r²=0.50 | 0.0% |
| 4,671,381 | T>G | 7.35e-10 | -0.173 | rs6107515 | r²=0.46 | 0.0% |

## Nota metodológica final
- Credible set formal exige modelo conjunta (SuSiE/FINEMAP) com genótipos;
  aqui reportamos MASSA POR CLUSTER de LD — suficiente para declarar que o
  sinal é um bloco haplotípico coeso, não um mosaico de falsos independentes.
- rs3747957 (índice Brain 2025): ver relatório QC; presente com p=9.7e-9.