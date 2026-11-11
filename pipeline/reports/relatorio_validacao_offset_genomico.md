# Validação do offset GRCh38→GRCh37 usado na colocalização de STX6
*`valida_offset_grch38_grch37.py` em 2026-11-11. Substitui a suposição de −30.864 bp por uma medição contra a API de assembly map do Ensembl.*

A colocalização de STX6 — o resultado mecanístico central do manuscrito — converte todas as posições eQTL de GRCh38 para GRCh37 somando uma constante de −30.864 bp. Até 2026-10-04 essa constante nunca tinha sido confrontada com uma fonte independente, e o docstring do script ainda afirmava que fora "validado por Ensembl MAP", o que era falso.

## 1. O offset é constante no span de STX6?

Amostra de **41** pontos no span de STX6 em GRCh37 (chr1:180.941.861–180.992.047). convertidos pela API de assembly map do Ensembl.

| offset observado (bp) | pontos |
|---|---|
| +30,864 | 41 ← **offset do código** |

## 2. As variantes da colocalização mudam de posição?

Para cada posição eQTL de STX6 (amostra espaçada dos 5 datasets), comparamos o casamento com o GWAS sob a constante de −30.864 bp contra o casamento confirmado pelo Ensembl.

| dataset | eQTL de STX6 | amostradas | casadas pela constante | casadas pelo Ensembl |
|---|---|---|---|---|
| QTD000051_stx6 | 1167 | 13 | 3 | 3 |
| QTD000075_stx6 | 887 | 13 | 4 | 4 |
| QTD000166_stx6 | 606 | 13 | 6 | 6 |
| QTD000176_stx6 | 582 | 13 | 11 | 11 |
| QTD000434_stx6 | 711 | 13 | 6 | 6 |

## 3. Veredito

**O offset −30.864 bp está correto e é constante.** O Ensembl devolve exatamente esse delta em 41/41 pontos amostrados do span de STX6, e as variantes que casam com o GWAS sob a constante são exatamente as mesmas sob o mapeamento verificado (30/30).

**Não há impacto no resultado da colocalização.** A limitação registrada no manuscrito é real como prudência editorial, mas hoje ela subdeclara: a conversão pode ser descrita como *verificada contra a API de assembly map do Ensembl*, o que é mais forte do que "não validada" e mais honesto do que o antigo "validado por Ensembl MAP" (que afirmava uma validação de cadeia inexistente).

## 4. Nota metodológica

Uma implementação própria de parser de chain file foi escrita e descartada antes desta versão: ela mapeava b38 180.972.725 para b37 30.545.990, enquanto o Ensembl devolve 180.941.861 — um erro constante de 150.426.735 bp, compatível com corrupção no acúmulo dos blocos do chain. A API de assembly map do Ensembl é usada como fonte primária por ser autoritativa e verificável.
