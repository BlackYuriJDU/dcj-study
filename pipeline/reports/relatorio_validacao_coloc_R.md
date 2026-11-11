# Validação da colocalização contra o pacote R `coloc`

*`crosscheck_coloc_R.py` em 2026-11-11. Compara a implementação única do projeto (`pipeline/lib/coloc_abf.py`) com a saída do `coloc.abf` real (R 4.5, coloc 5.2.3, executado em 2026-08-29) guardada em `tests/fixtures/`.*

## O que aconteceu (registro da correção de 2026-11)

A versão deste relatório publicada até a v0.4 concluía que "a máquina está correta — acordo a 6 decimais" e que STX6 não colocalizava (H3 = 0,994997). As duas conclusões estavam erradas:

1. **A parte (A) não era uma validação.** Ela reimplementava em R a mesma fórmula do Python, com o mesmo erro: o termo de H4 multiplicado por p1·p2·p12 em vez de p12. Duas cópias da mesma equação errada concordam entre si; a concordância não dizia nada sobre a equação.
2. **A parte (B) já mostrava o resultado certo.** O `coloc.abf` do pacote oficial devolvia PP.H4 = 0,9796 para os mesmos 390 pares. A diferença foi explicada como "dependência de prior". Não era: com o mesmo prior de efeito do código antigo (W = 0,0016) e apenas a combinação corrigida, PP.H4 = 0,9788.
3. **Havia mais dois erros na mesma função**: o fator de Bayes usava √r onde Wakefield usa √(1−r), e W = 0,04² era rotulado "coloc padrão" (o padrão do `coloc` é variância 0,04 para caso-controle e 0,0225·sdY² para traço quantitativo). Mudam pouco o H4 desta região, mas eram erros de fórmula.

A fórmula antiga foi removida de todo o código. O que valida a implementação atual são (i) este oráculo, (ii) casos sintéticos conferidos contra a enumeração explícita das configurações causais e (iii) as âncoras numéricas da auditoria — todos em `tests/test_coloc_abf.py`.

| Configuração, entrada de 390 variantes da v0.4 | H2 | H3 | H4 |
|---|---|---|---|
| Publicada até a v0.4 (três erros; valores históricos) | 0,005002 | 0,994997 | 4,6×10⁻⁷ |
| Fórmula correta com o prior de efeito antigo (W = 0,0016) | 0.000164 | 0.021160 | 0.978676 |
| **Conforme o `coloc` (W = 0,04 / 0,0225)** | 7.04e-06 | 0.020388 | **0.979605** |
| `coloc.abf` 5.2.3 real (registrado em 2026-08-29) | 7×10⁻⁶ | 0,0204 | 0,9796 |

## Oráculo: `coloc.abf` real × implementação do projeto

- Entrada: `tests/fixtures/stx6_crosscheck_input_v04.tsv` (390 variantes).
- Referência: `tests/fixtures/stx6_crosscheck_colocR_output_v04.tsv`, coluna `PP4_colocR` (posterior por variante, `SNP.PP.H4`).
- **Diferença máxima, variante a variante: 1.55e-15.**
- Resumo do projeto: PP.H4 = 0.9796, PP.H3 = 0.0204, PP.H2 = 7.0e-06 — igual ao que o pacote registrou (0,9796 / 0,0204 / 7×10⁻⁶).
- O `coloc.abf` foi chamado com `type = "cc"` para o GWAS e `type = "quant"`, `sdY = 1` para o eQTL-meta, com β e variância de β fornecidos. Nessa forma, N, MAF e a fração de casos não entram no fator de Bayes.

## Entrada atual (harmonização por variante completa)

- `pipeline/data/stx6_crosscheck_input.tsv`: 390 variantes; PP.H4 = 0.979605, PP.H3 = 0.020388, PP.H2 = 7.04e-06.
- Mesmas 390 posições da entrada da v0.4. O efeito do eQTL-meta mudou em 1 variante(s) (chr1:180,975,614): sítio multialélico em que a versão anterior somava efeitos de alelos diferentes.
- `coloc.abf` real executado sobre a entrada atual (`pipeline/data/stx6_crosscheck_colocR_output.tsv`; coloc.abf: coloc 5.2.3, R version 4.5.3 (2026-03-11); PP.H3=0.020388 PP.H4=0.979605; nsnps=390): **diferença máxima por variante 1.55e-15**.

## Fine-mapping: mesmo fator de Bayes

`finemap_ld.py` usa a mesma função de fator de Bayes de Wakefield da colocalização (`coloc_abf.abf`, W = 0,04). A forma já estava correta nesse script; a troca pela função única não altera os resultados.

## Reprodução

- `python3 pipeline/scripts/crosscheck_coloc_R.py --so-oraculo` (sem GWAS e sem R)
- `python3 -m pytest tests/test_coloc_abf.py -q`
- Opcional, com R e `coloc` instalados: `Rscript pipeline/scripts/coloc_crosscheck_R.R`
