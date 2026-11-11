# Correções de novembro de 2026

Registro do que mudou no repositório em resposta à auditoria de 2026-11
(`auditoria/DOSSIE_AUDITORIA_DCJ_STUDY.md`). Cada item traz o identificador do
dossiê, o que estava errado, o valor antes e depois, e o teste ou comando que
comprova. Base: commit `e0765256` (preprint v0.4). Resultado: preprint v0.5.

## O que muda para quem leu a v0.4

**A v0.4 afirmava que o sinal de STX6 no GWAS de sCJD não podia ser atribuído a
uma variante causal compartilhada com o eQTL cerebral de STX6 (H3 ≈ 0,995). A
afirmação estava errada e vinha de um erro de fórmula no nosso código.** Com a
fórmula correta e os mesmos dados, a probabilidade a posteriori de variante
compartilhada é 0,98 com os priors padrão do `coloc` (0,83 com p12 = 10⁻⁶).

O que o repositório afirma agora: suporte a uma variante causal compartilhada
entre risco de sCJD e expressão de STX6 no cérebro, sob o pressuposto de um
único sinal causal por traço, com a sensibilidade a p12 declarada. O que ele
**não** afirma: mediação. Colocalização não mostra que a expressão de STX6
esteja no caminho causal entre a variante e a doença. E o resultado não é novo:
reproduz o que Jones et al. (2020) e o estudo multiômico de 2025 já publicaram.

Nenhuma outra conclusão do projeto se inverte. As demais correções mudam
números, universos de teste, redação e infraestrutura.

## Números que mudaram

| Resultado | Antes (v0.4) | Depois (v0.5) | Achado |
|---|---|---|---|
| STX6, meta-eQTL × GWAS: PP.H4 | 4,6×10⁻⁷ | **0,9796** | G1 |
| idem: PP.H3 | 0,994997 | 0,0204 | G1 |
| idem: PP.H2 | 0,005002 | 7,0×10⁻⁶ | G1, G3, G4 |
| PP.H4 com p12 = 5×10⁻⁶ / 10⁻⁶ | não reportado | 0,9600 / 0,8277 | G7 |
| cenário p1 = p2 = 10⁻⁶, p12 = 10⁻⁷ | H3 = 0,665 ("conservador") | H4 = 0,9998; não é conservador nem usado como evidência | G7 |
| GTEx DLPFC sozinho: PP.H4 | 0,000 | 0,7015 | G1 |
| GTEx cerebelo sozinho: PP.H4 | 0,000 | 0,5672 | G1 |
| sQTL, 6 testes em genes vizinhos: maior PP.H4 | ≈ 0 em todos | 0,60 (0,13 com p12 = 10⁻⁶); demais ≤ 0,19 | G10 |
| Prior de efeito W | 0,0016 nos dois traços | 0,04 (GWAS) e 0,0225 (eQTL) | G4 |
| n GTEx DLPFC / cerebelo | 285 / 272 | 268 / 263 | G8 |
| Total de amostras da meta-eQTL | ~2.182 | 2.156 | G8 |
| Meta-eQTL em chr1:180.975.614 (β) | −0,1853 (alelos G>T e G>C somados) | −0,1807 (só G>C, cinco conjuntos) | G6 |
| Variantes com Q de Cochran p<0,05 / I²>50% | 37 / 52 (dossiê) | 36 / 51 | G6 |
| Menor p do eQTL-meta na região | 6,60×10⁻⁴⁷ em 180.949.780 (manuscrito) | 3,44×10⁻⁴⁷ em 180.951.692 | novo (N3) |
| λGC no relatório de QC | 1,0588 | 1,0587 (mesma constante do outro relatório) | D3 |
| Hits p<5×10⁻⁸ por locus | "30 dos 41 em PRNP" | 16 PRNP, 23 STX6, 2 GAL3ST1 | novo (N1) |
| miRNA, Welch sem ajuste, FDR<0,05 | 84 (BH sobre 939) | **89** (BH sobre 865 estimáveis); 84 no universo antigo | M2 |
| miRNA, idem ∩ detectáveis | 69 | 73; 69 no universo antigo | M2 |
| miRNA, modelo ajustado | 1 | 1 (nos dois universos) | M1, M2 |
| p dos 74 miRNAs constantes no modelo ajustado | 0,0203 (ruído numérico) | sem p: não estimável | M1 |
| Idade sozinha / RIN / grupo / grupo+idade / grupo+RIN / grupo+sexo / completo | 113 / 91 / 88 / 4 / 26 / 85 / 1 | 121 / 93 / 92 / 4 / 28 / 94 / 1 (universo antigo ao lado) | M2 |
| "A idade explica ~97% do colapso" | afirmado | retirado | M4 |
| Cérebro, FDR<0,05 | 437 (314↑/123↓), BH sobre 800 | **441** (315↑/126↓), BH sobre 742 estimáveis; 437 no universo antigo | E2 |
| "447 significativos" (validação cruzada) | citado | não reproduzido; o valor para 742 genes é 441 | E2 |
| 184 DEGs do artigo | "réplica exata" | 184 com magnitude antes do FDR; 183 com FDR no universo inteiro (CCL4 na fronteira) | E1 |
| Sangue × cérebro, miR-16-5p: esperado / p | 18,8 / 1,00 | 5,5 / 0,50 | I1 |
| idem, miR-93-5p | 6,3 / 0,98 | 2,5 / 0,51 | I1 |
| idem, let-7i-5p | 4,3 / 0,996 | 1,3 / 0,81 | I1 |
| Alvos validados (miR-16-5p / miR-93-5p) | 48 / 16 (com 1 de rato e 1 de camundongo) | 47 / 15 | I2 |
| Testes no FDR da integração | 5 (dois com p = 1 sem alvo) | 3 | I3 |
| Checagem do offset no Ensembl | 9/9 pontos | 41/41 pontos (reexecutada) | — |
| Suíte de testes | 37 | 201 | T1, T2 |

Reproduzidos sem mudança, com o GWAS real: 6.314.492 variantes, 0 linhas
malformadas, 41 hits p<5×10⁻⁸, lead de STX6 chr1:180.961.245 (p = 7,51×10⁻⁹),
massa de fine-mapping do STX6 de 90,5% (painéis ALL e EUR), 35 clusters no
clumping, concordância de direção 348/390 (89%), 390 variantes casadas.

## Por achado

### Genética e colocalização

**G1 — Fórmula de H4 (corrigido).** O termo de H4 era `p1·p2·p12·Σ`; o correto
é `p12·Σ`. Implementação única em `pipeline/lib/coloc_abf.py`, em espaço
logarítmico, conforme `coloc::combine.abf`. As cópias em `coloc_stx6_eqtl.py`,
`coloc_meta_stx6.py`, `coloc_sqtl_stx6.py` e `coloc_crosscheck_R.R` foram
removidas. Comprovação: `tests/test_coloc_abf.py` (enumeração explícita das
configurações causais; com a fórmula antiga, 19 dos 36 testes falham).

**G2 — "Validação" sem independência (corrigido).** A parte (A) do script R
copiava a fórmula errada. Foi removida; o script agora só roda o `coloc.abf`
real. Dois oráculos versionados: a saída de 2026-08-29 sobre a entrada da v0.4
(`tests/fixtures/`) e a saída de 2026-11-11 (R 4.5.3, coloc 5.2.3) sobre a
entrada atual (`pipeline/data/stx6_crosscheck_colocR_output.tsv`). Diferença
máxima por variante: 1,6×10⁻¹⁵ nos dois. `relatorio_validacao_coloc_R.md` foi
reescrito e diz o que aconteceu. Comprovação:
`python pipeline/scripts/crosscheck_coloc_R.py --so-oraculo`;
`test_oraculo_snp_pp_h4_bate_com_coloc_r`;
`test_segundo_oraculo_coloc_r_sobre_a_entrada_atual`.

**G3 — Fator de Bayes de Wakefield (corrigido).** `√r` → `√(1−r)`. Uma função
única (`coloc_abf.log_abf`) serve à colocalização e ao fine-mapping
(`finemap_ld.py` já tinha a forma certa; a massa de 90,5% não muda).
Efeito isolado no H4 da entrada da v0.4: 0,978809 → 0,978676.

**G4 — Prior de efeito W (corrigido).** `0,04²` → variância 0,04 para o GWAS
caso-controle e 0,0225·sdY² para o eQTL, como no `coloc`. sdY = 1 (expressão
normalizada); a estimativa de sdY a partir dos dados fica entre 0,67 e 1,00 por
conjunto e está no relatório. Sensibilidade a W (0,0016 a 0,1): H4 de 0,9787 a
0,9796. A causa exata do `0,04²` no histórico do arquivo não foi investigada.

**G5 — Hipóteses mal descritas (corrigido).** H2 é "associação só com o traço
2", não "duas variantes distintas". As descrições estão em
`coloc_abf.HIPOTESES`; "H2+H3+H4" saiu de scripts, relatórios e manuscrito; o
rótulo `H4/(H3+H4)` agora é a conta `H4/(H3+H4)` (`coloc_abf.h4_condicional`).

**G6 — Harmonização por posição (corrigido).** `pipeline/lib/harmoniza.py`:
chave (posição, REF, ALT), indels aparados, linhas repetidas por rsID
deduplicadas, meta só de alelos idênticos, GWAS casado pelos dois alelos com
inversão de sinal, descartes contados. Com a lógica antiga, 18 das 651 posições
da meta misturavam alelos (reproduzido em
`test_logica_antiga_misturava_18_posicoes`). Agora: 674 variantes em 650
posições; 390 casadas com o GWAS (as mesmas posições); 0 inversões de sinal, 0
incompatíveis, 0 só na fita complementar. O H4 não muda (0,979605).

**G7 — Sensibilidade a priors (corrigido).** Os scripts emitem a tabela de
sensibilidade a p12 com p1 e p2 fixos e a tabela de sensibilidade a W. O
cenário antes chamado "conservador" fica registrado com a explicação de que a
razão p12/(p1·p2) sobe de 10³ para 10⁵ e favorece H4.

**G8 — Tamanhos amostrais do GTEx (corrigido).** 285/272 → 268/263 (metadata
r8 do eQTL Catalogue, commit `7a5b65b`). Os scripts leem os n de
`pipeline/data/proveniencia_qtl.tsv`. Os n não entram no fator de Bayes.

**G9 — Interpretação (corrigido).** Relatórios, §4.7 do manuscrito e READMEs
dizem: suporte a variante compartilhada sob o modelo de um sinal causal por
traço; sensibilidade a p12; robustez numérica; limites (LD não modelado,
doadores compartilhados entre os tecidos do GTEx, tecido); não é mediação; não
é resultado novo. O diagnóstico "quanto vem do prior" foi recalculado: razão
H4/H3 a priori ≈ 2,57 com 390 variantes equivalentes; a posteriori ≈ 48.

**G10 — sQTL (corrigido).** Recalculado com a implementação única; o limiar
"Bonferroni" sobre PP.H4 foi retirado. Continua valendo: nenhum cluster de
sQTL do próprio STX6 nos cinco conjuntos. Dos seis testes em genes vizinhos, o
maior PP.H4 é 0,60 (KIAA1614, GTEx cerebelo), com 21% de concordância de
direção e sem repetição nos outros conjuntos; não é tratado como
colocalização. Corrigida também a chave de cluster, que nunca era normalizada
(o identificador tem quatro campos, e o código esperava cinco).

### miRNA de sangue (GSE140069)

**M1 — OLS instável (corrigido).** `pipeline/lib/ols.py`: QR de Householder;
resposta constante, posto deficiente, falta de graus de liberdade e ajuste
perfeito viram "não estimável", sem p. Defeito reproduzido: com o solver
antigo, os 74 miRNAs constantes recebiam p = 0,0203 com float do Python e
≈ 9×10⁻⁶ com `numpy.float64` neste ambiente. Comprovação: `tests/test_ols.py`;
`test_ols_ajustado_bate_com_statsmodels` (os 865 testes).

**M2 — Universo do BH (corrigido).** Universo declarado = testes estimáveis
(865). As contagens do universo antigo (939) aparecem ao lado em todos os
relatórios. O resultado ajustado é 1 nos dois.

**M3 — Codificações frágeis (corrigido).** Grupo e sexo são validados contra
as categorias esperadas; idade e RIN têm de ser numéricos; identificadores
repetidos no FDR são erro; planilha e series matrix têm de concordar no grupo.
Comprovação: `test_categorias_inesperadas_levantam_erro`.

**M4 — "A idade explica ~97%" (corrigido).** Retirado do relatório de idade e
do manuscrito. O texto agora dá as contagens por modelo e diz que uma razão
entre contagens de significativos não é fração atribuível.

**M5 — miR-500a-3p (verificado e descrito).** Acima do piso em 27 de 57 casos e
em 46 de 48 controles: o β de −7,9 mede a taxa de detecção. Dito no relatório
e no manuscrito.

### Expressão cerebral (GSE160208)

**E1 — "Réplica exata" (corrigido).** A ordem dos filtros e o universo são
parâmetros de `analise_gse160208.py`. O texto diz o que foi reproduzido: a
contagem de 184, com o filtro de magnitude antes do FDR.

**E2 — 437 × 447 e genes degenerados (corrigido).** 58 genes sem variância no
córtex frontal são marcados como não estimáveis. 441 no universo de 742; 437
no de 800. O "447" não foi reproduzido em nenhuma configuração.

**E3 — Docstring (corrigido).**

### Integração sangue × cérebro

**I1 a I5 (corrigidos).** K restrito ao painel testável; espécie do gene-alvo
humana; miRNAs sem alvo fora do teste e do FDR; leitura com `csv`; "controle
negativo" sem teste e leitura de mecanismo retirados; poder declarado.
As âncoras do dossiê (cálculo antigo e "correto" com N = 800) são reproduzidas
em `tests/test_integracao.py`. O resultado continua nulo.

### Reprodutibilidade e engenharia

**R1 — GWAS (corrigido).** `pipeline/scripts/obter_gwas.py` baixa do GWAS
Catalog, retoma se a conexão cair e confere o MD5
(`071790d80ccae0b41adabc7e0eefbf53`). Os scripts que leem o arquivo passam por
`gwas_fonte.exigir_gwas()`, que explica como obtê-lo.

**R2 — Dependências (corrigido).** `requirements.txt` com versões fixadas,
única fonte, usada pelo `make install` e pelo CI.

**R3 — Makefile (corrigido).** `clean` não apaga dados versionados (teste:
`test_make_clean_nao_apaga_nada_versionado`); instalação em `.venv`; alvos
`dados-gwas`, `relatorios-offline`, `relatorios-gwas`, `relatorios-rede`,
`relatorios`.

**R4 — Proveniência do eQTL (corrigido).** `pipeline/data/proveniencia_qtl.tsv`
registra release, commit da metadata, janela, data e SHA-256 de cada cache.
Em 2026-11-10 os nove caches eram idênticos, byte a byte, ao que o eQTL
Catalogue servia para a mesma janela. Se a r8 é final ou prévia não foi
apurado; o registro de checksums cobre os dois casos.

**R5 — Arquivos soltos (corrigido).** Dez índices `.tbi` removidos da raiz.

**R6 — Regeneração do zero (feito; ver `auditoria/PROGRESSO.md`).** Clone limpo
em outro diretório, `make install`, `make test`, `make dados-gwas`,
`make relatorios`. O teste
`test_relatorio_versionado_e_o_que_o_script_gera` roda os scripts offline e
exige que a saída seja igual ao relatório versionado.

### Testes

**T1 (corrigido).** O teste que exigia "0.9950" foi substituído por testes que
calculam. **T2 (corrigido).** Colocalização, hipergeométrico, regressão
ajustada e contagens de expressão têm testes que recalculam contra referência
independente. **T3 (corrigido).** Os testes acoplados a "97%", "v0.4" e "v3"
foram substituídos por `tests/test_textos.py`, que recalcula cada número antes
de procurá-lo no texto. Os guardas de privacidade e de integridade de citação
não foram tocados.

### Documentação e manuscrito

**D1 (corrigido).** As promessas de "tudo reproduzível com a biblioteca
padrão" foram trocadas pela descrição das dependências e do dado não
versionado. **D2 (corrigido).** Manuscrito v0.5 com nota de correção no início;
PDF regenerado. **D3 (corrigido).** λGC = 1,0587 nos dois relatórios.

## Pendências, com o motivo

- **Identidade gene a gene com a lista S1b e q-valores do Qlucore (seção 2.4
  do dossiê).** O suplemento do artigo não está no repositório; a comparação
  de 2026-08-24 não foi refeita.
- **R6, parte que depende de terceiros.** `finemap_ld.py` e
  `valida_offset_grch38_grch37.py` dependem do Ensembl REST e do NCBI
  E-utilities no momento da execução; não há cópia versionada dessas respostas.
- **Revisão científica do §4.7.** O texto está coerente com os números; não foi
  lido por alguém de genética estatística.
- **Versão publicada.** O registro no Zenodo e o site publicado continuam com a
  conclusão da v0.4 até uma nova versão ser depositada.

## Problemas novos encontrados durante a correção

- **N1 (corrigido).** Trechos de relatórios tinham sido escritos à mão depois
  de gerados (seção 3b do QC, "reconciliação" e "perda neuronal" do GSE160208):
  regenerar os apagava. Passaram para os scripts. Um deles estava errado
  ("30 dos 41 hits em PRNP").
- **N2 (corrigido).** `gera_figuras_v2.py` tinha os valores do gráfico de
  floresta dos miRNAs digitados à mão, e eles não correspondiam aos dados
  (d bruto de miR-16-5p: −0,17 na figura, −1,04 nos dados). Passou a calcular.
- **N3 (corrigido).** O manuscrito atribuía o menor p do eQTL-meta a
  chr1:180.949.780 (6,60×10⁻⁴⁷); o menor é 3,44×10⁻⁴⁷, em chr1:180.951.692.
- **N4 (corrigido).** A tabela de λ por estrato do manuscrito trazia n em
  dobro em relação ao relatório.
- **N5 (corrigido no que cita resultados).** O site ainda mostrava "60 miRNAs
  significantes" (número da v1), "réplica exata" e a leitura "o sangue reflete
  a periferia". Atualizado em português e inglês; o site não foi reimplantado.
- **N6 (não tratado).** A referência 2 do manuscrito tem título e terceiro
  autor diferentes dos que o próprio repositório registra para o mesmo DOI
  (`validacao_cruzada_gse160208_artigo_original.md`). Precisa de conferência na
  fonte.
- **N7 (não tratado).** `README.en.md` ainda lista no mapa do repositório
  arquivos que não são versionados (`ARQUIVO_COMPLETO.md`, `colaboracao/`,
  `memory/`).
- **N8 (não tratado).** `valida_offset_grch38_grch37.py` ainda casa variantes
  só pela posição; a conclusão dele (offset constante) não depende disso.
