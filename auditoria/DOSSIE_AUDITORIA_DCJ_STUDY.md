# Dossiê de auditoria — `BlackYuriJDU/dcj-study`

Complemento do `PROMPT_CORRECAO_DCJ_STUDY.md`. Reúne tudo o que três auditorias independentes estabeleceram sobre o repositório, para que a sessão de correção comece pelas correções e não por uma nova auditoria.

- **Commit auditado:** `e0765256110b521220902873d607b4dbe08fd558` (119 arquivos rastreados; suíte: 37 aprovados, 0 pulados).
- **Data da consolidação:** 2026-11-10.
- **Fontes:** auditoria A (ChatGPT), auditoria B (Kimi K3.1 Max) e auditoria C (Claude, reexecução direta sobre um clone limpo, com recálculo próprio).

## Como ler o estado de cada item

| Marca | Significado |
|---|---|
| **Verificado ×3** | A, B e C chegaram ao mesmo resultado; C recalculou a partir dos arquivos do repositório. |
| **Verificado ×2** | Duas fontes concordam (indicadas). |
| **Só C** | Achado novo da reexecução C; A e B não o registraram. |
| **Não reverificado** | Afirmado por A ou B; C não conseguiu refazer (motivo indicado). Tratar como hipótese forte, não como fato. |
| **Divergência** | As fontes discordam no número; o valor de C é o que consta, com a diferença anotada. |

Limite da auditoria C: o arquivo do GWAS (`GCST90001389_buildGRCh37.tsv.gz`) não está no Git e a rede do ambiente bloqueou o FTP do EBI. Tudo o que depende desse arquivo (QC, clumping, fine-mapping, colocalização por tecido, sQTL) ficou **não reverificado** por C. A colocalização da meta-eQTL foi reverificada porque sua entrada de 390 variantes está versionada.

---

## 1. Resumo executivo

1. **A conclusão "STX6 não colocaliza" vem de um erro de fórmula.** O termo de H4 foi multiplicado por `p1·p2·p12` em vez de `p12`. Corrigindo só isso, com os mesmos dados e o mesmo prior de efeito, PP.H4 passa de 4,6×10⁻⁷ para 0,979. **Verificado ×3.**
2. **A "validação em R" não validava.** A parte (A) do script R copia a mesma fórmula. A parte (B), que chama o `coloc.abf` oficial, já devolvia H4 = 0,98, e a divergência foi explicada no relatório e no manuscrito como "dependência de prior". Não é: o prior de efeito idêntico dá 0,979. **Verificado ×3.**
3. **Há um segundo erro na mesma função** (`abf` usa `√r` onde Wakefield usa `√(1−r)`) e um provável terceiro (`W = 0,04²`, rotulado "coloc padrão", quando o padrão é variância 0,04). Mudam pouco o H4 desta região, mas são erros de fórmula. **Só C.**
4. **O H4 corrigido é estável** a efeito aleatório, alelos estritos, leave-one-out e exclusão do GTEx (0,975–0,980), e **sensível a `p12`** (0,83 com `p12 = 10⁻⁶`). **Verificado ×2 (B, C).**
5. **Os demais problemas não invertem conclusões**, mas são reais: universo do teste hipergeométrico, OLS instável em respostas constantes, harmonização de alelos por posição, tamanhos amostrais do GTEx desatualizados, redação causal ("idade explica 97%", "réplica exata"), dependências não declaradas, GWAS sem download automatizado, e testes que travam números em vez de recalculá-los.
6. **Colocalização não é mediação.** O texto corrigido não pode trocar um excesso por outro: H4 alto apoia variante causal compartilhada sob os pressupostos do modelo, e nada além disso.

---

## 2. Âncoras numéricas

Valores que a sessão de correção deve reproduzir antes de mexer em qualquer coisa. Todos calculados por C sobre o commit auditado, salvo indicação.

### 2.1 Colocalização STX6 (entrada: `pipeline/data/stx6_crosscheck_input.tsv`, 390 variantes)

| Configuração | H2 | H3 | H4 |
|---|---|---|---|
| **[A] Como o repositório calcula** (`√r`, W = 0,0016, H4 × `p1·p2·p12`) | 0,005002 | 0,994997 | 4,619×10⁻⁷ |
| [A] idem, priors "conservador" (10⁻⁶, 10⁻⁶, 10⁻⁷) | 0,334571 | 0,665429 | 3,089×10⁻⁹ |
| **[B] Só a combinação corrigida** (resto igual ao repositório) | 0,000106 | 0,021085 | **0,978809** |
| [C] Combinação + ABF de Wakefield correto, W = 0,0016 | 0,000164 | 0,021160 | 0,978676 |
| **[D] Tudo conforme o `coloc` oficial** (W = 0,04 para caso-controle; 0,0225 para quantitativo, sdY = 1) | 7,04×10⁻⁶ | 0,020388 | **0,979605** |

Sensibilidade a `p12` na configuração [D], com `p1 = p2 = 10⁻⁴`:

| `p12` | H3 | H4 |
|---|---|---|
| 10⁻⁵ (padrão do `coloc`) | 0,020388 | 0,979605 |
| 5×10⁻⁶ | 0,039962 | 0,960025 |
| 10⁻⁶ | 0,172263 | 0,827677 |
| "conservador" do repositório (10⁻⁶, 10⁻⁶, 10⁻⁷) | 0,000208 | 0,999785 |

Sensibilidade a W (igual nos dois traços, combinação oficial, priors padrão): H4 = 0,978676 (W = 0,0016), 0,979530 (0,01), 0,979584 (0,0225), 0,979601 (0,04), 0,979611 (0,1).

Outros números da mesma entrada:

- Concordância de direção GWAS × meta-eQTL: 348/390 (89%).
- Maior `SNP.PP.H4`: chr1:180.951.692 (0,1048); depois 180.956.015 (0,0672), 180.956.985 (0,0656) e o lead do GWAS 180.961.245 (0,0634). O conjunto de 95% tem 17 variantes. A variante de maior posterior **não** é o lead do GWAS.
- **Oráculo já versionado:** `pipeline/data/stx6_crosscheck_colocR_output.tsv` guarda o `SNP.PP.H4` que o `coloc.abf` 5.2.3 real produziu. A reimplementação de C em Python coincide com ele variante a variante (diferença máxima 1,1×10⁻¹⁶), e o resumo coincide com o que o relatório registra para o pacote (H4 = 0,9796; H3 = 0,0204; H2 = 7×10⁻⁶). Serve de teste de regressão sem precisar de R.

Divergências entre fontes, todas na quinta casa decimal ou abaixo: B reporta H4 = 0,97965 e H3 = 0,02035 onde C obtém 0,979605 e 0,020388; B reporta 0,8280 com `p12 = 10⁻⁶` onde C obtém 0,8277. Nenhuma altera leitura.

### 2.2 Meta-eQTL (caches `pipeline/data/eqtl_QTD*_stx6.tsv`)

| Item | Valor | Estado |
|---|---|---|
| Posições na meta (≥ 3 datasets) | 651 | Verificado ×3 |
| Casadas com o GWAS | 390 | Verificado ×3 |
| Posições com REF/ALT discordantes combinadas em silêncio | 18 de 651 | Verificado ×3 |
| Dessas, dentro das 390 | 1 (chr1:180.975.614, G>T em dois datasets e G>C em três) | Verificado ×3 |
| Heterogeneidade nas 390: Q com p < 0,05 | 37 | Verificado ×2 (B, C) |
| Heterogeneidade nas 390: I² > 50% | 52 | Verificado ×2 (B, C) |
| Lead 180.961.245: k = 5, I² = 48,4%, z de efeito fixo = −14,37 | | Verificado ×2 (B, C) |
| Meta recalculada × entrada versionada | diferença máxima 0 | Só C |

Posições com mais de um par de alelos **dentro** de um mesmo dataset (C): CommonMind 29, ROSMAP 18, BrainSeq 38, GTEx DLPFC 0, GTEx cerebelo 0. A auditoria A reportou 74, 39 e 95; C não reproduziu esses números com a definição "posições", e a diferença parece ser de contagem (linhas em vez de posições). Há ainda linhas duplicadas do mesmo alelo, uma por rsID (86, 69, 91, 24, 24), sempre com β e erro-padrão idênticos, e indels em todos os arquivos (181, 137, 221, 44, 45 variantes).

Sensibilidades do H4 na configuração [D] (C; valores de B entre parênteses quando diferem):

| Variação | H4 | n |
|---|---|---|
| Efeito fixo, lógica do repositório | 0,97960 | 390 |
| DerSimonian–Laird | 0,97547 (B: 0,97431) | 390 |
| Efeito fixo, alelos estritos | 0,97960 (B: 0,97968) | 390 |
| Efeito aleatório, alelos estritos | 0,97547 (B: 0,97434) | 390 |
| Leave-one-out, cinco fontes | 0,97940 a 0,97959 (B: 0,97944 a 0,97968) | 365 a 380 |
| Sem os dois GTEx | 0,97947 (B: 0,97954) | 363 |
| Só CommonMind / só ROSMAP / só BrainSeq | 0,9782 / 0,9790 / 0,9792 | 389 / 365 / 374 |
| Só GTEx DLPFC / só GTEx cerebelo | 0,7015 / 0,5672 | 371 / 380 |

### 2.3 Tamanhos amostrais (metadata do eQTL Catalogue, commit `7a5b65b` de 2026-10-30, `data_tables/dataset_metadata_r8.tsv`)

| Dataset | Projeto declara | r8 oficial | r6/r7 (GTEx v8) |
|---|---|---|---|
| CommonMind QTD000075 | 586 | 586 | 586 |
| ROSMAP QTD000434 | 560 | 560 | 560 |
| BrainSeq QTD000051 | 479 | 479 | 479 |
| GTEx DLPFC QTD000176 | **285** | **268** | 175 |
| GTEx cerebelo QTD000166 | **272** | **263** | 209 |
| Total | ~2.182 | **2.156** | |

Verificado ×2 (B, C). Os n não entram no ABF calculado a partir de β e erro-padrão, portanto não mudam o H4.

### 2.4 Expressão cerebral GSE160208 (córtex frontal: 14 sCJD, 10 controles; 800 genes)

| Item | Valor | Estado |
|---|---|---|
| Welch + BH sobre 800, FDR < 0,05 | 437 (314 ↑, 123 ↓) | Verificado ×3 |
| Student + BH sobre 800 | 410 | Só C |
| p < 0,05 e \|Δ\| > 1, sem exigir q | 184 | Verificado ×2 (B, C) |
| p < 0,05, q ≤ 0,06 no universo de 800 e \|Δ\| > 1 | 183 (Welch), 182 (Student) | Verificado ×2 (B, C) |
| Filtro \|Δ\| > 1 primeiro (185 genes), depois BH | 184 (Welch e Student) | Verificado ×2 (B, C) |
| Genes com p indefinido (variância nula) dentro dos 800 | 58 | Só C |
| Identidade gene a gene com a lista S1b oficial | conjunto idêntico | Não reverificado (só B; C não tinha o suplemento) |
| q-values oficiais do Qlucore | não reproduzidos | Não reverificado (só B) |

### 2.5 miRNA de sangue GSE140069 (57 sCJD, 48 controles; 939 miRNAs)

| Item | Valor | Estado |
|---|---|---|
| Amostras com idade e RIN completos | 105 | Verificado ×2 (B, C) |
| Idade média, casos × controles | 66,4 × 53,6 anos | Verificado ×2 |
| Sexo | 60 F, 45 M, nenhum ausente | Só C |
| Detectáveis (≥ 25% acima do piso) | 269 | Verificado ×2 |
| A: Welch, BH sobre 939 | 84 | Verificado ×3 |
| A′: A ∩ detectáveis | 69 | Verificado ×3 |
| B: OLS ajustado, BH sobre 939 | 1 (`hsa-miR-500a-3p`, β = −7,9358, p = 8,03×10⁻⁶, q = 0,00754) | Verificado ×2 (B, C com `numpy.linalg.lstsq`) |
| B no universo filtrado (269) | 5 (`miR-1301-3p`, `miR-28-3p`, `miR-500a-3p`, `miR-502-3p`, `miR-93-5p`) | Verificado ×2 |
| miRNAs constantes em todas as amostras | 74 | Verificado ×2 |
| B com BH só sobre os 865 estimáveis | 1 | Só C |
| Posto da matriz de desenho / número de condição | 5 / 471 | Só C |

### 2.6 Integração sangue × cérebro (painel de 800; 314 genes ↑)

| miRNA | Alvos (repositório) | Alvos no painel | ∩ ↑ | Esperado (repo → correto) | p (repo → correto) |
|---|---|---|---|---|---|
| `hsa-miR-16-5p` | 48 (47 hsa, 1 rno) | 13 | 6 | 18,8 → 5,1 | 1,00 → 0,40 |
| `hsa-miR-93-5p` | 16 (15 hsa, 1 mmu) | 6 | 3 | 6,3 → 2,4 | 0,98 → 0,44 |
| `hsa-let-7i-5p` | 11 | 3 | 1 | 4,3 → 1,2 | 0,996 → 0,78 |
| `hsa-miR-106b-3p` | 0 | 0 | 0 | não estimável | não estimável |
| `hsa-miR-500a-3p` | 0 | 0 | 0 | não estimável | não estimável |

Verificado ×2 (A apontou o erro; C quantificou). Os valores "correto" ainda incluem os dois alvos não humanos; refazer após filtrar a espécie do alvo.

### 2.7 GWAS (não reverificado por C)

Segundo B e os relatórios do repositório: 6.314.492 variantes, 0 linhas malformadas, MD5 conferido com o GWAS Catalog, 41 hits com p < 5×10⁻⁸, λGC ≈ 1,0588, melhor variante em STX6 chr1:180.961.245 com p = 7,51×10⁻⁹. O relatório de QC diz 1,0588 e o de validação em R diz 1,0587. A massa de fine-mapping de 90,5% não foi reexecutada por nenhuma das três auditorias.

---

## 3. Achados por domínio

Cada achado traz Erro, Porquê e Como corrigir, com prioridade e estado de verificação. As ondas do prompt referenciam estes identificadores.

### 3.1 Genética e colocalização

**G1 — Fórmula de H4 (Crítica · Verificado ×3)**
- **Erro:** `p1 * p2 * p12 * diag` no termo de H4, em quatro lugares: `coloc_stx6_eqtl.py:131`, `coloc_meta_stx6.py:116`, `coloc_sqtl_stx6.py:88` e `coloc_crosscheck_R.R:34`. Os termos H0, H1 e H3 também diferem da forma oficial por fatores `(1−p)`, sem efeito numérico relevante.
- **Porquê:** no `coloc`, `p12` já é a probabilidade a priori de uma variante ser causal para os dois traços. Multiplicar por `p1·p2` penaliza H4 em 10⁸ com os priors padrão. Referência, `coloc::combine.abf`: `lH1 = log(p1) + logsum(l1)`; `lH2 = log(p2) + logsum(l2)`; `lH3 = log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(l1+l2))`; `lH4 = log(p12) + logsum(l1+l2)`; `lH0 = 0`.
- **Como corrigir:** uma única implementação em espaço logarítmico, em `pipeline/lib/`, importada pelos três scripts Python. Nenhuma cópia local. Resultado esperado: linhas [B] e [D] da seção 2.1.

**G2 — Validação cruzada sem independência (Crítica · Verificado ×3)**
- **Erro:** a parte (A) de `coloc_crosscheck_R.R` reimplementa a fórmula errada e o relatório conclui "a máquina está correta — acordo a 6 decimais". A parte (B), com o pacote oficial, deu H4 = 0,9796 e foi tratada como efeito de prior.
- **Porquê:** duas implementações da mesma equação errada concordam entre si. Concordância só valida quando uma das pontas é a referência.
- **Como corrigir:** o teste de regressão passa a comparar contra o oráculo versionado (`stx6_crosscheck_colocR_output.tsv`) e contra casos sintéticos de resposta conhecida. A parte (A) do script R é removida ou reescrita para chamar `coloc:::combine.abf`. `relatorio_validacao_coloc_R.md` é reescrito: o que ele registra hoje como validação é a evidência do erro.

**G3 — Fator de Bayes de Wakefield (Alta · Só C)**
- **Erro:** `coloc_stx6_eqtl.py:115-118` calcula `sqrt(r) * exp(...)` com `r = W/(V+W)`. O correto é `sqrt(1−r)`, ou seja, `sqrt(V/(V+W))`. A mesma forma errada está em `abf_w` do script R.
- **Porquê:** o próprio repositório tem a forma correta em `finemap_ld.py:96-101`. As duas funções chamadas `abf` dão valores diferentes para a mesma variante. O erro pesa mais onde os erros-padrão variam muito entre variantes; aqui muda o H4 de 0,978809 para 0,978676.
- **Como corrigir:** uma função única de `log ABF`, `0.5 * (log(1 − r) + r * z²)`, usada por colocalização e fine-mapping. Conferir que o fine-mapping não muda.

**G4 — Prior de efeito W (Média · Só C para a causa; B notou a diferença)**
- **Erro:** `W_ABF = 0.04 ** 2` com o comentário "coloc padrão".
- **Porquê:** o padrão do `coloc` é desvio-padrão 0,2 para caso-controle (variância 0,04) e 0,15·sdY para traço quantitativo (0,0225 com sdY = 1). O valor 0,04 já é a variância; o código o elevou ao quadrado de novo. Hipótese mais provável, a confirmar pelo histórico do arquivo.
- **Como corrigir:** priors de efeito por tipo de traço, com os valores do `coloc`, e a suposição sdY = 1 para o eQTL declarada e justificada (ou sdY estimado). Reportar a sensibilidade a W da seção 2.1.

**G5 — Hipóteses mal descritas (Alta · Verificado ×3)**
- **Erro:** texto afirma que "H2 não é hipótese nula: é duas variantes causais distintas", soma "H2+H3+H4" como "ambos os sinais reais" e rotula `H4/(H3+H4)` um valor calculado como `H4/H3`. Locais: `coloc_meta_stx6.py:124-126` e `136-145`, `coloc_stx6_eqtl.py:183` e `193`, `relatorio_validacao_coloc_R.md`, manuscrito §4.7 (linhas 171-182) e limitação da linha 267.
- **Porquê:** H2 é associação só com o traço 2. "Ambos associados" é H3 + H4.
- **Como corrigir:** descrições corretas das cinco hipóteses em código, relatórios e manuscrito; o rótulo e o cálculo passam a ser `H4/(H3+H4)`.

**G6 — Harmonização de alelos por posição (Alta · Verificado ×3)**
- **Erro:** `d[pos38] = (ref, alt, beta, se)` em `coloc_meta_stx6.py:77` e `crosscheck_coloc_R.py:42`. A última linha de cada posição vence, e `meta_analise` combina β de datasets sem conferir se REF/ALT coincidem.
- **Porquê:** em sítio multialélico, somam-se efeitos de alelos diferentes. Números na seção 2.2. O impacto medido no H4 desta região é nulo, mas a meta-análise está formalmente errada em 18 posições.
- **Como corrigir:** chave `(cromossomo, posição, REF, ALT)`, com indels normalizados; deduplicar linhas repetidas por rsID; combinar apenas alelos idênticos; o casamento com o GWAS confere os dois alelos e inverte o sinal quando necessário; registrar quantas variantes foram descartadas e por quê.

**G7 — Sensibilidade a priors (Média · Verificado ×2)**
- **Erro:** o cenário "conservador" (10⁻⁶, 10⁻⁶, 10⁻⁷) é apresentado como mais rigoroso.
- **Porquê:** a razão `p12/(p1·p2)` sobe de 10³ para 10⁵, o que favorece H4. A sensibilidade que importa, `p12` com `p1` e `p2` fixos, não é reportada.
- **Como corrigir:** tabela de sensibilidade a `p12` (seção 2.1) nos relatórios e no manuscrito; renomear ou retirar o cenário "conservador", explicando a razão de priors.

**G8 — Tamanhos amostrais do GTEx (Média · Verificado ×2)**
- **Erro:** 285 e 272 em `coloc_stx6_eqtl.py:52`, `coloc_meta_stx6.py:29-30`, `coloc_sqtl_stx6.py`, relatórios de colocalização e manuscrito (linhas 33 e 237, "~2.182").
- **Porquê:** a metadata r8 lista 268 e 263; o total é 2.156.
- **Como corrigir:** ler os n de um arquivo de proveniência versionado em vez de constantes digitadas; atualizar textos.

**G9 — Interpretação (Alta · Verificado ×2)**
- **Erro:** risco de trocar "não colocaliza" por "mediado por expressão de STX6".
- **Porquê:** colocalização testa variante causal compartilhada sob a suposição de um único sinal causal por traço. Não testa mediação. Os dois tecidos do GTEx compartilham doadores e não são coortes independentes. Sozinhos, dão H4 de 0,70 e 0,57. Dentro de um bloco em LD quase perfeito, parte da separação entre H3 e H4 vem do prior: se as 390 variantes fossem equivalentes, a razão a priori H4/H3 seria `p12/(p1·p2·(n−1))` ≈ 2,6, e a razão a posteriori observada é ≈ 48. Esse diagnóstico é de C e deve ser recalculado antes de entrar no texto.
- **Como corrigir:** redação que afirme o suporte a variante compartilhada com seus pressupostos, a sensibilidade a `p12`, a robustez numérica e os limites (sinal único, LD, sobreposição de doadores, tecido). Lembrar que a associação de STX6 com expressão em sCJD já é trabalho publicado (referência 4 do manuscrito); o resultado corrigido é uma reprodução, não um achado novo.

**G10 — Colocalização com sQTL (Média · Verificado ×2 para a fórmula; resto só C)**
- **Erro:** `coloc_sqtl_stx6.py` usa a mesma fórmula; o relatório aplica "Bonferroni" a um limiar de PP.H4 (0,5/6).
- **Porquê:** os PP.H4 ≈ 0 dos seis testes em genes vizinhos estão contaminados pelo erro G1. Bonferroni é correção de p-valores, não de probabilidades a posteriori. Os arquivos `cc` só contêm pares significativos, e ausência neles não é ausência de sQTL.
- **Como corrigir:** recalcular com a implementação única; retirar o limiar "Bonferroni"; manter a ressalva de poder que o manuscrito já tem.

### 3.2 miRNA de sangue (GSE140069)

**M1 — OLS instável em respostas constantes (Alta · Verificado ×2, mecanismo detalhado por C)**
- **Erro:** `resolver` em `analise_gse140069.py:91-120` usa equações normais com Gauss–Jordan e não detecta resposta constante.
- **Porquê:** na execução exata do repositório, os 74 miRNAs constantes recebem p = 0,0203, quando o teste não é estimável. O valor nasce de ruído de arredondamento: os mesmos números passados como `numpy.float64` (o que provavelmente troca o algoritmo de soma do Python) dão p = 5,5×10⁻⁷, o bastante para criar dezenas de falsos positivos após BH. B observou 87 resultados sob outra variante equivalente. O resultado publicado (1 miRNA) sobrevive porque esta execução específica caiu do lado inofensivo.
- **Como corrigir:** solver por QR ou SVD; resposta constante ou posto deficiente marcam o teste como não estimável (sem p); testes de borda com dados degenerados; conferência contra `numpy`/`statsmodels`.

**M2 — Universo do BH (Média · Só C)**
- **Erro:** os 74 não estimáveis contam entre os 939 testes.
- **Porquê:** universo de testes inclui hipóteses que não podem ser rejeitadas. Com 865, o resultado continua 1.
- **Como corrigir:** definir e registrar o universo antes do BH; reportar as duas contagens.

**M3 — Codificações frágeis (Média · Verificado ×2)**
- **Erro:** todo grupo diferente de `Control` vira caso (linha 215); sexo ausente vira 0 (linha 220); idade e RIN passam por `.isdigit()`; `fdr_bh` indexa por nome.
- **Porquê:** nenhum desses atalhos falha com os dados atuais (grupos `sCJD`/`Control`, sexo completo, 939 nomes únicos), mas todos falhariam em silêncio com outros dados.
- **Como corrigir:** validar categorias e levantar erro no inesperado; tratar ausência de forma explícita.

**M4 — "A idade explica ~97%" (Alta, redação · Verificado ×2)**
- **Erro:** manuscrito (linhas 39 e 257) e `relatorio_gse140069_idade.md` atribuem à idade ~97% do colapso.
- **Porquê:** a conta é `(88 − 4)/(88 − 1) = 96,6%`, uma comparação de contagens de significativos entre dois modelos. Não mede proporção de efeito. VIF baixo indica que o modelo é estimável, não que a atribuição é causal. Note-se ainda que 88 vem do OLS só com grupo, e o 84 do título vem do Welch.
- **Como corrigir:** redação descritiva com as contagens, os modelos e os universos. O teste `test_confundimento_idade_atribuido` exige a string "97%" e precisa mudar junto.

**M5 — Ponto a inspecionar (Só C, não é achado)**
O único sobrevivente, `hsa-miR-500a-3p`, tem β = −7,94 em log2 após ajuste. Um efeito dessa magnitude sugere miRNA do tipo detectado/não detectado. Olhar a distribuição por grupo antes de descrevê-lo no texto.

### 3.3 Expressão cerebral (GSE160208)

**E1 — "Réplica exata" (Média · Verificado ×2)**
- **Erro:** `relatorio_gse160208.md:81` ("critério EXATO"), `validacao_cruzada_gse160208_artigo_original.md:7` ("réplica ponto a ponto") e `README.en.md:16` ("exact replication").
- **Porquê:** a contagem 184 depende da ordem dos filtros e do universo do FDR (seção 2.4). Os q-values oficiais não foram reproduzidos.
- **Como corrigir:** afirmar que a lista de 184 genes foi reproduzida como conjunto aplicando o filtro de magnitude antes do BH, e que p e log2FC são próximos dos oficiais. O código deve expor a ordem dos filtros como parâmetro documentado.

**E2 — Inconsistência 437 × 447 e genes degenerados (Baixa · Verificado ×2)**
- **Erro:** `validacao_cruzada_gse160208_artigo_original.md:83` fala em "447 significativos"; o relatório principal diz 437. 58 genes sem variância ficam no universo de 800.
- **Como corrigir:** uma contagem, com o universo declarado; marcar os degenerados.

**E3 — Docstring desatualizada (Baixa · Verificado ×2)**
`analise_gse160208.py:15` ainda diz que o teste t "virá em versão futura".

### 3.4 Integração sangue × cérebro

**I1 — Universo do teste hipergeométrico (Alta · Verificado ×2)**
- **Erro:** `analise_integrada_cerebro_sangue.py:120-121` usa `len(tgt)` (todos os alvos do miRNA no genoma) com N = 800.
- **Porquê:** K tem de ser o número de alvos que pertencem ao painel de 800. Valores na seção 2.6. A conclusão nula se mantém, mas o "esperado" e os p estão errados.
- **Como corrigir:** restringir alvos ao painel medido antes do teste.

**I2 — Espécie do alvo (Média · Verificado ×2)**
- **Erro:** o filtro confere `Species (miRNA)` e ignora `Species (Target Gene)`.
- **Como corrigir:** exigir `hsa` nas duas colunas.

**I3 — Testes não estimáveis no BH (Média · Só C)**
- **Erro:** dois dos cinco miRNAs não têm nenhum alvo "Functional MTI" e entram no BH com p = 1.
- **Como corrigir:** reportá-los como não testáveis e ajustar o número de testes.

**I4 — Leitura do CSV (Baixa · Só C)**
`split(",")` em vez do módulo `csv`. Hoje dá o mesmo resultado (48, 16, 11); o arquivo tem linhas duplicadas.

**I5 — Redação (Média · Verificado ×2)**
O "controle negativo" só conta interseções, sem teste. A "leitura honesta" fala em mecanismo a partir de um resultado nulo com poder baixo (painel de 800 genes, 3 a 13 alvos por miRNA).

### 3.5 Reprodutibilidade e engenharia

**R1 — GWAS sem obtenção automatizada (Alta · Verificado ×3)**
- **Erro:** o arquivo está no `.gitignore` e sete scripts o abrem direto (`qc_gwas_gcst90001389.py`, `clumping_descoberta.py`, `coloc_stx6_eqtl.py`, `finemap_ld.py`, `finemap_stx6.py`, `gera_figuras_v2.py`, `valida_offset_grch38_grch37.py`).
- **Como corrigir:** script de obtenção com URL do GWAS Catalog e conferência de MD5; mensagem clara quando o arquivo falta.

**R2 — Dependências (Média · Verificado ×3)**
- **Erro:** `requirements.txt` lista matplotlib, pysam, pytest e scipy com limites mínimos. O código importa também `openpyxl` e `reportlab`; os testes usam `numpy` e `statsmodels`. O CI instala outra lista, sem `pysam` nem `reportlab`.
- **Como corrigir:** versões fixadas; uma fonte só, usada pelo CI.

**R3 — Makefile (Média · Só C)**
- **Erro:** `clean` executa `rm -f pipeline/data/eqtl_*.tsv`, que apaga caches **versionados** cuja origem é um servidor remoto. `install` usa `--break-system-packages`.
- **Como corrigir:** `clean` não toca em dados versionados; instalação em ambiente virtual.

**R4 — Proveniência do eQTL (Média · Verificado ×2)**
- **Erro:** os relatórios citam apenas "eQTL Catalogue r8".
- **Como corrigir:** registrar release, data de obtenção, commit da metadata e checksum de cada cache. B afirma que a r8 ainda é prévia mutável com versão final prevista para dezembro de 2026; C não confirmou isso (o repositório de metadata tem `dataset_metadata_r8.tsv` e `dataset_metadata_r8_beta.tsv`). O registro de proveniência resolve a dúvida nos dois casos.

**R5 — Arquivos soltos (Baixa · Só C)**
Dez índices `.tbi` na raiz do repositório.

**R6 — Regeneração do zero (Alta · Verificado ×2)**
Nenhuma auditoria regenerou o projeto inteiro a partir de um clone limpo. Componentes isolados, sim.

### 3.6 Testes

**T1 — Testes que travam o resultado errado (Crítica · Verificado ×3)**
- **Erro:** `test_coloc_meta_tem_h2_e_validacao` exige "0.9950" no relatório e "0.994997" na validação.
- **Porquê:** um teste que confere se o número antigo continua escrito impede a correção e não detecta fórmula errada.
- **Como corrigir:** substituir por testes que calculam: casos sintéticos de resposta conhecida, o oráculo versionado, e as âncoras da seção 2.1.

**T2 — Consistência por string (Alta · Verificado ×3)**
A maior parte de `tests/test_consistencia_docs.py` confere presença de texto. Vale como guarda de redação e de privacidade; não vale como validação científica. Todo resultado central precisa de um teste capaz de falhar quando a implementação está errada.

**T3 — Testes acoplados a textos que vão mudar**
`test_confundimento_idade_atribuido` (exige "97%"), `test_preprint_v04_e_assinado` (exige "v0.4"), `test_relatorio_gse140069_e_v3`.

**Não mexer:** os testes de privacidade (`test_nome_do_caso_ausente_do_repo`, `test_nada_rastreado_esta_no_gitignore`, `test_arquivos_privados_nao_rastreados`, `test_sem_figura_legada_com_nome_de_pessoa`, `test_doi_contaminado_nao_e_anunciado`) e os de integridade de citação estão corretos e protegem contra incidentes que já aconteceram.

### 3.7 Documentação e manuscrito

**D1 — Promessas de reprodutibilidade (Média · Verificado ×3)**
`README.en.md:5` ("Everything here is reproducible from standard-library Python"), `README.en.md:26` e manuscrito linha 259 ("standard-library statistical core"). Há dependências externas e um dado não versionado.

**D2 — Trechos do manuscrito afetados (`preprint/manuscrito_preprint.md`)**
Linha 33 (resumo, ~2.182 amostras); 39 (97%); 105 (§3.5, W = 0,04² e "sensibilidade 10× mais conservadora"); 160 (184 DEGs, "identical count"); 162-190 (§4.7 inteiro); 236-241 (discussão); 257 (idade); 259 e 282 (stdlib e validação "against R anchors"); 267 (limitação "cannot split H₃ from H₄"). O PDF versionado (`preprint/manuscrito_preprint.pdf`) é gerado por `gera_preprint_pdf.py` e precisa ser regenerado.

**D3 — λGC** 1,0587 em um relatório e 1,0588 em outro. Padronizar o arredondamento.

---

## 4. O que não é problema

Para a sessão de correção não gastar tempo aqui:

- O **offset constante GRCh38→GRCh37** foi verificado pelo próprio projeto para a janela de STX6. Fica como está; liftOver por chain file é melhoria futura.
- O **núcleo estatístico** (`welch`, `tcdf_p`, `fdr_bh`, `cohen_d`) confere com SciPy nos casos não degenerados.
- A **associação regional** GWAS × expressão de STX6 é forte (lead com p = 7,60×10⁻⁴⁷ na meta-eQTL, concordância de direção 89%). O que muda é a conclusão sobre variante compartilhada.
- As **simulações**, o **caso de referência** e o **site** não foram auditados e estão fora do escopo da correção, salvo números que citem resultados corrigidos.
- A **limpeza de privacidade** de outubro e seus guardas.

---

## 5. Engenharia de prompt: o que foi confirmado para o Fable 5.2

Conferido na documentação da Anthropic em 2026-11-10.

| Conhecimento do material original | Veredicto | Observação |
|---|---|---|
| Clareza direta e justificativa das instruções | Confirmado | A documentação recomenda dar o motivo junto do pedido; o modelo generaliza a partir da explicação. |
| Tags XML para separar instruções, contexto e dados | Confirmado | Nomes descritivos e consistentes; aninhar quando há hierarquia. |
| Documentos longos antes da consulta | Confirmado | Para entradas de 20 mil tokens ou mais. |
| "Até 30% de melhoria" | Confirmado com ressalva | "até 30% em testes", sobretudo com vários documentos. Não é garantia. |
| Few-shot com 3 a 5 exemplos | Confirmado | Exemplos relevantes, diversos e marcados. |
| Pré-preenchimento da resposta | **Desatualizado** | A partir dos modelos 4.6, prefill na última mensagem do assistente devolve erro 400. |
| Orçamento de tokens de raciocínio e pensamento visível | **Desatualizado para o Fable 5** | Só raciocínio adaptativo, saída de pensamento resumida, sem orçamentos. |
| Pedir que o modelo exponha o raciocínio | **Contraindicado** | Pode acionar a categoria de recusa `reasoning_extraction`. Pedir evidências, cálculos e resumo das ações é permitido. |
| Cache de prompt (1.024 tokens, TTL de 5 min, 4 pontos) | Não reverificado | Sem efeito neste uso: no Claude Code o cache é gerido pela ferramenta. |
| DSPy e otimização programática | Válido em geral | Não se aplica a uma execução única. |
| Sistemas multiagentes | Confirmado, com ajuste | O Fable 5 delega com mais facilidade; a orientação é dizer quando delegar e usar verificadores de contexto novo. |

Pontos específicos do Fable 5 usados no desenho do prompt:

- **Esforço:** cinco níveis (`low`, `medium`, `high`, `xhigh`, `max`). `high` é o padrão e o ponto de partida recomendado para o 5.2; `xhigh` onde houver ganho medido.
- **Instruções curtas bastam.** Prompts muito prescritivos, herdados de modelos anteriores, podem piorar o resultado. O prompt descreve objetivo, limites e critérios de aceite, e deixa a implementação para o modelo.
- **Limites explícitos.** O modelo pode tomar ações não pedidas; o prompt diz o que pode e o que não pode.
- **Relatos de progresso ancorados em evidência.** Instruir a conferir cada afirmação contra um resultado de ferramenta da sessão.
- **Pausas só quando necessárias:** ação destrutiva ou irreversível, mudança real de escopo, ou informação que só o usuário tem.
- **Verificação por subagente de contexto novo** supera a autocrítica.
- **Notas de memória** em arquivo ajudam em execuções longas; Git serve de registro de estado.
- **Classificadores de segurança:** o Fable 5 tem salvaguardas para conteúdo de biologia e ciências da vida (métodos de laboratório, mecanismos moleculares) que podem disparar em trabalho benigno. Este projeto é estatística e engenharia de software sobre dados públicos, mas, se uma etapa devolver recusa, a saída é rodar aquela etapa com outro modelo (a documentação cita o Claude Opus 4.8 como destino de fallback).

Correção de percurso das auditorias: B registrou no início que a documentação listava o Fable 5.1 e não o 5.2, e depois confirmou o 5.2. O modelo existe, com identificador `claude-fable-5-2`.

---

## 6. Repositórios semelhantes

Os seis citados por A existem (conferido em 2026-11-10). Servem de referência de método; nenhum código deve ser copiado.

| Repositório | O que é | Licença | O que aproveitar |
|---|---|---|---|
| `ssri-m/Systems-Biology-of-Creutzfeldt-Jakob-Disease` | GSE124571, limma, GO/KEGG/Reactome, STRING, estruturas | MIT | Organização por etapas com resultados intermediários |
| `HIGORMALz/creutzfeldt-jakob-mirna-ml` | Mesmo GSE140069; Random Forest, AUC ≈ 0,81 em validação cruzada de 5 partições (0,78 com 4 miRNAs) | Sem licença | Contraexemplo: classificador sem ajuste por idade na coorte em que a idade confunde |
| `arjunbabusb/CJD-Prion-Bioinformatics` | KEGG e redes STRING | Sem licença | Pouco |
| `htlin222/robust-lit-review` | Revisão sistemática com skills do Claude Code, validação de DOI, auditoria PRISMA | MIT | Rastreabilidade de fontes e auditoria de afirmações |
| `htlin222/meta-pipe` | Meta-análise em 9 etapas (protocolo a QA), revisão por agentes | Selo "Academic & Non-Commercial" | Portões entre etapas e revisor independente |
| `omertarikcicek/Clinical-Research-Data-Pipeline` | Pipeline clínico demonstrativo (dados simulados segundo A; C não confirmou) | Arquivo de licença presente, tipo não conferido | Validação de variáveis na camada de dados |

---

## 7. Fora do escopo desta execução

**Segunda doença.** B recomenda a ordem: corrigir o `dcj-study`, extrair um modelo independente de doença, e só então escolher outra. C concorda: abrir outro projeto agora propagaria os erros de estrutura. Candidatas levantadas por A, **não verificadas** por C: esclerose lateral amiotrófica (GSE234297, 96 casos esporádicos e 48 controles; GSE315586, 82,7 GB), Parkinson (GSE148434, GSE243639) e Alzheimer (campo muito explorado). Qualquer escolha exige conferir os conjuntos de dados e levar em conta o hardware disponível.

---

## Apêndice A — Verificação rápida das âncoras de colocalização

Só biblioteca padrão. Rodar na raiz do repositório:

```bash
python3 verifica_ancoras_coloc.py pipeline/data/stx6_crosscheck_input.tsv pipeline/data/stx6_crosscheck_colocR_output.tsv
```

```python
#!/usr/bin/env python3
"""Reproduz as ancoras da secao 2.1 do dossie (coloc::combine.abf em Python)."""
import csv, math, sys

def logsum(xs):
    m = max(xs)
    return m + math.log(sum(math.exp(x - m) for x in xs))

def logdiff(a, b):
    m = max(a, b)
    return m + math.log(math.exp(a - m) - math.exp(b - m))

def labf(beta, se, w, wakefield=True):
    r = w / (w + se * se)
    return 0.5 * (math.log(1 - r if wakefield else r) + r * (beta / se) ** 2)

def combina_oficial(l1, l2, p1, p2, p12):
    ls = [a + b for a, b in zip(l1, l2)]
    lh = [0.0, math.log(p1) + logsum(l1), math.log(p2) + logsum(l2),
          math.log(p1) + math.log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(ls)),
          math.log(p12) + logsum(ls)]
    d = logsum(lh)
    return [math.exp(x - d) for x in lh]

def combina_repo(l1, l2, p1, p2, p12):
    a1, a2 = [math.exp(x) for x in l1], [math.exp(x) for x in l2]
    sg, se, dg = sum(a1), sum(a2), sum(x * y for x, y in zip(a1, a2))
    h = [(1 - p1) * (1 - p2), p1 * (1 - p2) * sg, (1 - p1) * p2 * se,
         p1 * p2 * (1 - p12) * (sg * se - dg), p1 * p2 * p12 * dg]
    return [x / sum(h) for x in h]

rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
g = [(float(r["beta_gwas"]), float(r["se_gwas"])) for r in rows]
e = [(float(r["beta_eqtl_meta"]), float(r["se_eqtl_meta"])) for r in rows]
PAD = (1e-4, 1e-4, 1e-5)

def mostra(rot, h):
    print(f"{rot:44s} H2={h[2]:.6g} H3={h[3]:.6f} H4={h[4]:.6g}")

rg = [labf(b, s, 0.0016, False) for b, s in g]
re_ = [labf(b, s, 0.0016, False) for b, s in e]
mostra("[A] como o repositorio", combina_repo(rg, re_, *PAD))
mostra("[B] so a combinacao corrigida", combina_oficial(rg, re_, *PAD))
cg = [labf(b, s, 0.0016) for b, s in g]
ce = [labf(b, s, 0.0016) for b, s in e]
mostra("[C] + ABF de Wakefield correto", combina_oficial(cg, ce, *PAD))
og = [labf(b, s, 0.04) for b, s in g]
oe = [labf(b, s, 0.0225) for b, s in e]
mostra("[D] oficial (W 0.04 / 0.0225)", combina_oficial(og, oe, *PAD))
for p12 in (5e-6, 1e-6):
    mostra(f"[D] p12={p12:g}", combina_oficial(og, oe, 1e-4, 1e-4, p12))

if len(sys.argv) > 2:  # oraculo: SNP.PP.H4 do coloc.abf real, ja versionado
    ref = [float(r["PP4_colocR"]) for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")]
    ls = [a + b for a, b in zip(og, oe)]
    d = logsum(ls)
    print("max |SNP.PP.H4 - coloc.abf R| =", max(abs(math.exp(x - d) - y) for x, y in zip(ls, ref)))
```

Saída esperada:

```text
[A] como o repositorio                       H2=0.0050023 H3=0.994997 H4=4.61894e-07
[B] so a combinacao corrigida                H2=0.000106015 H3=0.021085 H4=0.978809
[C] + ABF de Wakefield correto               H2=0.000163796 H3=0.021160 H4=0.978676
[D] oficial (W 0.04 / 0.0225)                H2=7.04298e-06 H3=0.020388 H4=0.979605
[D] p12=5e-06                                H2=1.38044e-05 H3=0.039962 H4=0.960025
[D] p12=1e-06                                H2=5.95068e-05 H3=0.172263 H4=0.827677
max |SNP.PP.H4 - coloc.abf R| = 1.1102230246251565e-16
```

## Apêndice B — Fontes externas

- `coloc`, funções `combine.abf` e `approx.bf.estimates`: https://rdrr.io/cran/coloc/src/R/claudia.R
- eQTL Catalogue, metadata por release: https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources (`data_tables/dataset_metadata_r8.tsv`, commit `7a5b65b`)
- Guia de prompting do Claude Fable 5: https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/prompting-claude-fable-5
- Boas práticas gerais de prompting: https://platform.claude.com/docs/en/build-with-claude/prompt-engineering/claude-prompting-best-practices
- Níveis de esforço: https://platform.claude.com/docs/en/build-with-claude/effort
