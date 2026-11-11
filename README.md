# dcj-study — Projeto de Pesquisa Simulada

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22164909.svg)](https://doi.org/10.5281/zenodo.22164909)

> Simulação de pesquisa sobre **Doença de Creutzfeldt-Jakob esporádica (DCJ)**,
> centrada no caso fictício **"Caso Referência"**, usando apenas dados públicos reais
> e anonimizados que já existem.

> **Correção de novembro de 2026.** Versões anteriores deste repositório e do preprint (até a
> v0.4) afirmavam que o sinal de STX6 no GWAS não colocalizava com o eQTL cerebral de STX6. A
> afirmação vinha de um erro de fórmula no nosso próprio código. A análise corrigida, conferida
> variante a variante contra o pacote R `coloc`, dá suporte a uma variante causal compartilhada
> (probabilidade a posteriori 0,98 com os priors padrão; 0,83 com p12 = 10⁻⁶), o que reproduz
> trabalho já publicado. Colocalização não é mediação. Tudo o que mudou, com o valor antes e
> depois, está em [`CORRECOES.md`](CORRECOES.md).

## Missão
Obter, administrar, gerir e gerar valor a partir de dados priônicos públicos —
curados e analisados da melhor forma possível — para produzir materiais úteis
a laboratórios reais de pesquisa em príons.

## Mapa do projeto
| Caminho | Conteúdo | Status |
|---|---|---|
| [`MEMORIA.md`](MEMORIA.md) | Diário de sessão — contexto e pendências | vivo |
| `caso_referencia/dossie_clinico.md` | Dossiê clínico simulado (sCJD MM1) validado por fontes | ✅ v0.2 |
| `caso_referencia/linha_do_tempo.csv` · `exames_simulados.csv` · `fontes.md` | Dados estruturados + citações | ✅ |
| `research/estado_da_arte_dcj.md` | Estado da arte 2024–2025 completo com URLs | ✅ |
| `research/datasets_publicos.md` | Catálogo GEO (curadoria direta via E-utilities) | ✅ |
| `research/catalogo_datasets_prionicas_CJD.md` | Catálogo completo verificado: WGS/ENA, GWAS, proteômica PRIDE + o que NÃO existe | ✅ |
| `pipeline/data/GSE160208_series_matrix.txt.gz` | Dataset REAL baixado (47 amostras sCJD/CT, Copenhagen) | ✅ |
| `pipeline/scripts/tavily_search.sh` | Busca padrão via API Tavily | ✅ |
| `pipeline/scripts/analise_caso_referencia.py` → `reports/relatorio_caso_referencia.md` | Validação + critérios CDC/NPDPSC sobre o caso simulado | ✅ |
| `pipeline/scripts/analise_gse160208.py` → `reports/relatorio_gse160208.md` | Análise REAL: covariáveis + Welch/FDR + estratificação MM1 vs. VV2 | ✅ v3 |
| `pipeline/scripts/ponte_caso_referencia.py` → `reports/relatorio_ponte_caso_referencia.md` | Tabela-ponte: cada achado do caso de referência × evidência na coorte real | ✅ |
| `pipeline/scripts/analise_gse140069.py` → `reports/relatorio_gse140069.md` | Cross-modal REAL: miRNA SANGUE (57 sCJD vs. 48 CT) — 89 significativos sem ajuste entre 865 testes estimáveis (84 no universo antigo de 939) → 1 após ajuste idade/sexo/RIN; direção ↓ e núcleo nominal preservados | ✅ v4 |
| `pipeline/scripts/gera_figuras.py` → `reports/figuras/` | Volcanos ×2, heatmap top-25, timeline do caso de referência | ✅ 4 figuras |
| `colaboracao/centros_alvo.md` | NPDPSC, UCL, vigilância BR — contatos e estratégia | ✅ |
| `colaboracao/carta_projeto.md` | Carta de apresentação (rascunho — não enviar sem revisão) | ✍️ |
| `memory/` | Aprendizados (decisões, erros, padrões) | vivo |

> **Nota de privacidade (2026-08-29)**: `MEMORIA.md`, `memory/`, `colaboracao/` e
> `ARQUIVO_COMPLETO.md` continuam a existir localmente, mas deixaram de ser
> rastreados neste repo público — a face pública do projeto é o código,
> relatórios, preprint e este README.

## Como reproduzir
```bash
make install             # ambiente virtual .venv com as versões fixadas em requirements.txt
make test                # suíte de testes; nenhum teste pode ser pulado
make dados-gwas          # baixa o GWAS GCST90001389 (188 MB) e confere o MD5
make relatorios-offline  # expressão cerebral, miRNA, sangue × cérebro (só dados versionados)
make relatorios-gwas     # QC do GWAS, clumping, colocalização (precisa do arquivo do GWAS)
make relatorios-rede     # fine-mapping com LD e checagem de offset (Ensembl e NCBI)
```
O núcleo estatístico (`pipeline/lib/`) usa só a biblioteca padrão; ler a planilha, gerar figuras e
o PDF e rodar os testes exigem os pacotes de `requirements.txt`. Cada resultado central é
recalculado nos testes contra uma referência independente (`coloc.abf`, SciPy, statsmodels).

## Destaques científicos
- **STX6**: a colocalização entre o sinal do GWAS e o eQTL cerebral de STX6 (meta de 5 conjuntos, 2.156 amostras) dá PP.H4 = 0,98 com os priors padrão do `coloc` e 0,83 com p12 = 10⁻⁶, sob o modelo de uma variante causal por traço. Reproduz resultado publicado; não demonstra mediação.
- **Tratamento**: frente mais promissora = redução do alvo PrP (ASO ION717, NCT06153966; siRNA Broad, NCT07444580); quinacrina/doxiciclina/pentosano: capítulos encerrados.
- **Dados reais**: GSE160208 mostra forte ativação de genes de neuroinflamação na sCJD (MYD88, TLR2, C1QA, CSF1 com FDR<10⁻⁴); 441 de 742 genes testáveis com FDR<0,05.

## Regras invioláveis
1. **Sem pacientes reais ou novos.** Somente datasets públicos, anonimizados.
2. **Honestidade científica**: distinguir evidência estabelecida de experimento preliminar.
3. Tudo documentado em PT-BR, rastreável à fonte original.
