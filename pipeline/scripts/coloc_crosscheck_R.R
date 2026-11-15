#!/usr/bin/env Rscript
# coloc_crosscheck_R.R — Roda o coloc.abf do pacote R `coloc` (a referência do
# campo) sobre a região STX6 harmonizada e grava o SNP.PP.H4 por variante.
#
# Entrada : pipeline/data/stx6_crosscheck_input.tsv (crosscheck_coloc_R.py)
# Saída   : pipeline/data/stx6_crosscheck_colocR_output.tsv
#
# Este script NÃO reimplementa nenhuma equação. Até a v0.4 ele tinha uma
# "parte (A)" que copiava a fórmula do Python (com o mesmo erro) e concluía
# que a implementação estava validada; essa parte foi removida na correção de
# 2026-11. A comparação com a implementação do projeto é feita em Python
# (crosscheck_coloc_R.py e tests/test_coloc_abf.py), contra a saída deste
# script.
#
# Uso (na raiz do repositório): Rscript pipeline/scripts/coloc_crosscheck_R.R
# Requer: R >= 4 e o pacote `coloc` (>= 5.2).
# Executado na correção de 2026-11 com R 4.5.3 e coloc 5.2.3 (conda-forge).
# Os testes usam dois oráculos: a saída de 2026-08-29 sobre a entrada da v0.4
# (congelada em tests/fixtures/) e a saída deste script sobre a entrada atual.

suppressMessages(library(coloc))

d <- read.delim("pipeline/data/stx6_crosscheck_input.tsv", check.names = FALSE)
cat("variantes:", nrow(d), "\n")

N_GWAS <- 17679
S_GWAS <- 4110 / N_GWAS
P1 <- 1e-4; P2 <- 1e-4; P12 <- 1e-5      # padrão do coloc

# Com beta e varbeta fornecidos, o fator de Bayes só depende deles e do prior
# de efeito (sd 0,2 para "cc"; 0,15*sdY para "quant"). N, s e MAF são exigidos
# pela interface, mas não entram na conta.
p_m <- pmin(d$maf, 1 - d$maf)
n_eff <- 1 / (2 * d$se_eqtl_meta^2 * p_m * (1 - p_m))
snp <- paste0("chr1:", d$pos_b37, ":", d$ref, ":", d$alt)

D1 <- list(snp = snp, position = d$pos_b37,
           beta = d$beta_gwas, varbeta = d$se_gwas^2,
           type = "cc", s = S_GWAS, N = N_GWAS)
D2 <- list(snp = snp, position = d$pos_b37,
           beta = d$beta_eqtl_meta, varbeta = d$se_eqtl_meta^2,
           type = "quant", sdY = 1, N = n_eff, MAF = p_m)

res <- coloc.abf(D1, D2, p1 = P1, p2 = P2, p12 = P12)
s <- res$summary
cat(sprintf("nsnps=%d · PP.H0=%.3g H1=%.3g H2=%.3g H3=%.6f H4=%.6f\n",
            s["nsnps"], s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"],
            s["PP.H3.abf"], s["PP.H4.abf"]))

# coloc.abf pode reordenar as variantes: alinhar pela identificação
r <- res$results
pp4 <- r$SNP.PP.H4[match(snp, r$snp)]
saida <- "pipeline/data/stx6_crosscheck_colocR_output.tsv"
writeLines(sprintf("# coloc.abf: coloc %s, %s; PP.H3=%.6f PP.H4=%.6f; nsnps=%d",
                   as.character(packageVersion("coloc")), R.version.string,
                   s["PP.H3.abf"], s["PP.H4.abf"], s["nsnps"]), saida)
suppressWarnings(write.table(
  data.frame(pos_b37 = d$pos_b37, PP4_colocR = pp4), saida,
  sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE))
cat("[ok] saida: pipeline/data/stx6_crosscheck_colocR_output.tsv\n")
