# Diagnóstico λ_GC — GCST90001389
*`finemap_ld.py` em 2026-11-11. λ = mediana(χ²)/0.454936.*

- λ global (amostra sistemática de 631,449): **1.0587**

| Estrato de MAF | n | λ do estrato |
|---|---|---|
| MAF<0.05 | 162,618 | 1.0617 |
| 0.05–0.25 | 266,332 | 1.0579 |
| 0.25–0.45 | 164,460 | 1.0547 |
| >0.45 | 38,039 | 1.0703 |

## Leitura honesta
- Gradiente de λ entre estratos: **0.0156**.
- Estratificação populacional clássica infla MAIS os alelos comuns;
  gradiente pequeno (<0.02) sugere inflação majoritariamente poligênica/
  residual, não estratificação grave. Gradiente grande (>0.05) pede PCA.
- Conclusão para o manuscrito: λ global 1.059 é limítrofe-saudável;
  declaramos correção por genomic control nas inferências primárias.