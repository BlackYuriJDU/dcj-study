# dcj-study — alvos de reprodução. Tudo roda dentro do ambiente virtual .venv.
#
#   make install        cria .venv e instala requirements.txt (versões fixadas)
#   make test           suíte de testes (nenhum teste pode ser pulado)
#   make dados-gwas     baixa o GWAS (188 MB) do GWAS Catalog e confere o MD5
#   make relatorios     regenera TODOS os relatórios e figuras das análises
#                       (= relatorios-offline + relatorios-gwas + relatorios-rede)
#   make pdf            regenera o PDF do preprint
#   make clean          remove só derivados não versionados (nunca dados)

PY := .venv/bin/python
S  := pipeline/scripts

.PHONY: install test dados-gwas conferir-dados relatorios relatorios-offline \
        relatorios-gwas relatorios-rede figuras pdf simulacoes clean

install:
	python3 -m venv .venv
	.venv/bin/pip install -r requirements.txt

test:
	$(PY) -m pytest tests/ -q -rs

dados-gwas:
	$(PY) $(S)/obter_gwas.py

conferir-dados:
	$(PY) $(S)/obter_gwas.py --conferir
	$(PY) $(S)/registra_proveniencia_qtl.py

# --- Só dados versionados: não precisa de rede nem do GWAS -------------------
relatorios-offline:
	$(PY) $(S)/analise_gse160208.py
	$(PY) $(S)/analise_gse140069.py
	$(PY) $(S)/analise_gse140069_idade.py
	$(PY) $(S)/analise_integrada_cerebro_sangue.py
	$(PY) $(S)/gera_figuras.py

# --- Precisa do arquivo do GWAS (make dados-gwas); sem rede ------------------
relatorios-gwas:
	$(PY) $(S)/obter_gwas.py --conferir
	$(PY) $(S)/qc_gwas_gcst90001389.py
	$(PY) $(S)/clumping_descoberta.py
	$(PY) $(S)/coloc_stx6_eqtl.py
	$(PY) $(S)/coloc_meta_stx6.py
	$(PY) $(S)/coloc_sqtl_stx6.py
	$(PY) $(S)/gera_figuras_v2.py

# --- Precisa do GWAS E de rede (Ensembl REST, NCBI E-utilities) --------------
# finemap_ld.py reescreve relatorio_finemap_loci.md e relatorio_lambda_gc.md;
# crosscheck_coloc_R.py vem por último porque cita a massa do fine-mapping.
relatorios-rede:
	$(PY) $(S)/finemap_ld.py
	LD_POP=1000GENOMES:phase_3:EUR $(PY) $(S)/finemap_ld.py
	$(PY) $(S)/valida_offset_grch38_grch37.py
	$(PY) $(S)/crosscheck_coloc_R.py

relatorios: relatorios-offline relatorios-gwas relatorios-rede

pdf:
	$(PY) $(S)/gera_preprint_pdf.py

simulacoes:
	$(PY) $(S)/simulacao_prion.py
	$(PY) $(S)/varredura_blindagem.py
	$(PY) $(S)/simulacao_calibrada.py

# `clean` NÃO apaga nada versionado. Até 2026-11 este alvo fazia
# `rm -f pipeline/data/eqtl_*.tsv`, o que destruía caches versionados cuja
# origem é um servidor remoto. O GWAS baixado também fica: para removê-lo,
# apague pipeline/data/GCST90001389_buildGRCh37.tsv.gz à mão.
clean:
	find . -name __pycache__ -type d -not -path './.venv/*' -prune -exec rm -rf {} +
	rm -rf .pytest_cache
	rm -f *.tbi pipeline/data/*.part pipeline/data/offset_ensembl_cache.json
