# -*- coding: utf-8 -*-
"""
test_reprodutibilidade.py — Infraestrutura de reprodução (auditoria 2026-11:
R1, R2, R3, R5 e D3).

- o GWAS não versionado tem obtenção automatizada e mensagem clara se faltar;
- a leitura do GWAS é em fluxo e não sobrescreve variantes;
- dependências: uma fonte só, com versões fixadas, usada pelo CI;
- `make clean` não apaga nada versionado.
"""
import ast
import gzip
import hashlib
import re
import subprocess
import sys
from pathlib import Path

import pytest
from scipy.stats import chi2

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
import gwas_fonte as GF  # noqa: E402

CAB = ("chromosome\tbase_pair_location\tother_allele\teffect_allele\t"
       "effect_allele_frequency\todds_ratio\tci_lower\tci_upper\tp_value\t"
       "beta\tstandard_error\n")
LINHAS = [
    "1\t100\tA\tG\t0.30\t1.1\t1.0\t1.2\t0.01\t0.10\t0.04\n",
    "1\t100\tA\tT\t0.05\t0.8\t0.6\t1.1\t0.20\t-0.20\t0.15\n",   # mesmo sítio
    "1\t250\tC\tT\t0.40\t1.0\t0.9\t1.1\t0.90\t0.00\t0.03\n",
    "11\t100\tA\tG\t0.30\t1.1\t1.0\t1.2\t0.01\t0.10\t0.04\n",    # chr 11 ≠ chr 1
    "2\t100\tG\tC\t0.10\t1.3\t1.1\t1.5\t0.001\t0.26\t0.08\n",
]


@pytest.fixture
def gwas_sintetico(tmp_path, monkeypatch):
    p = tmp_path / "gwas.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(CAB + "".join(LINHAS))
    monkeypatch.setattr(GF, "GWAS_GZ", p)
    return p


# ------------------------------------------------------------------ R1: GWAS
def test_gwas_ausente_da_mensagem_com_o_comando(tmp_path, monkeypatch):
    monkeypatch.setattr(GF, "GWAS_GZ", tmp_path / "nao_existe.tsv.gz")
    with pytest.raises(SystemExit) as e:
        GF.exigir_gwas()
    msg = str(e.value)
    assert "obter_gwas.py" in msg and "GCST90001389" in msg


def test_carregar_regiao_em_fluxo_sem_sobrescrever(gwas_sintetico):
    r = GF.carregar_regiao("1", 50, 200)
    assert [(x["pos"], x["outro"], x["efeito"]) for x in r] == \
        [(100, "A", "G"), (100, "A", "T")]          # os dois alelos do sítio
    assert r[0]["beta"] == 0.10 and r[0]["se"] == 0.04 and r[0]["eaf"] == 0.30
    assert GF.carregar_regiao("1", 0, 1000)[-1]["pos"] == 250
    assert len(GF.carregar_regiao(2, 0, 1000)) == 1  # cromossomo como inteiro
    assert GF.carregar_regiao("1", 300, 400) == []


def test_md5_em_fluxo(gwas_sintetico):
    assert GF.md5_arquivo(gwas_sintetico, bloco=7) == \
        hashlib.md5(gwas_sintetico.read_bytes()).hexdigest()


def test_md5_e_url_fixados():
    assert GF.GWAS_MD5 == "071790d80ccae0b41adabc7e0eefbf53"
    assert GF.GWAS_URL.startswith("https://ftp.ebi.ac.uk/pub/databases/gwas/")
    assert GF.GWAS_URL.endswith("GCST90001389_buildGRCh37.tsv.gz")


def test_scripts_nao_abrem_o_gwas_direto():
    """Todo script que lê o GWAS passa por gwas_fonte (mensagem clara e MD5)."""
    for arq in (BASE / "pipeline" / "scripts").glob("*.py"):
        txt = arq.read_text(encoding="utf-8")
        if arq.name == "obter_gwas.py":
            continue
        assert "GCST90001389_buildGRCh37" not in txt, \
            f"{arq.name} referencia o arquivo do GWAS sem passar por gwas_fonte"


def test_gwas_continua_fora_do_versionamento():
    rastreados = subprocess.run(["git", "ls-files"], cwd=BASE, text=True,
                                capture_output=True).stdout.splitlines()
    assert not any("GCST90001389_buildGRCh37.tsv.gz" in r for r in rastreados)


# ------------------------------------------------------------- D3: lambda_GC
def test_lambda_gc_constante_unica():
    assert GF.MEDIANA_CHI2_1GL == pytest.approx(float(chi2.ppf(0.5, 1)), abs=5e-7)
    valores = sorted(float(chi2.ppf(q / 1001, 1)) for q in range(1, 1001))
    assert GF.lambda_gc(valores) == pytest.approx(1.0, abs=5e-3)
    assert GF.lambda_gc([v * 1.0587 for v in valores]) == pytest.approx(
        1.0587, abs=6e-3)
    for nome in ("qc_gwas_gcst90001389.py", "finemap_ld.py"):
        txt = (BASE / "pipeline" / "scripts" / nome).read_text(encoding="utf-8")
        assert "MEDIANA_CHI2_1GL" in txt
        assert not re.search(r"/\s*0\.4549\b", txt), f"{nome}: constante solta"


def test_lambda_gc_igual_nos_dois_relatorios():
    qc = (BASE / "pipeline/reports/relatorio_qc_gwas_gcst90001389.md").read_text(
        encoding="utf-8")
    lam = (BASE / "pipeline/reports/relatorio_lambda_gc.md").read_text(
        encoding="utf-8")
    a = re.search(r"lambda_GC = ([\d.]+)", qc).group(1)
    b = re.search(r"λ global[^*]*\*\*([\d.]+)\*\*", lam).group(1)
    assert a == b == "1.0587"


# ------------------------------------------------------------ R2: dependências
def _imports_externos():
    stdlib = set(sys.stdlib_module_names)
    locais = {p.stem for d in ("pipeline/lib", "pipeline/scripts", "tests")
              for p in (BASE / d).glob("*.py")}
    achados = set()
    for d in ("pipeline/lib", "pipeline/scripts", "tests"):
        for arq in (BASE / d).glob("*.py"):
            for no in ast.walk(ast.parse(arq.read_text(encoding="utf-8"))):
                if isinstance(no, ast.Import):
                    achados |= {a.name.split(".")[0] for a in no.names}
                elif isinstance(no, ast.ImportFrom) and no.module and no.level == 0:
                    achados.add(no.module.split(".")[0])
    return achados - stdlib - locais


def test_requirements_fixa_tudo_o_que_o_codigo_importa():
    linhas = [l.split("#")[0].strip() for l in
              (BASE / "requirements.txt").read_text(encoding="utf-8").splitlines()]
    linhas = [l for l in linhas if l]
    fixados = {}
    for l in linhas:
        m = re.fullmatch(r"([A-Za-z0-9_.-]+)==([0-9][0-9A-Za-z.]*)", l)
        assert m, f"dependência sem versão fixada (use ==): {l!r}"
        fixados[m.group(1).lower()] = m.group(2)
    externos = {m.lower() for m in _imports_externos()}
    assert externos <= set(fixados), \
        f"importados mas ausentes do requirements.txt: {sorted(externos - set(fixados))}"


def test_ci_usa_a_mesma_fonte():
    ci = (BASE / ".github/workflows/ci.yml").read_text(encoding="utf-8")
    assert "pip install -r requirements.txt" in ci
    assert not re.search(r"pip install (?!-r requirements\.txt)", ci), \
        "CI instala pacotes fora do requirements.txt"
    mk = (BASE / "Makefile").read_text(encoding="utf-8")
    assert "pip install -r requirements.txt" in mk
    assert "--break-system-packages" not in mk


# ------------------------------------------------------------- R3 / R5: Makefile
def test_make_clean_nao_apaga_nada_versionado():
    mk = (BASE / "Makefile").read_text(encoding="utf-8")
    receita = mk.split("\nclean:\n", 1)[1]
    assert "eqtl_" not in receita and "sqtl_" not in receita
    r = subprocess.run(["make", "-n", "clean"], cwd=BASE, text=True,
                       capture_output=True)
    assert r.returncode == 0, r.stderr
    rastreados = subprocess.run(["git", "ls-files"], cwd=BASE, text=True,
                                capture_output=True).stdout.splitlines()
    alvos = re.findall(r"rm -r?f ([^\n]+)", r.stdout)
    padroes = [p for linha in alvos for p in linha.split()]
    assert padroes, "receita de clean sem rm (o teste perdeu o alvo)"
    import fnmatch
    for pad in padroes:
        atingidos = [f for f in rastreados if fnmatch.fnmatch(f, pad)
                     or f.startswith(pad.rstrip("/") + "/")]
        assert not atingidos, f"`make clean` apagaria versionados: {atingidos[:5]}"
    assert not [f for f in rastreados if "__pycache__" in f]


def test_sem_indices_tbi_versionados():
    rastreados = subprocess.run(["git", "ls-files"], cwd=BASE, text=True,
                                capture_output=True).stdout.splitlines()
    assert not [f for f in rastreados if f.endswith(".tbi")]


def test_makefile_tem_alvo_que_regenera_todos_os_relatorios():
    """Cada relatório com script gerador é coberto por `make relatorios`."""
    r = subprocess.run(["make", "-n", "relatorios"], cwd=BASE, text=True,
                       capture_output=True)
    assert r.returncode == 0, r.stderr
    chamados = set(re.findall(r"pipeline/scripts/(\w+\.py)", r.stdout))
    esperados = {"analise_gse160208.py", "analise_gse140069.py",
                 "analise_gse140069_idade.py",
                 "analise_integrada_cerebro_sangue.py", "gera_figuras.py",
                 "qc_gwas_gcst90001389.py", "clumping_descoberta.py",
                 "coloc_stx6_eqtl.py", "coloc_meta_stx6.py",
                 "coloc_sqtl_stx6.py", "gera_figuras_v2.py", "finemap_ld.py",
                 "valida_offset_grch38_grch37.py", "crosscheck_coloc_R.py"}
    assert esperados <= chamados, sorted(esperados - chamados)
    for s in chamados:
        assert (BASE / "pipeline" / "scripts" / s).exists(), s


# ----------------------------------- relatórios versionados × scripts (R6)
def _sem_data(txt):
    return re.sub(r"\b20\d\d-\d\d-\d\d\b", "<data>", txt)


@pytest.mark.parametrize("modulo,relatorio", [
    ("analise_gse160208", "relatorio_gse160208.md"),
    ("analise_gse140069", "relatorio_gse140069.md"),
    ("analise_gse140069_idade", "relatorio_gse140069_idade.md"),
    ("analise_integrada_cerebro_sangue", "relatorio_integracao_cerebro_sangue.md"),
])
def test_relatorio_versionado_e_o_que_o_script_gera(modulo, relatorio, tmp_path,
                                                    monkeypatch, capsys):
    """Roda o script e compara a saída, linha a linha, com o relatório
    versionado (ignorando só a data). Pega relatório editado à mão depois de
    gerado e relatório esquecido depois de mudar o script — os dois casos
    aconteceram neste repositório."""
    import importlib
    sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
    mod = importlib.import_module(modulo)
    monkeypatch.setattr(mod, "REPORTS", tmp_path)
    mod.main()
    capsys.readouterr()
    gerado = _sem_data((tmp_path / relatorio).read_text(encoding="utf-8"))
    versionado = _sem_data((BASE / "pipeline" / "reports" / relatorio).read_text(
        encoding="utf-8"))
    assert gerado == versionado, (
        f"{relatorio} está defasado em relação a {modulo}.py — "
        "rode `make relatorios-offline`")


def test_relatorio_de_validacao_e_o_que_o_script_gera(tmp_path, monkeypatch, capsys):
    """Idem para a validação da colocalização (modo sem GWAS e sem R)."""
    sys.path.insert(0, str(BASE / "pipeline" / "scripts"))
    import crosscheck_coloc_R as X
    monkeypatch.setattr(X, "REPORTS_SAIDA", tmp_path, raising=False)
    monkeypatch.setattr(sys, "argv", ["crosscheck_coloc_R.py", "--so-oraculo"])
    X.main()
    capsys.readouterr()
    nome = "relatorio_validacao_coloc_R.md"
    assert _sem_data((tmp_path / nome).read_text(encoding="utf-8")) == _sem_data(
        (BASE / "pipeline" / "reports" / nome).read_text(encoding="utf-8"))
