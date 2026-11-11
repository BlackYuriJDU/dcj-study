# -*- coding: utf-8 -*-
"""
test_coloc_abf.py — Testes que CALCULAM a colocalização.

Três camadas, todas capazes de falhar quando a implementação está errada:

1. Casos sintéticos de resposta conhecida: as posteriores são conferidas contra
   a enumeração explícita de todas as configurações causais (definição do
   modelo de Giambartolomei 2014), não contra uma segunda cópia da fórmula.
2. Oráculo: `tests/fixtures/stx6_crosscheck_*_v04.tsv` guarda a entrada de 390
   variantes e o `SNP.PP.H4` que o `coloc.abf` 5.2.3 real (R) devolveu para ela.
3. Âncoras numéricas da auditoria de 2026-11 (dossiê, seção 2.1).

Histórico: até a v0.4 o termo de H4 era multiplicado por p1·p2·p12 (em vez de
p12) e o fator de Bayes usava sqrt(r) (em vez de sqrt(1−r)). O teste antigo
apenas conferia se "0.9950" continuava escrito no relatório, e por isso travava
o resultado errado. Qualquer um dos dois erros derruba os testes abaixo.
"""
import csv
import itertools
import math
import random
import re
import sys
from pathlib import Path

import pytest

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "pipeline" / "lib"))
import coloc_abf as C  # noqa: E402

FIX = BASE / "tests" / "fixtures"
PADRAO = (1e-4, 1e-4, 1e-5)


# ------------------------------------------------------------------ utilitários
def enumeracao_explicita(l1, l2, p1, p2, p12):
    """Posteriores H0..H4 pela DEFINIÇÃO do modelo: soma sobre configurações.

    No máximo uma variante causal por traço. Cada configuração tem prior
    p1 (só traço 1 na variante i), p2 (só traço 2 na variante j), p1·p2
    (traço 1 em i e traço 2 em j, i≠j) ou p12 (os dois na mesma variante i),
    e verossimilhança relativa igual ao produto dos fatores de Bayes.
    """
    n = len(l1)
    b1 = [math.exp(x) for x in l1]
    b2 = [math.exp(x) for x in l2]
    h0 = 1.0
    h1 = sum(p1 * b1[i] for i in range(n))
    h2 = sum(p2 * b2[j] for j in range(n))
    h3 = sum(p1 * p2 * b1[i] * b2[j]
             for i, j in itertools.product(range(n), repeat=2) if i != j)
    h4 = sum(p12 * b1[i] * b2[i] for i in range(n))
    s = h0 + h1 + h2 + h3 + h4
    return [h / s for h in (h0, h1, h2, h3, h4)]


def ler_oraculo():
    with (FIX / "stx6_crosscheck_input_v04.tsv").open() as fh:
        ent = list(csv.DictReader(fh, delimiter="\t"))
    with (FIX / "stx6_crosscheck_colocR_output_v04.tsv").open() as fh:
        sai = list(csv.DictReader(fh, delimiter="\t"))
    assert [r["pos_b37"] for r in ent] == [r["pos_b37"] for r in sai]
    gwas = [(float(r["beta_gwas"]), float(r["se_gwas"])) for r in ent]
    eqtl = [(float(r["beta_eqtl_meta"]), float(r["se_eqtl_meta"])) for r in ent]
    pos = [int(r["pos_b37"]) for r in ent]
    pp4_r = [float(r["PP4_colocR"]) for r in sai]
    return pos, gwas, eqtl, pp4_r


@pytest.fixture(scope="module")
def oraculo():
    return ler_oraculo()


@pytest.fixture(scope="module")
def res_oficial(oraculo):
    _, gwas, eqtl, _ = oraculo
    return C.coloc(gwas, eqtl, tipo1="cc", tipo2="quant", sdY2=1.0)


# ------------------------------------------------------- 1. fator de Bayes
def test_log_abf_e_o_de_wakefield():
    """ABF = sqrt(V/(V+W)) · exp(z²/2 · W/(V+W)). Conferido em forma fechada."""
    for beta, se, w in [(0.2, 0.05, 0.04), (-0.07, 0.012, 0.0225),
                        (0.0, 0.3, 0.04), (1.3, 0.9, 0.0016)]:
        v = se * se
        esperado = math.sqrt(v / (v + w)) * math.exp(
            (beta / se) ** 2 / 2 * w / (v + w))
        assert math.exp(C.log_abf(beta, se, w)) == pytest.approx(esperado,
                                                                 rel=1e-12)


def test_abf_sem_efeito_nao_favorece_associacao():
    """Com beta = 0, o fator de Bayes tem de ser < 1 (evidência contra).
    A forma errada sqrt(r) também dá < 1, mas com o valor trocado: com
    V = W o correto é sqrt(1/2); com V = W/3 o correto é 1/2 (o errado dá
    sqrt(3)/2). Este é o caso que separa sqrt(1−r) de sqrt(r)."""
    w = 0.04
    assert math.exp(C.log_abf(0.0, math.sqrt(w), w)) == pytest.approx(
        math.sqrt(0.5), rel=1e-12)
    assert math.exp(C.log_abf(0.0, math.sqrt(w / 3), w)) == pytest.approx(
        0.5, rel=1e-12)


def test_abf_de_estudo_sem_informacao_tende_a_um():
    """Erro-padrão enorme = dado não informativo: ABF → 1 (log ABF → 0).
    Na forma errada, sqrt(r) → 0, o que 'provaria' a hipótese nula sem dado."""
    assert C.log_abf(0.5, 1e4, 0.04) == pytest.approx(0.0, abs=1e-8)


def test_priors_de_efeito_sao_os_do_pacote_coloc():
    """coloc: desvio-padrão 0,2 (caso-controle) e 0,15·sdY (quantitativo)."""
    assert C.prior_w("cc") == pytest.approx(0.04)
    assert C.prior_w("quant") == pytest.approx(0.0225)
    assert C.prior_w("quant", sdY=2.0) == pytest.approx(0.09)
    with pytest.raises(ValueError):
        C.prior_w("binario")


def test_finemap_usa_o_mesmo_fator_de_bayes():
    """Fine-mapping e colocalização compartilham uma única função de ABF."""
    txt = (BASE / "pipeline" / "scripts" / "finemap_ld.py").read_text(
        encoding="utf-8")
    assert "from coloc_abf import" in txt, "finemap_ld.py não importa a lib única"
    for arq in (BASE / "pipeline" / "scripts").glob("*.py"):
        assert not re.search(r"^def abf\(", arq.read_text(encoding="utf-8"),
                             re.M), f"{arq.name} mantém cópia local de abf()"


# ------------------------------------------- 2. combinação: resposta conhecida
def test_uma_variante_forma_fechada():
    """n = 1: H3 é impossível e H4 = p12·B1·B2 / (1 + p1·B1 + p2·B2 + p12·B1·B2)."""
    l1, l2 = [math.log(50.0)], [math.log(20.0)]
    p1, p2, p12 = PADRAO
    pp = C.combinar(l1, l2, p1, p2, p12)
    den = 1 + p1 * 50 + p2 * 20 + p12 * 1000
    assert pp[3] == 0.0
    assert pp[4] == pytest.approx(p12 * 1000 / den, rel=1e-12)
    assert pp[1] == pytest.approx(p1 * 50 / den, rel=1e-12)
    assert pp[2] == pytest.approx(p2 * 20 / den, rel=1e-12)
    assert sum(pp) == pytest.approx(1.0, abs=1e-12)


def test_combinacao_bate_com_enumeracao_explicita():
    """Posteriores == soma sobre todas as configurações causais (definição)."""
    rng = random.Random(20261110)
    for n in (2, 3, 7, 25):
        for _ in range(20):
            l1 = [rng.uniform(-3, 18) for _ in range(n)]
            l2 = [rng.uniform(-3, 18) for _ in range(n)]
            for priors in (PADRAO, (1e-4, 1e-4, 1e-6), (1e-3, 1e-5, 5e-6)):
                obtido = C.combinar(l1, l2, *priors)
                esperado = enumeracao_explicita(l1, l2, *priors)
                for a, b in zip(obtido, esperado):
                    assert a == pytest.approx(b, rel=1e-9, abs=1e-15)


def test_sinal_compartilhado_da_h4():
    """Mesmo causal nos dois traços, vizinhos sem sinal → H4 domina."""
    n = 50
    l1 = [-1.0] * n
    l2 = [-1.0] * n
    l1[10] = l2[10] = 30.0
    pp = C.combinar(l1, l2, *PADRAO)
    assert pp[4] > 0.99


def test_sinais_em_variantes_distintas_dao_h3():
    """Causais em variantes diferentes (sem LD) → H3 domina, H4 ≈ 0."""
    n = 50
    l1 = [-1.0] * n
    l2 = [-1.0] * n
    l1[10] = 30.0
    l2[40] = 30.0
    pp = C.combinar(l1, l2, *PADRAO)
    assert pp[3] > 0.99
    assert pp[4] < 1e-6


def test_sinal_em_um_so_traco_da_h1_ou_h2():
    n = 50
    forte = [-1.0] * n
    forte[5] = 30.0
    fraco = [-3.0] * n        # evidência contra associação em toda variante
    assert C.combinar(forte, fraco, *PADRAO)[1] > 0.99      # H1: só traço 1
    assert C.combinar(fraco, forte, *PADRAO)[2] > 0.99      # H2: só traço 2


def test_sem_sinal_da_h0():
    n = 50
    pp = C.combinar([-1.0] * n, [-1.0] * n, *PADRAO)
    assert pp[0] > 0.99


def test_p12_entra_sozinho_no_prior_de_h4():
    """Com um só causal compartilhado, a razão H4/H3 a posteriori tem de ser
    p12/(p1·p2) vezes a razão de verossimilhanças. Multiplicar H4 por
    p1·p2·p12 (erro da v0.4) encolhe H4 em 10⁸ e inverte este caso."""
    n = 10
    l1 = [0.0] * n
    l2 = [0.0] * n            # todas as variantes equivalentes, ABF = 1
    p1, p2, p12 = PADRAO
    pp = C.combinar(l1, l2, p1, p2, p12)
    # n configurações de H4 (prior p12 cada) e n(n−1) de H3 (prior p1·p2 cada)
    assert pp[4] / pp[3] == pytest.approx(p12 / (p1 * p2 * (n - 1)), rel=1e-9)


def test_espaco_logaritmico_nao_estoura():
    """z = 60 nos dois traços: exp() em espaço linear estouraria o float."""
    l1 = [C.log_abf(6.0, 0.1, 0.04)] + [0.0] * 5      # log ABF ≈ 1439
    l2 = [C.log_abf(6.0, 0.1, 0.0225)] + [0.0] * 5
    assert l1[0] > 709                                 # exp() estouraria
    pp = C.combinar(l1, l2, *PADRAO)
    assert all(math.isfinite(x) for x in pp)
    assert pp[4] == pytest.approx(1.0, abs=1e-6)


def test_snp_pp_h4_soma_um_e_segue_o_produto():
    l1 = [1.0, 5.0, 2.0]
    l2 = [0.5, 4.0, 6.0]
    pp = C.snp_pp_h4(l1, l2)
    assert sum(pp) == pytest.approx(1.0, abs=1e-12)
    bruto = [math.exp(a + b) for a, b in zip(l1, l2)]
    for a, b in zip(pp, bruto):
        assert a == pytest.approx(b / sum(bruto), rel=1e-12)


def test_h4_condicional_e_h4_sobre_h3_mais_h4():
    """O rótulo 'H4/(H3+H4)' tem de ser a conta H4/(H3+H4), não H4/H3."""
    pp = [0.0, 0.0, 0.0, 0.2, 0.6]
    assert C.h4_condicional(pp) == pytest.approx(0.75)


def test_conjunto_credivel():
    pp = [0.5, 0.3, 0.15, 0.04, 0.01]
    assert C.conjunto_credivel(pp, 0.95) == [0, 1, 2]
    assert C.conjunto_credivel(pp, 0.96) == [0, 1, 2, 3]


def test_tamanhos_diferentes_sao_erro():
    with pytest.raises(ValueError):
        C.combinar([0.0, 1.0], [0.0], *PADRAO)
    with pytest.raises(ValueError):
        C.combinar([], [], *PADRAO)


# ------------------------------------------------ 3. oráculo: coloc.abf real
def test_oraculo_snp_pp_h4_bate_com_coloc_r(oraculo, res_oficial):
    """SNP.PP.H4 variante a variante == saída do coloc.abf 5.2.3 (R)."""
    _, _, _, pp4_r = oraculo
    assert len(pp4_r) == 390
    dif = max(abs(a - b) for a, b in zip(res_oficial["snp_pp_h4"], pp4_r))
    assert dif < 1e-12, f"diferença máxima {dif:.3g} contra o coloc.abf real"


def test_oraculo_resumo_bate_com_coloc_r(res_oficial):
    """Resumo registrado para o pacote: H4 = 0,9796; H3 = 0,0204; H2 = 7e-6."""
    pp = res_oficial["pp"]
    assert pp[4] == pytest.approx(0.9796, abs=5e-5)
    assert pp[3] == pytest.approx(0.0204, abs=5e-5)
    assert pp[2] == pytest.approx(7e-6, abs=5e-7)
    assert pp[0] < 1e-40 and pp[1] < 1e-40


# --------------------------------------- 4. âncoras da auditoria (seção 2.1)
def test_ancora_oficial_stx6(res_oficial):
    pp = res_oficial["pp"]
    assert pp[2] == pytest.approx(7.04298e-06, rel=1e-5)
    assert pp[3] == pytest.approx(0.020388, abs=5e-7)
    assert pp[4] == pytest.approx(0.979605, abs=5e-7)
    assert C.h4_condicional(pp) == pytest.approx(0.979612, abs=5e-7)


@pytest.mark.parametrize("priors,h3,h4", [
    ((1e-4, 1e-4, 1e-5), 0.020388, 0.979605),
    ((1e-4, 1e-4, 5e-6), 0.039962, 0.960025),
    ((1e-4, 1e-4, 1e-6), 0.172263, 0.827677),
    ((1e-6, 1e-6, 1e-7), 0.000208, 0.999785),
])
def test_ancora_sensibilidade_p12(oraculo, priors, h3, h4):
    _, gwas, eqtl, _ = oraculo
    pp = C.coloc(gwas, eqtl, "cc", "quant", priors=priors)["pp"]
    assert pp[3] == pytest.approx(h3, abs=5e-7)
    assert pp[4] == pytest.approx(h4, abs=5e-7)


@pytest.mark.parametrize("w,h4", [
    (0.0016, 0.978676), (0.01, 0.979530), (0.0225, 0.979584),
    (0.04, 0.979601), (0.1, 0.979611),
])
def test_ancora_sensibilidade_w(oraculo, w, h4):
    _, gwas, eqtl, _ = oraculo
    pp = C.coloc(gwas, eqtl, w1=w, w2=w)["pp"]
    assert pp[4] == pytest.approx(h4, abs=5e-7)


def test_tabelas_de_sensibilidade_da_lib(oraculo):
    """As tabelas que os scripts emitem saem da lib e batem com as âncoras."""
    _, gwas, eqtl, _ = oraculo
    t12 = dict((p12, pp) for p12, pp in
               C.sensibilidade_p12(gwas, eqtl, "cc", "quant"))
    assert set(t12) >= {1e-5, 5e-6, 1e-6}
    assert t12[1e-6][4] == pytest.approx(0.827677, abs=5e-7)
    tw = dict(C.sensibilidade_w(gwas, eqtl))
    assert tw[0.04][4] == pytest.approx(0.979601, abs=5e-7)


def test_ancora_variantes_de_maior_posterior(oraculo, res_oficial):
    """A variante de maior SNP.PP.H4 NÃO é o lead do GWAS; conjunto de 95%
    tem 17 variantes."""
    pos, _, _, _ = oraculo
    pp = res_oficial["snp_pp_h4"]
    ordem = sorted(range(len(pp)), key=lambda i: -pp[i])
    assert [pos[i] for i in ordem[:4]] == [180951692, 180956015, 180956985,
                                           180961245]
    assert pp[ordem[0]] == pytest.approx(0.1048, abs=5e-5)
    assert pp[ordem[3]] == pytest.approx(0.0634, abs=5e-5)
    assert len(C.conjunto_credivel(pp, 0.95)) == 17


def test_ancora_concordancia_de_direcao(oraculo):
    _, gwas, eqtl, _ = oraculo
    conc = sum(1 for (bg, _), (be, _) in zip(gwas, eqtl) if bg * be > 0)
    assert (conc, len(gwas)) == (348, 390)


def test_diagnostico_prior_vs_posterior_da_razao_h4_h3(res_oficial):
    """G9: se as 390 variantes fossem equivalentes, a razão a priori H4/H3
    seria p12/(p1·p2·(n−1)) ≈ 2,57; a razão a posteriori observada é ≈ 48."""
    n = res_oficial["n"]
    p1, p2, p12 = PADRAO
    assert C.razao_h4_h3_priori(n, p1, p2, p12) == pytest.approx(2.5707, abs=5e-4)
    pp = res_oficial["pp"]
    assert pp[4] / pp[3] == pytest.approx(48.05, abs=0.05)


# ------------------------------- 5. uma única implementação, usada por todos
SCRIPTS_COLOC = ("coloc_stx6_eqtl.py", "coloc_meta_stx6.py",
                 "coloc_sqtl_stx6.py", "crosscheck_coloc_R.py")


def test_scripts_usam_a_implementacao_unica():
    """Nenhum script redefine a combinação de hipóteses nem o fator de Bayes."""
    for nome in SCRIPTS_COLOC:
        txt = (BASE / "pipeline" / "scripts" / nome).read_text(encoding="utf-8")
        assert "coloc_abf" in txt, f"{nome} não importa pipeline/lib/coloc_abf"
        assert not re.search(r"^\s*def (coloc|abf|h_post)\(", txt, re.M), \
            f"{nome} mantém cópia local de coloc()/abf()"


def test_nenhuma_copia_da_formula_antiga():
    """A combinação antiga (H4 × p1·p2·p12) e o ABF com sqrt(r) não podem
    reaparecer em código executável do pipeline (Python ou R)."""
    padroes = [r"p1\s*\*\s*p2\s*\*\s*p12",              # H4 × p1·p2·p12
               r"sqrt\(r\)\s*\*\s*(math\.)?exp",         # sqrt(r) em vez de sqrt(1−r)
               r"0\.04\s*(\*\*|\^)\s*2"]                 # W = 0,04² ("coloc padrão")
    arquivos = (list((BASE / "pipeline" / "scripts").glob("*.py"))
                + list((BASE / "pipeline" / "scripts").glob("*.R"))
                + list((BASE / "pipeline" / "lib").glob("*.py")))
    achados = []
    for arq in arquivos:
        txt = arq.read_text(encoding="utf-8")
        for pad in padroes:
            if re.search(pad, txt):
                achados.append((arq.name, pad))
    assert not achados, f"fórmula antiga presente: {achados}"


def test_descricao_das_hipoteses():
    """H2 é 'associação só com o traço 2' — não 'duas variantes distintas'."""
    h = C.HIPOTESES
    assert set(h) == {"H0", "H1", "H2", "H3", "H4"}
    assert "só" in h["H1"] and "só" in h["H2"]
    assert "distint" in h["H3"]
    assert "compartilhad" in h["H4"]
    assert "distint" not in h["H2"]
