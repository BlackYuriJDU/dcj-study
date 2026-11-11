# -*- coding: utf-8 -*-
"""test_consistencia_docs.py — turning the project's consistency rule into a CI gate.

Project rule: after correcting a script, regenerate the derived reports. A stale
report plus a new script is a silent inconsistency. This file turns that rule —
and, since 2026-10-04, a set of scientific-integrity invariants — into gates that
fail the build.

Text files only (tracked in git); large datasets stay out of CI by design.
"""
import json
import re
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
REPORTS = BASE / "pipeline" / "reports"
FIGS = REPORTS / "figuras"
SCRIPTS = BASE / "pipeline" / "scripts"


def ler(rel: str) -> str:
    p = BASE / rel
    assert p.exists(), f"arquivo ausente: {rel}"
    return p.read_text(encoding="utf-8")


# ---------------------------------------------------------------- GSE140069
def test_relatorio_gse140069_e_v3():
    txt = ler("pipeline/reports/relatorio_gse140069.md")
    assert "v3" in txt.splitlines()[0] + txt.splitlines()[1], \
        "relatório GSE140069 não é v3"
    script = ler("pipeline/scripts/analise_gse140069.py")
    assert "v3" in script[:1200], "script analise_gse140069.py não declara v3"
    for numero in ("84", "**1**", "0.048"):
        assert numero in txt, f"número-chave v3 ausente no relatório: {numero}"


# ---------------------------------------------------------------- coloc STX6
# O teste que havia aqui (test_coloc_meta_tem_h2_e_validacao) exigia as strings
# "0.9950" e "0.994997" nos relatórios: travava o resultado de uma fórmula
# errada (auditoria 2026-11, achado T1). Foi substituído por testes que
# CALCULAM a colocalização e exigem que os relatórios citem o valor calculado:
# tests/test_coloc_abf.py (sintéticos, oráculo do coloc.abf real, âncoras) e
# tests/test_coloc_pipeline.py (harmonização → meta → colocalização).
def test_relatorios_de_coloc_nao_citam_o_resultado_da_formula_antiga():
    for rel in ("pipeline/reports/relatorio_coloc_meta_stx6.md",
                "pipeline/reports/relatorio_coloc_stx6.md",
                "pipeline/reports/relatorio_coloc_sqtl_stx6.md"):
        txt = ler(rel)
        assert "0.9950" not in txt and "0.6654" not in txt, \
            f"{rel} ainda cita H3 da fórmula antiga"
        assert "Bonferroni: limiar H4" not in txt


def test_validacao_eur_existe():
    txt = ler("pipeline/reports/relatorio_finemap_loci_EUR.md")
    assert "phase_3:EUR" in txt
    assert "90.5%" in txt, "massa STX6 (90.5%) ausente no run EUR"


def test_offset_genomico_verificado_sem_exagero():
    """O offset GRCh38→GRCh37 de STX6 foi medido em 2026-10-04 contra a API de
    assembly map do Ensembl (relatorio_validacao_offset_genomico.md).

    Isso torna FALSOS os dois extremos que o texto já afirmou em versões
    diferentes: "não foi validado" (docstring de 2026-10-03) e "validado por
    Ensembl MAP" (manuscrito, que sugeria uma validação de cadeia inexistente —
    o parser de chain file próprio foi escrito, provado defeituoso e
    descartado). O que é verdade hoje é estreito: verificação por assembly map,
    sem cadeia. O gate prende os dois lados.
    """
    alvos = [f"pipeline/scripts/{n}" for n in
             ("coloc_stx6_eqtl.py", "coloc_meta_stx6.py",
              "coloc_sqtl_stx6.py", "crosscheck_coloc_R.py")]

    # (1) não afirmar validação por chain/liftOver — nunca foi feita
    for rel in alvos:
        txt = ler(rel)
        for padrao in (r"validad\w*.{0,60}(chain|liftOver)",
                       r"(chain|liftOver).{0,60}validad"):
            assert not re.search(padrao, txt, re.I), \
                f"{rel} afirma validação por chain/liftOver, que não ocorreu"

    # (2) não continuar chamando de suposição não validada — deixou de ser verdade
    for rel in alvos:
        txt = ler(rel)
        assert not re.search(r"n[aã]o (foi )?validad|n[aã]o validada",
                             txt, re.I), \
            f"{rel} ainda chama o offset de não validado"

    # (3) o manuscrito não pode superdeclarar: assembly map não é liftOver
    ms = ler("preprint/manuscrito_preprint.md")
    assert not re.search(r"validad\w*.{0,60}(chain|liftOver)", ms, re.I), \
        "manuscrito afirma validação por chain/liftOver, que não ocorreu"

    # (4) a medição existe e registra o veredito
    rel = ler("pipeline/reports/relatorio_validacao_offset_genomico.md")
    for marca in ("+30,864", "30/30"):
        assert marca in rel, f"marcador da validação ausente no relatório: {marca}"
    # todos os pontos amostrados concordam (n/n), qualquer que seja o n usado
    # na execução (o padrão do script é 41; a versão de 2026-10-04 usou 9)
    m = re.search(r"em (\d+)/(\d+) pontos amostrados", rel)
    assert m and m.group(1) == m.group(2) and int(m.group(1)) >= 9, \
        "relatório do offset não registra concordância em todos os pontos"

    # (5) o manuscrito reflete o estado verificado
    assert "verified this offset against the Ensembl" in ms, \
        "manuscrito não registra a verificação do offset"


def test_cache_ensembl_consistente_com_a_constante():
    """O cache em disco guarda pares b37→b38 do Ensembl; o script regenerate o
    relatório a partir dele.

    Bug real de 2026-10-04: o cache continha pares com delta +30.865 (um +1 a
    mais). Como §1 calcula `p38 - p37` direto, um re-run teria reescrito o
    veredito do relatório de "+30.864, correto" para "+30.865, incorreto" sem
    erro, sem aviso e sem rede — o ensaio inteiro passaria. E o par errado era
    plausível à vista: 180941860→180972725 é o mapeamento *correto* de
    180941861, então a inspeção visual não pegaria.

    O cache é descartável (o script o reconstrói); se um dia for commitado, tem
    de concordar com a constante.
    """
    p = BASE / "pipeline" / "data" / "offset_ensembl_cache.json"
    if not p.exists():
        return                      # cache opcional; o script reconstrói
    d = json.loads(p.read_text(encoding="utf-8"))
    for k, v in d.items():
        assert v is not None, f"cache Ensembl com falha em b37={k}"
        delta = v - int(k)
        assert delta == 30_864, (
            f"cache Ensembl inconsistente com OFFSET_CONSTANTE: b37 {k} -> "
            f"b38 {v} (delta {delta}, esperado 30864). Um re-run trocaria o "
            f"veredito do relatório em silêncio.")


# ---------------------------------------------------------------- figuras
def test_figuras_principais_existem():
    for png in ("volcano_gse160208.png", "volcano_gse140069.png",
                "heatmap_top_genes.png", "timeline_caso_referencia.png",
                "manhattan_gwas.png", "forest_mirnas.png",
                "coloc_stx6_regional.png"):
        assert (FIGS / png).exists(), f"figura ausente: {png}"


def test_sem_figura_legada_com_nome_de_pessoa():
    """A linha do tempo de um caso real foi substituída por uma versão genérica.

    O guard não hardcoda o nome da pessoa (seria reintroduzi-la no repo): ele
    exige que a figura saneada exista e seja a ÚNICA da série timeline, de modo
    que a versão pessoal não possa voltar sem que este teste caia.
    """
    series = sorted(FIGS.glob("timeline_*.png"))
    assert [p.name for p in series] == ["timeline_caso_referencia.png"], (
        f"série timeline deve conter só a figura saneada; encontrada: "
        f"{[p.name for p in series]}")


# ---------------------------------------------------------------- privacidade
# O termo é montado em duas partes de propósito. Escrever o nome inteiro aqui
# publicaria o identificador no próprio guard — e o repositório é público. Com o
# termo em partes, o arquivo passa a varrer a si mesmo, e não precisa de
# autoexclusão (que antes era um furo: qualquer arquivo poderia se declarar
# isento nomeando-se).
_TERMO = "li" + "to"

# Fila única, resolvida por caminho — não por conteúdo. 2026-10-05: o montador
# do dossiê privado saiu do versionamento (ver .gitignore), então a fila está
# vazia. Ela existe para o gate não precisar de autoexclusão: um arquivo nunca
# se declara isento pelo próprio nome, e qualquer exclusão futura fica escrita
# aqui, visível em review.
PENDENTES: set[str] = set()


def test_nome_do_caso_ausente_do_repo():
    """O nome do caso real não pode reaparecer em arquivo versionado.

    Motivo histórico: a lista de renomeação de 2026-10-04 cobriu manuscrito,
    CITATION.cff, LICENSE e READMEs, mas não varreu os relatórios do caso de
    referência nem as linhas de proveniência, que citavam scripts já renomeados.
    O guard é por substantivo, para pegar o nome e o slug antigo.
    """
    # Lookarounds em vez de \b: \b trata "_" como caractere de PALAVRA, então a
    # forma que o nome assume em nome de arquivo (carta_<termo>.md) escapava do
    # gate. Delimitar por qualquer caractere não-alfanumérico pega "_", "-" e "/".
    exato = re.compile(rf"(?<![0-9A-Za-z]){_TERMO}(?![0-9A-Za-z])"
                       rf"|(?<![0-9A-Za-z])dcj ?- ?{_TERMO}(?![0-9A-Za-z])",
                       re.I)
    # Pré-filtro barato: todo caso do padrão exato contém o termo, então isto é
    # um superconjunto — não pode haver falso negativo. Feito com lower()+`in`
    # (C puro) e não com re.I, que faz case-folding caractere a caractere. O
    # dataset de miRTarBase são 25 MB versionados e domina este gate.
    extensoes = (".md", ".cff", ".py", ".sh", ".yml", ".yaml", ".txt",
                 ".csv", ".R", ".html", ".js", ".css", ".svg", ".json",
                 ".gitignore")
    alvos = subprocess.run(
        ["git", "ls-files"], capture_output=True, text=True,
        encoding="utf-8", cwd=BASE).stdout.splitlines()
    # Só texto. Binários (.pdf, .gz, .xlsx, .png, .tbi) ficam de fora: o título
    # do PDF vai em metadata comprimida e se verifica descompactando; o resto
    # não tem onde esconder texto legível.
    #
    # Leitura em blocos, não read_text(). O gargalo é um arquivo só:
    # pipeline/data/mirtarbase_MTI.csv, 25,50 MB — 622x o maior arquivo de
    # texto do repo (Litho Foundation Website/app.js, 40 KB) e 49x todos os
    # outros arquivos de texto somados (81 arquivos, 513 KB). Slurpá-lo levava
    # o pico do gate a 403 MB (tracemalloc: 128 MB no read_text, +276 MB no
    # .lower() sobre a string UCS-2 de 50 MB) e 728 ms. Em blocos de 1 MB o
    # pico cai para 19 MB e o tempo para 293 ms — melhor nos dois eixos, e o
    # pré-filtro `in` roda sobre a string em cache em vez de realocar o
    # arquivo inteiro a cada execução da suíte.
    #
    # A cauda de sobreposição é o que torna a troca segura: qualquer casamento
    # do padrão cabe em ~12 caracteres (dcj ?- ?TERMO, mais os lookarounds), e
    # `janela` sempre carrega os últimos 64 da leitura anterior, então um
    # padrão cortado na divisa entre dois blocos é remontado inteiro antes do
    # regex. Verificado com blocos de 64 B sobre 13 casos sintéticos,
    # incluindo termo pisado na divisa: as três estratégias (slurp, linha a
    # linha, bloco) concordam em todos.
    #
    # Registro honesto: a versão anterior deste comentário afirmava que a
    # versão slurp "caía com MemoryError dentro da suíte completa". Não
    # reproduzi em 15 execuções da suíte nesta máquina (36 passed em todas).
    # A troca se sustenta pelo ganho medido de memória e tempo, não pelo
    # crash — que pode ter vindo de pressão de outra máquina.
    TAMANHO_BLOCO = 1 << 20  # 1 MB
    SOBRA = 64              # > 2x o comprimento máximo de um casamento
    vazamentos = []
    for rel in alvos:
        if rel in PENDENTES or not rel.endswith(extensoes):
            continue
        try:
            with (BASE / rel).open(encoding="utf-8") as fh:
                cauda = ""
                while bloco := fh.read(TAMANHO_BLOCO):
                    janela = cauda + bloco
                    if _TERMO in janela.lower() and exato.search(janela):
                        vazamentos.append(rel)
                        break
                    cauda = janela[-SOBRA:]
        except (UnicodeDecodeError, OSError):
            continue
    assert not vazamentos, (
        f"nome do caso real reapareceu em: {vazamentos}")


def test_nada_rastreado_esta_no_gitignore():
    """Um arquivo não pode estar versionado E ignorado ao mesmo tempo.

    A intuição é "gitignore protege", e aí está o problema: se o arquivo já
    está rastreado, o .gitignore não age sobre ele e o repo segue publicando
    o conteúdo. Foi assim com o montador do dossiê privado — entrou no
    .gitignore em 2026-10-05 e continuou no índice até `git rm --cached`.

    Por isso o gate compara as duas listas em vez de confiar no .gitignore.
    """
    ignorados = subprocess.run(
        ["git", "ls-files", "--ignored", "--exclude-standard", "--cached"],
        capture_output=True, text=True, encoding="utf-8", cwd=BASE).stdout
    vazando = sorted(p for p in ignorados.splitlines() if p.strip())
    assert not vazando, (
        f"versionados mas ignorados (o .gitignore não protege estes): {vazando}")


# ---------------------------------------------------------------- citável
def test_citation_tem_doi_e_autor():
    txt = ler("CITATION.cff")
    assert "10.5281/zenodo.23173672" in txt, "DOI ausente no CITATION.cff"
    assert "Araújo" in txt and "Arthur" in txt


def test_readme_badge_doi():
    for readme in ("README.md", "README.en.md"):
        txt = ler(readme)
        assert "zenodo.22164909" in txt, f"badge DOI ausente em {readme}"
    assert not re.search(r"60 miRNAs sig\b", ler("README.md")), \
        "número v1 (60 miRNAs sig) voltou ao README"


def test_doi_contaminado_nao_e_anunciado():
    """10.5281/zenodo.22164910 é a versão do Zenodo de 2026-08-29, anterior à
    limpeza de privacidade. O arquivo que ela serve tem 17 arquivos com o
    identificador do caso real, entre eles o dossiê clínico e o contexto do
    paciente.

    A causa foi de encadeamento, não de conteúdo: em 2026-10-04 a release foi
    re-cortada sob um tag novo, o que não fez o Zenodo mintir versão nova, e os
    badges continuaram anunciando o registro antigo. O repositório inteiro
    estava limpo enquanto as superfícies citáveis mandavam o leitor para o
    arquivo sujo — o vazamento que importa é o que continua resolvível.

    O gate é sobre o DOI errado, não sobre o certo: vigiar apenas que o DOI
    válido está presente não impede uma linha desatenta de apontar para o
    antigo, que é exatamente o que aconteceu.
    """
    contaminado = "22164910"
    for rel in ("CITATION.cff", "README.md", "README.en.md",
                "preprint/manuscrito_preprint.md",
                "Litho Foundation Website/index.html",
                "Litho Foundation Website/dados.html"):
        assert contaminado not in ler(rel), (
            f"{rel} anuncia o DOI contaminado (10.5281/zenodo.{contaminado}) "
            f"— resolve para o arquivo anterior à limpeza de privacidade")


# ------------------------------------------------- integridade científica
# Gates adicionados em 2026-10-04 após auditoria que encontrou alegações de
# ineditismo falsas e números atribuídos à variante errada. Cada regra abaixo
# existiu porque o erro correspondente já aconteceu uma vez.

def test_preprint_nao_afirma_ineditismo():
    """rs3747957 e a mediação de expressão por STX6 já foram publicadas [1,4].
    O manuscrito não pode reapresentá-las como descobertas próprias."""
    txt = ler("preprint/manuscrito_preprint.md")
    for proibido in ("previously unreported", "unremarked for five years",
                     "sleep in deposited data"):
        assert proibido not in txt, \
            f"alegação de ineditismo reintroduzida: {proibido!r}"
    # a ressalva do locus chr16 precisa continuar presente, em qualquer redação
    assert re.search(r"not a new locus|as an observation", txt), \
        "ressalva de que chr16 não é locus novo foi removida"


def test_preprint_cita_primeiro_autor_coreto():
    """Jones E é o primeiro autor; Mead S é o último dos 65 (autor sênior)."""
    txt = ler("preprint/manuscrito_preprint.md")
    assert "Jones E" in txt, "citação primária correta (Jones E) ausente"
    assert "Mead S" not in txt, "manuscrito voltou a citar 'Mead S' como primeiro autor"


def test_preprint_referencia_4_completa():
    txt = ler("preprint/manuscrito_preprint.md")
    assert "10.1093/brain/awaf032" in txt, \
        "referência 4 (Brain 2025) sem DOI verificado"


def test_p_eqtl_atribuido_a_variante_certa():
    """O lead (chr1:180,961,245) tem p = 7.60e-47. O valor 6.60e-47 pertence a
    outra variante (180,949,780, mínimo regional). Atribuir o menor ao lead é
    um erro que já ocorreu e foi corrigido uma vez."""
    txt = ler("preprint/manuscrito_preprint.md")
    assert "7.60×10⁻⁴⁷" in txt, "p do eQTL do lead (7.60e-47) ausente"
    # toda menção a 6.60e-47 precisa ser qualificada como mínimo regional, seja
    # antes ("the regional minimum ... reaches p = 6.60e-47") ou depois
    # ("a different variant, at chr1:180,949,780, with p = 6.60e-47")
    qualificadores = ("180,949,780", "regional minimum", "mínimo regional",
                      "different variant", "outra variante")
    for m in re.finditer(r"6\.60×10⁻⁴⁷", txt):
        janela = txt[max(0, m.start() - 220): m.end() + 220]
        assert any(q in janela for q in qualificadores), \
            f"6.60e-47 sem qualificação de mínimo regional: ...{janela[:120]}..."


def test_pvals_de_headline_preservados():
    txt = ler("preprint/manuscrito_preprint.md")
    for pv in ("1.62×10⁻¹⁵", "6.18×10⁻¹⁰", "7.51×10⁻⁹"):
        assert pv in txt, f"p-valor de headline ausente: {pv}"


def test_sqtl_tem_ressalva_de_poder():
    txt = ler("preprint/manuscrito_preprint.md")
    assert "not evidence of absence" in txt, \
        "ressalva de poder do sQTL removida (zero detecção ≠ ausência)"
    assert "exist in any cohort" not in txt, \
        "sQTL: 'exist' (prova negativa) reintroduzido"


def test_preprint_v04_e_assinado():
    txt = ler("preprint/manuscrito_preprint.md")
    assert "v0.4" in txt, "preprint não está em v0.4"
    assert "Arthur Araújo" in txt
    assert "84 FDR-significant miRNAs" in txt, \
        "wording FDR (correção 29/08) reverteu"


def test_confundimento_idade_atribuido():
    """Fase 2 (2026-10-04): o colapso 84→1 é atribuído a idade ~97% (88→4),
    RIN 62 (88→26), sexo 3 (88→85), com VIF < 1.4. O manuscrito não pode voltar
    a dizer "fully explains" (superdeclara) nem omitir a atribuição medida."""
    txt = ler("preprint/manuscrito_preprint.md")
    assert "fully explains" not in txt, \
        "\"fully explains\" reintroduzido: a idade responde por ~97%, não 100%"
    assert "97%" in txt, "atribuição medida (~97%) ausente do manuscrito"
    assert "variance-inflation" in txt, \
        "VIF não citado: sem ele a atribuição não é identificável"

    rel = ler("pipeline/reports/relatorio_gse140069_idade.md")
    for numero in ("113", "91", "88", "26", "1.37", "1.32"):
        assert numero in rel, f"número-chave do teste de idade ausente: {numero}"
    assert "não aditivas" in rel, \
        "contribuições de covariáveis apresentadas como aditivas (somariam >100%)"


# ---------------------------------------------------------------- scripts
def test_oraculos_estatisticos_disponiveis():
    """test_stats_core.py faz pytest.importorskip("scipy.stats") — sem scipy a
    regressão de welch/tcdf_p/fdr_bh/cohen_d contra as implementações de
    referência é silenciosamente pulada. Este gate falha alto em vez de deixar
    o build verde sem validar o núcleo estatístico."""
    faltando = []
    for mod, pkg in (("scipy.stats", "scipy"), ("statsmodels.stats.multitest",
                                                "statsmodels"), ("numpy", "numpy")):
        try:
            import importlib
            importlib.import_module(mod)
        except ImportError:
            faltando.append(pkg)
    assert not faltando, (
        f"pacotes ausentes: {', '.join(faltando)} — "
        "pip install numpy scipy statsmodels")


def test_finemap_parametriza_populacao():
    txt = ler("pipeline/scripts/finemap_ld.py")
    assert "LD_POP" in txt, "finemap_ld.py perdeu a parametrização de população"
    assert "POP_TAG" in txt, "cache de LD não é chaveado por população"


# ---------------------------------------------------------------- privacidade
def test_arquivos_privados_nao_rastreados():
    """Nenhum arquivo de caso real pode voltar ao índice do git."""
    import subprocess
    r = subprocess.run(["git", "ls-files"], cwd=BASE, capture_output=True, text=True)
    rastreados = set(r.stdout.split("\n"))
    for proibido in ("pasted.txt", "pasted.png", "MEMORIA.md",
                     "ARQUIVO_COMPLETO.md", "research/caso_real_contexto.md",
                     "research/ecossistema_ciencia_aberta_mapa.md",
                     "research/verificacao_ensaios_ion717_prism_2026-08-26.md"):
        assert proibido not in rastreados, \
            f"arquivo privado rastreado de novo: {proibido}"