/* Lithos Foundation — app.js (i18n PT-BR/EN + interações, 3 páginas) */
(function () {
  "use strict";

  var PAGE = document.documentElement.getAttribute("data-page") || "index";

  var I18N = {
    pt: {
      _title_index: "Lithos Foundation — Laboratório aberto de doenças priônicas",
      _title_processo: "Processo — Lithos Foundation",
      _title_dados: "Dados & Relatórios — Lithos Foundation",

      nav_inicio: "Início", nav_processo: "Processo", nav_dados: "Dados",
      nav_roadmap: "A seguir", nav_transp: "Transparência",

      hero_status: "PRiSM (NCT07444580) recrutando — única via anti-PRNP aberta a novos participantes",
      hero_eyebrow: "LABORATÓRIO ABERTO · DOENÇAS PRIÔNICAS · BRASIL",
      hero_h1: "Pedra sobre <em>pedra.</em>",
      hero_lead: "A Lithos Foundation é um laboratório de ciência aberta dedicado às doenças priônicas — incluída a Doença de Creutzfeldt-Jakob. Nascemos da urgência de um caso real e público no Brasil e trabalhamos do único jeito que a janela clínica permite: verificar dados que já existem, replicar o que foi publicado e transformar evidência em caminho navegável para famílias e médicos.",
      hero_cta1: "Ver os dados & relatórios", hero_cta2: "Como trabalhamos",
      stat1: "datasets reais analisados", stat2: "loci GWAS replicados", stat3: "r das magnitudes vs. artigo original", stat4: "miRNAs: sem ajuste → ajustado",
      hero_logo_nota: "marca provisória — a pedra lapidada é o registro que fica; a faceta laranja, a diferença que procuramos",

      op_eyebrow: "O QUE FAZEMOS · SEMPRE",
      op_h2: "Operação contínua",
      op_sub: "O que roda toda semana, sem exceção — porque o cenário muda em dias e a informação atrasada não ajuda ninguém. O método completo está na página <a href=\"processo.html\">Processo</a>.",
      op1_t: "VIGILÂNCIA DE ENSAIOS",
      op1_p: "Semanal: ClinicalTrials.gov + comunicados da CJD Foundation, CureFFI e Prion Alliance — historicamente mais rápidos que o registro oficial.",
      op2_t: "INTELIGÊNCIA CIENTÍFICA",
      op2_p: "Varredura PubMed de terapia priônica dos últimos 30 dias, a cada sessão de trabalho.",
      op3_t: "DOSSIÊS PARA MÉDICOS",
      op3_p: "Material informativo pronto (imagem, LCR, RT-QuIC, MRC-PDRS) para a decisão dos assistentes — nunca conselho clínico direto.",
      op4_t: "CURADORIA DE DADOS",
      op4_p: "GEO · ENA · PRIDE · GWAS Catalog: o que existe, o que falta e o que não existe — com checksums MD5. Veja a página <a href=\"dados.html\">Dados</a>.",
      op5_t: "REDE BRASIL",
      op5_p: "RT-QuIC (HC-FMUSP · A.C. Camargo · UFRJ/LAB3DDI), centros de referência por região, pesquisadores-chave e canal oficial do MS.",
      op6_t: "GUIA DE FAMÍLIAS",
      op6_p: "PT-BR claro: o que é a doença, o que é um ensaio, o que é promessa — sem jargão e sem esperança vendida.",

      sit_eyebrow: "SITUAÇÃO ATUAL · VERIFICADO 29·08·2026",
      sit_h2: "Ensaios e vias, agora",
      sit_sub: "O cenário muda em dias — este quadro reflete a última verificação viva dos registros oficiais e comunicados das organizações de pacientes.",
      sit1_b: "RECRUTANDO · VIA PRIORITÁRIA",
      sit1_p: "Fase 1, <span class=\"destaque\">dose única intratecal</span> com distribuição pan-cerebral. Em camundongos: −49% de PrP e <span class=\"destaque\">+64% de sobrevida mesmo pós-sintomas</span>. Cinco sites nos EUA; custos de viagem apoiados por doações à Prion Alliance.",
      sit1_n: "fase 1 avalia segurança — participação é pesquisa, não tratamento garantido",
      sit2_b: "FECHADO A NOVOS PARTICIPANTES",
      sit2_p: "“Active, not recruiting” desde 26·08 (regimes 1–2 completos; regime 3 encerrou a inclusão). <span class=\"destaque\">Leitura primária estimada para fev/2027</span> — vigilância semanal para eventual reabertura.",
      sit3_b: "AINDA NÃO RECRUTA",
      sit3_p: "Medicamento reposicionado (anti-HIV), via oral. Xuanwu Hospital / Capital Medical University, Pequim. <span class=\"destaque\">Sem via prática</span> para paciente brasileiro no momento.",
      sit4_t: "Acesso expandido / compassivo",
      sit4_b: "NENHUM PROGRAMA PÚBLICO",
      sit4_p: "Nenhum sponsor anunciou programa público; a droga do PRiSM <span class=\"destaque\">não sai dos EUA</span>. A via realista continua sendo a inclusão direta no ensaio, pelo neurologista assistente.",
      gargalo: "<b>O gargalo decisivo:</b> RT-QuIC+ (ou mutação PRNP documentada) é o passaporte diagnóstico; MRC-PDRS ≥ 15 é o filtro de elegibilidade. Dossiê pronto em <b>dias — não semanas</b> — é o que encurta o caminho.",

      road_eyebrow: "OBJETIVOS · PRÓXIMO ESTRATO",
      road_h2: "O que vem a seguir",
      road_sub: "A ordem importa: publicar antes de falar, constituir antes de captar, presença antes de promessa.",
      r1_t: "Preprint & DOI",
      r1_p: "Manuscrito em inglês no bioRxiv, arquivo completo no Zenodo (DOI citável) e repositório público com o código de análise.",
      r2_t: "Constituição formal",
      r2_p: "Associação sem fins lucrativos: estatuto, CNPJ, CEBAS e registro de marca no INPI — a fundação sai do papel.",
      r3_t: "Presença própria",
      r3_p: "Reserva dos domínios verificados livres em 28·08 (lithosfoundation.org e afins) e endereço de contato institucional.",
      r4_t: "Plataforma de famílias",
      r4_p: "Navegação em PT-BR: do primeiro sintoma ao ensaio, com checklists de exames, perguntas para o médico e fontes verificadas.",
      r5_t: "Advocacy",
      r5_p: "Ministério da Saúde, ANVISA e Congresso: diagnóstico RT-QuIC acessível no SUS e ensaios clínicos de doenças priônicas no Brasil.",

      tr_eyebrow: "TRANSPARÊNCIA",
      tr_h2: "Regras de casa",
      tr1_t: "SÓ DADOS PÚBLICOS",
      tr1_p: "Nenhum dado de paciente individual. Saúde é dado sensível (LGPD): a trilha real usa somente datasets públicos e anonimizados que já existiam.",
      tr2_t: "EVIDÊNCIA COM NÍVEL",
      tr2_p: "Pré-clínico ≠ humano. Réplica ≠ descoberta. Fase 1 = segurança, não promessa. Cada afirmação carrega o seu nível de evidência.",
      tr3_t: "RASTREÁVEL",
      tr3_p: "Toda afirmação aponta para fonte primária com URL e data de acesso. Dados com checksum MD5; código versionado no repositório.",
      aviso: "ESTA PÁGINA É INFORMATIVA E NÃO SUBSTITUI DECISÃO MÉDICA. Nada aqui é conselho individual; todo material clínico é entregue aos médicos assistentes, que decidem.",

      proc_eyebrow: "PROCESSO · METODOLOGIA ABERTA",
      proc_h1: "Como cada número nasce, é posto à prova e vira registro",
      proc_lead: "Nada aqui começa numa intuição: começa num repositório público. O que garante a confiança não é o resultado — é o caminho que o resultado percorreu até aqui. Este é esse caminho, de ponta a ponta.",

      fluxo_eyebrow: "O FLUXO",
      fluxo_h2: "Seis estágios, uma regra: nenhum atalho",
      fluxo_aria: "Fluxo em seis estágios: fontes públicas, aquisição com checksum, análise estatística, réplica contra o artigo original, auditoria adversarial, preprint e engajamento; vigilância clínica semanal alimenta os dossiês.",
      fl1: "FONTES PÚBLICAS", fl2: "AQUISIÇÃO + MD5", fl3: "ANÁLISE ESTATÍSTICA",
      fl4: "RÉPLICA × ARTIGO", fl5: "AUDITORIA ADVERSARIAL", fl6: "PREPRINT & ENGAJAMENTO",
      fl7: "VIGILÂNCIA CLÍNICA SEMANAL → DOSSIÊS",
      fluxo_cap: "O estágio vermelho é deliberado: todo resultado é atacado por nós mesmos antes de ver alguém de fora.",

      fontes_eyebrow: "ESTÁGIOS 1–2 · FONTES & AQUISIÇÃO",
      fontes_h2: "De onde vem cada byte",
      f1_t: "REPOSITÓRIOS DE DADOS",
      f1_p: "GEO e ENA (transcriptomas, WGS), PRIDE (proteômica), GWAS Catalog (sumários estatísticos). Download direto da fonte primária — nunca de espelho ou repositório de terceiros.",
      f2_t: "REGISTROS DE ENSAIOS",
      f2_p: "ClinicalTrials.gov verificado a cada sessão, cruzado com os comunicados da CJD Foundation, CureFFI e Prion Alliance — que historicamente antecipam o registro oficial.",
      f3_t: "LITERATURA & SUPLEMENTOS",
      f3_p: "PubMed (janela de 30 dias) e os suplementos oficiais dos autores: quando o GEO vem sem tabela (caso GSE140069), a análise usa o arquivo oficial distribuído com a publicação.",
      aq_t: "Regras de aquisição",
      aq1: "Somente dados públicos e anonimizados — saúde é dado sensível (LGPD); nenhum dado individual entra no laboratório.",
      aq2: "Checksum MD5 por arquivo baixado; dados grandes ficam fora do git, mas são re-baixáveis com o checksum documentado (ex.: GWAS de 197 MB).",
      aq3: "Catálogo inclui os vazios: registrar o que NÃO existe (ex.: ausência de certos cohortes abertos) evita que outros procurem o que não há.",
      aq4: "Data de acesso anotada em toda fonte — datasets mudam; o registro fixa o estado usado pela análise.",
      an_t: "Regras de análise",
      an1: "Covariáveis explícitas em todo modelo: sexo, códon 129 e subtipo no cérebro; idade, sexo e RIN no sangue (OLS).",
      an2: "Testes clássicos auditáveis — Welch com correção BH-FDR; QC de GWAS com λGC reportado (1.059 em 6,31 milhões de variantes).",
      an3: "Anotação funcional via fontes oficiais (Ensembl GRCh37 REST + NCBI), nunca por memória de modelo.",
      an4: "Resultado negativo é resultado: a integração cérebro×sangue (miRTarBase 10.0 × DEGs) foi negativa e está documentada.",

      reg_eyebrow: "O REGISTRO · 24—28·08·2026",
      reg_h2: "O processo em ação, camada por camada",
      reg_sub: "Cada camada abaixo é datada, rastreada à fonte primária e verificável no repositório. Nada aqui se apaga: o registro é a nossa matéria-prima.",
      reg_como_t: "COMO LER",
      reg_como_p: "O registro é organizado em formações geológicas: o que foi colocado primeiro sustenta o que veio depois. Todas as afirmações apontam para fonte primária, com URL e data de acesso.",
      reg_leg1: "dados & verificação", reg_leg2: "clínica & urgência",

      e1h: "<b>FORMAÇÃO 1 · DADOS REAIS</b><time>24—26·08·2026</time>",
      r11: "<time>24·08</time><b>Fundação do laboratório.</b> Duas trilhas: “caso referência” (simulação didática sCJD MM1, dossiê validado pelos critérios CDC/NPDPSC) e trilha real — somente dados públicos e anonimizados. Pacote completo distribuído: 36 documentos, 7 scripts, checksums MD5.",
      r12: "<time>24·08</time><b>GSE160208</b> (47 amostras sCJD/controles, Copenhagen): análise completa — 441 de 742 genes testáveis significantes (Welch + BH-FDR; 437 de 800 na contagem anterior, que incluía 58 genes sem variância); neuroinflamação massiva (MYD88, TLR2, C1QA, CSF1 — FDR < 10⁻⁴).",
      r13: "<time>24·08</time><b>Catálogos verificados:</b> GEO, ENA (WGS PRJEB57720), GWAS Catalog (GCST90001389), proteômica PRIDE (PXD050656) — incluindo o que <b>não</b> existe.",
      e2h: "<b>FORMAÇÃO 2 · VERIFICAÇÃO INDEPENDENTE</b><time>24·08·2026</time>",
      r21: "<time>24·08</time><b>Contagem do artigo do GSE160208 reproduzida:</b> 184 DEGs com o filtro de magnitude antes do FDR (183 com o FDR no painel inteiro); top-10 de genes na mesma ordem, r = 1.000. Não é réplica exata do procedimento dos autores.",
      r22: "<time>24·08</time><b>GSE140069</b> (miRNA no sangue, 57 sCJD × 48 controles): 89 miRNAs significantes sem ajuste e 1 após ajuste por idade, sexo e RIN; os 4 achados do artigo mantêm direção e p nominal < 0,05 (let-7i-5p, miR-16-5p, miR-93-5p, miR-106b-3p); AUC 0.788 (sCJD×CT) e 0.924 (sCJD×AD).",
      r23: "<time>24·08</time><b>GWAS GCST90001389</b> (197 MB): QC de 6.314.492 variantes (0 malformadas, λGC = 1.059) e <b>réplica independente 3/3 dos loci publicados</b> — PRNP p=1,6e-15 · STX6 p=7,5e-9 · GAL3ST1 p=6,2e-10 — com fine-mapping de rs3747957 (mesma direção do Brain 2025).",
      r24: "<time>24·08</time><b>Auditoria cética adversarial:</b> máquina estatística validada (Welch/BH corretos a 1e-13); fragilidades de assinaturas publicadas, documentadas com honestidade.",
      e3h: "<b>FORMAÇÃO 3 · HIPÓTESES & SIMULAÇÃO</b><time>24·08·2026</time>",
      r31: "<time>24·08</time><b>Hipótese geradora “alfândega intercelular”</b> (checkpoint biofísico seletivo nas pontes entre células) submetida ao Institut Pasteur — grupo Chiara Zurzolo — com 3 previsões testáveis.",
      r32: "<time>24·08</time><b>Blindagem G127V:</b> varredura encontrou limiar de percolação em ~41% de cobertura para colapsar o espalhamento priônico.",
      r33: "<time>24·08</time><b>Simulação calibrada por literatura</b> (Weibull estocástico): 3 validações — incluindo sobrevida MM1 de 4,3 meses contra 4–5 publicado.",
      e4h: "<b>FORMAÇÃO 4 · CLÍNICA & REDE</b><time>28·08·2026</time>",
      r41: "<time>28·08</time><b>Cenário de ensaios mapeado e repriorizado:</b> PRiSM (NCT07444580) recrutando → via prioritária; ION717 fechado a novos desde 26·08; efavirenz (NCT07482085) ainda sem recrutamento.",
      r42: "<time>28·08</time><b>Logística internacional:</b> visto B-2 médico, cronograma de deslocamento do PRiSM, contatos dos 5 sites e do sponsor (Broad Institute).",
      r43: "<time>28·08</time><b>Dossiê médico modelo:</b> RM (DWI/FLAIR), LCR (14-3-3/t-tau), RT-QuIC e escala MRC-PDRS — pronto para a decisão dos médicos assistentes.",
      r44: "<time>28·08</time><b>Rede Brasil consolidada:</b> RT-QuIC (HC-FMUSP · A.C. Camargo · UFRJ/LAB3DDI), centros de referência por região e canal oficial do Ministério da Saúde.",
      e5h: "<b>FORMAÇÃO 5 · INSTITUCIONAL</b><time>28·08·2026</time>",
      r51: "<time>28·08</time><b>Engajamento científico:</b> cartas formais ao Prion Alliance (Broad Institute) e ao HC-FMUSP (Neurologia Cognitiva), com o pacote completo do laboratório.",
      r52: "<time>28·08</time><b>Guia de famílias em PT-BR</b> publicado — sem jargão, sem promessa.",
      r53: "<time>28·08</time><b>Batismo: Lithos Foundation</b> — homenagem mineral à origem da missão; terreno verificado (sem colisões de marca, domínios livres).",

      gar_eyebrow: "GARANTIAS DO PROCESSO",
      gar_h2: "Por que acreditar neste registro",
      g1_t: "MD5 EM TUDO",
      g1_p: "Cada arquivo baixado carrega checksum documentado; o pacote completo do laboratório (36 documentos + 7 scripts) foi distribuído com verificação de integridade.",
      g2_t: "CÓDIGO VERSIONADO",
      g2_p: "Os 7 scripts de análise vivem no repositório com histórico de commits — cada figura e cada número são regeneráveis a partir do dado bruto.",
      g3_t: "AUDITORIA ANTES DE PUBLICAR",
      g3_p: "Crítica estatística adversarial interna obrigatória: erros de desenho identificados foram corrigidos antes da distribuição — e documentados quando não foram evitados.",
      g4_t: "NADA DE PACIENTE INDIVIDUAL",
      g4_p: "O caso real que motiva a missão nunca entra no banco de dados: a trilha clínica orienta famílias e médicos com informação pública, sem tocar em prontuário.",

      dados_eyebrow: "DADOS · RELATÓRIOS · VERIFICAÇÃO",
      dados_h1: "Os dados, abertos de ponta a ponta",
      dados_lead: "Tudo nesta página vem de repositórios públicos que qualquer pessoa pode baixar hoje. Cada gráfico é uma leitura-sumário de relatórios completos, com código versionado e checksum MD5 — e cada réplica é um serviço de verificação ao campo, não uma descoberta nossa.",

      cat_eyebrow: "CATÁLOGO",
      cat_h2: "Os datasets do laboratório",
      cat_sub: "Analisados a fundo ou catalogados como próximos na fila — incluindo o que decidimos não perseguir, e por quê.",
      th1: "ACESSO", th2: "CONTEÚDO", th3: "ESCALA", th4: "SITUAÇÃO",
      d1: "Transcriptoma cerebral (sCJD × controles), Copenhagen",
      d2: "miRNA no sangue (sCJD × CT × Alzheimer)",
      d3: "GWAS de sCJD — sumários estatísticos completos",
      d4: "WGS (sequenciamento completo) — ENA",
      d5: "Proteômica de príons — PRIDE",
      st_analisado: "ANALISADO", st_cat: "CATALOGADO",

      rel_eyebrow: "RELATÓRIOS VISUAIS",
      rel_h2: "Gráficos, com os números por trás",
      p1_t: "GWAS · réplica independente 3/3",
      p1_aria: "Três cromossomos com os loci replicados nas posições reais: PRNP no cromossomo 20 a 4,67 megabases (p igual a 1,6e-15); STX6 no cromossomo 1 a 180,96 megabases (p igual a 7,5e-9); GAL3ST1 no cromossomo 22 a 30,95 megabases (p igual a 6,2e-10).",
      p1_badge: "réplica 3/3 dos loci",
      p1_gws: "significância genômica (−log10 p = 7.3)",
      p1_note: "posições dos marcadores em coordenadas reais do genoma (GRCh37); alturas proporcionais a −log10(p) · anotação Ensembl + NCBI",
      p2_t: "Transcriptoma cerebral",
      p2_aria: "Barra: 441 genes significantes de 742 testáveis, 59,4 por cento; contagem de 184 genes do artigo original reproduzida, com r igual a 1.000",
      p2_linha1: "neuroinflamação: os 4 hubs acima do limiar em toda a coorte",
      p2_linha2: "47 amostras · Copenhagen · MM1 vs. VV2 estratificado",
      p2_note: "duas contagens, com o universo declarado: 441 DEGs pelo nosso critério (BH-FDR sobre 742 genes testáveis; 437 se contados os 800 do painel) e 184 pelo critério do artigo original, reproduzido com o filtro de magnitude antes do FDR",
      p3_t: "miRNA no sangue",
      p3_aria: "Dois medidores semicirculares: AUC 0,788 para sCJD contra controles; AUC 0,924 para sCJD contra Alzheimer; quatro de quatro achados do artigo replicados",
      p3_linha: "89 miRNAs sem ajuste → 1 com ajuste (idade+sexo+RIN)",
      p3_note: "honestidade da casa: com covariáveis (idade+sexo+RIN), 1 miRNA sobrevive ao FDR — a assinatura final do artigo é de 3 miRNAs por qPCR, e a fragilidade está documentada no nosso relatório",
      p5_t: "Integração cérebro × sangue",
      p5_badge: "RESULTADO NEGATIVO · DOCUMENTADO",
      p5_p: "Cruzamos os miRNAs do sangue com os genes mais expressos no cérebro (GSE160208), usando alvos funcionais validados: <b>nenhum enriquecimento detectado.</b> O teste tem pouco poder — só 3 a 13 alvos de cada miRNA estão no painel — e o resultado nulo não apoia nem contraria uma relação entre sangue e cérebro.",
      p5_note: "negativos se publicam: este resultado integra o §4.5 do preprint, com o poder do teste declarado",
      p4_t: "Simulação calibrada",
      p4_aria: "Barras de dano esperado: sem intervenção 100 por cento; alfândega intercelular 50,2 por cento; muro total 16,3 por cento; limiar de blindagem em 41 por cento",
      p4_svg1: "sem intervenção · 100%",
      p4_svg2: "alfândega 80/5 · 50,2%",
      p4_svg5: "muro total · 16,3%",
      p4_svg3: "↑ limiar de blindagem G127V: ~41% de cobertura colapsa o espalhamento (percolação p_c = 0,593)",
      p4_chip: "hipótese “alfândega intercelular” → Institut Pasteur",
      p4_chip2: "3 previsões testáveis",
      p4_note: "modelo gerador de hipóteses — não é previsão clínica; a síntese não existia na literatura (verificado por busca amplo antes do envio)",
      dados_nota: "FIGURAS COMPLETAS NO REPOSITÓRIO: vulcões (cérebro e sangue), heatmap dos top-25 genes e linha do tempo do caso de referência — com código e checksums MD5. Esta página mostra leituras-sumário; o detalhe estatístico integral vive nos relatórios.",

      foot_nota: "Nome provisório. Nasce da urgência de um caso real e público no Brasil, anunciado pela família em 21·08·2026 — e de uma pergunta simples: por que a informação certa chega tarde demais?",
      foot_meta: "registro vivo · atualizado em 29·08·2026<br>contato: em breve, junto da constituição formal<br>feito em PT-BR · tudo verificável"
    },

    en: {
      _title_index: "Lithos Foundation — Open laboratory for prion diseases",
      _title_processo: "Process — Lithos Foundation",
      _title_dados: "Data & Reports — Lithos Foundation",

      nav_inicio: "Home", nav_processo: "Process", nav_dados: "Data",
      nav_roadmap: "Next", nav_transp: "Transparency",

      hero_status: "PRiSM (NCT07444580) recruiting — the only anti-PRNP route open to new participants",
      hero_eyebrow: "OPEN LABORATORY · PRION DISEASES · BRAZIL",
      hero_h1: "Stone upon <em>stone.</em>",
      hero_lead: "Lithos Foundation is an open-science laboratory dedicated to prion diseases — including Creutzfeldt-Jakob Disease. We were born from the urgency of a real, public case in Brazil, and we work the only way the clinical window allows: verifying data that already exists, replicating what has been published, and turning evidence into a navigable path for families and physicians.",
      hero_cta1: "See the data & reports", hero_cta2: "How we work",
      stat1: "real datasets analyzed", stat2: "GWAS loci replicated", stat3: "r of effect sizes vs. the original paper", stat4: "miRNAs: unadjusted → adjusted",
      hero_logo_nota: "provisional mark — the cut stone is the record that remains; the orange facet is the difference we look for",

      op_eyebrow: "WE DO · ALWAYS",
      op_h2: "Continuous operation",
      op_sub: "What runs every single week — because the landscape changes in days, and late information helps no one. The full method lives on the <a href=\"processo.html\">Process</a> page.",
      op1_t: "TRIAL WATCH",
      op1_p: "Weekly: ClinicalTrials.gov + statements from the CJD Foundation, CureFFI and Prion Alliance — historically faster than the official registry.",
      op2_t: "SCIENTIFIC INTELLIGENCE",
      op2_p: "PubMed sweep of prion therapy from the last 30 days, every working session.",
      op3_t: "PHYSICIAN DOSSIERS",
      op3_p: "Ready-to-use informational material (imaging, CSF, RT-QuIC, MRC-PDRS) for the attending physicians' decision — never direct clinical advice.",
      op4_t: "DATA CURATION",
      op4_p: "GEO · ENA · PRIDE · GWAS Catalog: what exists, what is missing and what does not exist — with MD5 checksums. See the <a href=\"dados.html\">Data</a> page.",
      op5_t: "BRAZIL NETWORK",
      op5_p: "RT-QuIC (HC-FMUSP · A.C. Camargo · UFRJ/LAB3DDI), referral centers by region, key researchers and the Ministry of Health's official channel.",
      op6_t: "FAMILY GUIDE",
      op6_p: "Plain PT-BR: what the disease is, what a trial is, what a promise is — no jargon, no sold hope.",

      sit_eyebrow: "CURRENT STATUS · VERIFIED AUG 29, 2026",
      sit_h2: "Trials and routes, right now",
      sit_sub: "The landscape changes in days — this board reflects the latest live check of official registries and patient-organization statements.",
      sit1_b: "RECRUITING · PRIORITY ROUTE",
      sit1_p: "Phase 1, <span class=\"destaque\">single intrathecal dose</span> with pan-brain distribution. In mice: −49% PrP and <span class=\"destaque\">+64% survival even post-symptom onset</span>. Five US sites; travel costs supported by donations to Prion Alliance.",
      sit1_n: "phase 1 measures safety — participation is research, not guaranteed treatment",
      sit2_b: "CLOSED TO NEW PARTICIPANTS",
      sit2_p: "\"Active, not recruiting\" since Aug 26 (regimens 1–2 complete; regimen 3 finished enrollment). <span class=\"destaque\">Primary readout estimated for Feb 2027</span> — weekly watch for a possible reopening.",
      sit3_b: "NOT YET RECRUITING",
      sit3_p: "Repurposed drug (anti-HIV), oral. Xuanwu Hospital / Capital Medical University, Beijing. <span class=\"destaque\">No practical route</span> for a Brazilian patient right now.",
      sit4_t: "Expanded / compassionate access",
      sit4_b: "NO PUBLIC PROGRAM",
      sit4_p: "No sponsor has announced a public program; the PRiSM drug <span class=\"destaque\">does not leave the US</span>. The realistic route remains direct enrollment, led by the attending neurologist.",
      gargalo: "<b>The decisive bottleneck:</b> RT-QuIC+ (or a documented PRNP mutation) is the diagnostic passport; MRC-PDRS ≥ 15 is the eligibility filter. A dossier ready in <b>days — not weeks</b> — is what shortens the path.",

      road_eyebrow: "GOALS · THE UPCOMING STRATUM",
      road_h2: "What comes next",
      road_sub: "Order matters: publish before speaking, incorporate before fundraising, presence before promises.",
      r1_t: "Preprint & DOI",
      r1_p: "English manuscript on bioRxiv, full archive on Zenodo (citable DOI) and a public repository with the analysis code.",
      r2_t: "Formal incorporation",
      r2_p: "Nonprofit association: bylaws, CNPJ, CEBAS and trademark registration with INPI — the foundation leaves the drawing board.",
      r3_t: "Own presence",
      r3_p: "Reserve of the domains verified as available on Aug 28 (lithosfoundation.org and kin) plus an institutional contact address.",
      r4_t: "Family platform",
      r4_p: "PT-BR navigation: from first symptom to trial, with exam checklists, questions for the physician and verified sources.",
      r5_t: "Advocacy",
      r5_p: "Ministry of Health, ANVISA and Congress: accessible RT-QuIC diagnosis in the public system and prion-disease clinical trials in Brazil.",

      tr_eyebrow: "TRANSPARENCY",
      tr_h2: "House rules",
      tr1_t: "PUBLIC DATA ONLY",
      tr1_p: "No individual patient data. Health is sensitive data (LGPD/Brazilian GDPR): the real track uses only public, anonymized datasets that already existed.",
      tr2_t: "EVIDENCE WITH A LEVEL",
      tr2_p: "Preclinical ≠ human. Replication ≠ discovery. Phase 1 = safety, not promise. Every claim carries its own level of evidence.",
      tr3_t: "TRACEABLE",
      tr3_p: "Every claim points to a primary source with URL and access date. Data with MD5 checksums; code versioned in the repository.",
      aviso: "THIS PAGE IS INFORMATIONAL AND DOES NOT REPLACE MEDICAL DECISION-MAKING. Nothing here is individual advice; all clinical material is delivered to the attending physicians, who decide.",

      proc_eyebrow: "PROCESS · OPEN METHODOLOGY",
      proc_h1: "How every number is born, put to the test, and becomes record",
      proc_lead: "Nothing here starts with a hunch: it starts in a public repository. What earns trust is not the result — it is the path the result traveled to get here. This is that path, end to end.",

      fluxo_eyebrow: "THE FLOW",
      fluxo_h2: "Six stages, one rule: no shortcuts",
      fluxo_aria: "Six-stage flow: public sources, acquisition with checksum, statistical analysis, replication against the original paper, adversarial audit, preprint and engagement; weekly clinical surveillance feeds the dossiers.",
      fl1: "PUBLIC SOURCES", fl2: "ACQUISITION + MD5", fl3: "STATISTICAL ANALYSIS",
      fl4: "REPLICATION × PAPER", fl5: "ADVERSARIAL AUDIT", fl6: "PREPRINT & ENGAGEMENT",
      fl7: "WEEKLY CLINICAL SURVEILLANCE → DOSSIERS",
      fluxo_cap: "The red stage is deliberate: every result is attacked by us before it ever sees an outsider.",

      fontes_eyebrow: "STAGES 1–2 · SOURCES & ACQUISITION",
      fontes_h2: "Where every byte comes from",
      f1_t: "DATA REPOSITORIES",
      f1_p: "GEO and ENA (transcriptomes, WGS), PRIDE (proteomics), GWAS Catalog (summary statistics). Downloaded straight from the primary source — never from a mirror or a third party.",
      f2_t: "TRIAL REGISTRIES",
      f2_p: "ClinicalTrials.gov checked every session, cross-read with statements from the CJD Foundation, CureFFI and Prion Alliance — historically ahead of the official registry.",
      f3_t: "LITERATURE & SUPPLEMENTS",
      f3_p: "PubMed (30-day window) and the authors' official supplements: when GEO ships without a data table (the GSE140069 case), the analysis uses the official file distributed with the paper.",
      aq_t: "Acquisition rules",
      aq1: "Public, anonymized data only — health is sensitive data (LGPD); no individual data ever enters the laboratory.",
      aq2: "MD5 checksum per downloaded file; large datasets stay out of git but are re-downloadable with their documented checksum (e.g., the 197 MB GWAS).",
      aq3: "The catalog includes the empties: recording what does NOT exist (e.g., missing open cohorts) stops others from hunting for what isn't there.",
      aq4: "Access date noted on every source — datasets change; the registry freezes the state the analysis used.",
      an_t: "Analysis rules",
      an1: "Covariates explicit in every model: sex, codon 129 and subtype in brain; age, sex and RIN in blood (OLS).",
      an2: "Classic auditable tests — Welch with BH-FDR correction; GWAS QC with reported λGC (1.059 across 6.31 million variants).",
      an3: "Functional annotation via official sources (Ensembl GRCh37 REST + NCBI), never from model memory.",
      an4: "A negative result is a result: the brain×blood integration (miRTarBase 10.0 × DEGs) came out negative and is documented.",

      reg_eyebrow: "THE REGISTRY · AUG 24—28, 2026",
      reg_h2: "The process in action, layer by layer",
      reg_sub: "Every layer below is dated, traced to its primary source and verifiable in the repository. Nothing here is erased: the registry is our raw material.",
      reg_como_t: "HOW TO READ",
      reg_como_p: "The registry is organized as geological formations: what was laid down first supports what came after. Every claim points to a primary source, with URL and access date.",
      reg_leg1: "data & verification", reg_leg2: "clinical & urgency",

      e1h: "<b>FORMATION 1 · REAL DATA</b><time>Aug 24—26, 2026</time>",
      r11: "<time>Aug 24</time><b>Laboratory founded.</b> Two tracks: a “reference case” (didactic sCJD MM1 simulation, dossier validated against CDC/NPDPSC criteria) and the real track — public, anonymized data only. Full package distributed: 36 documents, 7 scripts, MD5 checksums.",
      r12: "<time>Aug 24</time><b>GSE160208</b> (47 sCJD/control samples, Copenhagen): full analysis — 441 of 742 testable genes significant (Welch + BH-FDR; 437 of 800 in the earlier count, which included 58 genes with no variance); massive neuroinflammation (MYD88, TLR2, C1QA, CSF1 — FDR < 10⁻⁴).",
      r13: "<time>Aug 24</time><b>Catalogs verified:</b> GEO, ENA (WGS PRJEB57720), GWAS Catalog (GCST90001389), PRIDE proteomics (PXD050656) — including what does <b>not</b> exist.",
      e2h: "<b>FORMATION 2 · INDEPENDENT VERIFICATION</b><time>Aug 24, 2026</time>",
      r21: "<time>Aug 24</time><b>Count of the GSE160208 paper reproduced:</b> 184 DEGs with the magnitude filter applied before the FDR (183 with the FDR over the whole panel); top-10 gene list in the same order, r = 1.000. Not an exact replication of the authors' procedure.",
      r22: "<time>Aug 24</time><b>GSE140069</b> (blood miRNA, 57 sCJD × 48 controls): 89 significant miRNAs unadjusted and 1 after adjustment for age, sex and RIN; the paper's 4 findings keep their direction and nominal p < 0.05 (let-7i-5p, miR-16-5p, miR-93-5p, miR-106b-3p); AUC 0.788 (sCJD×CT) and 0.924 (sCJD×AD).",
      r23: "<time>Aug 24</time><b>GWAS GCST90001389</b> (197 MB): QC of 6,314,492 variants (0 malformed, λGC = 1.059) and an <b>independent 3/3 replication of the published loci</b> — PRNP p=1.6e-15 · STX6 p=7.5e-9 · GAL3ST1 p=6.2e-10 — with fine-mapping of rs3747957 (same direction as Brain 2025).",
      r24: "<time>Aug 24</time><b>Adversarial skeptical audit:</b> statistical machinery validated (Welch/BH correct to 1e-13); fragilities of published signatures documented with honesty.",
      e3h: "<b>FORMATION 3 · HYPOTHESES & SIMULATION</b><time>Aug 24, 2026</time>",
      r31: "<time>Aug 24</time><b>Generative hypothesis \"intercellular customs checkpoint\"</b> (selective biophysical checkpoint at intercellular bridges) submitted to Institut Pasteur — Chiara Zurzolo's group — with 3 testable predictions.",
      r32: "<time>Aug 24</time><b>G127V shielding:</b> a sweep found a percolation threshold at ~41% coverage that collapses prion spread.",
      r33: "<time>Aug 24</time><b>Literature-calibrated simulation</b> (stochastic Weibull): 3 validations — including MM1 survival of 4.3 months against 4–5 published.",
      e4h: "<b>FORMATION 4 · CLINICAL & NETWORK</b><time>Aug 28, 2026</time>",
      r41: "<time>Aug 28</time><b>Trial landscape mapped and reprioritized:</b> PRiSM (NCT07444580) recruiting → priority route; ION717 closed to new participants since Aug 26; efavirenz (NCT07482085) not yet recruiting.",
      r42: "<time>Aug 28</time><b>International logistics:</b> B-2 medical visa, PRiSM travel schedule, contacts for all 5 sites and the sponsor (Broad Institute).",
      r43: "<time>Aug 28</time><b>Model clinical dossier:</b> MRI (DWI/FLAIR), CSF (14-3-3/t-tau), RT-QuIC and the MRC-PDRS scale — ready for the attending physicians' decision.",
      r44: "<time>Aug 28</time><b>Brazil network consolidated:</b> RT-QuIC (HC-FMUSP · A.C. Camargo · UFRJ/LAB3DDI), referral centers by region and the Ministry of Health's official channel.",
      e5h: "<b>FORMATION 5 · INSTITUTIONAL</b><time>Aug 28, 2026</time>",
      r51: "<time>Aug 28</time><b>Scientific engagement:</b> formal letters to Prion Alliance (Broad Institute) and HC-FMUSP (Cognitive Neurology), with the laboratory's complete package.",
      r52: "<time>Aug 28</time><b>Family guide in PT-BR</b> published — no jargon, no promises.",
      r53: "<time>Aug 28</time><b>Naming: Lithos Foundation</b> — a mineral tribute to the mission's origin; terrain verified (no brand collisions, domains available).",

      gar_eyebrow: "PROCESS GUARANTEES",
      gar_h2: "Why trust this registry",
      g1_t: "MD5 ON EVERYTHING",
      g1_p: "Every downloaded file carries a documented checksum; the laboratory's complete package (36 documents + 7 scripts) was distributed with integrity verification.",
      g2_t: "VERSIONED CODE",
      g2_p: "The 7 analysis scripts live in the repository with commit history — every figure and every number is regenerable from the raw data.",
      g3_t: "AUDIT BEFORE PUBLISHING",
      g3_p: "Mandatory internal adversarial statistical review: identified design flaws were fixed before distribution — and documented when they weren't avoided.",
      g4_t: "NO INDIVIDUAL PATIENT DATA",
      g4_p: "The real case that motivates the mission never enters the database: the clinical track guides families and physicians with public information, never touching a medical record.",

      dados_eyebrow: "DATA · REPORTS · VERIFICATION",
      dados_h1: "The data, open end to end",
      dados_lead: "Everything on this page comes from public repositories anyone can download today. Each chart is a summary read-out of full reports, with versioned code and MD5 checksums — and every replication is a verification service to the field, not our own discovery.",

      cat_eyebrow: "CATALOG",
      cat_h2: "The laboratory's datasets",
      cat_sub: "Fully analyzed or cataloged as next in line — including what we chose not to chase, and why.",
      th1: "ACCESSION", th2: "CONTENT", th3: "SCALE", th4: "STATUS",
      d1: "Brain transcriptome (sCJD × controls), Copenhagen",
      d2: "Blood miRNA (sCJD × CT × Alzheimer's)",
      d3: "sCJD GWAS — full summary statistics",
      d4: "WGS (whole-genome sequencing) — ENA",
      d5: "Prion proteomics — PRIDE",
      st_analisado: "ANALYZED", st_cat: "CATALOGED",

      rel_eyebrow: "VISUAL REPORTS",
      rel_h2: "Charts, with the numbers behind them",
      p1_t: "GWAS · independent 3/3 replication",
      p1_aria: "Three chromosomes with the replicated loci at their real positions: PRNP on chromosome 20 at 4.67 megabases (p = 1.6e-15); STX6 on chromosome 1 at 180.96 megabases (p = 7.5e-9); GAL3ST1 on chromosome 22 at 30.95 megabases (p = 6.2e-10).",
      p1_badge: "3/3 published loci replicated",
      p1_gws: "genome-wide significance (−log10 p = 7.3)",
      p1_note: "marker positions at real genome coordinates (GRCh37); heights proportional to −log10(p) · Ensembl + NCBI annotation",
      p2_t: "Brain transcriptome",
      p2_aria: "Bar: 441 significant genes out of 742 testable, 59.4 percent; the original paper's count of 184 genes reproduced, with r equal to 1.000",
      p2_linha1: "neuroinflammation: all 4 hubs above threshold across the cohort",
      p2_linha2: "47 samples · Copenhagen · MM1 vs. VV2 stratified",
      p2_note: "two counts, with the test universe stated: 441 DEGs by our criterion (BH-FDR over 742 testable genes; 437 if all 800 panel genes are counted) and 184 by the original paper's criterion, reproduced with the magnitude filter applied before the FDR",
      p3_t: "Blood miRNA",
      p3_aria: "Two semicircular gauges: AUC 0.788 for sCJD versus controls; AUC 0.924 for sCJD versus Alzheimer's; four of four paper findings replicated",
      p3_linha: "89 miRNAs unadjusted → 1 adjusted (age+sex+RIN)",
      p3_note: "house honesty: with covariates (age+sex+RIN), 1 miRNA survives FDR — the paper's final signature is 3 miRNAs by qPCR, and the fragility is documented in our report",
      p5_t: "Brain × blood integration",
      p5_badge: "NEGATIVE RESULT · DOCUMENTED",
      p5_p: "We crossed the blood miRNAs with the genes up-regulated in brain (GSE160208), using validated functional targets: <b>no enrichment detected.</b> The test has little power — only 3 to 13 targets of each miRNA are on the panel — and the null result neither supports nor argues against a blood–brain relationship.",
      p5_note: "negatives get published: this result is part of the preprint's §4.5, with the power of the test stated",
      p4_t: "Calibrated simulation",
      p4_aria: "Bars of expected damage: no intervention 100 percent; intercellular customs 50.2 percent; total wall 16.3 percent; shielding threshold at 41 percent",
      p4_svg1: "no intervention · 100%",
      p4_svg2: "customs 80/5 · 50.2%",
      p4_svg5: "total wall · 16.3%",
      p4_svg3: "↑ G127V shielding threshold: ~41% coverage collapses the spread (percolation p_c = 0.593)",
      p4_chip: "\"intercellular customs\" hypothesis → Institut Pasteur",
      p4_chip2: "3 testable predictions",
      p4_note: "a hypothesis-generating model — not a clinical prediction; the synthesis did not exist in the literature (verified by broad search before submission)",
      dados_nota: "FULL FIGURES IN THE REPOSITORY: volcanoes (brain and blood), top-25 gene heatmap and the reference-case timeline — with code and MD5 checksums. This page shows summary read-outs; the full statistical detail lives in the reports.",

      foot_nota: "A provisional name. It is born from the urgency of a real, public case in Brazil, announced by the family on Aug 21, 2026 — and from a simple question: why does the right information arrive too late?",
      foot_meta: "living registry · updated Aug 29, 2026<br>contact: coming soon, with formal incorporation<br>bilingual PT/EN · everything verifiable"
    }
  };

  /* ================= aplicação do idioma ================= */
  function applyLang(lang) {
    var dict = I18N[lang] || I18N.pt;
    document.documentElement.lang = lang === "pt" ? "pt-BR" : "en";
    var title = dict["_title_" + PAGE];
    if (title) document.title = title;

    var nodes = document.querySelectorAll("[data-i18n]");
    for (var i = 0; i < nodes.length; i++) {
      var k = nodes[i].getAttribute("data-i18n");
      if (dict[k] !== undefined) nodes[i].textContent = dict[k];
    }
    var nodesH = document.querySelectorAll("[data-i18n-html]");
    for (var j = 0; j < nodesH.length; j++) {
      var kh = nodesH[j].getAttribute("data-i18n-html");
      if (dict[kh] !== undefined) nodesH[j].innerHTML = dict[kh];
    }
    var nodesA = document.querySelectorAll("[data-i18n-aria]");
    for (var m = 0; m < nodesA.length; m++) {
      var ka = nodesA[m].getAttribute("data-i18n-aria");
      if (dict[ka] !== undefined) nodesA[m].setAttribute("aria-label", dict[ka]);
    }

    var btns = document.querySelectorAll(".lang-sw button");
    for (var b = 0; b < btns.length; b++) {
      var ativo = btns[b].getAttribute("data-lang") === lang;
      btns[b].classList.toggle("ativo", ativo);
      btns[b].setAttribute("aria-pressed", ativo ? "true" : "false");
    }
    try { localStorage.setItem("lithos-lang", lang); } catch (e) {}
  }

  var langSalvo = null;
  try { langSalvo = localStorage.getItem("lithos-lang"); } catch (e) {}
  var langInicial = langSalvo || ((navigator.language || "pt").toLowerCase().indexOf("pt") === 0 ? "pt" : "en");
  applyLang(langInicial);

  var sw = document.querySelectorAll(".lang-sw button");
  for (var s = 0; s < sw.length; s++) {
    sw[s].addEventListener("click", function () { applyLang(this.getAttribute("data-lang")); });
  }

  /* ================= menu móvel ================= */
  var hamb = document.getElementById("hamb");
  var mob = document.getElementById("navmob");
  hamb.addEventListener("click", function () {
    var aberto = mob.classList.toggle("aberta");
    hamb.setAttribute("aria-expanded", aberto ? "true" : "false");
    hamb.setAttribute("aria-label", aberto ? "Fechar menu" : "Abrir menu");
  });
  mob.addEventListener("click", function (e) {
    if (e.target.tagName === "A") { mob.classList.remove("aberta"); hamb.setAttribute("aria-expanded", "false"); }
  });

  /* ================= revelação ao rolar ================= */
  var reduz = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var els = document.querySelectorAll(".reveal");
  if (reduz || !("IntersectionObserver" in window)) {
    for (var i = 0; i < els.length; i++) els[i].classList.add("ok");
    return;
  }
  var io = new IntersectionObserver(function (entries) {
    entries.forEach(function (en) {
      if (en.isIntersecting) { en.target.classList.add("ok"); io.unobserve(en.target); }
    });
  }, { threshold: 0.12 });
  for (var k = 0; k < els.length; k++) io.observe(els[k]);
})();
