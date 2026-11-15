#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gera_preprint_pdf.py — Markdown → PDF limpo para submissão (medRxiv/bioRxiv).

Por que existe: o PDF de 2026-08-30 levou caixas ■ para expoentes/superscritos
(⁻ ⁴ ⁵...) e subscritos (₂ ₃ ₄) — fora do Latin-1 das fontes base do PDF.
Este script sanitiza esses glifos ANTES de renderizar:
  10⁻¹⁵ → 10^-15      H₂/H₃/H₄ → H2/H3/H4      π_i → pi_i
Greekglyphs (λ β ρ × → ≈ ≥ √ ∝ ↑ ↓ — –) passam intactos: o ReportLab os
mapeia à fonte Symbol, que os contém (verificado na extração do PDF).

Entrada:  preprint/manuscrito_preprint.md
Saídas:   preprint/manuscrito_preprint.pdf (canonical)
"""
import re
import unicodedata
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (Paragraph, SimpleDocTemplate, Spacer, Table,
                                TableStyle)

BASE = Path(__file__).resolve().parents[2]
SRC = BASE / "preprint" / "manuscrito_preprint.md"
DST = BASE / "preprint" / "manuscrito_preprint.pdf"

SUP = dict(zip("⁻⁰¹²³⁴⁵⁶⁷⁸⁹", "-0123456789"))
SUB = dict(zip("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
SUB_LETRAS = {"ᵢ": "_i", "ᵣ": "_r", "ᵤ": "_u", "ₑ": "_e", "ₒ": "_o",
              "ₓ": "_x", "ₐ": "_a", "ₕ": "_h", "ₖ": "_k", "ₗ": "_l",
              "ₘ": "_m", "ₙ": "_n", "ₚ": "_p", "ₛ": "_s", "ₜ": "_t",
              "ⱼ": "_j"}


def sanear(texto: str) -> str:
    texto = re.sub(r"([⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+)",
                   lambda m: "^" + "".join(SUP[c] for c in m.group(1)), texto)
    texto = re.sub(r"([₀-₉ᵢᵣᵤₑₒₓₐₕₖₗₘₙₚₛₜⱼ]+)",
                   lambda m: "".join(SUB_LETRAS.get(c) or SUB[c]
                                     for c in m.group(1)),
                   texto)
    # Latin Extended (ex.: Areškevičiūtė → Areskeviciute): fora do Latin-1 =
    # caixa no PDF; NFKD remove diacríticos com segurança.
    return "".join(
        ("".join(b for b in unicodedata.normalize("NFKD", c) if ord(b) < 128)
         if 0x100 <= ord(c) <= 0x24F else c)
        for c in texto)


estilos = getSampleStyleSheet()
H1 = ParagraphStyle("h1", parent=estilos["Heading1"], fontSize=15, spaceAfter=8)
H2 = ParagraphStyle("h2", parent=estilos["Heading2"], fontSize=13,
                    spaceBefore=12, spaceAfter=6)
H3 = ParagraphStyle("h3", parent=estilos["Heading3"], fontSize=11.5,
                    spaceBefore=8, spaceAfter=4)
BODY = ParagraphStyle("body", parent=estilos["BodyText"], fontSize=9.5,
                      leading=13, alignment=4)
MONO = ParagraphStyle("mono", parent=estilos["Code"], fontSize=7.5,
                      leading=9.5)
# Bloco de citação do Markdown ("> ..."): usado na nota de correção da v0.5.
NOTA = ParagraphStyle("nota", parent=BODY, leftIndent=10, rightIndent=6,
                      backColor=colors.whitesmoke, borderPadding=3)


def inline(t: str) -> str:
    t = sanear(t)
    t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"\*(.+?)\*", r"<i>\1</i>", t)
    t = re.sub(r"`(.+?)`", r"<font face='Courier'>\1</font>", t)
    return t


def para_tabela(linhas):
    cel = [[Paragraph(inline(c.strip()), MONO) for c in l.split("|")[1:-1]]
           for l in linhas if l.strip("|").strip()]
    cel = [c for c in cel if c and not all(re.match(r"^:?-+:?$", p.text) for p in c)]
    if not cel:
        return None
    ncol = len(cel[0])
    largura = (A4[0] - 4 * cm) / ncol
    t = Table(cel, colWidths=[largura] * ncol)
    t.setStyle(TableStyle([
        ("GRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
    ]))
    return t


def main():
    m = re.search(r"Preprint draft (v[0-9.]+)", SRC.read_text(encoding="utf-8"))
    versao = m.group(1) if m else ""
    doc = SimpleDocTemplate(str(DST), pagesize=A4,
                            leftMargin=2 * cm, rightMargin=2 * cm,
                            topMargin=2 * cm, bottomMargin=2 * cm,
                            title=f"dcj-study preprint {versao}")
    story, pbuf, tabbuf = [], [], []
    estado = {"nota": False}

    def flush_para():
        if pbuf:
            texto = " ".join(pbuf).strip()
            if texto:
                story.append(Paragraph(inline(texto),
                                       NOTA if estado["nota"] else BODY))
                story.append(Spacer(1, 4))
            pbuf.clear()

    def flush_tab():
        if tabbuf:
            t = para_tabela(tabbuf)
            if t:
                story.append(t)
                story.append(Spacer(1, 6))
            tabbuf.clear()

    for raw in SRC.read_text(encoding="utf-8").splitlines():
        l = raw.rstrip()
        if l.startswith("|") and l.endswith("|"):
            flush_para()
            tabbuf.append(l)
            continue
        flush_tab()
        # bloco de citação: tira o marcador e muda o estilo do parágrafo
        em_nota = l.startswith(">")
        if em_nota != estado["nota"]:
            flush_para()
            estado["nota"] = em_nota
        if em_nota:
            l = l[1:].strip()
        # item de lista: cada item é um parágrafo próprio, com marcador
        if l.startswith("- "):
            flush_para()
            pbuf.append("• " + l[2:])
            continue
        if l.startswith("### "):
            flush_para(); story.append(Paragraph(inline(l[4:]), H3))
        elif l.startswith("## "):
            flush_para(); story.append(Paragraph(inline(l[3:]), H2))
        elif l.startswith("# "):
            flush_para(); story.append(Paragraph(inline(l[2:]), H1))
        elif l.strip() in ("---", ""):
            flush_para()
        else:
            pbuf.append(l)
    flush_para(); flush_tab()
    doc.build(story)
    print(f"[ok] {DST}")


if __name__ == "__main__":
    main()
