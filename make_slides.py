"""
Generate DeSurv PowerPoint presentation.
10-minute talk for statisticians familiar with genomics but not NMF/semi-supervised NMF.
Minimum font size: 16pt throughout.
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.enum.text import PP_ALIGN
from pptx.dml.color import RGBColor
from lxml import etree
from pptx.oxml.ns import qn

# ── Palette ──────────────────────────────────────────────────────────────────
NAVY   = RGBColor(0x1B, 0x35, 0x5E)
TEAL   = RGBColor(0x00, 0x7A, 0x7C)
ORANGE = RGBColor(0xE0, 0x7B, 0x1C)
LGRAY  = RGBColor(0xF2, 0xF4, 0xF7)
MGRAY  = RGBColor(0xCC, 0xD0, 0xD8)
DKGRAY = RGBColor(0x44, 0x4C, 0x5A)
WHITE  = RGBColor(0xFF, 0xFF, 0xFF)
RED    = RGBColor(0xC0, 0x39, 0x2B)
GOLD   = RGBColor(0xF5, 0xC5, 0x18)

MIN_FONT = 16   # global minimum


def fs(size):
    """Return font size, enforcing minimum."""
    return max(size, MIN_FONT)


# ── Low-level helpers ────────────────────────────────────────────────────────

def add_rect(slide, left, top, width, height, fill_color, line_color=None):
    shape = slide.shapes.add_shape(1, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if line_color:
        shape.line.color.rgb = line_color
    else:
        shape.line.fill.background()
    return shape


def add_text_box(slide, text, left, top, width, height,
                 font_size=18, bold=False, italic=False,
                 color=DKGRAY, align=PP_ALIGN.LEFT,
                 wrap=True, font_name="Calibri"):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = wrap
    para = tf.paragraphs[0]
    para.alignment = align
    run = para.add_run()
    run.text = text
    run.font.size = Pt(fs(font_size))
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    run.font.name = font_name
    return txBox


def add_bullet_box(slide, bullets, left, top, width, height,
                   font_size=17, color=DKGRAY, indent_levels=None,
                   font_name="Calibri"):
    """Bullet list with proper PPTX bullet formatting (not symbol characters).
    indent_levels: 0 = top-level bullet, 1 = sub-bullet."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    if indent_levels is None:
        indent_levels = [0] * len(bullets)

    # EMU constants for indentation
    MAR_L_BASE = 342900   # ~0.375 inch per indent level
    HANG_EMU   = 228600   # ~0.25 inch hanging indent

    for i, (text, level) in enumerate(zip(bullets, indent_levels)):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()

        pPr = para._p.get_or_add_pPr()

        # Left margin and hanging indent (creates bullet indent)
        mar_l = MAR_L_BASE * (level + 1)
        pPr.set('marL', str(mar_l))
        pPr.set('indent', str(-HANG_EMU))

        # Remove any buNone that python-pptx may have inserted
        for bn in pPr.findall(qn('a:buNone')):
            pPr.remove(bn)

        # Real PPTX bullet character
        bullet_char = '\u2022' if level == 0 else '\u2013'   # • or –
        buChar = etree.SubElement(pPr, qn('a:buChar'))
        buChar.set('char', bullet_char)

        # Space before top-level bullets (skip first)
        if level == 0 and i > 0:
            spcBef = etree.SubElement(pPr, qn('a:spcBef'))
            spcPts = etree.SubElement(spcBef, qn('a:spcPts'))
            spcPts.set('val', '120')

        # Split on \n: first part goes into the run, subsequent parts follow
        # an <a:br/> line-break element so they stay in the same bullet
        parts = text.split('\n')
        run = para.add_run()
        run.text = parts[0].strip()
        run.font.size = Pt(fs(font_size - 1 * level))
        run.font.color.rgb = color
        run.font.name = font_name

        for part in parts[1:]:
            # Line break within the paragraph (no new bullet)
            br = etree.SubElement(para._p, qn('a:br'))
            rPr_src = run._r.find(qn('a:rPr'))
            if rPr_src is not None:
                import copy
                br_rPr = copy.deepcopy(rPr_src)
                br.append(br_rPr)
            # Continuation run (inherits same formatting)
            cont_run = para.add_run()
            cont_run.text = part.strip()
            cont_run.font.size = Pt(fs(font_size - 1 * level))
            cont_run.font.color.rgb = color
            cont_run.font.name = font_name

    return txBox


# ── Composite layout helpers ─────────────────────────────────────────────────

def add_header_bar(slide, title, subtitle=None):
    W = Inches(13.33)
    H = Inches(1.2) if subtitle else Inches(0.9)
    add_rect(slide, 0, 0, W, H, NAVY)
    add_text_box(slide, title,
                 Inches(0.35), Inches(0.08), Inches(12.6), Inches(0.55),
                 font_size=26, bold=True, color=WHITE)
    if subtitle:
        add_text_box(slide, subtitle,
                     Inches(0.35), Inches(0.68), Inches(12.6), Inches(0.42),
                     font_size=16, italic=True, color=MGRAY)


def setup_layout_footer(prs, text="Young et al.  |  DeSurv  |  2026"):
    """Inject a footer text box into the blank slide layout XML so every slide inherits it."""
    from pptx.oxml.ns import nsmap as _nsmap

    layout = prs.slide_layouts[6]   # Blank layout used by all our slides

    # Dimensions in EMU
    LEFT  = int(Inches(0.2))
    TOP   = int(Inches(7.1))
    CX    = int(Inches(13.0))
    CY    = int(Inches(0.38))
    # Colour as hex string
    r, g, b = MGRAY
    clr_hex = f"{r:02X}{g:02X}{b:02X}"

    sp_xml = f"""\
<p:sp xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"
      xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main"
      xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <p:nvSpPr>
    <p:cNvPr id="9001" name="SlideFooter"/>
    <p:cNvSpPr txBox="1"><a:spLocks noGrp="1"/></p:cNvSpPr>
    <p:nvPr/>
  </p:nvSpPr>
  <p:spPr>
    <a:xfrm><a:off x="{LEFT}" y="{TOP}"/><a:ext cx="{CX}" cy="{CY}"/></a:xfrm>
    <a:prstGeom prst="rect"><a:avLst/></a:prstGeom>
    <a:noFill/>
  </p:spPr>
  <p:txBody>
    <a:bodyPr wrap="none" lIns="0" tIns="0" rIns="0" bIns="0"/>
    <a:lstStyle/>
    <a:p>
      <a:pPr algn="ctr"/>
      <a:r>
        <a:rPr lang="en-US" sz="1600" dirty="0">
          <a:solidFill><a:srgbClr val="{clr_hex}"/></a:solidFill>
          <a:latin typeface="Calibri"/>
        </a:rPr>
        <a:t>{text}</a:t>
      </a:r>
    </a:p>
  </p:txBody>
</p:sp>"""

    sp_elem = etree.fromstring(sp_xml)
    layout.shapes._spTree.append(sp_elem)


def slide_num(slide, n):
    add_text_box(slide, str(n), Inches(12.8), Inches(7.1),
                 Inches(0.4), Inches(0.3), font_size=16, color=MGRAY)


def blank_slide(prs):
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = WHITE
    return slide


def make_prs():
    prs = Presentation()
    prs.slide_width  = Inches(13.33)
    prs.slide_height = Inches(7.5)
    return prs


# ── Single-cell callout (reused in two slides) ───────────────────────────────

def add_sc_callout(slide, left, top, width, height):
    """Gold-bordered callout box for the single-cell / bulk RNA argument."""
    add_rect(slide, left, top, width, height, RGBColor(0xFF, 0xF8, 0xE1), GOLD)
    add_text_box(slide, "Why not just use single-cell data?",
                 left + Inches(0.12), top + Inches(0.1),
                 width - Inches(0.24), Inches(0.36),
                 font_size=17, bold=True, color=RGBColor(0x7B, 0x5B, 0x00))
    add_text_box(slide,
                 "scRNA-seq & spatial transcriptomics now resolve individual cell populations — "
                 "but bulk expression remains the primary source of large, clinically annotated "
                 "cohorts with the sample sizes needed for stable survival modeling.\n"
                 "Hundreds of patients with mature follow-up are available in bulk, rarely in sc.",
                 left + Inches(0.12), top + Inches(0.52),
                 width - Inches(0.24), height - Inches(0.58),
                 font_size=17, color=DKGRAY)


# =============================================================================
# SLIDES
# =============================================================================

def slide_title(prs):
    slide = blank_slide(prs)
    add_rect(slide, 0, 0, Inches(13.33), Inches(7.5), NAVY)
    add_rect(slide, Inches(0.5), Inches(2.0), Inches(12.33), Pt(3), ORANGE)

    add_text_box(slide,
                 "Survival-Guided Matrix Factorization\nIdentifies Reproducible Prognostic Programs\nin Pancreatic Cancer",
                 Inches(0.5), Inches(2.12), Inches(12.33), Inches(2.5),
                 font_size=32, bold=True, color=WHITE)

    add_text_box(slide, "Amber M. Young",
                 Inches(0.5), Inches(4.72), Inches(7.0), Inches(0.48),
                 font_size=22, color=MGRAY)
    add_text_box(slide, "Department of Biostatistics  |  University of North Carolina at Chapel Hill",
                 Inches(0.5), Inches(5.24), Inches(10.0), Inches(0.4),
                 font_size=18, color=MGRAY)
    add_text_box(slide, "R package: github.com/ayoung31/DeSurv",
                 Inches(0.5), Inches(6.55), Inches(8.0), Inches(0.38),
                 font_size=16, italic=True, color=RGBColor(0x99, 0xBB, 0xFF))
    slide_num(slide, 1)
    return slide


def slide_motivation(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "The Problem: Variance ≠ Prognosis",
                   "Standard NMF maximizes explained variance — but variance-dominant signals are often not prognostic")

    # Left: problem bullets
    add_bullet_box(slide,
                   ["Bulk tumor expression is a mixture: malignant cells,\n  fibroblasts, immune infiltrate",
                    "Standard NMF minimizes reconstruction error\n  → captures highest-variance signals",
                    "Highest-variance signals ≈ tissue composition / tumor purity\n  → often outcome-neutral",
                    "Standard workflow: discover factors unsupervised,\n  then evaluate clinical relevance retrospectively",
                    "Misalignment between discovery objective and\n  evaluation criterion is a structural cost"],
                   Inches(0.35), Inches(1.35), Inches(6.1), Inches(4.0),
                   font_size=17, color=DKGRAY)

    # PDAC example box
    add_rect(slide, Inches(0.35), Inches(5.45), Inches(6.1), Inches(1.2), LGRAY, ORANGE)
    add_text_box(slide, "In PDAC specifically:",
                 Inches(0.5), Inches(5.52), Inches(5.8), Inches(0.3),
                 font_size=17, bold=True, color=ORANGE)
    add_text_box(slide,
                 "Exocrine cell signal dominates variance. Reflects tumor purity — not cancer biology. Prognostically neutral.",
                 Inches(0.5), Inches(5.86), Inches(5.9), Inches(0.7),
                 font_size=17, color=DKGRAY)

    # Key tension box
    add_rect(slide, Inches(6.6), Inches(1.35), Inches(6.35), Inches(2.0), LGRAY, NAVY)
    add_text_box(slide, "Objective mismatch:",
                 Inches(6.75), Inches(1.43), Inches(6.1), Inches(0.32),
                 font_size=17, bold=True, color=NAVY)
    add_text_box(slide,
                 "Minimizes:    ‖ X − WH ‖²   (reconstruction error)\n\nEvaluates:    concordance with survival (C-index)",
                 Inches(6.75), Inches(1.82), Inches(6.1), Inches(1.3),
                 font_size=17, color=DKGRAY)

    # Single-cell callout (prominent)
    add_sc_callout(slide, Inches(6.6), Inches(3.5), Inches(6.35), Inches(3.15))


    slide_num(slide, 2)
    return slide


def slide_nmf_background(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "NMF Background",
                   "Nonnegative Matrix Factorization — additive, interpretable decomposition")

    # Equation
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(12.6), Inches(1.1), LGRAY, MGRAY)
    add_text_box(slide, "X  ≈  W · H",
                 Inches(0.5), Inches(1.38), Inches(5.0), Inches(0.65),
                 font_size=30, bold=True, color=NAVY)
    add_text_box(slide, "(p × n)   =   (p × k)  ×  (k × n)     W, H ≥ 0",
                 Inches(5.2), Inches(1.52), Inches(7.5), Inches(0.45),
                 font_size=18, color=DKGRAY)
    add_text_box(slide, "p = genes,  n = samples,  k = latent factors",
                 Inches(0.5), Inches(2.08), Inches(7.0), Inches(0.3),
                 font_size=16, italic=True, color=DKGRAY)

    col_w = Inches(4.0)
    col_h = Inches(4.35)
    top   = Inches(2.6)
    gap   = Inches(0.24)

    # W box
    add_rect(slide, Inches(0.35), top, col_w, col_h, LGRAY, TEAL)
    add_text_box(slide, "W — Gene Programs  (p × k)",
                 Inches(0.5), top + Inches(0.1), col_w - Inches(0.2), Inches(0.38),
                 font_size=17, bold=True, color=TEAL)
    add_bullet_box(slide,
                   ["W[g, j] = weight of gene g in program j",
                    "Shared across all samples",
                    "Nonneg → interpret as membership / loading",
                    "In DeSurv: survival gradient acts on W"],
                   Inches(0.5), top + Inches(0.55), col_w - Inches(0.2), col_h - Inches(0.6),
                   font_size=16, color=DKGRAY)

    # H box
    left_h = Inches(0.35) + col_w + gap
    add_rect(slide, left_h, top, col_w, col_h, LGRAY, ORANGE)
    add_text_box(slide, "H — Sample Loadings  (k × n)",
                 left_h + Inches(0.15), top + Inches(0.1), col_w - Inches(0.2), Inches(0.38),
                 font_size=17, bold=True, color=ORANGE)
    add_bullet_box(slide,
                   ["H[j, i] = activity of program j in sample i",
                    "Nonneg → interpret as mixture coefficient",
                    "Sample-specific",
                    "In DeSurv: NO survival gradient\n  → preserves mixture interpretation"],
                   left_h + Inches(0.15), top + Inches(0.55), col_w - Inches(0.2), col_h - Inches(0.6),
                   font_size=16, color=DKGRAY)

    # Z box
    left_z = left_h + col_w + gap
    add_rect(slide, left_z, top, col_w, col_h, LGRAY, NAVY)
    add_text_box(slide, "Z = Wᵀ X — Factor Scores  (k × n)",
                 left_z + Inches(0.15), top + Inches(0.1), col_w - Inches(0.2), Inches(0.38),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["Z[j, i] = expression of program j in sample i",
                    "Covariates in the Cox model",
                    "New patients: Z_new = W̃ᵀ X_new\n  (projection, no re-fitting)"],
                   left_z + Inches(0.15), top + Inches(0.55), col_w - Inches(0.2), col_h - Inches(0.6),
                   font_size=16, color=DKGRAY)


    slide_num(slide, 3)
    return slide


def slide_desurv_model(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "DeSurv: Survival-Supervised Deconvolution",
                   "Joint NMF + Cox proportional hazards — α controls supervision strength")

    # Objective
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(12.6), Inches(1.3), LGRAY, NAVY)
    add_text_box(slide, "Joint objective:",
                 Inches(0.5), Inches(1.38), Inches(3.0), Inches(0.3),
                 font_size=16, bold=True, color=NAVY)
    add_text_box(slide,
                 "ℒ(W, H, β)  =  (1 − α) · ℒ_NMF(W, H)  −  α · ℒ_Cox(W, β)",
                 Inches(0.5), Inches(1.72), Inches(12.2), Inches(0.65),
                 font_size=22, bold=True, color=NAVY)
    add_text_box(slide,
                 "α ∈ [0, 1)   |   α = 0 → standard NMF   |   larger α → more survival supervision",
                 Inches(0.5), Inches(2.38), Inches(12.2), Inches(0.3),
                 font_size=16, italic=True, color=TEAL)

    col_w = Inches(4.0)
    col_h = Inches(3.75)
    top   = Inches(2.8)
    gap   = Inches(0.24)

    # Box 1: where supervision enters
    add_rect(slide, Inches(0.35), top, col_w, col_h, LGRAY, TEAL)
    add_text_box(slide, "1. Where supervision enters",
                 Inches(0.5), top + Inches(0.08), col_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=TEAL)
    add_bullet_box(slide,
                   ["Cox likelihood is a function of W and β",
                    "Survival gradient acts on W only",
                    "H gets no survival gradient\n  → preserves mixture interpretation",
                    "W is directed toward outcome-relevant structure"],
                   Inches(0.5), top + Inches(0.5), col_w - Inches(0.2), col_h - Inches(0.55),
                   font_size=16, color=DKGRAY)

    # Box 2: model selection
    left2 = Inches(0.35) + col_w + gap
    add_rect(slide, left2, top, col_w, col_h, LGRAY, ORANGE)
    add_text_box(slide, "2. Model selection",
                 left2 + Inches(0.15), top + Inches(0.08), col_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=ORANGE)
    add_bullet_box(slide,
                   ["Bayesian optimization over (k, α, λ, ξ, n_top)",
                    "Criterion: cross-validated C-index",
                    "1-SE rule → smallest k within 1 SE of maximum",
                    "PDAC result: k = 3,  α ≈ 0.37",
                    "All selection within training data"],
                   left2 + Inches(0.15), top + Inches(0.5), col_w - Inches(0.2), col_h - Inches(0.55),
                   font_size=16, color=DKGRAY)

    # Box 3: out-of-sample
    left3 = left2 + col_w + gap
    add_rect(slide, left3, top, col_w, col_h, LGRAY, NAVY)
    add_text_box(slide, "3. Out-of-sample scoring",
                 left3 + Inches(0.15), top + Inches(0.08), col_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["W fixed at training time",
                    "New patients: Z_new = W̃ᵀ X_new",
                    "No re-fitting; no validation survival data needed",
                    "Risk score = Ẑβ  (linear predictor)"],
                   left3 + Inches(0.15), top + Inches(0.5), col_w - Inches(0.2), col_h - Inches(0.55),
                   font_size=16, color=DKGRAY)


    slide_num(slide, 4)
    return slide


def slide_simulations(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Simulation: DeSurv Wins When Variance ≠ Prognosis",
                   "100 replicates | p = 3,000 genes, n = 200 samples, true k = 3 | both methods tuned by BO over CV C-index")

    # Design strip
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(12.6), Inches(0.9), LGRAY, MGRAY)
    add_text_box(slide,
                 "Primary scenario: prognostic programs explain LOW variance vs. outcome-neutral background genes.  "
                 "β = (2, 0, 0)  →  only factor 1 drives survival.",
                 Inches(0.5), Inches(1.38), Inches(12.2), Inches(0.72),
                 font_size=17, color=DKGRAY)

    box_w = Inches(4.0)
    box_h = Inches(4.55)
    top   = Inches(2.3)
    gap   = Inches(0.2)

    # C-index
    add_rect(slide, Inches(0.35), top, box_w, box_h, LGRAY, TEAL)
    add_text_box(slide, "C-index (test set)",
                 Inches(0.5), top + Inches(0.1), box_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=TEAL)
    add_text_box(slide,
                 "DeSurv median:  0.80+\nNMF median:       ~0.62\nΔ ≈ 0.18\nPaired Wilcoxon P < 0.001",
                 Inches(0.5), top + Inches(0.55), box_w - Inches(0.2), Inches(1.1),
                 font_size=17, color=DKGRAY)
    add_rect(slide, Inches(0.5), top + Inches(1.75), box_w - Inches(0.3), Inches(2.5), WHITE, TEAL)
    add_text_box(slide, "[Fig. 2A — C-index box plots]",
                 Inches(0.7), top + Inches(2.85), box_w - Inches(0.6), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Precision
    left2 = Inches(0.35) + box_w + gap
    add_rect(slide, left2, top, box_w, box_h, LGRAY, ORANGE)
    add_text_box(slide, "Gene Program Precision",
                 left2 + Inches(0.15), top + Inches(0.1), box_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=ORANGE)
    add_text_box(slide,
                 "Fraction of top-weighted genes in learned factor\noverlapping the true prognostic program\n\nNMF: near zero  (variance-driven, misses signal)\nDeSurv: substantially higher",
                 left2 + Inches(0.15), top + Inches(0.55), box_w - Inches(0.2), Inches(1.3),
                 font_size=17, color=DKGRAY)
    add_rect(slide, left2 + Inches(0.15), top + Inches(1.95), box_w - Inches(0.3), Inches(2.3), WHITE, ORANGE)
    add_text_box(slide, "[Fig. 2B — Precision box plots]",
                 left2 + Inches(0.3), top + Inches(2.95), box_w - Inches(0.6), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Scaling behavior
    left3 = left2 + box_w + gap
    add_rect(slide, left3, top, box_w, box_h, LGRAY, NAVY)
    add_text_box(slide, "Advantage Scales with Signal Gap",
                 left3 + Inches(0.15), top + Inches(0.1), box_w - Inches(0.2), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["Primary (variance ≠ prognosis): Δ ≈ 0.18",
                    "Mixed (partial overlap): Δ attenuated",
                    "Null (β = 0): both ≈ 0.5  (no spurious gains)",
                    "k selection: DeSurv → k = 3 reliably;\n  NMF under-selects (k ≈ 2)"],
                   left3 + Inches(0.15), top + Inches(0.55), box_w - Inches(0.2), box_h - Inches(0.6),
                   font_size=16, color=DKGRAY)


    slide_num(slide, 5)
    return slide


def slide_model_selection(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Model Selection in PDAC: BO + 1-SE Rule",
                   "Training cohorts: TCGA-PAAD + CPTAC-3")

    # Left: problem
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(5.9), Inches(5.8), LGRAY, MGRAY)
    add_text_box(slide, "Standard NMF rank selection in PDAC:",
                 Inches(0.5), Inches(1.38), Inches(5.65), Inches(0.35),
                 font_size=17, bold=True, color=RED)
    add_bullet_box(slide,
                   ["Reconstruction residuals: no clear elbow",
                    "Cophenetic correlation: ambiguous",
                    "Silhouette width: inconsistent",
                    "→ No single k clearly indicated"],
                   Inches(0.5), Inches(1.82), Inches(5.65), Inches(2.0),
                   font_size=17, color=DKGRAY)
    add_rect(slide, Inches(0.5), Inches(3.95), Inches(5.6), Inches(2.85), WHITE, MGRAY)
    add_text_box(slide, "[SI Fig. S4 — NMF heuristics (all inconsistent)]",
                 Inches(1.1), Inches(5.2), Inches(4.3), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Arrow
    add_text_box(slide, "→", Inches(6.0), Inches(3.5),
                 Inches(0.6), Inches(0.6), font_size=30, bold=True, color=TEAL,
                 align=PP_ALIGN.CENTER)

    # Right: solution
    add_rect(slide, Inches(6.5), Inches(1.3), Inches(6.45), Inches(5.8), LGRAY, TEAL)
    add_text_box(slide, "DeSurv: BO over CV C-index",
                 Inches(6.65), Inches(1.38), Inches(6.2), Inches(0.35),
                 font_size=17, bold=True, color=TEAL)
    add_bullet_box(slide,
                   ["Bayesian optimization explores (k, α, λ, …)",
                    "Objective = cross-validated C-index\n  (matches the clinical goal)",
                    "GP-predicted mean C-index over k × α grid",
                    "1-SE rule: smallest k within 1 SE of max"],
                   Inches(6.65), Inches(1.82), Inches(6.2), Inches(2.0),
                   font_size=17, color=DKGRAY)

    add_rect(slide, Inches(6.65), Inches(4.0), Inches(6.1), Inches(0.9), NAVY)
    add_text_box(slide,
                 "PDAC result:   k = 3,   α ≈ 0.37   (CV C-index = 0.63)",
                 Inches(6.8), Inches(4.12), Inches(5.8), Inches(0.65),
                 font_size=18, bold=True, color=WHITE)

    add_rect(slide, Inches(6.65), Inches(5.05), Inches(6.1), Inches(1.8), WHITE, TEAL)
    add_text_box(slide, "[Fig. 2D — GP C-index heatmap over k × α]",
                 Inches(7.5), Inches(5.85), Inches(4.3), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)


    slide_num(slide, 6)
    return slide


def slide_pdac_factors(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "PDAC Factor Structure: DeSurv vs. Standard NMF  (both k = 3)",
                   "Gene program identity measured by Spearman correlation to established PDAC reference signatures")

    # DeSurv
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(6.1), Inches(5.85), LGRAY, TEAL)
    add_text_box(slide, "DeSurv  (k = 3,  α = 0.37)",
                 Inches(0.5), Inches(1.38), Inches(5.9), Inches(0.35),
                 font_size=17, bold=True, color=TEAL)
    add_bullet_box(slide,
                   ["D1  — classical tumor + restCAF/iCAF stroma",
                    "Tumor-stroma coupling in one factor",
                    "MOST prognostic factor (largest ΔCox)",
                    "D2  — proCAF / activated immune-stromal",
                    "D3  — basal-like tumor programs",
                    "No exocrine/acinar factor\n  → purity signal suppressed by survival gradient"],
                   Inches(0.5), Inches(1.82), Inches(5.9), Inches(5.1),
                   font_size=17, color=DKGRAY,
                   indent_levels=[0, 1, 1, 0, 0, 0])

    # NMF
    add_rect(slide, Inches(6.6), Inches(1.3), Inches(6.35), Inches(5.85), LGRAY, RED)
    add_text_box(slide, "Standard NMF  (k = 3)",
                 Inches(6.75), Inches(1.38), Inches(6.1), Inches(0.35),
                 font_size=17, bold=True, color=RED)
    add_bullet_box(slide,
                   ["N1  — classical tumor programs",
                    "N2  — exocrine / acinar programs",
                    "High variance; reflects tumor purity",
                    "Prognostically NEUTRAL",
                    "1 of 3 factors devoted to purity signal",
                    "N3  — mixed microenvironment",
                    "Aggregates immune + fibroblast + ECM\n  → does not separate CAF subtypes",
                    "No direct counterpart to DeSurv D1"],
                   Inches(6.75), Inches(1.82), Inches(6.1), Inches(5.1),
                   font_size=17, color=DKGRAY,
                   indent_levels=[0, 0, 1, 1, 1, 0, 0, 0])


    slide_num(slide, 7)
    return slide


def slide_variance_survival(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Variance vs. Survival Contribution in PDAC",
                   "Semi-partial R² (expression variance)  vs.  ΔCox log-likelihood (survival contribution per factor)")

    # Figure placeholder
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(7.6), Inches(5.85), LGRAY, MGRAY)
    add_text_box(slide,
                 "[Fig. 3C — scatter plot\n x: variance explained\n y: ΔCox partial log-likelihood\n\n Each point = one factor\n DeSurv (circles)  vs  NMF (triangles)]",
                 Inches(1.7), Inches(3.6), Inches(4.5), Inches(1.8),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Findings
    add_rect(slide, Inches(8.15), Inches(1.3), Inches(4.85), Inches(5.85), LGRAY, NAVY)
    add_text_box(slide, "Key findings:",
                 Inches(8.3), Inches(1.38), Inches(4.6), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["NMF highest-variance factor (N2, exocrine):",
                    "~38% expression variance",
                    "ΔCox ≈ 0  (no survival contribution)",
                    "DeSurv most prognostic factor (D1):",
                    "Moderate variance ~18%",
                    "ΔCox >> 0  (all survival signal here)",
                    "All NMF factors: negligible ΔCox\n  despite spanning range of variance",
                    "→ Variance and prognosis diverge in PDAC"],
                   Inches(8.3), Inches(1.82), Inches(4.6), Inches(5.1),
                   font_size=16, color=DKGRAY,
                   indent_levels=[0, 1, 1, 0, 1, 1, 0, 0])


    slide_num(slide, 8)
    return slide


def slide_validation(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "External Validation: 5 Independent PDAC Cohorts",
                   "Scored by projection only — no retraining, no access to validation survival data during training")

    # Forest plot placeholder
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(7.9), Inches(5.85), LGRAY, MGRAY)
    add_text_box(slide,
                 "[Fig. 4A — forest plot\n Per-cohort hazard ratios\n\n DeSurv (D1–D3) vs NMF (N1–N3)\n 5 validation cohorts\n + pooled estimate (black diamond)]",
                 Inches(2.0), Inches(3.6), Inches(4.3), Inches(1.8),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Cohorts
    add_rect(slide, Inches(8.45), Inches(1.3), Inches(4.55), Inches(2.15), LGRAY, NAVY)
    add_text_box(slide, "Validation cohorts:",
                 Inches(8.6), Inches(1.38), Inches(4.3), Inches(0.32),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["Puleo (array)",
                    "Moffitt (GEO array)",
                    "PACA-AU (RNA-seq)",
                    "PACA-AU (array)",
                    "Dijk (array)"],
                   Inches(8.6), Inches(1.78), Inches(4.3), Inches(1.55),
                   font_size=16, color=DKGRAY)

    # Pooled result
    add_rect(slide, Inches(8.45), Inches(3.6), Inches(4.55), Inches(1.3), NAVY)
    add_text_box(slide, "Pooled HR per SD  (DeSurv D1):",
                 Inches(8.6), Inches(3.7), Inches(4.3), Inches(0.35),
                 font_size=17, bold=True, color=WHITE)
    add_text_box(slide,
                 "HR = 1.45  (95% CI 1.29–1.63)\nP < 0.001   (stratified Cox)",
                 Inches(8.6), Inches(4.08), Inches(4.3), Inches(0.72),
                 font_size=18, bold=True, color=ORANGE)

    # NMF comparison
    add_rect(slide, Inches(8.45), Inches(5.05), Inches(4.55), Inches(2.1), LGRAY, RED)
    add_text_box(slide, "Standard NMF at same rank (k = 3):",
                 Inches(8.6), Inches(5.13), Inches(4.3), Inches(0.35),
                 font_size=17, bold=True, color=RED)
    add_text_box(slide,
                 "Weaker pooled effect\nHR per SD ~1.2–1.3\nKM: HR 1.59  (95% CI 1.03–2.45)",
                 Inches(8.6), Inches(5.55), Inches(4.3), Inches(1.4),
                 font_size=17, color=DKGRAY)


    slide_num(slide, 9)
    return slide


def slide_km_curves(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Kaplan–Meier Curves: Pooled External Validation",
                   "Risk groups dichotomized at cross-validated cutpoint learned from training data")

    # DeSurv KM
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(6.1), Inches(5.85), LGRAY, TEAL)
    add_text_box(slide, "DeSurv",
                 Inches(0.5), Inches(1.38), Inches(5.9), Inches(0.35),
                 font_size=18, bold=True, color=TEAL)
    add_rect(slide, Inches(0.5), Inches(1.82), Inches(5.8), Inches(4.95), WHITE, TEAL)
    add_text_box(slide,
                 "[Fig. 4B — KM curves\n high vs low risk\n pooled validation\n\n HR = 2.17  (95% CI 1.68–2.80)\n P < 0.001]",
                 Inches(1.3), Inches(3.8), Inches(3.5), Inches(1.6),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # NMF KM
    add_rect(slide, Inches(6.6), Inches(1.3), Inches(6.35), Inches(5.85), LGRAY, RED)
    add_text_box(slide, "Standard NMF  (k = 3)",
                 Inches(6.75), Inches(1.38), Inches(6.1), Inches(0.35),
                 font_size=18, bold=True, color=RED)
    add_rect(slide, Inches(6.75), Inches(1.82), Inches(6.1), Inches(4.95), WHITE, RED)
    add_text_box(slide,
                 "[Fig. 4C — KM curves\n high vs low risk\n pooled validation\n\n HR = 1.59  (95% CI 1.03–2.45)\n P = 0.036]",
                 Inches(8.1), Inches(3.8), Inches(3.5), Inches(1.6),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    add_rect(slide, Inches(0.35), Inches(6.95), Inches(12.65), Inches(0.32), NAVY)
    add_text_box(slide, "DeSurv: clearer risk stratification with better-separated survival curves",
                 Inches(0.5), Inches(6.97), Inches(12.4), Inches(0.28),
                 font_size=17, bold=True, color=WHITE, align=PP_ALIGN.CENTER)


    slide_num(slide, 10)
    return slide


def slide_robustness(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "K = 3 Is Robust and Provides Independent Prognostic Value",
                   "Sensitivity analysis: K ∈ {2,3,5,7,9} × α ∈ {0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95}  (35 combinations)")

    # Left: table summary
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(7.1), Inches(5.85), LGRAY, MGRAY)
    add_text_box(slide, "Significance in external validation (adjusted for PurIST, DeCAF):",
                 Inches(0.5), Inches(1.38), Inches(6.85), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)

    for i, (k, n_sig, note, bc) in enumerate([
        ("K = 3", "4 / 7 α values", "Signal independent of known classifiers", TEAL),
        ("K = 7", "2 / 7 α values", "K=7 content largely recoverable from PurIST/DeCAF", RED),
    ]):
        top_r = Inches(1.88) + i * Inches(2.2)
        add_rect(slide, Inches(0.5), top_r, Inches(6.7), Inches(2.0), WHITE, bc)
        add_text_box(slide, k,
                     Inches(0.65), top_r + Inches(0.12), Inches(1.2), Inches(0.38),
                     font_size=18, bold=True, color=bc)
        add_text_box(slide, f"Significant at {n_sig}",
                     Inches(0.65), top_r + Inches(0.58), Inches(6.2), Inches(0.35),
                     font_size=17, color=DKGRAY)
        add_text_box(slide, note,
                     Inches(0.65), top_r + Inches(1.0), Inches(6.2), Inches(0.8),
                     font_size=17, italic=True, color=DKGRAY)

    add_rect(slide, Inches(0.35), Inches(6.35), Inches(7.1), Inches(0.8), NAVY)
    add_text_box(slide,
                 "DeSurv K = 3 provides prognostic value INDEPENDENT of PurIST / DeCAF classifiers",
                 Inches(0.5), Inches(6.43), Inches(6.9), Inches(0.65),
                 font_size=17, bold=True, color=WHITE)

    # Right: why k=3 + biology
    add_rect(slide, Inches(7.6), Inches(1.3), Inches(5.35), Inches(3.35), LGRAY, NAVY)
    add_text_box(slide, "Why k = 3 works:",
                 Inches(7.75), Inches(1.38), Inches(5.1), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["Supervision creates broad concordance plateau:\n  k = 3–12 within 1 SE in PDAC",
                    "1-SE rule selects most parsimonious k on plateau",
                    "Standard NMF: concordance increases steadily\n  → needs k = 7 to approach DeSurv"],
                   Inches(7.75), Inches(1.82), Inches(5.1), Inches(2.65),
                   font_size=16, color=DKGRAY)

    add_rect(slide, Inches(7.6), Inches(4.8), Inches(5.35), Inches(2.35), LGRAY, TEAL)
    add_text_box(slide, "Biology of D1:",
                 Inches(7.75), Inches(4.88), Inches(5.1), Inches(0.35),
                 font_size=17, bold=True, color=TEAL)
    add_bullet_box(slide,
                   ["Classical tumor + restCAF (iCAF) stroma in one factor",
                    "Tumor-intrinsic and CAF subtypes jointly drive PDAC prognosis",
                    "Not recoverable from existing classifiers alone"],
                   Inches(7.75), Inches(5.3), Inches(5.1), Inches(1.7),
                   font_size=16, color=DKGRAY)


    slide_num(slide, 11)
    return slide


def slide_conclusions(prs):
    slide = blank_slide(prs)
    add_rect(slide, 0, 0, Inches(13.33), Inches(7.5), NAVY)
    add_rect(slide, 0, 0, Inches(13.33), Inches(1.1), TEAL)
    add_text_box(slide, "Summary",
                 Inches(0.35), Inches(0.12), Inches(12.6), Inches(0.75),
                 font_size=28, bold=True, color=WHITE)

    col_w = Inches(3.95)
    col_h = Inches(2.85)
    top   = Inches(1.25)
    gap   = Inches(0.18)
    DIM   = RGBColor(0x1E, 0x3D, 0x6A)
    LT    = RGBColor(0xCC, 0xDD, 0xFF)

    for i, (title, bc, bullets) in enumerate([
        ("Methodological", TEAL, [
            "DeSurv = NMF + Cox, joint optimization",
            "Survival gradient on W, not H\n  → mixture interpretation preserved",
            "Projection scoring: Z = W̃ᵀ X_new",
            "α = 0 recovers standard NMF",
        ]),
        ("Empirical", ORANGE, [
            "Simulations: advantage when variance ≠ prognosis;\n  zero advantage under null",
            "Pooled validation HR per SD = 1.45\n  (5 cohorts, P < 0.001)",
            "k = 3 independent of PurIST / DeCAF",
        ]),
        ("Biological", WHITE, [
            "D1 = classical tumor + restCAF stroma coupling",
            "NMF N2 = exocrine purity signal (not prognostic)",
            "Tumor purity dominates variance across cancer types\n  → mismatch is likely the norm",
        ]),
    ]):
        left = Inches(0.35) + i * (col_w + gap)
        add_rect(slide, left, top, col_w, col_h, DIM, bc)
        add_text_box(slide, title,
                     left + Inches(0.15), top + Inches(0.1), col_w - Inches(0.25), Inches(0.35),
                     font_size=17, bold=True, color=bc)
        add_bullet_box(slide, bullets,
                       left + Inches(0.15), top + Inches(0.52), col_w - Inches(0.25), col_h - Inches(0.58),
                       font_size=16, color=LT)

    # Central insight strip
    add_rect(slide, Inches(0.35), Inches(4.25), Inches(12.65), Inches(0.9), TEAL)
    add_text_box(slide,
                 "Aligning the factorization objective with the clinical question—during discovery, not just evaluation—\nyields more transportable prognostic signatures with fewer, biologically interpretable factors.",
                 Inches(0.5), Inches(4.32), Inches(12.4), Inches(0.82),
                 font_size=17, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Single-cell callout (prominent, gold)
    add_sc_callout(slide, Inches(0.35), Inches(5.3), Inches(12.65), Inches(1.85))

    add_text_box(slide, "R package: github.com/ayoung31/DeSurv  |  Paper code: github.com/ayoung31/DeSurv-paper",
                 Inches(0.35), Inches(7.22), Inches(12.65), Inches(0.28),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    slide_num(slide, 12)
    return slide


# =============================================================================
# BACKUP SLIDES
# =============================================================================

def slide_backup_divider(prs):
    slide = blank_slide(prs)
    add_rect(slide, 0, 0, Inches(13.33), Inches(7.5), MGRAY)
    add_text_box(slide, "BACKUP SLIDES",
                 Inches(1.0), Inches(2.5), Inches(11.33), Inches(1.1),
                 font_size=40, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
    add_text_box(slide, "Supplementary material and additional details",
                 Inches(1.0), Inches(3.7), Inches(11.33), Inches(0.55),
                 font_size=22, color=DKGRAY, align=PP_ALIGN.CENTER)
    return slide


def slide_backup_algorithm(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Algorithm Details",
                   "Alternating optimization — converges to stationary point under mild conditions")

    add_bullet_box(slide,
                   ["Initialize W, H using consensus NMF across multiple random starts",
                    "Alternating updates:",
                    "Update H: fix W, β → minimize ℒ_NMF  (ℓ₂ penalty)",
                    "Update W: fix H, β → minimize joint objective  (survival gradient here)",
                    "Update β: fix W, H → maximize ℒ_Cox  (elastic-net penalized Cox)",
                    "Repeat until convergence (relative objective change < ε)",
                    "Multiple initializations → select best by training C-index",
                    "W truncated to BO-selected n_top = 270 genes per factor (PDAC)",
                    "Convergence: stationary point guaranteed (SI Appendix)"],
                   Inches(0.5), Inches(1.42), Inches(12.3), Inches(5.5),
                   font_size=17, color=DKGRAY,
                   indent_levels=[0, 0, 1, 1, 1, 1, 0, 0, 0])


    slide_num(slide, "B1")
    return slide


def slide_backup_null_mixed(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Null and Mixed Simulation Scenarios",
                   "Verifying DeSurv gains reflect genuine signal, not overfitting")

    # Null
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(6.1), Inches(5.85), LGRAY, MGRAY)
    add_text_box(slide, "Null scenario  (β = 0)",
                 Inches(0.5), Inches(1.38), Inches(5.9), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["No true survival signal in any factor",
                    "Both DeSurv and NMF: C-index ≈ 0.50",
                    "DeSurv does NOT impose spurious associations",
                    "α tuned to zero → reduces to standard NMF",
                    "Method adapts to data structure"],
                   Inches(0.5), Inches(1.82), Inches(5.9), Inches(2.5),
                   font_size=17, color=DKGRAY)
    add_rect(slide, Inches(0.5), Inches(4.45), Inches(5.8), Inches(2.45), WHITE, MGRAY)
    add_text_box(slide, "[SI Fig. S2 — null scenario box plots]",
                 Inches(1.5), Inches(5.55), Inches(3.8), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    # Mixed
    add_rect(slide, Inches(6.6), Inches(1.3), Inches(6.35), Inches(5.85), LGRAY, ORANGE)
    add_text_box(slide, "Mixed scenario  (partial overlap)",
                 Inches(6.75), Inches(1.38), Inches(6.1), Inches(0.35),
                 font_size=17, bold=True, color=ORANGE)
    add_bullet_box(slide,
                   ["Survival driven by factor-1 markers + 150 background genes",
                    "Background genes load on all factors\n  → partial overlap with variance-dominant programs",
                    "DeSurv advantage attenuated vs. primary scenario",
                    "Δ C-index smaller but still present"],
                   Inches(6.75), Inches(1.82), Inches(6.1), Inches(2.5),
                   font_size=17, color=DKGRAY)
    add_rect(slide, Inches(6.75), Inches(4.45), Inches(6.1), Inches(2.45), WHITE, ORANGE)
    add_text_box(slide, "[SI Fig. S2 — mixed scenario box plots]",
                 Inches(7.7), Inches(5.55), Inches(4.0), Inches(0.3),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)


    slide_num(slide, "B2")
    return slide


def slide_backup_ksens(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "K-Sensitivity: Adjusted P-values Across K × α Grid",
                   "Stratified Cox in 5 external validation cohorts, adjusted for PurIST and DeCAF classifiers")

    # Simple text-based table description
    add_rect(slide, Inches(0.35), Inches(1.3), Inches(12.6), Inches(5.55), LGRAY, MGRAY)

    add_text_box(slide,
                 "Rows: K ∈ {2, 3, 5, 7, 9}    Columns: α ∈ {0, 0.25, 0.35, 0.55, 0.75, 0.85, 0.95}",
                 Inches(0.5), Inches(1.38), Inches(12.2), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)

    add_bullet_box(slide,
                   ["K = 2: not significant at any α",
                    "K = 3: significant (P < 0.05) at 4 / 7 α values after classifier adjustment",
                    "Unadjusted: 6 / 7",
                    "K = 5: significant at 1–2 / 7 α values",
                    "K = 7: significant at 2 / 7 α values after adjustment  (unadjusted: 6 / 7)",
                    "K = 7 content largely recoverable from PurIST/DeCAF alone",
                    "K = 9: not significant after adjustment",
                    "K = 3 at α = 0.55 (production model): P < 0.05 unadjusted;  independent of classifiers",
                    "Full tables in SI Appendix Tables S1–S3"],
                   Inches(0.5), Inches(1.88), Inches(12.2), Inches(4.7),
                   font_size=17, color=DKGRAY,
                   indent_levels=[0, 0, 1, 0, 0, 1, 0, 0, 0])


    slide_num(slide, "B3")
    return slide


def slide_backup_datasets(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "Dataset Details",
                   "Training and external validation cohorts")

    rows = [
        ("Role", "Cohort", "Platform", "n", "Source"),
        ("Training", "TCGA-PAAD", "RNA-seq", "~150", "GDC / TCGA"),
        ("Training", "CPTAC-3", "RNA-seq", "~140", "GDC / CPTAC"),
        ("Validation", "Puleo et al.", "Microarray", "~309", "ArrayExpress E-MTAB-6134"),
        ("Validation", "Moffitt et al.", "Microarray", "~145", "GEO GSE71729"),
        ("Validation", "PACA-AU (array)", "Microarray", "~90", "ICGC / EGA"),
        ("Validation", "PACA-AU (seq)", "RNA-seq", "~90", "ICGC / EGA"),
        ("Validation", "Dijk et al.", "Microarray", "~64", "ArrayExpress E-MTAB-6830"),
    ]

    col_ws = [Inches(1.6), Inches(2.2), Inches(1.8), Inches(1.0), Inches(5.7)]
    row_h = Inches(0.62)
    left0 = Inches(0.45)
    top0  = Inches(1.38)

    for i, row in enumerate(rows):
        is_hdr = (i == 0)
        bg = NAVY if is_hdr else (LGRAY if i % 2 == 0 else WHITE)
        fc = WHITE if is_hdr else DKGRAY
        x = left0
        for cell, cw in zip(row, col_ws):
            add_rect(slide, x, top0 + i * row_h, cw, row_h, bg, MGRAY)
            add_text_box(slide, cell,
                         x + Inches(0.06), top0 + i * row_h + Inches(0.12),
                         cw - Inches(0.12), row_h - Inches(0.2),
                         font_size=16, bold=is_hdr, color=fc)
            x += cw

    add_text_box(slide,
                 "Preprocessing: log₂(TPM+1)  →  within-sample rank transform (preserves nonneg).   "
                 "Gene selection: top 1,970 genes by mean expression + variance (intersection across training cohorts).",
                 Inches(0.35), Inches(6.7), Inches(12.6), Inches(0.55),
                 font_size=16, italic=True, color=DKGRAY)


    slide_num(slide, "B4")
    return slide


def slide_backup_w_corr(prs):
    slide = blank_slide(prs)
    add_header_bar(slide, "W-Matrix Correspondence Between Methods",
                   "Spearman rank correlation of gene loadings (all genes) between DeSurv and NMF factors")

    add_rect(slide, Inches(0.35), Inches(1.3), Inches(6.6), Inches(5.85), LGRAY, MGRAY)
    add_text_box(slide,
                 "[Fig. 3D — 3×3 heatmap\n Spearman correlation\n of W columns between\n DeSurv (rows) and NMF (cols)]",
                 Inches(1.5), Inches(3.8), Inches(3.8), Inches(1.5),
                 font_size=16, italic=True, color=MGRAY, align=PP_ALIGN.CENTER)

    add_rect(slide, Inches(7.15), Inches(1.3), Inches(5.8), Inches(5.85), LGRAY, NAVY)
    add_text_box(slide, "Interpretation:",
                 Inches(7.3), Inches(1.38), Inches(5.55), Inches(0.35),
                 font_size=17, bold=True, color=NAVY)
    add_bullet_box(slide,
                   ["N1 (classical) ↔ D3 (basal-like): strong positive correlation",
                    "N3 (mixed microenvironment) ↔ D2 (proCAF): moderate correlation",
                    "D1 (most prognostic): weak, uniform correlation with ALL NMF factors",
                    "→ D1's tumor-stroma coupling does not exist in the NMF solution",
                    "N2 (exocrine): no strong DeSurv counterpart\n  → exocrine signal suppressed by DeSurv",
                    "Two methods share some biology but differ on:\n  exocrine suppression  &  tumor-stroma coupling"],
                   Inches(7.3), Inches(1.82), Inches(5.55), Inches(5.1),
                   font_size=16, color=DKGRAY,
                   indent_levels=[0, 0, 0, 1, 0, 0])


    slide_num(slide, "B5")
    return slide


# =============================================================================
# ASSEMBLE
# =============================================================================

def build_presentation(out_path):
    prs = make_prs()
    setup_layout_footer(prs)          # once — inherited by all slides via layout
    slide_title(prs)
    slide_motivation(prs)
    slide_nmf_background(prs)
    slide_desurv_model(prs)
    slide_simulations(prs)
    slide_model_selection(prs)
    slide_pdac_factors(prs)
    slide_variance_survival(prs)
    slide_validation(prs)
    slide_km_curves(prs)
    slide_robustness(prs)
    slide_conclusions(prs)

    slide_backup_divider(prs)
    slide_backup_algorithm(prs)
    slide_backup_null_mixed(prs)
    slide_backup_ksens(prs)
    slide_backup_datasets(prs)
    slide_backup_w_corr(prs)

    prs.save(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    out = r"C:\Users\amyou\Documents\GitHub\DeSurv-paper\DeSurv_presentation.pptx"
    build_presentation(out)
