# 72: datasheet for the score panels from 71 (3D, 4A, 4B, 4C, 4D) plus the NFAT gene-removal
# robustness from 70. One sheet per panel: what is plotted (per-cell values, bars as Excel formulas), the statistics
# (with the alternative gene sets alongside), and the genes in each score.
# Usage (from .): python3 04_cell_level/72_score_panels_workbook.py <out.xlsx>
import csv, sys, math, statistics as stt
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter as CL
R = "output/cell_level/data/score_panels/"; A = "output/cell_level/data/nfat_audit/"
rd = lambda f: list(csv.DictReader(open(f)))
B = Font(bold=True); T = Font(bold=True, size=13); G = Font(color="808080", italic=True)
HDR = PatternFill("solid", fgColor="DDDDDD"); MAIN = PatternFill("solid", fgColor="E2EFDA"); ALT = PatternFill("solid", fgColor="F2F2F2")
GL = {"Young_Ctrl": "3 mo CTRL", "Aged_Ctrl": "24 mo CTRL", "Young_KO": "3 mo P1KO", "Aged_KO": "24 mo P1KO"}
G4 = ["Young_Ctrl", "Aged_Ctrl", "Young_KO", "Aged_KO"]
def ptext(p):
    p = float(p)
    if p >= 0.995: return "p = 1.00"
    if p >= 0.1: return f"p = {p:.2f}"
    if p >= 0.001: return f"p = {p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"p = {m} × 10^{int(e)}"
def hdr(ws, cols, fill=HDR, c0=1):
    r = ws.max_row + 1
    for j, h in enumerate(cols): x = ws.cell(r, c0+j, h); x.font = B; x.fill = fill; x.alignment = Alignment(wrap_text=True, vertical="top")
    return r
def lines(ws, title, ls):
    ws.append([title]); ws.cell(ws.max_row, 1).font = T
    for l in ls: ws.append([l])
    ws.append([])
def widths(ws, w):
    for j, x in enumerate(w, 1): ws.column_dimensions[CL(j)].width = x
wb = Workbook(); wb.remove(wb.active)

# ---- correlation panels
cs = rd(R+"correlations_paper_vs_alternative.csv"); pc = rd(R+"correlation_panels_per_cell.csv")
INFO = {"3D": ("Piezo1 vs synthetic phenotype score", "x = synthetic phenotype score, y = Piezo1 expression"),
        "4B": ("Piezo1 vs NFAT pathway score", "x = Piezo1 expression, y = NFAT pathway score"),
        "4C": ("NFAT pathway score vs synthetic phenotype score", "x = synthetic phenotype score, y = NFAT pathway score")}
for pn, (t, ax) in INFO.items():
    ws = wb.create_sheet(f"Fig {pn}")
    lines(ws, f"Fig {pn}: {t}", [f"Plot: {ax}; one linear fit with 95% band per age (green 3 mo, red 24 mo); Spearman rho per age. CTRL SMCs only.",
                                               "Per-cell values from column J (x and y, by age). Rows labelled 'alternative' use the alternative gene sets, for comparison."])
    r = hdr(ws, ["Version", "Age", "SMCs (n)", "Spearman rho", "p value", "p, as written", "Star", "Genes"], MAIN)
    for s in [s for s in cs if s["panel"] == pn]:
        r += 1
        for j, v in enumerate([s["version"], s["age"], int(s["n"]), float(s["rho"]), float(s["p"]), ptext(s["p"]), s["star"], s["genes"]], 1):
            c = ws.cell(r, j, v); c.fill = MAIN if s["version"] == "paper" else ALT
        ws.cell(r, 4).number_format = "0.000"; ws.cell(r, 5).number_format = "0.0E+00" if float(s["p"]) < 1e-3 else "0.000"
    ws.append([]); ws.append(["Star levels: * p < 0.05, ** < 0.01, *** < 0.001 (as in the Fig 3 and 4 correlation panels)."]); ws.cell(ws.max_row, 1).font = G
    hr = 4; C0 = 10
    for j, h in enumerate(["3 mo: x", "3 mo: y", "24 mo: x", "24 mo: y", "3 mo: mouse", "24 mo: mouse"]):
        c = ws.cell(hr, C0+j, h); c.font = B; c.fill = HDR
    for j, a in enumerate(["Young", "Aged"]):
        rows = [p for p in pc if p["panel"] == pn and p["age"] == a]
        for i, p in enumerate(rows, hr+1):
            ws.cell(i, C0+2*j, float(p["x_paper"])); ws.cell(i, C0+2*j+1, float(p["y_paper"])); ws.cell(i, C0+4+j, p["sample"]).font = G
    widths(ws, [13, 7, 9, 10, 10, 16, 6, 60, 3, 11, 11, 11, 11, 12, 12])

# ---- 4A
a4 = rd(R+"fig4A_paper_vs_alternative.csv")
ws = wb.create_sheet("Fig 4A")
lines(ws, "Fig 4A: Spearman rho of Piezo1 with each pathway score, CTRL SMCs (469 young, 477 aged)",
      ["Plot: one dot per pathway and age (green 3 mo, red 24 mo); star from the p value. Rows marked 'yes' differ from the alternative gene sets:",
       "NFAT Pathway (no Nfat5), NFAT-Associated Remodeling (no Fas), Synthetic Phenotype (SI Methods set)."])
r = hdr(ws, ["Pathway", "Differs from alternative?", "3 mo rho", "3 mo p", "3 mo star", "24 mo rho", "24 mo p", "24 mo star", "Alternative set: 3 mo rho (star)", "Alternative set: 24 mo rho (star)", "Genes"], MAIN)
sub = {}
for s in a4:
    if s["version"] == "alternative": sub[(s["pathway"].replace("Responsive", "Associated"), s["age"])] = s
for nm in dict.fromkeys(s["pathway"] for s in a4 if s["version"] == "paper"):
    y = [s for s in a4 if s["version"] == "paper" and s["pathway"] == nm and s["age"] == "Young"][0]
    g = [s for s in a4 if s["version"] == "paper" and s["pathway"] == nm and s["age"] == "Aged"][0]
    sy, sg = sub[(nm, "Young")], sub[(nm, "Aged")]; r += 1
    vals = [nm, "yes" if y["differs_from_alternative"] == "TRUE" else "", float(y["rho"]), float(y["p"]), y["star"], float(g["rho"]), float(g["p"]), g["star"],
            f'{float(sy["rho"]):+.3f} ({sy["star"]})', f'{float(sg["rho"]):+.3f} ({sg["star"]})', y["genes"]]
    for j, v in enumerate(vals, 1): ws.cell(r, j, v)
    for c in (3, 6): ws.cell(r, c).number_format = "0.000"
    for c in (4, 7): ws.cell(r, c).number_format = "0.0E+00"
widths(ws, [26, 9, 9, 10, 7, 9, 10, 7, 16, 16, 70])

# ---- 4D superplot
D = R+"fig4D_superplot_inputs/"
cells = rd(D+"per_cell_values.csv"); pa = rd(D+"per_animal_values.csv"); st = rd(R+"fig4D_stats_paper_vs_alternative.csv"); bars0 = rd(R+"fig4D_bars_paper_vs_alternative.csv")
ws = wb.create_sheet("Fig 4D")
lines(ws, "Fig 4D: NFAT-associated remodeling score per SMC, 5 genes (Wnt5a, Mmp2, Mmp9, Vcam1, Timp1; Fas removed)",
      ["Plot: bar = mean across SMCs; error bar = SEM; small dots = SMCs (columns N–Q); large symbols = per-mouse means. Group order 3 mo CTRL, 24 mo CTRL, 3 mo P1KO, 24 mo P1KO."])
ws.append([]); ws.append(["BARS (green = what to plot)"]); ws.cell(ws.max_row, 1).font = B
r0 = hdr(ws, ["Group", "SMCs (n)", "BAR: mean", "ERROR BAR: SEM", "Mouse 1", "Mouse 1 mean", "Mouse 1 n", "Mouse 2", "Mouse 2 mean", "Mouse 2 n", "6 genes (with Fas): mean", "6 genes: SEM"], MAIN)
hr = 4; C0 = 14; r = r0
for j, g in enumerate(G4):
    v = [(c["sample"], float(c["value"])) for c in cells if c["group"] == g]; v.sort()
    ws.cell(hr, C0+j, f"{GL[g]}: per-cell values").font = B; ws.cell(hr, C0+j).fill = HDR
    for i, (_, val) in enumerate(v, hr+1): ws.cell(i, C0+j, val)
    col = CL(C0+j); rng = f"{col}{hr+1}:{col}{hr+len(v)}"; r += 1
    row = [GL[g], f"=COUNT({rng})", f"=AVERAGE({rng})", f"=STDEV({rng})/SQRT(COUNT({rng}))"]
    ms = sorted([a for a in pa if a["group"] == g], key=lambda a: a["sample"])
    for a in ms: row += [a["sample"], float(a["animal_mean"]), int(a["n_cells"])]
    row += [None]*(10-len(row)); b0 = [b for b in bars0 if b["group"] == g][0]; row += [float(b0["mean_6genes"]), float(b0["sem_6genes"])]
    for k, val in enumerate(row, 1): ws.cell(r, k, val)
    for c in (3, 4, 6, 9, 11, 12): ws.cell(r, c).number_format = "0.0000"
    assert abs(stt.mean([x for _, x in v]) - float(b0["mean_5genes"])) < 1e-9
ws.append([]); ws.append(["BRACKETS: two-sided Mann–Whitney U across SMCs, Bonferroni ×4 (one p per comparison)"]); ws.cell(ws.max_row, 1).font = B
r = hdr(ws, ["Comparison", "SMCs (n)", "p value", "p, as written", "Star", "6 genes (with Fas): p", "6 genes: star"], MAIN)
for s in st:
    r += 1
    for k, val in enumerate([f'{GL[s["group1"]]} vs {GL[s["group2"]]}', f'{s["n1"]} vs {s["n2"]}', float(s["p_bonf"]), ptext(s["p_bonf"]), s["sig"],
                             float(s["p_bonf_6genes"]), s["sig_6genes"]], 1): ws.cell(r, k, val)
    for c in (3, 6): ws.cell(r, c).number_format = "0.0E+00"
widths(ws, [16, 9, 11, 11, 12, 11, 8, 12, 11, 8, 13, 11, 3, 14, 14, 14, 14])

# ---- robustness (SI table candidate)
rb = rd(A+"nfat_activity_gene_removal_tests.csv")
ws = wb.create_sheet("NFAT score robustness")
lines(ws, "NFAT-associated remodeling score: robustness of the 4D result to removing genes",
      ["Fold change of group means and Bonferroni-adjusted (×4) Mann–Whitney p across SMCs, the same test as Fig 4D."])
r = hdr(ws, ["Genes in score", "Comparison", "Fold change", "p (Bonferroni ×4)", "p, as written", "Star"], MAIN)
for s in rb:
    r += 1
    c1, c2 = s["comparison"].split(" vs ")
    for k, val in enumerate([s["set"], f"{GL[c1]} vs {GL[c2]}", float(s["fold"]), float(s["p_bonf"]), ptext(s["p_bonf"]), s["sig"]], 1): ws.cell(r, k, val)
    ws.cell(r, 3).number_format = "0.00"; ws.cell(r, 4).number_format = "0.0E+00"
widths(ws, [44, 26, 10, 14, 18, 7])

ws = wb.create_sheet("How to read", 0)
for t, b in [("Score panels Fig 3D, 4A, 4B, 4C, 4D", True), ("", False),
             ("One sheet per panel. Green rows/headers = the plotted values. Grey = the alternative gene sets, for comparison.", False), ("", False),
             ("Gene sets", True),
             ("  Synthetic phenotype score (3D, 4A, 4C): Vim, Col1a2, Dcn, Igfbp6, Efemp1 (SI Methods; also Fig 3F).", False),
             ("     Alternative set: Klf4, Vim, Mmp2, Mmp9, Timp1 (shares 3 genes with the NFAT score in 4D).", False),
             ("  NFAT pathway score (4A, 4B, 4C): 12 genes; Nfat5 (tonicity-regulated, calcineurin-independent) is not included. Alternative: 13 genes with Nfat5.", False),
             ("  NFAT-associated remodeling score (4A, 4D): Wnt5a, Mmp2, Mmp9, Vcam1, Timp1. Alternative: the same plus Fas (Fasl is detected in almost no SMCs).", False),
             ("", False),
             ("Scores: per-cell mean of log-normalized expression (zeros included). Correlations: CTRL SMCs (469 young, 477 aged). 4D: all visceral SMCs.", False),
             ("Built by 04_cell_level/72_score_panels_workbook.py from the CSVs written by 70 and 71.", False)]:
    ws.append([t]); ws.cell(ws.max_row, 1).font = B if b else Font()
ws.column_dimensions["A"].width = 150
wb.move_sheet("Fig 4A", offset=-2)
wb.save(sys.argv[1]); print("saved", sys.argv[1], wb.sheetnames)
