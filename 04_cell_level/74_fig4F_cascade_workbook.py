# 74: datasheet for Fig 4F (pseudotime cascade) from 73. Sheets: the plotted LOESS curves (with the 6-gene
# NFAT curve alongside), the 205 per-cell pseudotime values and module scores, the cascade order (where each curve
# crosses Z = 0), the with/without-Fas comparison, and the genes.
# Usage (from .): python3 04_cell_level/74_fig4F_cascade_workbook.py <out.xlsx>
import csv, sys
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter as CL
R = "output/cell_level/data/score_panels/fig4F/"
rd = lambda f: list(csv.DictReader(open(f)))
B = Font(bold=True); T = Font(bold=True, size=13); G = Font(color="808080", italic=True)
HDR = PatternFill("solid", fgColor="DDDDDD"); MAIN = PatternFill("solid", fgColor="E2EFDA"); ALT = PatternFill("solid", fgColor="F2F2F2")
COL = {"Contractile": "2CA02C", "Cyto Remodeling": "FF7F00", "Piezo1": "1F77B4", "NFAT-assoc. remodeling": "9467BD", "NFAT Effectors": "9467BD", "Synthetic": "E41A1C"}
def hdr(ws, cols, fill=HDR, c0=1, r=None):
    r = r or ws.max_row + 1
    for j, h in enumerate(cols): x = ws.cell(r, c0+j, h); x.font = B; x.fill = fill; x.alignment = Alignment(wrap_text=True, vertical="top")
    return r
def lines(ws, title, ls):
    ws.append([title]); ws.cell(ws.max_row, 1).font = T
    for l in ls: ws.append([l])
    ws.append([])
def widths(ws, w):
    for j, x in enumerate(w, 1): ws.column_dimensions[CL(j)].width = x
def ptext(p):
    p = float(p)
    if p >= 0.1: return f"p = {p:.2f}"
    if p >= 0.001: return f"p = {p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"p = {m} × 10^{int(e)}"

order = rd(R+"fig4F_cascade_order_paper_vs_alternative.csv"); lo = rd(R+"fig4F_loess_grid.csv")
pc = rd(R+"fig4F_per_cell_pseudotime_and_modules.csv"); chk = rd(R+"fig4F_nfat_module_with_without_Fas.csv")
MODS = ["Contractile", "Cyto Remodeling", "Piezo1", "NFAT-assoc. remodeling", "Synthetic"]
wb = Workbook(); wb.remove(wb.active)

ws = wb.create_sheet("How to read")
lines(ws, "Fig 4F: gene-module dynamics along pseudotime in Piezo1+ CTRL SMCs (n = 205)", [
    "The purple module is the NFAT-associated remodeling score (Mmp2, Mmp9, Vcam1, Timp1, Wnt5a). For comparison, the 6-gene module with Fas ('NFAT Effectors')",
    "is computed on the same cells and pseudotime.",
    "The pseudotime values are those exported by pipeline module Fig4I (Fig4I_data.xlsx).",
    "",
    "Sheet 'Plotted curves': the 200-point LOESS curves that are drawn (Z, lower and upper band = fit ± 1 SE, clipped at ±2 / ±2.5 as in the panel), one block per",
    "    module, plus the local aged fraction that sets line opacity (0.25 + 0.75 × aged fraction). The 6-gene NFAT curve (with Fas) is the last block, grey.",
    "Sheet 'Per cell': the 205 cells with pseudotime (blended, Slingshot, kNN aged fraction), age, mouse, cluster, and each module's mean log-normalized",
    "    expression and Z-score (Z across the 205 cells; Piezo1 Z is the residual of Piezo1 on log10 UMI, as in the pipeline).",
    "Sheet 'Cascade order': where each curve first crosses Z = 0 (down for Contractile, up for the rest), where it peaks, and its Spearman rho with pseudotime,",
    "    for both NFAT module versions.",
    "Sheet 'NFAT with-without Fas': how close the two purple curves are.",
    "Sheet 'Genes': the genes in each module.",
    "",
    "How the curves are made: module score = mean log-normalized expression of the genes (zeros included) → Z across the 205 cells → LOESS on blended pseudotime",
    "(span 0.45) → predicted on a 200-point grid trimmed 2% at each end. Blended pseudotime = rank of (0.6 × Slingshot pseudotime + 0.4 × kNN aged fraction, k = 15).",
    "Script: 04_cell_level/73_fig4F_cascade.R."])
widths(ws, [160])

# ---- plotted curves
ws = wb.create_sheet("Plotted curves")
lines(ws, "LOESS curves as drawn. Columns: pseudotime, aged fraction, then Z / lower / upper per module.", [
    "Line colour per module: Contractile #2CA02C, Cyto Remodeling #FF7F00, Piezo1 #1F77B4, NFAT-assoc. remodeling #9467BD, Synthetic #E41A1C. Opacity = 0.25 + 0.75 × aged fraction."])
cols = ["Pseudotime (grid)", "Aged fraction (local)"]
for m in MODS: cols += [f"{m}: Z", f"{m}: lower", f"{m}: upper"]
cols += ["NFAT Effectors, 6 genes with Fas: Z", "6 genes: lower", "6 genes: upper"]
r0 = hdr(ws, cols)
for j, m in enumerate(MODS):
    for k in range(3): ws.cell(r0, 3+3*j+k).fill = PatternFill("solid", fgColor=COL[m]); ws.cell(r0, 3+3*j+k).font = Font(bold=True, color="FFFFFF")
for k in range(3): ws.cell(r0, 3+3*len(MODS)+k).fill = ALT
rec = {m: [d for d in lo if d["version"] == "paper" and d["module"] == m] for m in MODS}
sub6 = [d for d in lo if d["version"] == "alternative" and d["module"] == "NFAT Effectors"]
n = len(rec["Contractile"]); assert n == 200 and all(len(v) == n for v in rec.values()) and len(sub6) == n
for i in range(n):
    row = [float(rec["Contractile"][i]["pt"]), float(rec["Contractile"][i]["age_frac"])]
    for m in MODS: row += [float(rec[m][i][k]) for k in ("z", "z_lo", "z_hi")]
    row += [float(sub6[i][k]) for k in ("z", "z_lo", "z_hi")]
    ws.append(row)
for rr in range(r0+1, r0+1+n):
    for c in range(1, len(cols)+1): ws.cell(rr, c).number_format = "0.000"
widths(ws, [12, 12] + [11]*(3*len(MODS)) + [14, 11, 11]); ws.freeze_panes = ws.cell(r0+1, 3)

# ---- per cell
ws = wb.create_sheet("Per cell")
lines(ws, "The 205 Piezo1+ CTRL SMCs: pseudotime and module scores (5-gene NFAT module in green; 6-gene NFAT module in grey).", [
    "Z = (score − mean) / SD across these 205 cells. Piezo1 Z = standardized residual of Piezo1 log-expression on log10(total UMI)."])
keys = list(pc[0].keys()); r0 = hdr(ws, keys)
for j, k in enumerate(keys, 1):
    if k.startswith("NFAT-assoc"): ws.cell(r0, j).fill = MAIN
    if k.startswith("NFAT Effectors"): ws.cell(r0, j).fill = ALT
for d in pc:
    ws.append([(float(d[k]) if k not in ("cell", "sample", "age", "cluster") else d[k]) for k in keys])
for rr in range(r0+1, r0+1+len(pc)):
    for c in range(4, len(keys)+1): ws.cell(rr, c).number_format = "0.000"
widths(ws, [36, 12, 7] + [12]*(len(keys)-3)); ws.freeze_panes = ws.cell(r0+1, 2)

# ---- cascade order
ws = wb.create_sheet("Cascade order")
lines(ws, "Where each curve crosses Z = 0 and peaks (pseudotime 0–1), and its rank correlation with pseudotime", [
    "Crossing = first grid point where the curve goes below 0 (Contractile) or above 0 (others), interpolated. Peak = position of the extreme Z.",
    "The order of the crossings is the same in both versions: Contractile falls first; Cyto Remodeling and Synthetic cross next, together; the NFAT module and Piezo1 cross later and rise to the end."])
r0 = hdr(ws, ["Version", "Module", "Crosses Z = 0 at", "Peak at", "Z at start", "Z at end", "Spearman rho vs pseudotime", "p", "p, as written", "Genes"])
for d in order:
    ws.append([d["version"], d["module"], float(d["pt_cross_zero"]), float(d["pt_peak"]), float(d["z_start"]), float(d["z_end"]),
               float(d["spearman_rho_vs_pseudotime"]), float(d["p"]), ptext(d["p"]), d["genes"]])
    rr = ws.max_row
    for c in range(1, 11): ws.cell(rr, c).fill = MAIN if d["version"] == "paper" else ALT
    for c in (3, 4, 5, 6, 7): ws.cell(rr, c).number_format = "0.000"
    ws.cell(rr, 8).number_format = "0.0E+00"
widths(ws, [14, 22, 14, 10, 10, 10, 14, 10, 16, 60])

# ---- with/without Fas
ws = wb.create_sheet("NFAT with-without Fas")
lines(ws, "How different is the NFAT module curve without Fas?", ["Both computed on the same 205 cells and pseudotime."])
r0 = hdr(ws, ["Metric", "Value"])
for d in chk: ws.append([d["metric"], float(d["value"])]); ws.cell(ws.max_row, 2).number_format = "0.000"
widths(ws, [48, 12])

# ---- genes
ws = wb.create_sheet("Genes")
lines(ws, "Genes in each Fig 4F module", ["The Fig 4F 'Synthetic' module is a cascade module (transcription-factor / dedifferentiation markers) and is NOT the synthetic phenotype score of Fig 3D/3F/4C."])
r0 = hdr(ws, ["Version", "Module", "Genes", "n genes"])
for d in order: ws.append([d["version"], d["module"], d["genes"], len(d["genes"].split(", "))]);
for rr in range(r0+1, ws.max_row+1):
    for c in range(1, 5): ws.cell(rr, c).fill = MAIN if ws.cell(rr, 1).value == "paper" else ALT
widths(ws, [14, 22, 60, 8])

wb.save(sys.argv[1]); print("saved", sys.argv[1])
