# PAPER PANEL: Fig S3B (Excel workbook behind the panel)
# 64: one shareable Excel workbook behind Fig S3B ("only Piezo1 rises
# with age"). Raw per-cell UMI for every expressing SMC, mean/SD/SEM as Excel formulas over those
# cells, a simple unmatched log2FC with a delta-method interval the reader can recompute by hand,
# and the 200 bootstrap draws from 46 with the plotted median/interval as formulas over them.
# Usage: python3 64_S3B_expressing_cells_workbook.py <out.xlsx>
import csv, sys, math, statistics as st
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
D = "output/cell_level/channels_aging/data/"
out = sys.argv[1]
res = list(csv.DictReader(open(D+"conditional_positive_channels.csv")))
CH = [r["gene"] for r in res]; R = {r["gene"]: r for r in res}
cells = list(csv.DictReader(open(D+"expressing_cells_counts_long.csv")))
draws = list(csv.DictReader(open(D+"conditional_positive_draws.csv")))
GR = [("Young_Ctrl", 3), ("Aged_Ctrl", 24)]
bold = Font(bold=True); hdr = PatternFill("solid", fgColor="DDDDDD")
def header(ws, cols):
    ws.append(cols)
    for c in ws[1]: c.font = bold; c.fill = hdr; c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"
def widths(ws, w):
    for i, x in enumerate(w, 1): ws.column_dimensions[get_column_letter(i)].width = x
wb = Workbook()

# Raw cells, sorted so each channel x group is one contiguous block
cells.sort(key=lambda c: (CH.index(c["gene"]), c["group"] != "Young_Ctrl", c["animal"], c["cell_barcode"]))
raw = wb.active; raw.title = "Raw cells"
header(raw, ["Channel","Group","Age (months)","Mouse","Cell barcode","UMI of channel (plotted, S3B)","Total UMI of cell"])
block = {}
for c in cells:
    r = raw.max_row + 1; key = (c["gene"], c["group"])
    block.setdefault(key, [r, r]); block[key][1] = r
    raw.append([c["gene"], c["group"], 3 if c["group"]=="Young_Ctrl" else 24, c["animal"], c["cell_barcode"], int(c["counts"]), int(c["total_UMI"])])
widths(raw, [10,12,8,13,40,14,12])

# Bootstrap draws, one contiguous block per channel
bs = wb.create_sheet("Bootstrap draws")
header(bs, ["Channel","Draw","Young expressing cells in matched set (n)","Aged expressing cells in matched set (n)","Mean UMI, young","Mean UMI, aged","log2FC (aged / young)"])
dblock = {}
for d in draws:
    r = bs.max_row + 1; dblock.setdefault(d["gene"], [r, r]); dblock[d["gene"]][1] = r
    bs.append([d["gene"], int(d["draw"]), int(d["n_young_expressing"]), int(d["n_aged_expressing"]),
               float(d["mean_young"]), float(d["mean_aged"]), float(d["log2FC"]) if d["log2FC"] not in ("", "NA") else None])
widths(bs, [10,7,14,14,11,11,12])

# Summary: per group mean +/- SEM (formulas over Raw cells)
sm = wb.create_sheet("Mean and SEM", 0)
header(sm, ["Channel","Age (months)","Expressing SMCs (n)","Mean UMI per expressing SMC (line in S3B)","SD","SEM = SD/√n","Mice"])
srow = {}; check = []
for g in CH:
    for grp, age in GR:
        a, b = block[(g, grp)]; rng = f"'Raw cells'!F{a}:F{b}"; r = sm.max_row + 1; srow[(g, grp)] = r
        mice = sorted({c["animal"] for c in cells if c["gene"]==g and c["group"]==grp})
        sm.append([g, age, f"=COUNT({rng})", f"=AVERAGE({rng})", f"=STDEV({rng})", f"=E{r}/SQRT(C{r})", ", ".join(mice)])
        for col in "DEF": sm[f"{col}{r}"].number_format = "0.000"
        v = [int(c["counts"]) for c in cells if c["gene"]==g and c["group"]==grp]
        check.append((g, grp, len(v), st.mean(v), st.stdev(v)/math.sqrt(len(v))))
widths(sm, [10,8,11,16,9,11,26])

# Fold change: unmatched (hand-checkable) vs the plotted matched bootstrap
fc = wb.create_sheet("Fold change and tests", 1)
header(fc, ["Channel","Simple log2FC, unmatched = log2(mean aged / mean young)","SE of log2FC (delta method)",
            "Simple 95% CI, lower (±1.96 SE)","Simple 95% CI, upper",
            "PLOTTED in S3B: matched bootstrap median log2FC","PLOTTED: 95% interval, lower (2.5th percentile of draws)","PLOTTED: 95% interval, upper (97.5th percentile)",
            "Background genes: median log2FC","Background 95% range, lower","Background 95% range, upper","Background genes (n)","Empirical p vs background","Star in S3B"])
for g in CH:
    y, a = srow[(g,"Young_Ctrl")], srow[(g,"Aged_Ctrl")]; r = fc.max_row + 1; d0, d1 = dblock[g]
    dr = f"'Bootstrap draws'!G{d0}:G{d1}"; s = R[g]; p = float(s["p_emp"])
    fc.append([g, f"=LOG('Mean and SEM'!D{a}/'Mean and SEM'!D{y},2)",
               f"=SQRT(('Mean and SEM'!F{a}/'Mean and SEM'!D{a})^2+('Mean and SEM'!F{y}/'Mean and SEM'!D{y})^2)/LN(2)",
               f"=B{r}-1.96*C{r}", f"=B{r}+1.96*C{r}",
               f"=MEDIAN({dr})", f"=PERCENTILE({dr},0.025)", f"=PERCENTILE({dr},0.975)",
               float(s["bg_median"]), float(s["bg_lo"]), float(s["bg_hi"]), int(s["n_bg"]), p,
               ("*" if p < 0.05 else "")])
    for col in "BCDEFGHIJK": fc[f"{col}{r}"].number_format = "0.000"
    fc[f"M{r}"].number_format = "0.000"
widths(fc, [10,16,12,11,11,15,15,15,13,12,12,10,11,8])

# Per mouse
pm = wb.create_sheet("Per mouse", 2)
header(pm, ["Channel","Age (months)","Mouse","Expressing SMCs (n)","Mean UMI per expressing SMC","SEM"])
for g in CH:
    for grp, age in GR:
        a, b = block[(g, grp)]
        for m in sorted({c["animal"] for c in cells if c["gene"]==g and c["group"]==grp}):
            r = pm.max_row + 1; rg = f"'Raw cells'!F{a}:F{b}"; mr = f"'Raw cells'!D{a}:D{b}"
            pm.append([g, age, m, f'=COUNTIF({mr},"{m}")', f'=AVERAGEIF({mr},"{m}",{rg})',
                       f'=SQRT((SUMPRODUCT(({mr}="{m}")*{rg}^2)-D{r}*E{r}^2)/(D{r}-1))/SQRT(D{r})'])
            pm[f"E{r}"].number_format = pm[f"F{r}"].number_format = "0.000"
widths(pm, [10,8,13,11,14,9])

# Notes
nt = wb.create_sheet("Notes", 0)
notes = [
 ("Figure S3B: UMI per expressing SMC, young (3 mo) vs aged (24 mo) CTRL", True),
 ("Data: every CTRL SMC with ≥1 UMI of the channel (sheet 'Raw cells'). The violins in S3B plot these cells; the line is the mean.", False),
 ("", False),
 ("Sheet 'Mean and SEM': mean ± SEM across expressing cells, computed by formula from 'Raw cells'. SEM = SD/√n.", False),
 ("Sheet 'Fold change and tests', columns B–E: a simple log2 fold change of the two means, with SE by the delta method,", False),
 ("   SE = √[(SEM_aged/mean_aged)² + (SEM_young/mean_young)²] / ln 2, and a 95% CI of ±1.96 SE. This can be checked by hand.", False),
 ("Columns F–H: the values PLOTTED in Figure S3B. In each of 200 bootstrap draws, young and aged SMCs were resampled with", False),
 ("   replacement and each young cell matched to an aged cell with total UMI within 10% (aged cells are slightly deeper);", False),
 ("   log2FC was computed among expressing cells of the matched set. Bar = median of draws; error bar = 2.5th–97.5th percentile.", False),
 ("   The draws are in sheet 'Bootstrap draws', and F–H are formulas over them.", False),
 ("Columns I–M: each channel's median compared with 120 genes of matched detection rate (±2 percentage points) processed", False),
 ("   identically. Empirical p = 2 × fraction of background genes at least as extreme. * in S3B = p < 0.05.", False),
 ("", False),
 ("Why the S3B interval is wide (Piezo1 +0.08 to +0.86):", True),
 ("   Piezo1 is sparse. ~100 expressing cells per age; 73% of young and 67% of aged expressing cells have exactly 1 UMI.", False),
 ("   The aged mean is pulled up by a tail of a few cells (7, 10, 18 UMI), so it is uncertain: aged SEM 0.23 vs young 0.08.", False),
 ("   The simple delta-method CI (columns D–E) is about as wide, so the width is a property of the data, not of the method.", False),
 ("", False),
 ("Unit of analysis: cells. With 2 mice per age, the statistics describe the sampled cells; per-mouse means are in 'Per mouse'.", False),
 ("Source: scripts 46 (bootstrap, background) and 48 (per-cell export), 04_cell_level/; workbook built by 64.", False),
]
for t, b in notes:
    nt.append([t]); nt.cell(nt.max_row, 1).font = Font(bold=b)
nt.column_dimensions["A"].width = 130
wb.save(out)
for row in check: print("%-8s %-11s n=%3d mean=%.3f sem=%.3f" % row)
