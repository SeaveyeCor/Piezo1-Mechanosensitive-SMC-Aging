# PAPER PANEL: Fig S3E (Excel workbook behind the panel)
# 53c: one shareable Excel workbook that reproduces the "Channel Detection in
# Contractile SMCs" plot (Fig S3E right). Counts come straight from the per-cell export; % and SEM
# are Excel formulas so the reader can see how each bar is computed. Brackets from 52/52b.
# Usage: python3 53c_channel_detection_workbook.py <out.xlsx> [--with-piezo1 | --s3d]
# --s3d: Fig S3D instead, Piezo1 in all SMCs, tests from 62 (same method as S3E).
import csv, sys, math
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
D = "output/cell_level/data/"
out = sys.argv[1]; withP1 = "--with-piezo1" in sys.argv; s3d = "--s3d" in sys.argv
CH = ["Piezo1"] if s3d else (["Piezo1"] if withP1 else []) + ["Piezo2","Trpc1","Trpc3","Tmem63a","Tmem63b","Pkd1","Trpm7","Pkd2"]
GR = [("Young_Ctrl",3,"CTRL"),("Young_KO",3,"P1KO"),("Aged_Ctrl",24,"CTRL"),("Aged_KO",24,"P1KO")]
cells = [r for r in csv.DictReader(open(D+"per_cell_channel_counts_SMC.csv")) if s3d or r["compartment"].startswith("contractile")]
stats = {}
if s3d:
    for r in csv.DictReader(open(D+"piezo1_pct_S3D_tests.csv")): stats[(r["gene"], 3 if r["age"]=="Young" else 24)] = r
else:
    for f, age in [("channels_contractile_compartment_younggeno.csv",3),("channels_contractile_compartment_agedgeno.csv",24)]:
        for r in csv.DictReader(open(D+f)): stats[(r["gene"],age)] = r
bold = Font(bold=True); hdr = PatternFill("solid", fgColor="DDDDDD")
def header(ws, cols):
    ws.append(cols)
    for c in ws[1]: c.font = bold; c.fill = hdr; c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "A2"
def widths(ws, w):
    for i, x in enumerate(w, 1): ws.column_dimensions[get_column_letter(i)].width = x
wb = Workbook()
ws = wb.active; ws.title = "Plot data"
header(ws, ["Channel","Age (months)","Genotype","SMCs (n)","SMCs with transcript detected (n)","% SMCs detecting (bar)","SEM (error bar, %)"])
check = []
for g in CH:
    for grp, age, geno in GR:
        sub = [c for c in cells if c["group"]==grp]; k = sum(int(c[g])>0 for c in sub); n = len(sub)
        r = ws.max_row + 1
        ws.append([g, age, geno, n, k, f"=E{r}/D{r}*100", f"=SQRT((E{r}/D{r})*(1-E{r}/D{r})/D{r})*100"])
        ws[f"F{r}"].number_format = ws[f"G{r}"].number_format = "0.00"
        check.append((g, grp, 100*k/n, 100*math.sqrt((k/n)*(1-k/n)/n)))
widths(ws, [11,8,10,9,16,14,12])
ws2 = wb.create_sheet("Per mouse")
header(ws2, ["Channel","Age (months)","Genotype","Mouse","SMCs (n)","SMCs with transcript detected (n)","% SMCs detecting"])
for g in CH:
    for grp, age, geno in GR:
        for m in sorted({c["animal"] for c in cells if c["group"]==grp}):
            sub = [c for c in cells if c["animal"]==m]; k = sum(int(c[g])>0 for c in sub); r = ws2.max_row + 1
            ws2.append([g, age, geno, m, len(sub), k, f"=F{r}/E{r}*100"]); ws2[f"G{r}"].number_format = "0.00"
widths(ws2, [11,8,10,13,9,16,12])
ws3 = wb.create_sheet("Statistics")
# ONE p per comparison, and the star comes from it. S3E: logistic p Bonferroni-corrected
# for the channels tested at each age; S3D: one gene, so the logistic p as is.
NADJ = 1 if s3d else len(CH)
header(ws3, ["Channel","Age (months)","Comparison","Matched difference, P1KO − CTRL (percentage points)","95% interval, lower","95% interval, upper",
             "p value (logistic, adjusted for log genes detected" + (")" if s3d else f"; Bonferroni ×{NADJ} channels)"),"Bracket label"]
             + ([] if s3d else ["(Unadjusted p, reference only)"]))
for g in CH:
    for age in (3,24):
        s = stats[(g,age)]; lo, hi, p0 = float(s["lo"]), float(s["hi"]), float(s["logit_p"]); p = min(1.0, p0*NADJ)
        ws3.append([g, age, f"P1KO vs CTRL, {age} mo", round(float(s["matched_diff_pp"]),2), round(lo,2), round(hi,2), p,
                    "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "ns"] + ([] if s3d else [p0]))
        if s3d: assert ws3.cell(ws3.max_row, 8).value == s["label"]
for row in ws3.iter_rows(min_row=2, min_col=7, max_col=9):
    for c in row:
        if isinstance(c.value, float): c.number_format = "0.0E+00" if c.value < 0.001 else "0.000"
        if c.column == 9: c.font = Font(color="808080")
widths(ws3, [11,8,20,20,11,11,18,10,14])
ws4 = wb.create_sheet("Notes")
notes = [
 ("Figure S3D: % of SMCs with Piezo1 detected" if s3d else "Channel Detection in Contractile SMCs" + (" (version including Piezo1, the deleted gene)" if withP1 else " (Figure S3E, right)")),
 "",
 (f"Cells: all bona fide SMCs ({len(cells):,}). Mice:" if s3d else f"Cells: contractile SMCs, i.e. bona fide SMCs outside the synthetic subcluster (Cl5) ({len(cells):,}). Mice:") + ", ".join(f" {age} mo {geno} n = {len({c['animal'] for c in cells if c['group'] == grp})}" for grp, age, geno in GR) + ".",
 "Detected: at least 1 UMI for that channel in the cell.",
 "Bar: % of SMCs detecting the channel, all cells of the group pooled (Plot data, column F).",
 "Error bar: SEM across SMCs, sqrt(p(1 − p)/n) x 100 (Plot data, column G).",
 "Per mouse: the same % for each mouse separately.",
 "Brackets (Statistics): P1KO vs CTRL within each age. (1) Logistic regression, detected ~ genotype + log(genes detected per cell). (2) Bootstrap: cells of both genotypes resampled with replacement and matched 1:1 on genes detected per cell (within 10%); median and 95% interval of 300 draws. One p value per comparison: the logistic p" + ("" if s3d else ", Bonferroni-corrected for the 8 channels tested at each age (p × 8, max 1)") + "; stars * <0.05, ** <0.01, *** <0.001. The bootstrap interval is shown as the effect size and does not decide the star. Same test in Figure S3D and Figure S3E.",
 "Statistics treat cells as independent observations. With 1-2 mice per group, they describe the sampled cells and do not support inference at the level of the mouse.",
]
for t in notes: ws4.append([t])
ws4["A1"].font = Font(bold=True, size=12); ws4.column_dimensions["A"].width = 140
for row in ws4.iter_rows():
    for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
wb.move_sheet("Notes", offset=-3)
wb.save(out)
# check against the values the plot was drawn from
ref = ({("Piezo1",r["grp"]):(float(r["pct"]),float(r["se"])) for r in csv.DictReader(open(D+"piezo1_pct_S3D_by_group.csv"))} if s3d else
       {(r["gene"],r["grp"]):(float(r["pct_pooled"]),float(r["se"])) for r in csv.DictReader(open(D+"channels_contractile_compartment_by_group.csv"))})
bad = [x for x in check if x[:2] in ref and (abs(ref[x[:2]][0]-x[2])>1e-6 or abs(ref[x[:2]][1]-x[3])>1e-6)]
print(f"{out}: {len(check)} bars, {sum(x[:2] in ref for x in check)} checked against plot data, {len(bad)} mismatches")
