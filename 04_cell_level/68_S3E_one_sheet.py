# PAPER PANEL: Fig S3E (Excel workbook behind the panel)
# 68: ONE easy-to-read sheet for the S3E compensation plot (contractile SMCs, Piezo1 excluded):
# counts, % and SEM (formulas), the depth-matched difference with its 95% interval, and ONE p value per
# comparison: the depth-adjusted logistic p, Bonferroni-adjusted for the 8 channels tested at each age.
# The unadjusted p is kept in a clearly labeled last column. Usage: python3 68_S3E_one_sheet.py <out.xlsx>
import csv, sys
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
D = "output/cell_level/data/"
CH = ["Piezo2","Trpc1","Trpc3","Tmem63a","Tmem63b","Pkd1","Trpm7","Pkd2"]  # Piezo1 (deleted gene) excluded
grp = {(r["gene"], r["grp"]): r for r in csv.DictReader(open(D+"channels_contractile_compartment_by_group.csv"))}
st = {}
for f, age in [("channels_contractile_compartment_younggeno.csv", 3), ("channels_contractile_compartment_agedgeno.csv", 24)]:
    for r in csv.DictReader(open(D+f)): st[(r["gene"], age)] = r
wb = Workbook(); ws = wb.active; ws.title = "S3E channel detection"
H = ["Channel","Age (months)","CTRL SMCs (n)","CTRL detecting (n)","CTRL % detecting","CTRL SEM (%)",
     "P1KO SMCs (n)","P1KO detecting (n)","P1KO % detecting","P1KO SEM (%)",
     "Difference, P1KO − CTRL (pp, depth-matched)","95% interval","p value","Significant?","(Unadjusted p, for reference only)"]
ws.append(["Figure S3E: detection of mechanically gated ion channels in contractile SMCs, P1KO vs CTRL within age. No channel differs significantly at either age."])
ws["A1"].font = Font(bold=True, size=12)
ws.append(["p value = logistic regression adjusted for sequencing depth (log genes detected per cell), Bonferroni-corrected for the 8 channels tested at each age (unadjusted p × 8, max 1). Bars in the plot = % detecting ± SEM, SEM = √(p(1−p)/n). Cells are the unit of analysis. Piezo1, the deleted gene, is not included."])
ws.append([]); ws.append(H); hr = ws.max_row
for c in ws[hr]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="DDDDDD"); c.alignment = Alignment(wrap_text=True, vertical="top")
for age, (gc, gk) in [(3, ("Young_Ctrl","Young_KO")), (24, ("Aged_Ctrl","Aged_KO"))]:
    for g in CH:
        c, k, s = grp[(g, gc)], grp[(g, gk)], st[(g, age)]
        nc, nk = int(c["n_pooled"]), int(k["n_pooled"])
        dc, dk = round(float(c["pct_pooled"])*nc/100), round(float(k["pct_pooled"])*nk/100)
        p = float(s["logit_p"]); padj = min(1.0, p*len(CH)); r = ws.max_row + 1
        ws.append([g, age, nc, dc, f"=D{r}/C{r}*100", f"=SQRT((D{r}/C{r})*(1-D{r}/C{r})/C{r})*100",
                   nk, dk, f"=H{r}/G{r}*100", f"=SQRT((H{r}/G{r})*(1-H{r}/G{r})/G{r})*100",
                   round(float(s["matched_diff_pp"]), 1), f'{float(s["lo"]):+.1f} to {float(s["hi"]):+.1f}',
                   padj, "Yes" if padj < 0.05 else "No (ns)", p])
        for col in "EFIJ": ws[f"{col}{r}"].number_format = "0.0"
        ws[f"M{r}"].number_format = ws[f"O{r}"].number_format = "0.00"
        ws[f"O{r}"].font = Font(color="808080")
        assert abs(100*dc/nc - float(c["pct_pooled"])) < 1e-6 and abs(100*dk/nk - float(k["pct_pooled"])) < 1e-6
        print(f"{g:8s} {age:2d}  {100*dc/nc:5.1f}% vs {100*dk/nk:5.1f}%  p_adj={padj:.2f}  (raw {p:.3f})")
for j, w in enumerate([10,8,9,10,10,9,9,10,10,9,14,15,9,11,13], 1): ws.column_dimensions[get_column_letter(j)].width = w
ws.freeze_panes = f"B{hr+1}"
wb.save(sys.argv[1])
