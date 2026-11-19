# ==============================================================================
# GENE SETS FOR SMC ANALYSIS PIPELINE
# ==============================================================================
# Centralized gene set definitions used across all figure modules
# ==============================================================================

#' Get all defined gene sets
#' @return Named list of gene vectors
get_gene_sets <- function() {
  list(
    # ============================================================================
    # SMC PHENOTYPE MARKERS
    # ============================================================================

    # Core contractile markers (for phenotype classification)
    contractile_core = c("Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn"),

    # Extended contractile markers (for detailed analysis)
    contractile_extended = c(
      "Myh11", "Acta2", "Tagln", "Cnn1", "Mylk", "Myocd", "Smtn",
      "Lmod1", "Myl6", "Myl9", "Cald1", "Tpm1", "Tpm2", "Des",
      "Actg2", "Myl12a", "Myl12b"
    ),

    # Contractile score gene set (comprehensive)
    contractile_score = c(
      # Core contractile
      "Acta2", "Myh11", "Myl9",
      # Cytoskeletal/structural
      "Cnn1", "Cald1", "Tagln", "Tpm1", "Tpm2",
      # Regulatory/modulatory
      "Calm1", "Calm2", "Calm3", "Mylk", "Rock1", "Rock2",
      "Prkca", "Prkcb", "Prkcd", "Dapk3", "Ppp1r14a"
    ),

    # Core synthetic markers (for phenotype classification)
    # Literature-supported panel:
    #   Col1a1: Yap 2021 (fibroblast-like vSMC), Schneider 2023 (colon), Wolfson 2025 (synthetic gene)
    #   Fn1: Yap 2021 (fibroblast-like vSMC), Wolfson 2025 (synthetic gene)
    #   Spp1: Yap 2021 (fibroblast-like vSMC marker)
    #   Klf4: Yap 2021 (master dedifferentiation switch), Wolfson 2025 (synthetic gene, log2FC=2.11)
    #   Vim: Schneider 2023 (fibroblast marker), Wolfson 2025 ("synthetic SMC marker")
    synthetic_core = c("Col1a1", "Fn1", "Spp1", "Klf4", "Vim"),

    # Extended synthetic markers
    synthetic_extended = c(
      "Col1a1", "Fn1", "Spp1", "Klf4", "Vim",
      "Col3a1", "Vcam1", "Tnc", "S100a4", "Mmp2", "Timp1", "Col1a2", "Sparc"
    ),

    # Synthetic score gene set (DEPRECATED - replaced by synthetic_phenotype)
    # Too heterogeneous: includes proliferation genes that dilute Piezo1/NFAT signal
    synthetic_score = c(
      "Klf4", "Klf5", "Mki67", "Pcna", "Top2a",
      "Col1a1", "Col1a2", "Col3a1", "Tgfb1", "Fn1", "Vcam1"
    ),

    # Synthetic phenotype score
    # "Phenotype Switch" score for Fig3G (Piezo1 vs Synthetic) and Fig4G (NFAT vs Synthetic)
    # Fibroblast-like markers of contractile-to-synthetic transition.
    # ZERO overlap with NFAT Activity score (Wnt5a, Mmp2, Mmp9, Vcam1, Fas, Timp1),
    # enabling clean, non-autocorrelated correlation analyses.
    #
    # Validation (bona fide SMCs, n=4,886):
    #   Piezo1 corr: Young rho=+0.051*, Aged rho=+0.182*** (strengthens with age)
    #   NFAT Pathway corr: Young ns, Aged rho=+0.232*** (strengthens with age)
    #   NFAT Activity corr: Young rho=+0.272***, Aged rho=+0.542*** (zero overlap)
    #   Aging increase: p=5.2e-26***  |  KO rescue: p=2.6e-26***
    #
    # Literature support:
    #   - Vim: Mesenchymal/synthetic marker; intermediate filament switch
    #          desmin→vimentin defines SMC dedifferentiation (Rensen 2007, Wolfson 2025)
    #   - Col1a2: Fibrillar collagen; hallmark of fibroblast-like/fibromyocyte
    #             SMC state (Wirka 2019, Yap 2021)
    #   - Dcn: Decorin; fibroblast-type proteoglycan, upregulated in modulated
    #          SMCs (Pan 2020)
    #   - Igfbp6: IGF binding protein 6; senescence-associated secretory marker,
    #             age-dependent upregulation in SMCs
    #   - Efemp1: Fibulin-3; matricellular protein associated with fibrosis
    #             and cellular senescence
    synthetic_phenotype = c(
      "Vim",      # Vimentin - mesenchymal/synthetic intermediate filament
      "Col1a2",   # Collagen I alpha-2 - fibroblast-like ECM production
      "Dcn",      # Decorin - fibroblast-type proteoglycan
      "Igfbp6",   # IGFBP6 - senescence-associated secretory marker
      "Efemp1"    # Fibulin-3 - matricellular, age/fibrosis-associated
    ),

    # ============================================================================
    # CONTROL CELL TYPE MARKERS
    # ============================================================================

    # PDGFRa+ cell markers
    pdgfra_markers = c("Pdgfra", "Kcnn3"),

    # ICC markers
    icc_markers = c("Kit", "Ano1"),

    # ============================================================================
    # PATHWAY GENE SETS
    # ============================================================================

    # Dedifferentiation markers (validated Piezo1-correlated, stronger with aging)
    # Young rho=+0.113**, Aged rho=+0.144***
    dedifferentiation = c(
      "Klf4",     # Master dedifferentiation TF
      "Vim",      # Vimentin - synthetic/mesenchymal marker
      "S100a4",   # Fibroblast marker, SMC dedifferentiation
      "Cd34",     # Progenitor marker
      "Cd44",     # Adhesion/migration
      "Ly6a"      # Sca-1, stem cell antigen
    ),

    # MMP/TIMP remodeling axis (for Fig4G: NFAT vs MMP correlation)
    # NFAT vs MMP: Young rho=+0.275***, Aged rho=+0.415***
    mmp_axis = c(
      "Mmp2",     # Matrix metalloproteinase 2
      "Mmp9",     # Matrix metalloproteinase 9
      "Timp1",    # Tissue inhibitor of MMPs 1
      "Timp2",    # Tissue inhibitor of MMPs 2
      "Timp3"     # Tissue inhibitor of MMPs 3
    ),

    # NFAT pathway (transcription factors and regulators)
    nfat_pathway = c(
      "Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4", "Nfat5",
      "Ppp3ca", "Ppp3cb", "Ppp3cc", "Rcan1", "Rcan2", "Rcan3",
      "Calm1", "Calm3"
    ),

    # ============================================================================
    # NFAT-RESPONSIVE REMODELING GENE SET (DATA-VALIDATED in cluster-based SMCs)
    # ============================================================================
    # These 6 genes are part of the NFAT-associated remodeling program.
    # They show: UP with aging, DOWN with Piezo1 KO (rescue).
    #
    # IMPORTANT — EVIDENCE LEVELS:
    #   Mmp2  — DIRECT NFATc3 target (ChIP + luciferase, Zhang et al. 2025)
    #   Mmp9  — DIRECT NFATc3 target (ChIP + luciferase, Zhang et al. 2025)
    #   Vcam1 — NFAT-family responsive (NFAT1 promoter binding, not NFATc3-specific)
    #   Fas   — NFAT-associated (NFAT regulates FasL; Fas co-regulated in context)
    #   Timp1 — Co-regulated with MMP2/9 (no direct NFAT promoter evidence)
    #   Wnt5a — UPSTREAM AMPLIFIER, not downstream target. Wnt5a activates
    #           NFAT via non-canonical Wnt/Ca2+ pathway (Dejmek 2006, Dodge 2012).
    #           No evidence for NFATc3 transcribing Wnt5a.
    #
    # Display label "NFAT Activity" (in the paper: NFAT-associated remodeling score),
    # the downstream readout, as opposed to the upstream "NFAT Pathway" machinery score.
    # The variable name nfatc3_remodeling is used by the downstream scripts.
    #
    # SELECTION: Unbiased screen of 26,309 genes identified 691
    # showing NFAT-like pattern. These 6 are NFAT-associated genes that
    # ALSO show the expected pattern in our cluster-based SMC data:
    #   - >50% increase with aging (Aged_Ctrl vs Young_Ctrl)
    #   - >30% decrease with KO (Aged_KO vs Aged_Ctrl)
    #   - Both effects p<0.05
    #
    # REMOVED: Tnf (goes UP with KO, opposite of expected)
    #          Rcan2 (no aging effect in Myh11+, weak in cluster-based)
    # ADDED:   Fas (apoptosis receptor, +103% aging, -51% rescue)
    #          Timp1 (MMP inhibitor, +1369% aging, -77% rescue)

    # NFAT Activity score (DATA-VALIDATED)
    # Variable name kept as nfatc3_remodeling for backward compatibility
    # Display label: "NFAT Activity" (not "NFATc3 Targets",
    # not "NFAT-Responsive Remodeling")
    nfatc3_remodeling = c(
      "Wnt5a",    # Wnt signaling (UPSTREAM AMPLIFIER) - +127% aging, -82% rescue
      "Mmp2",     # Matrix metalloproteinase 2 (DIRECT target) - +179% aging, -71% rescue
      "Mmp9",     # Matrix metalloproteinase 9 (DIRECT target) - +427% aging, -93% rescue
      "Vcam1",    # Vascular cell adhesion (NFAT-responsive) - +123% aging, -37% rescue
      "Fas",      # Apoptosis receptor (NFAT-associated) - +103% aging, -41% rescue
      "Timp1"     # Tissue inhibitor of MMPs (co-regulated) - +1369% aging, -92% rescue
    ),

    # DEPRECATED: Original literature-based set (kept for reference)
    # Some genes did NOT show expected pattern in our data
    nfatc3_remodeling_original = c(
      "Wnt5a", "Mmp2", "Mmp9", "Vcam1",
      "Tnf",   # WARNING: Goes UP with KO (+227% in cluster-based)
      "Rcan2"  # WARNING: No aging effect (-2% in Myh11+)
    ),

    # Immediate early genes (calcium-responsive transcription factors)
    # Pattern: DOWN with aging, UP with KO (restored calcium responsiveness)
    immediate_early_genes = c(
      "Fos",      # AP-1 component - rapid calcium response
      "Fosb",     # AP-1 component
      "Jun",      # AP-1 component
      "Junb",     # AP-1 component
      "Egr1"      # Early growth response 1 - calcium-induced
    ),

    # Data-driven aging signature (broader set, not all NFAT-validated)
    # Genes upregulated with aging, rescued by Piezo1 KO
    aging_rescue_signature = c(
      # Inflammation/complement
      "C3",       # Complement C3 - inflammation/innate immunity
      "C4b",      # Complement C4b - inflammation
      "Vcam1",    # Vascular cell adhesion molecule 1
      # ECM/fibrosis markers
      "Vim",      # Vimentin - synthetic phenotype marker
      "Lum",      # Lumican - ECM proteoglycan
      "Mgp",      # Matrix Gla protein - calcification
      "Efemp1",   # Fibulin-like ECM protein
      # Growth/signaling
      "Igfbp6",   # IGF binding protein 6
      "Wnt5a",    # Wnt signaling
      # Matrix remodeling
      "Mmp2",     # Matrix metalloproteinase 2
      # Other
      "Thsd4",    # Thrombospondin domain 4
      "Gpm6a",    # Glycoprotein M6A
      "Cadm2",    # Cell adhesion molecule 2
      "Upk3b",    # Uroplakin 3B
      "Cst3"      # Cystatin C
    ),

    # nfat_targets: alias gene list used by some downstream scripts
    nfat_targets = c(
      "C3", "Igfbp6", "Vim", "Thsd4", "Gpm6a", "Cadm2", "Upk3b",
      "C4b", "Mgp", "Cst3", "Lum", "Efemp1", "Wnt5a", "Mmp2", "Vcam1"
    ),

    # ============================================================================
    # JCI INSIGHT MECHANOTRANSDUCTION MARKERS
    # ============================================================================
    # Reference: JCI Insight 2025;10(20):e188669
    # Pathological mechanotransduction genes in HISMCs
    # Pattern: UP with aging, DOWN with Piezo1 KO (rescue)

    # Core significant genes (p<0.001 for rescue, strongest effects)
    jci_mechano_core = c(
      "Wnt5a",    # Wnt signaling - 2.27x aging, 0.18x rescue (p=1.01e-16)
      "Il6st",    # IL6 signal transducer (gp130) - 1.22x aging, 0.79x rescue (p=5.64e-05)
      "Bmp6",     # BMP6 - 1.78x aging, 0.53x rescue (p=1.05e-04)
      "Fgf5"      # FGF5 - 3.58x aging, 0.05x rescue (p=7.93e-04)
    ),

    # Extended significant genes - 8-gene rescue panel (all p<0.05 for rescue)
    # These genes show consistent pattern: UP with aging, DOWN with Piezo1 KO
    jci_mechano_extended = c(
      "Wnt5a",    # Differentiation - 2.27x aging, 0.18x rescue (p=1.01e-16 ***)
      "Il6st",    # IL6-STAT3 - 1.22x aging, 0.79x rescue (p=5.64e-05 ***)
      "Bmp6",     # Differentiation - 1.78x aging, 0.53x rescue (p=1.05e-04 ***)
      "Fgf5",     # Growth Factors - 3.58x aging, 0.05x rescue (p=7.93e-04 ***)
      "Lif",      # IL6-STAT3 - 6.98x aging, 0.26x rescue (p=8.34e-03 **)
      "Pdgfb",    # Growth Factors - 1.13x aging, 0.60x rescue (p=1.97e-03 **)
      "Fgf18",    # Growth Factors - 2.63x aging, 0.18x rescue (p=8.73e-03 **)
      "Bmp2"      # Differentiation - 4.40x aging, 0.39x rescue (p=3.70e-02 *)
    ),

    # IL6-STAT3 pathway (complete)
    jci_il6_stat3 = c(
      "Il6",      # Interleukin 6
      "Il11",     # Interleukin 11
      "Il6st",    # IL6 signal transducer (gp130) - shared receptor
      "Lif",      # Leukemia inhibitory factor
      "Stat3"     # Signal transducer and activator of transcription 3
    ),

    # Growth factors from JCI paper
    jci_growth_factors = c(
      "Ereg",     # Epiregulin - EGFR ligand
      "Areg",     # Amphiregulin - EGFR ligand
      "Hbegf",    # Heparin-binding EGF-like growth factor
      "Fgf5",     # Fibroblast growth factor 5
      "Fgf7",     # Fibroblast growth factor 7
      "Fgf18",    # Fibroblast growth factor 18
      "Nrg1",     # Neuregulin 1
      "Pdgfa",    # PDGF-A
      "Pdgfb",    # PDGF-B
      "Lif"       # Leukemia inhibitory factor
    ),

    # Cytokines from JCI paper
    jci_cytokines = c(
      "Il6",      # Interleukin 6
      "Il11",     # Interleukin 11
      "Csf3",     # G-CSF
      "Clcf1"     # Cardiotrophin-like cytokine factor 1
    ),

    # Differentiation regulators from JCI paper
    jci_differentiation = c(
      "Wnt5a",    # Wnt5a
      "Bmp2",     # BMP2
      "Bmp6",     # BMP6
      "Dkk1",     # Dickkopf 1 (Wnt inhibitor)
      "Tgfb1",    # TGF-beta 1
      "Jag1",     # Jagged 1 (Notch ligand)
      "Dll4",     # Delta-like 4 (Notch ligand)
      "Inhbb",    # Inhibin beta B
      "Inhba"     # Inhibin/Activin beta A
    ),

    # Signaling mediators from JCI paper
    jci_signaling = c(
      "Acvr1",    # Activin receptor type 1
      "Stat3",    # STAT3
      "Fgfr1",    # FGF receptor 1
      "Il6st",    # IL6 signal transducer
      "Egfr",     # EGF receptor
      "Erbb4",    # ERBB4/HER4
      "Irf1",     # Interferon regulatory factor 1
      "Egr1"      # Early growth response 1
    ),

    # TGF-beta / stress pathway
    tgfbeta_pathway = c(
      "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Tgfbr3",
      "Smad2", "Smad3", "Smad4", "Smad7",
      "Junb", "Jun", "Fos", "Fosb"  # stress response TFs
    ),

    # Calcium signaling (SOCE/IP3)
    calcium_pathway = c(
      "Calm1", "Calm2", "Calm3", "Camk2a", "Camk2b",
      "Atp2a1", "Atp2a2", "Atp2b1", "Atp2b4",
      "Itpr1", "Itpr2", "Itpr3", "Ryr1", "Ryr2",
      "Stim1", "Stim2", "Orai1", "Orai2"
    ),

    # Focal adhesion / integrins
    focal_adhesion = c(
      "Itgav", "Itga5", "Itgb1", "Itgb3", "Itgb5",
      "Pxn", "Ptk2", "Vcl", "Tln1", "Tln2",
      "Ilk", "Fermt2", "Zyx", "Vasp"
    ),

    # YAP/TAZ mechanotransduction
    yap_taz_pathway = c(
      "Yap1", "Wwtr1", "Tead1", "Tead2", "Tead3", "Tead4",
      "Lats1", "Lats2", "Stk3", "Stk4", "Mob1a", "Mob1b",
      "Ctgf", "Cyr61", "Ankrd1"
    ),

    # MAPK/PI3K/Akt pathway
    mapk_pi3k_akt = c(
      "Mapk1", "Mapk3", "Mapk8", "Mapk14",
      "Pik3ca", "Pik3cb", "Pik3r1",
      "Akt1", "Akt2", "Mtor", "Rps6kb1",
      "Raf1", "Braf", "Map2k1", "Map2k2"
    ),

    # Cellular senescence
    cellular_senescence = c(
      "Cdkn1a", "Cdkn2a", "Cdkn2b", "Trp53", "Rb1",
      "Ccnd1", "Cdk4", "Cdk6", "E2f1",
      "Serpine1", "Il6", "Cxcl1", "Ccl2"
    ),

    # Rho pathway
    rho_pathway = c(
      "Rhoa", "Rhob", "Rhoc", "Rac1", "Cdc42",
      "Rock1", "Rock2", "Arhgef1", "Arhgef2",
      "Limk1", "Limk2", "Cofl1", "Pfn1"
    ),

    # Mechanosensing (YAP/TAZ integrated)
    mechanosensing = c(
      "Piezo1", "Piezo2", "Trpc1", "Trpc6",
      "Cav1", "Yap1", "Wwtr1", "Tead1"
    ),

    # ============================================================================
    # ECM GENE SETS (data-validated)
    # ============================================================================
    # Validation criteria: >20% UP with aging, >20% DOWN with KO, both p<0.05
    # in cluster-based SMCs (5,681 cells)
    #
    # Structure: no overlap with the NFAT scores (Mmp2, Mmp9 and Timp1 are
    # only in nfatc3_remodeling).

    # ECM Structural — UNIFIED ECM SCORE (5 genes; no overlap with NFAT scores)
    # Bgn, Col6a1, Dcn, Sparc carry the structural/glycoprotein signal;
    # Timp2 is the only ECM protease-axis gene with no NFAT-target overlap.
    ecm_structural = c(
      "Bgn",      # +50% aging, -21% rescue (proteoglycan)
      "Col6a1",   # +29% aging, -31% rescue (collagen VI)
      "Dcn",      # +54% aging, -28% rescue (proteoglycan)
      "Sparc",    # +26% aging, -42% rescue (matricellular)
      "Timp2"     # +23% aging, -40% rescue (MMP inhibitor; not in nfatc3_remodeling)
    ),

    # ---- DEPRECATED gene sets (kept for backward compatibility) ----
    # Use `ecm_structural` instead. These names are retained so existing
    # analyses that loaded older pipeline state continue to function.

    # Not used — replaced by `ecm_structural`.
    # 3/4 of these genes (Mmp2, Mmp9, Timp1) overlap with nfatc3_remodeling.
    ecm_remodeling_validated = c(
      "Mmp2",     # +179% aging, -71% rescue   [also in nfatc3_remodeling]
      "Mmp9",     # +427% aging, -93% rescue   [also in nfatc3_remodeling]
      "Timp1",    # +1369% aging, -92% rescue  [also in nfatc3_remodeling]
      "Timp2"     # +23% aging, -40% rescue
    ),

    # ECM Remodeling - Original (for reference)
    ecm_remodeling = c(
      "Mmp2", "Mmp3", "Mmp9", "Mmp14",
      "Timp1", "Timp2", "Timp3",
      "Adamts1", "Adamts4", "Adamts5"
    ),

    # Not used — replaced by `ecm_structural` (which adds Timp2).
    ecm_organization_validated = c(
      "Col6a1",   # +29% aging, -31% rescue
      "Sparc",    # +26% aging, -42% rescue
      "Bgn",      # +50% aging, -21% rescue
      "Dcn"       # +54% aging, -28% rescue
    ),

    # ECM Organization - Original (for reference)
    ecm_organization = c(
      "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col6a1",
      "Fn1", "Sparc", "Bgn", "Dcn", "Vcan"
    ),

    # Collagens - DATA-VALIDATED (only 2 show pattern)
    ecm_collagens_validated = c(
      "Col6a1",   # +29% aging, -31% rescue
      "Col6a2"    # +91% aging, -42% rescue
    ),

    # Collagens - Full set (for reference)
    ecm_collagens = c(
      "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
      "Col5a1", "Col5a2", "Col6a1", "Col6a2", "Col6a3",
      "Col8a1", "Col11a1", "Col12a1", "Col14a1"
    ),

    # ECM glycoproteins - DATA-VALIDATED (4 genes)
    ecm_glycoproteins_validated = c(
      "Sparc",    # +26% aging, -42% rescue
      "Bgn",      # +50% aging, -21% rescue
      "Dcn",      # +54% aging, -28% rescue
      "Fbln5"     # +79% aging, -55% rescue
    ),

    # ECM glycoproteins - Full set (for reference)
    ecm_glycoproteins = c(
      "Fn1", "Tnc", "Thbs1", "Thbs2", "Sparc", "Postn",
      "Vcan", "Bgn", "Dcn", "Lum", "Fbln1", "Fbln2", "Fbln5"
    ),

    # ECM catabolism (degradation focused)
    ecm_catabolism = c(
      "Mmp2", "Mmp3", "Mmp9", "Mmp14", "Mmp12", "Mmp13",
      "Adamts1", "Adamts4", "Adamts5", "Adamts9",
      "Ctsk", "Ctsl", "Ctsb",  # Cathepsins
      "Plau", "Plaur"  # uPA system
    ),

    # ============================================================================
    # MUSCLE CONTRACTION GO GENES
    # ============================================================================

    muscle_contraction = c(
      # Core contractile apparatus
      "Myh11", "Myh9", "Myl9", "Myl6", "Mylk", "Mylk3",
      "Acta2", "Actg2", "Actc1",
      "Tpm1", "Tpm2", "Tpm4",
      "Tagln", "Tagln2",
      "Cnn1", "Cnn2", "Cnn3",
      "Smtn", "Lmod1",
      # Thin filament regulation
      "Tnnc1", "Tnnc2", "Tnni1", "Tnni2", "Tnnt1", "Tnnt2",
      # Calcium handling
      "Calm1", "Calm2", "Calm3", "Atp2a2", "Ryr2",
      # Z-disc/cytoskeletal
      "Des", "Actn1", "Actn2", "Actn4", "Csrp1", "Csrp2"
    ),

    # ============================================================================
    # CASCADE MODULE GENE SETS (for Fig4I pseudotime cascade)
    # ============================================================================

    # Cytoskeletal remodeling module (intermediate between contractile and synthetic)
    cascade_cyto_remodeling = c("Ezr", "Syne2", "Lmna", "Vim"),

    # ============================================================================
    # CYTOSKELETAL REMODELING SCORE
    # ============================================================================
    # Aging-up / KO-rescue cytoskeletal remodeling genes. Same data-validation
    # criteria as ecm_structural and nfatc3_remodeling: each gene shows the
    # canonical aging-up / KO-rescue pattern in cluster-based SMCs.
    #
    # Starting set is the cascade_cyto_remodeling gene panel (Ezr, Syne2, Lmna,
    # Vim) which was already validated as the cytoskeletal step in the Fig4I
    # pseudotime cascade. Vim overlaps with synthetic_phenotype (Vim, Col1a2,
    # Dcn, Igfbp6, Efemp1) — this is intentional (it's a mesenchymal/cytoskeletal
    # marker that drives both phenotypes), but flagged in the score's Methods text.
    #
    # Used by FigS4H_cyto_remodeling.R as the canonical 4-group cyto panel,
    # parallel to Fig4E (NFAT-Resp Remodeling) and FigS3K (ECM Structural).
    cyto_remodeling = c(
      "Ezr",      # Ezrin — actin-membrane crosslinker, mesenchymal
      "Syne2",    # Nesprin-2 — nuclear envelope linker, cytoskeletal coupling
      "Lmna",     # Lamin A/C — nuclear lamina, mechanotransduction
      "Vim"       # Vimentin — intermediate filament, mesenchymal/synthetic marker
    ),

    # NFAT effectors (downstream remodeling — same as nfatc3_remodeling)
    cascade_nfat_effectors = c("Wnt5a", "Mmp2", "Mmp9", "Vcam1", "Fas", "Timp1"),

    # Synthetic transcription/dedifferentiation markers
    cascade_synthetic = c("Klf4", "Klf5", "Lgals3", "Cdkn1a", "Icam1", "Vcan"),

    # ============================================================================
    # P1+ vs P1- COMPARISON GENE SETS (for Fig4J)
    # ============================================================================

    # Synthetic markers panel (broad)
    p1_comparison_synthetic = c(
      "Klf4", "Klf5", "Vim", "Fn1", "Col1a1", "Spp1",
      "Lgals3", "Cdkn1a", "Icam1", "Vcan", "S100a4"
    ),

    # Cytoskeletal remodeling panel
    p1_comparison_cyto = c(
      "Ezr", "Vim", "Myl12b", "Add3", "Arpc1b", "Syne2",
      "Tagln2", "Lmna", "Flna", "Vcl", "Actn1"
    ),

    # NFAT / Ca2+ signaling panel
    p1_comparison_nfat_ca = c(
      "Nfatc1", "Nfatc2", "Nfatc4", "Ppp3ca", "Calm1",
      "Itpr2", "Itpr3", "Trpm7", "Rcan1"
    ),

    # ECM / Remodeling panel
    p1_comparison_ecm = c(
      "Mmp2", "Mmp9", "Timp1", "Timp2", "Col6a1",
      "Sparc", "Bgn", "Dcn", "Loxl1"
    ),

    # ============================================================================
    # PATHWAY CORRELATION GENE SETS (for FigS3G/H)
    # ============================================================================

    pathway_correlation_genes = c(
      # Positive correlation expected
      "Nfatc1", "Klf4", "Fosb", "Itpr2", "Tgfb1", "Junb",
      "Ppp3ca", "Rcan1", "Itpr1", "Cdkn1a", "Jun", "Pxn",
      "Tgfbr2", "Nfatc3", "Itgav", "Smad3", "Itga5", "Smad7",
      "Rock1", "Eln", "Lats2", "Stim1", "Itgb3", "Mapk14",
      "Smad2", "Fos", "Adamts1", "Orai1", "Nfatc2", "Tgfbr1",
      "Mapk1", "Mapk8", "Mtor", "Fbln5", "Ptk2", "Wwtr1",
      "Adamts4", "Stim2", "Trp53", "Col4a1", "Akt1", "Pik3ca",
      "Lats1", "Il6", "Yap1", "Pcna", "Smad4", "Cxcl1", "Rock2",
      # Synthetic classification panel genes
      "Col1a1", "Fn1", "Spp1", "Vim",
      # Negative correlation expected
      "Cnn1", "Myh11", "Acta2", "Tagln"
    )
  )
}

#' Get a specific gene set by name
#' @param set_name Character string of gene set name
#' @return Character vector of genes
get_gene_set <- function(set_name) {
  sets <- get_gene_sets()
  if (!set_name %in% names(sets)) {
    stop(sprintf("Gene set '%s' not found. Available: %s",
                 set_name, paste(names(sets), collapse = ", ")))
  }
  return(sets[[set_name]])
}

#' Filter gene set to available genes in expression matrix
#' @param genes Character vector of genes
#' @param available_genes Character vector of genes in dataset
#' @return List with available and missing genes
filter_available_genes <- function(genes, available_genes) {
  available <- genes[genes %in% available_genes]
  missing <- genes[!genes %in% available_genes]

  list(
    available = available,
    missing = missing,
    n_available = length(available),
    n_missing = length(missing),
    pct_available = 100 * length(available) / length(genes)
  )
}
