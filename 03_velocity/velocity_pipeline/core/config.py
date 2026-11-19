"""
Velocity pipeline configuration — group order, colors, valid stat pairs, gene sets.
"""
from pathlib import Path

# Data paths
ROOT = Path('data/velocity')
H5AD_DYNAMICAL = ROOT / 'velocity_output' / 'velocity_smc.h5ad'
H5AD_STOCHASTIC = ROOT / 'velocity_output' / 'velocity_smc_stochastic.h5ad'
OUTPUT_DIR = ROOT / 'velocity_pipeline' / 'output'

# Groups
GROUP_ORDER = ['Young_Ctrl', 'Aged_Ctrl', 'Young_KO', 'Aged_KO']
# 4-group bar labels: age in months only — bar color encodes genotype
# (Ctrl grey, P1KO pink) plus a legend, so the labels stay uncluttered.
GROUP_LABELS = {'Young_Ctrl': '3',
                'Aged_Ctrl':  '24',
                'Young_KO':   '3',
                'Aged_KO':    '24'}

# Manuscript-consistent colors — match smc-analysis pipeline get_colors():
# Ctrl bars (Young and Aged) all use the same grey, KO bars (Young and Aged)
# all use the same pink. Genotype is the only thing color encodes; age is on
# the x-axis label.
GROUP_COLORS = {'Young_Ctrl': '#7F7F7F', 'Aged_Ctrl': '#7F7F7F',
                'Young_KO':   '#ff7fb2', 'Aged_KO':   '#ff7fb2'}

# Pairwise comparisons. Each tuple = (group_a, group_b, label).
# Excludes the two diagonal cross-group comparisons that conflate aging and
# genotype simultaneously and have no clean biological interpretation:
#   - Young_Ctrl vs Aged_KO  (both age and genotype change)
#   - Young_KO   vs Aged_Ctrl (both age and genotype change)
# Included pairs cleanly isolate one variable each:
PAIRWISE_COMPARISONS = [
    ('Young_Ctrl', 'Aged_Ctrl', 'Aging (WT)'),       # age effect at WT
    ('Young_KO',   'Aged_KO',   'Aging (KO)'),       # age effect under KO
    ('Young_Ctrl', 'Young_KO',  'KO effect (young)'),# genotype effect at young
    ('Aged_Ctrl',  'Aged_KO',   'KO rescue (aged)'), # genotype effect at aged
]
EXCLUDED_PAIRS = [('Young_Ctrl', 'Aged_KO'), ('Young_KO', 'Aged_Ctrl')]

# Gene sets (from smc-analysis pipeline/core/gene_sets.R)
CONTRACTILE_CORE = ['Myh11', 'Acta2', 'Tagln', 'Cnn1', 'Mylk', 'Myocd', 'Smtn']

# Subset that survives the dynamical-mode kinetic fit (verified against
# velocity_smc.h5ad). Used for the dynamical CCI.
CONTRACTILE_DYNAMICAL = ['Myh11', 'Cnn1', 'Lmod1', 'Actg2']

# Genes that need stochastic mode for velocity recovery (sparse TFs).
CONTRACTILE_STOCHASTIC_ONLY = ['Myocd']  # add Srf, Mef2c here if recovered later
