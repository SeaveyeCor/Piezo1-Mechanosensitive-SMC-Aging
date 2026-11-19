#!/usr/bin/env python3
"""
Velocity pipeline orchestrator.

Generates the RNA velocity panels (Fig 4H-4J, Fig S4C).
Run from the velocity-analysis project root:

    cd data/velocity
    python -m velocity_pipeline.run_velocity_pipeline

Prerequisites:
  - velocity_output/velocity_smc.h5ad        (dynamical model)
  - velocity_output/velocity_smc_stochastic.h5ad  (stochastic model;
        produced by run_velocity_stochastic.py — required for Myocd
        and other low-abundance TFs that the dynamical model filters)
"""
from velocity_pipeline.figure_modules import Fig_velocity_summary


def main():
    print('=' * 70)
    print('VELOCITY PIPELINE — Piezo1 manuscript (Fig 4H-J)')
    print('=' * 70)

    print('\n── Velocity panels (Fig 4H-J) ──')
    Fig_velocity_summary.generate()


    print('\n' + '=' * 70)
    print('All velocity figures generated.')
    print('  Panels → velocity_pipeline/output/Fig4F*/')
    print('=' * 70)


if __name__ == '__main__':
    main()
