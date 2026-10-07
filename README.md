# MHW2CoralEcoSys

Primary demo: **[Physiology-constrained MaxEnt training and projection](Physiology_MaxEnt_demo.ipynb)**. [中文说明](README.zh-CN.md)

This executed notebook starts at the beginning of the study: real occurrence/background inputs → native MaxEnt training → spatial holdout evaluation → baseline and SSP projections → physiological constraint → figures. It is not a conservation-prioritization demo or website.

## Run

Install Python 3.12 and Java 17+ (`java` on PATH), then create and activate a virtual environment and run:

```bash
pip install -r requirements.txt
jupyter lab Physiology_MaxEnt_demo.ipynb
```

Open from the repository root and select **Run → Run All Cells**. Tables and figures are already saved for preview. New model, logs, predictions and PNG/SVG/PDF figures go to `outputs/physiology_maxent_demo/`. Original project sources stay read-only.

Automated checks execute the same notebook:

```bash
python -m unittest discover -s tests
python scripts/execute_notebook.py --notebook Physiology_MaxEnt_demo.ipynb
```

## Actual inputs and scope

- Seeded subset of 1,500 actual Acropora presence records and 10,000 target-group background records, retaining parent row identifiers. Original final seven predictors, FC=LQHP and RM=0.5 are preserved.
- One newly fitted native Java MaxEnt 3.4.4 model with fixed 5° block holdout. Full tuning and the parent ten-bootstrap ensemble are not rerun.
- 2,000 fixed global background locations sampled from actual GFDL-ESM4 baseline and SSP126/245/585 rasters. Static bathymetry/rugosity come from parent SWD; missing environment stays missing.
- Accepted coral constraint: `HSI_constrained = HSI_MaxEnt × Φ^0.5`, with asymmetric temperature performance × aragonite-saturation sigmoid. It is **post-projection multiplication**, not a physiological feature injected into training or a modified MaxEnt objective.
- Physiological coefficients are original structural assumptions, not newly estimated experimental parameters. The demo does not establish improved prediction from the constraint.
- Holdout AUC is presence-background discrimination, not true-absence accuracy. Logistic output is an HSI index; sample-location means are not global reef-area-weighted means. Range checks are not full MESS or future validation.

This demonstrates the requested starting module using actual data, not all four taxa, three climate members, full global rasters or the entire paper. [Lineage, hashes and parent scripts](analysis/modules/physiology_maxent_demo/) document the boundary.

## Layout

```text
Physiology_MaxEnt_demo.ipynb             Primary executed research demo
analysis/modules/physiology_maxent_demo/ Actual subset, engine, constraints, provenance
scripts/curate_maxent_demo.py            Re-curation from original project
scripts/execute_notebook.py             Notebook execution checks
tests/                                  Formula and regression checks
MHW2CoralEcoSys_demo.ipynb               Earlier optional conservation module
archive/                                Historical webpage, CLI and earlier guides
```

MaxEnt 3.4.4 and its license notices are bundled under `vendor/`, from the authors' [official source repository](https://github.com/mrmaxent/Maxent/tree/963d0114e05a39f92f832dfbe800c22442d9d8a6/ArchivedReleases/3.4.4). Scientific data retain original source reuse conditions.

[Target contribution: PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1). Before maintainer merge, use the contribution branch or contributor fork rather than target `main`.
