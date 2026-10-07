# MHW2CoralEcoSys

Primary demo: **[Small-sample MaxEnt training and full native-grid projection](Physiology_MaxEnt_demo.ipynb)**. [中文说明](README.zh-CN.md)

[Current full-grid notebook package](https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/tag/v5.0.0-full-grid-maxent-demo).

The release also contains `full_grid_results.zip`: twenty verified GeoTIFFs, the newly fitted lambdas and continuous maps. These are derived demo results, not the paper's full four-taxon ensemble.

This executed notebook starts at the beginning of the study: real occurrence/background inputs → native MaxEnt training → spatial holdout evaluation → baseline and SSP projections → physiological constraint → figures. It is not a conservation-prioritization demo or website.

## Run

Install Python 3.12 and Java 17+ (`java` on PATH), then create and activate a virtual environment and run:

```bash
pip install -r requirements.txt
python scripts/acquire_full_grid_inputs.py
jupyter lab Physiology_MaxEnt_demo.ipynb
```

Open from the repository root and select **Run → Run All Cells**. Tables and figures are already saved for preview. New model, logs, predictions and PNG/SVG/PDF figures go to `outputs/physiology_maxent_demo/`. Original project sources stay read-only.

The complete input archive is 624,654,011 bytes (about 596 MiB), separate from Git. The acquisition command verifies its archive and individual-file SHA-256 checksums. It contains 26 actual full source rasters, not downsampled grids. Obtain it from [full-grid-inputs-v1](https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/tag/full-grid-inputs-v1). Local re-curation is also available through `scripts/prepare_full_grid_inputs.py --project ORIGINAL_PROJECT_PATH`.

Automated checks execute the same notebook:

```bash
python -m unittest discover -s tests
python scripts/execute_notebook.py --notebook Physiology_MaxEnt_demo.ipynb
```

## Actual inputs and scope

- Seeded subset of 1,500 actual Acropora presence records and 10,000 target-group background records, retaining parent row identifiers. Original final seven predictors, FC=LQHP and RM=0.5 are preserved.
- One newly fitted native Java MaxEnt 3.4.4 model with fixed 5° block holdout. Full tuning and the parent ten-bootstrap ensemble are not rerun.
- The native training uses the original default 500-iteration limit. Successful execution does not independently establish optimizer convergence or biological validity; the actual training sample count and holdout AUC are reported separately.
- Complete 4338×1928 native grids for GFDL-ESM4 baseline and SSP126/245/585. Every complete predictor cell is sent to native `density.Project`; there is no spatial subsampling or interpolation of sampled predictions. Static bathymetry/rugosity use the original full rasters, exactly cropped to the aligned member extent. Missing environments stay NoData.
- Accepted coral constraint: `HSI_constrained = HSI_MaxEnt × Φ^0.5`, with asymmetric temperature performance × aragonite-saturation sigmoid. It is **post-projection multiplication**, not a physiological feature injected into training or a modified MaxEnt objective.
- Physiological coefficients are original structural assumptions, not newly estimated experimental parameters. The demo does not establish improved prediction from the constraint.
- Holdout AUC is presence-background discrimination, not true-absence accuracy. Logistic output is an HSI index. Full-grid means are cell-area-weighted over the common complete environmental domain, not restricted to coral-reef polygons. Range checks are not full MESS or future validation.

This implements small-sample training and full spatial projection for one taxon/member, not all four taxa, three climate members or the entire paper. The post-projection combination matches manuscript_第八版20260806 Eq. 1. All 24,381 values of the active parent Acropora physiology curve match exactly. [Formula audit](docs/mainline_physiology_audit.json) and [lineage](analysis/modules/physiology_maxent_demo/) document the sources. The earlier sampled-point notebook is historical material under `archive/point_projection_20261007/`.

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
