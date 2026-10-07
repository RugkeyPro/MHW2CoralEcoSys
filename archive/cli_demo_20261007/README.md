# MHW2CoralEcoSys

Research code, a real reef-cell demonstration dataset and a runnable scientific analysis for marine-heatwave-related coral-community suitability, microplastic exposure and conservation prioritization. [中文说明](README.zh-CN.md).

The delivery follows the scientific-code format of [Lake_Microplastics_Analysis_System](https://github.com/NKUHuLab/Lake_Microplastics_Analysis_System): analysis scripts, selected actual inputs, a command-line demo, calculated outputs, expected results and installation/reproduction instructions.

[Standalone scientific-demo package](https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/download/v2.0.0-scientific-demo/MHW2CoralEcoSys-scientific-demo.zip) · [Target-repository contribution](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1) · [Successful Linux scientific reproduction](https://github.com/RugkeyPro/MHW2CoralEcoSys/actions/runs/37599587090).

## 1. System requirements

- Python **3.11 or 3.12**; the validated local environment uses Python 3.12.
- NumPy 2.3.5, pandas 2.3.3 and Matplotlib 3.10.7, pinned in `requirements.txt`.
- A normal CPU workstation. The standalone dataset is 1.2 MB compressed and has 29,746 reef cells; no GPU, Node.js, R, database or map service is needed for this demo.
- Windows and Linux are checked through the local run and GitHub CI. Other operating systems are not independently tested.

## 2. Installation

```bash
git clone https://github.com/NKUHuLab/MHW2CoralEcoSys.git
cd MHW2CoralEcoSys
python -m venv .venv
```

Activate the environment:

```powershell
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

```bash
# Linux / macOS
source .venv/bin/activate
```

Then install:

```bash
pip install -r requirements.txt
```

Before the contribution is merged into the target repository, use [PR #1](https://github.com/NKUHuLab/MHW2CoralEcoSys/pull/1) or its contribution branch rather than the target's unchanged `main`.

## 3. Demo: actual reef-cell inputs to new analysis results

Run from the repository root:

```bash
python run_demo.py
```

The equivalent module entry is:

```bash
python analysis/modules/conservation_demo/code/run_all.py
```

This is an actual numerical analysis. It **does not use the frozen result tables as calculation inputs**. Starting with the bundled 29,746 real reef cells and nine climate-member support fields, it:

1. Validates coordinates, area units, support masks, MPEI fields and source hashes.
2. Computes habitat-only and exposure-aware equal-area reef priority selections with the original tie-handling rule.
3. Calculates support-weighted future MPEI, representation and reallocated area for three SSPs and four overlapping reporting regions.
4. Runs the original **1,000 joint climate-member / 5° spatial-block bootstrap replicates**, with seed **20260804**.
5. Derives percentile intervals and a new Table S5 from the newly calculated draws.
6. Compares the new point estimates, all 12,000 draw records, interval summaries and Table S5 to frozen accepted references, then creates scientific figures.

The seven calculation functions are copied unchanged from the accepted September 6 parent workflow. The portable wrapper supplies bundled paths and validation; it does not change support thresholds, score normalization, the priority rule or bootstrap design.

The complete validated Windows run took **29.9 seconds** on the development machine after installation. Runtime depends on the machine; installation is separate. Its 275,631 finite bootstrap values matched the recorded accepted values exactly. Full validation details are included under `docs/` and `examples/expected_outputs/`.

### Generated outputs

All new outputs go to `outputs/conservation_demo/`:

| Output | Content |
| --- | --- |
| `point_estimates.csv` | 12 newly calculated region × scenario estimates |
| `bootstrap_draws.csv.gz` | 12,000 newly calculated records, including legitimate missing outcomes |
| `summary.csv` | Estimates, bootstrap percentile intervals and effective replicate counts |
| `Table_S5_recomputed.csv` | Compact two-metric conservation table |
| `selected_reef_cells.csv.gz` | Actual candidate-cell selections and fractional reallocation, for each scenario |
| `conservation_comparison.png/.svg/.pdf` | Recomputed point estimates and 95% intervals |
| `reef_priority_reallocation.png` | Actual-cell selection changes under SSP5–8.5 |
| `run_checks.json` | Hash checks, record counts, numerical differences, effective draws, runtime and environment |

See [the module guide](analysis/modules/conservation_demo/README.md) for expected values and numerical comparison rules. A shortened operational check is available as `python run_demo.py --replicates 25 --no-plots`; it is labelled a shortened run and cannot reproduce the full published uncertainty intervals.

## 4. Research code and module organization

```text
run_demo.py
requirements.txt
analysis/modules/conservation_demo/
    code/conservation_core.py    Original numerical functions and constants
    code/run_all.py              Portable input-to-analysis-to-validation workflow
    data/reef_cells.csv.gz       Actual numeric/boolean reef-cell inputs
    reference/                  Accepted frozen results, used only for checks
    provenance.json             Original lineage, selected fields and SHA-256
research_workflow/              Unchanged larger parent analysis scripts
data/                           Dataset index and field dictionary
docs/                           Methods, validation and run instructions
examples/expected_outputs/      Captured results of the verified full demonstration
scripts/                        Source curation and scientific-package generation
tests/                          Input, selection, source-code and regression checks
archive/web_demo_20261006/       Earlier web presentation, retained as an archive
```

The parent scripts in `research_workflow/` include repaired-habitat preparation, period-matched MPEI conservation and future-risk analysis. Their original external scientific inputs and path assumptions remain explicit. They are source/provenance material and are not claimed to run against the small demo dataset. See [parent-workflow requirements](research_workflow/README.md).

## 5. Use on your own data

```bash
python run_demo.py --input /path/to/reef_cells.csv.gz --output /path/to/results
```

Use the [input dictionary](data/README.md). Keep one unique reef-cell record, reef area in km², MPA fractions in [0,1], member-level OR10/MESS support as 0/1/NaN, and MPEI in kg m⁻³. Preserve the global-reef MPEI normalization and the exact input row order when attempting frozen-draw reproduction. Candidate cells have complete member support and positive support fraction.

For genuinely different input data, the physical and selection checks still run, but frozen accepted-result comparisons are marked inapplicable. Altering `--seed` or `--replicates` is also recorded and cannot be presented as the original uncertainty calculation.

## 6. Validation and scientific scope

Run the unit/regression checks:

```bash
python -m unittest discover -s tests
```

The full demo separately verifies the source hash manifest, fractional area-budget closure, unchanged habitat-support capture, original point estimates, every compatible frozen bootstrap record and the accepted table. GitHub Actions runs the same scientific calculation on Linux and saves its generated outputs as an artifact. Read `run_checks.json`; process exit alone does not establish a successful scientific comparison.

This reproduces the **accepted conservation module from archived derived reef-cell inputs**, not the entire original ocean model, environmental preparation or MaxEnt training. The MPEI field is the shared 2045–2055 36-tracer mean. It is external modelled exposure, not internal dose or toxicological effect. The output compares prioritization schemes; it does not measure pollution removal or realized MPA success. The original global MPEI-change intervals include zero; GBR has 991 effective draws in the full 1,000-replicate design.

Underlying original project migration gaps and manuscript inconsistencies are not resolved by this module. All 29,746 archived positive-area cells are retained without result-dependent filtering. WKB geometry and unused component columns are omitted from the demo; exact numeric round-trip and source hashes are recorded.

See [methods and validation](docs/methods_and_validation.md), [dataset index](data/Dataset_index.csv), [module provenance](analysis/modules/conservation_demo/provenance.json) and [verified expected outputs](examples/expected_outputs/).

## Attribution and license

Software is under the [MIT License](LICENSE). The sample contains derived project inputs; original reef/MPA and other third-party datasets retain their original reuse conditions. This subset does not include their full raw boundary geometries. Cite this repository, commit and the accepted September 6 conservation lineage when reusing the demonstration, and consult the research team for the full study and upstream sources.
