# Equal-area reef conservation: runnable scientific example

Run from the repository root after `pip install -r requirements.txt`:

```bash
python analysis/modules/conservation_demo/code/run_all.py
```

`python run_demo.py` is the same entry. No webpage or JavaScript runtime is required.

## Calculation inputs

`data/reef_cells.csv.gz` contains 29,746 real reef cells and 25 fields. The numeric subset comes from the accepted September 6 `future_mpei_mean_2045_2055_reef_cells_20260903.parquet` and `repaired_member_support.parquet`.

All input rows and their original order are retained. Coordinates, area, MPA coverage, period-matched MPEI, original global min-max MPEI score, scenario support validity and all nine individual climate-member masks are included. WKB geometry and unused polymer/size component fields are omitted. Floating values were written with 17 significant digits and re-read using pandas `round_trip`; every retained column matched the original source exactly.

`reference/` holds the accepted point estimates, all 12,000 bootstrap records, the interval summary and compact Table S5. These files never provide numerical inputs to selection or resampling. They are used after new results are computed.

## Unchanged scientific rule

The equal-area budget is the area integral of existing MPA fractions within the candidate domain. Habitat-only priority ranks by climate-model support `P`. Exposure-aware priority ranks by `P + 0.0001 × E*`, retaining the support tier as the primary criterion. `E*` is the original global-reef min-max-normalized future MPEI, not a newly normalized regional score.

When the budget cuts a tied score group, all cells tied at that boundary receive the same fractional selection. The original `exact_budget_selection` function preserves the exact area budget. A higher selected-area MPEI under the secondary criterion describes exposure prioritization, not a reduction of pollution.

Three climate members are resampled with replacement. Five-degree spatial blocks are sampled with replacement independently within each reporting region. Each spatial multiplicity changes cell area, not the underlying cell concentration. The RNG seed is 20260804, default repetitions are 1,000, and the original region iteration and input order are retained to reproduce the frozen draw sequence.

## Expected results

The full run should produce 12 point estimates and 12,000 draw records. Under SSP1–2.6, global reallocated priority area is approximately **16.0193549244%**, and the selected-area mean-MPEI change is approximately **6.5484043952%**. Under SSP5–8.5 the corresponding global values are **2.1138870013%** and **2.3635300482%**.

Global, Caribbean and Southeast Asia have 1,000 effective replicates for the limiting metrics. GBR has 991 effective replicates because some resamples cannot define the required exposure-weighted denominator. Missing results are retained as missing and excluded only from the relevant metric's percentile calculation, exactly as in the original workflow.

Frozen keys, missingness and all compatible numeric fields are compared. Concentrations use absolute tolerance 1e-18 kg m⁻³; other metrics use absolute tolerance 1e-8 in their recorded units, with relative tolerance 1e-9. `run_checks.json` reports per-field maximum errors and actual finite comparison counts. These tolerances account for recorded decimal and floating-point differences; they are not applied to identifiers or source hashes.

## Outputs and configuration

Default output directory: repository `outputs/conservation_demo/`. The full run creates tables, selections, plots and an explicit validation report. Source files are read-only. `--replicates`, `--seed`, `--input`, `--output` and `--no-plots` are documented CLI controls.

The first N draw records can be compared to the original sequence for an unchanged input and seed with N≤1,000. Only the full unchanged design compares the final uncertainty summary and Table S5 to the accepted references. A different input or seed is clearly labelled a customized calculation; it is not silently declared a reproduction of the frozen study.

## Scientific limits

This example begins at derived reef-cell inputs. It does not regenerate the parent NetCDF/HSI/MESS fields or rerun the full ocean and SDM models. It retains the original 36-tracer 2045–2055 shared MPEI trajectory. The regions overlap; global and regional totals must not be added together. The point estimates and bootstrap intervals measure priority redistribution and external modelled exposure, not realized protection, internal dose or toxicological effect.
