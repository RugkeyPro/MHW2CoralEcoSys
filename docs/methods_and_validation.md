# Methods and verification boundary

## Notebook entry and execution

`MHW2CoralEcoSys_demo.ipynb` is the primary demo. Its eight code cells separately read inputs, calculate estimates, bootstrap, summarize, compare, render and save. It does not launch the previous root script or read expected output images to stand in for fresh calculations. Outputs are saved in `outputs/notebook_demo/` and embedded in the executed Notebook. `scripts/execute_notebook.py` supplies a temporary kernelspec pointing to the currently validated Python interpreter, executes all cells without allowing errors, checks full-design results and package versions, and saves the executed artifact. GitHub CI executes this same notebook.

## Accepted model lineage

The sample follows the September 6 repaired conservation analysis: August repaired HSI/MESS support and the mean 2045–2055 36-tracer MPEI field. The support masks belong to three CMIP6 members and SSP126/245/585. All scenarios share the same future plastic field. Four taxa are represented in the upstream OR10/MESS constraints; the demo does not retrain those SDMs.

The sample includes every archived positive-area reef cell, all original runtime numeric/boolean fields and all member support masks. Source order is unchanged, global-reef MPEI normalization is retained and numeric serialization was checked for exact round-trip. This is a column subset of derived inputs, not a subset selected because it gives a desired result.

## Prioritization

The area budget is `sum(reef_area × existing_MPA_fraction)` on each valid candidate domain. The two priorities rank by model support `P` and by `P + 0.0001 × E*`. The latter adds exposure ranking within support tiers. Exact-score ties at a budget boundary are selected proportionally; the same partial fraction is assigned to every tied cell.

The numerical core is copied from the original accepted calculation, not rewritten from final result tables. Source AST comparison tests verify all seven function bodies against the unchanged archived original. Independent tests check budget closure, partial tie selection and primary-support retention.

Selected-area MPEI is weighted by reef area and model support, excluding unsupported MPEI values from its numerator and denominator. Reallocation is based on the overlap of fractional selections. These are configuration/exposure quantities, not measured plastic removal or realized habitat protection.

## Bootstrap

The original RNG seed is 20260804. Three climate members are sampled with replacement. Each region's five-degree blocks are also sampled with replacement; cell area is weighted by the block multiplicity. Default repetitions are 1,000, producing 12,000 region/scenario records. Input row order, region iteration order and support-mask ordering are kept identical to the original design.

Percentile intervals are computed on finite outcomes for each metric. Original missing outcomes remain missing. Global, Caribbean and Southeast Asia have 1,000 limiting effective replicates; GBR has 991. The global MPEI-change intervals include zero. Point estimates and percentile intervals are different objects; the wrapper does not force point estimates inside intervals or replace missing draws with zero.

## Reproduction checks

The runtime calculation reads `reef_cells.csv.gz` and the core functions. Frozen point/draw/summary/Table-S5 references are opened only after newly computed outputs exist. The runner verifies bundled file hashes, input schemas, area closure and preserved support capture, then aligns reference records by their exact keys and checks missingness and numeric values.

Concentration comparisons use `atol=1e-18 kg m−3`; other metrics use `atol=1e-8` in their stored units, with `rtol=1e-9`. Identifiers, row coverage, flags and source hashes are not approximate. Reports list each metric's maximum absolute difference and the actual finite comparison counts.

Shortened/customized runs remain usable analyses but are clearly labelled. Only the unchanged complete 1,000-replicate run compares its final interval summary and Table S5 to the accepted full result. No PASS marker is written when a check fails; the process raises an error and saves `status: failed`.

## Verification versus scientific validity

Successful reproduction establishes that this input/algorithm/configuration package produces the archived conservation calculations. It does not independently validate original HSI/MESS models, source emissions, spatial reef/MPA records or upstream ocean physics. Annual endpoint MPEI means are not within-year concentration means; the shared future model trajectory is not a separate SSP-specific plastic-emissions projection.

The complete project has other environmental/SDM inputs, model-environment reconstruction needs and prior migration/manuscript issues. Those are not closed by this portable module. The previously built webpage remains in a dated archive and is not the current scientific execution entry.
