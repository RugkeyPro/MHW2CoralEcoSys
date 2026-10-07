# Original parent research workflow

These four scripts are copied unchanged from the accepted project. Byte hashes are recorded in the conservation-demo provenance manifest. They retain original scientific algorithms and original full-project path assumptions.

| Script | Role | Full upstream inputs needed |
| --- | --- | --- |
| `prepare_figureS6_repaired_20260821.py` | Repaired HSI/MESS member support and equal-area habitat configuration | Reef/MPA area grid, OR10 thresholds, three members' repaired HSI/MESS rasters, older static-exposure support |
| `process_future_mpei_conservation_20260903.py` | Map the future MPEI field, calculate priorities, bootstrap and summarize | Period-matched 36-tracer NetCDF and reef-cell inputs |
| `run_repaired_analysis_20260906.py` | Join repaired habitat with accepted 2045–2055 future MPEI and verify promotion | Both previous modules and the complete original project trees / comparison tables |
| `analyze_future_risk_period_matched_mpei_20260902.py` | Period-matched future habitat quality, MPEI burden/pressure and sensitivity | Historical/future HSI, three-member scenarios and period-matched MPEI fields |

The standalone `conservation_demo` begins after these modules have generated the archived reef-cell inputs. It copies their seven independent numerical functions unchanged and supplies all required numeric inputs. It does not call external NetCDF, raster or full project paths.

To execute the larger originals, supply the complete upstream files and update their explicit path variables to a verified local full-project layout. NumPy/pandas are not sufficient for those original IO stages: they additionally use rasterio, netCDF4, SciPy and related project modules. No installation of those optional dependencies makes the missing upstream data appear. They are not part of the standalone demo's default command.
