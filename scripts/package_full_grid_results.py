"""Package verified derived grids and maps without machine-specific runtime logs."""
import json
import zipfile
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.modules.physiology_maxent_demo.full_grid import digest

ROOT = Path(__file__).resolve().parents[1]
outputs = ROOT / 'outputs/physiology_maxent_demo'
checks = json.loads((outputs / 'run_checks.json').read_text())
audit = json.loads((ROOT / 'docs/full_grid_projection_validation.json').read_text())
if checks['status'] != 'passed' or not checks['full_spatial_projection'] or audit['status'] != 'passed':
    raise ValueError('Verified full outputs are required')
files = list((outputs / 'full_grid').glob('*/*.tif')) + list((outputs / 'full_grid').glob('*/projection_checks.json'))
if len([p for p in files if p.suffix == '.tif']) != 20:
    raise ValueError('All twenty full-grid products are required')
files += [outputs / 'run_checks.json', outputs / 'full_grid_summary.csv', outputs / 'model/Acropora.lambdas']
files += list(outputs.glob('full_grid_hsi_comparison.*'))
destination = ROOT / 'deliverables/full_grid_results.zip'
with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as archive:
    for path in files:
        archive.write(path, path.relative_to(outputs).as_posix())
    for filename in ['full_grid_projection_validation.json', 'mainline_physiology_audit.json']:
        archive.write(ROOT / 'docs' / filename, 'validation/' + filename)
with zipfile.ZipFile(destination) as archive:
    if archive.testzip() is not None:
        raise ValueError('Results ZIP CRC failed')
manifest = {'file': destination.name, 'bytes': destination.stat().st_size,
            'sha256': digest(destination), 'verified_native_grid_products': 20, 'zip_crc_test': 'passed'}
(ROOT / 'deliverables/full_grid_result_package_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps(manifest, indent=2))
