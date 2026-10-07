"""Package only the requested starting module, its actual inputs and executed notebook."""
import hashlib
import json
import zipfile
from pathlib import Path
import nbformat

ROOT = Path(__file__).resolve().parents[1]
notebook = nbformat.read(ROOT / 'Physiology_MaxEnt_demo.ipynb', as_version=4)
nbformat.validate(notebook)
cells = [cell for cell in notebook.cells if cell.cell_type == 'code']
if len(cells) != 7 or any(cell.execution_count is None for cell in cells):
    raise RuntimeError('All seven notebook stages must be executed')
if any(item.output_type == 'error' for cell in cells for item in cell.outputs):
    raise RuntimeError('Notebook execution errors')
checks = json.loads((ROOT / 'outputs/physiology_maxent_demo/run_checks.json').read_text())
audit = json.loads((ROOT / 'docs/maxent_input_validation.json').read_text())
if checks['status'] != 'passed' or audit['status'] != 'passed':
    raise RuntimeError('Execution and independent input checks are required')
module = ROOT / 'analysis/modules/physiology_maxent_demo'
paths = [ROOT / name for name in ['README.md', 'README.zh-CN.md', 'requirements.txt', 'LICENSE',
                                 '.gitattributes',
                                 'Physiology_MaxEnt_demo.ipynb', 'docs/maxent_input_validation.json',
                                 'scripts/execute_notebook.py', 'scripts/curate_maxent_demo.py',
                                 'scripts/create_maxent_notebook.py', 'scripts/verify_maxent_inputs.py',
                                 'scripts/package_maxent_demo.py', 'tests/test_physiology_maxent.py']]
paths += [path for path in module.rglob('*') if path.is_file() and '__pycache__' not in path.parts]
paths += [ROOT / 'examples/expected_outputs/maxent_run_checks.json',
          ROOT / 'examples/expected_outputs/maxent_scenario_summary.csv']
output = ROOT / 'deliverables'
output.mkdir(exist_ok=True)
destination = output / 'MHW2CoralEcoSys-physiology-maxent-demo.zip'
with zipfile.ZipFile(destination, 'w', zipfile.ZIP_DEFLATED) as archive:
    for path in paths:
        archive.write(path, 'MHW2CoralEcoSys/' + path.relative_to(ROOT).as_posix())
with zipfile.ZipFile(destination) as archive:
    if archive.testzip() is not None:
        raise ValueError('ZIP CRC test failed')
    manifest = json.loads((module / 'provenance.json').read_text(encoding='utf-8'))
    for record in manifest['packaged_files']:
        if hashlib.sha256(archive.read('MHW2CoralEcoSys/' + record['file'])).hexdigest() != record['sha256']:
            raise ValueError('Packaged source hash mismatch')
result = {'file': destination.name, 'bytes': destination.stat().st_size,
          'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
          'primary_entrypoint': 'Physiology_MaxEnt_demo.ipynb', 'zip_crc_test': 'passed',
          'native_maxent_training_executed': True, 'full_parent_reproduction': False}
(output / 'maxent_package_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
