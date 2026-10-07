"""Package immutable full-grid inputs separately from the small Git repository."""
import json
import zipfile
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.modules.physiology_maxent_demo.full_grid import digest

ROOT = Path(__file__).resolve().parents[1]
folder = ROOT / 'analysis/modules/physiology_maxent_demo/data/full_grid'
manifest = json.loads((folder / 'manifest.json').read_text(encoding='utf-8'))
archive_path = ROOT / 'deliverables/full_grid_inputs.zip'
with zipfile.ZipFile(archive_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as archive:
    archive.write(folder / 'manifest.json', 'manifest.json')
    for record in manifest['files']:
        path = folder / record['file']
        if digest(path) != record['sha256']:
            raise ValueError('Full source hash mismatch')
        archive.write(path, record['file'])
with zipfile.ZipFile(archive_path) as archive:
    if archive.testzip() is not None:
        raise ValueError('Full input ZIP CRC failed')
result = {'file': archive_path.name, 'bytes': archive_path.stat().st_size,
          'sha256': digest(archive_path), 'files': len(manifest['files']),
          'zip_crc_test': 'passed', 'spatial_subsampling': False}
(ROOT / 'deliverables/full_grid_input_package_manifest.json').write_text(json.dumps(result, indent=2) + '\n')
(ROOT / 'analysis/modules/physiology_maxent_demo/full_grid_download.json').write_text(
    json.dumps(dict(result, url='https://github.com/RugkeyPro/MHW2CoralEcoSys/releases/download/full-grid-inputs-v1/full_grid_inputs.zip'), indent=2) + '\n')
print(json.dumps(result, indent=2))
