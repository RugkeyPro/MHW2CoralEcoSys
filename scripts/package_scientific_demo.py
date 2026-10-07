"""Create a standalone scientific-code/data package, excluding the legacy website."""
import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
output = ROOT / 'deliverables'
output.mkdir(exist_ok=True)
paths = [ROOT / name for name in ('README.md','README.zh-CN.md','requirements.txt','run_demo.py','LICENSE','.gitattributes')]
for directory in ('analysis','research_workflow','data','docs','examples','scripts','tests'):
    paths.extend(path for path in (ROOT / directory).rglob('*') if path.is_file()
                 and '__pycache__' not in path.parts and path.suffix != '.pyc')
checks = json.loads((ROOT / 'outputs/conservation_demo/run_checks.json').read_text(encoding='utf-8'))
if checks.get('status') != 'passed' or not checks.get('full_frozen_design'):
    raise RuntimeError('A verified full scientific run is required before packaging')
archive_path = output / 'MHW2CoralEcoSys-scientific-demo.zip'
with zipfile.ZipFile(archive_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(set(paths)):
        archive.write(path, str(Path('MHW2CoralEcoSys') / path.relative_to(ROOT)))
with zipfile.ZipFile(archive_path) as archive:
    if archive.testzip() is not None:
        raise RuntimeError('Packaged ZIP CRC test failed')
    manifest = json.loads((ROOT / 'analysis/modules/conservation_demo/provenance.json').read_text(encoding='utf-8'))
    for record in manifest['files']:
        entry = f'MHW2CoralEcoSys/{record["file"]}'
        if hashlib.sha256(archive.read(entry)).hexdigest() != record['sha256']:
            raise RuntimeError(f'Packaged source hash mismatch: {entry}')
result = {'file': archive_path.name, 'sha256': hashlib.sha256(archive_path.read_bytes()).hexdigest(),
          'bytes': archive_path.stat().st_size, 'verified_full_analysis': True, 'zip_crc_test': 'passed',
          'legacy_website_included': False}
(output / 'package_manifest.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
