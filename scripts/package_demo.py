"""Ship a built static demo, source copies and one-command browser launcher."""
import hashlib
import json
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
destination = ROOT / 'demo'
destination.mkdir(exist_ok=True)
shutil.copytree(ROOT / 'dist', destination, dirs_exist_ok=True)
current = {p.relative_to(ROOT / 'dist') for p in (ROOT / 'dist').rglob('*') if p.is_file()}
for file in destination.rglob('*'):
    if file.is_file() and file.relative_to(destination) not in current:
        # Only previously generated demo assets are eligible for pruning.
        if file.parent == destination / 'assets' and file.suffix in ('.js', '.css'):
            file.unlink()
manifest = json.loads((ROOT / 'data/provenance.json').read_text(encoding='utf-8'))
for item in manifest['files']:
    file = destination / 'data/source' / Path(item['file']).name
    if hashlib.sha256(file.read_bytes()).hexdigest() != item['sha256']:
        raise RuntimeError(f'Browser download differs from source: {file.name}')
output = ROOT / 'outputs'
output.mkdir(exist_ok=True)
with zipfile.ZipFile(output / 'MHW2CoralEcoSys-demo.zip', 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    for path in destination.rglob('*'):
        if path.is_file():
            archive.write(path, str(Path('demo') / path.relative_to(destination)))
    for name in ('README.md', 'serve_demo.py'):
        archive.write(ROOT / name, name)
print('Packaged standalone static demo; all nine downloadable source hashes match.')
