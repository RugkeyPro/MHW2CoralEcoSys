"""Download and validate the explicitly documented complete input archive."""
import hashlib
import json
import tempfile
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
module = ROOT / 'analysis/modules/physiology_maxent_demo'
download = json.loads((module / 'full_grid_download.json').read_text())
manifest = json.loads((module / 'full_grid_sources.json').read_text(encoding='utf-8'))
destination = module / 'data/full_grid'
destination.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(prefix='mhw2coral-full-inputs-') as directory:
    archive_path = Path(directory) / 'full_grid_inputs.zip'
    print(f"Downloading complete source rasters: {download['bytes'] / 1024 ** 2:.1f} MiB", flush=True)
    value = hashlib.sha256()
    with urllib.request.urlopen(download['url'], timeout=60) as response, archive_path.open('wb') as stream:
        for block in iter(lambda: response.read(1024 * 1024), b''):
            value.update(block)
            stream.write(block)
    if value.hexdigest() != download['sha256']:
        raise ValueError('Downloaded archive SHA-256 mismatch')
    with zipfile.ZipFile(archive_path) as archive:
        # Only explicitly manifested paths are materialized; no unrestricted extractall.
        for record in manifest['files']:
            relative = Path(record['file'])
            if relative.is_absolute() or '..' in relative.parts:
                raise ValueError('Unsafe archive path')
            target = destination / relative
            payload = archive.read(relative.as_posix())
            if hashlib.sha256(payload).hexdigest() != record['sha256']:
                raise ValueError('Archived raster SHA-256 mismatch')
            if target.is_file() and hashlib.sha256(target.read_bytes()).hexdigest() != record['sha256']:
                raise ValueError(f'Refusing to replace a different existing file: {relative}')
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.is_file():
                target.write_bytes(payload)
    (destination / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(f"Verified {len(manifest['files'])} complete input rasters")
