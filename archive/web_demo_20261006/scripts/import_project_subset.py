"""Copy selected, unchanged project results and record their byte provenance."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = {
    'figure4_future_period_matched_mpei_20260902/source': [
        'figure4_period_matched_mpei_change_ensemble_20260902.csv',
        'figure4_period_matched_mpei_change_by_member_20260902.csv',
        'figure4_quality_reduction_ensemble_20260902.csv',
        'figure4_quality_reduction_by_member_20260902.csv',
        *[f'{s}_2050_ensemble_mean_quality_reduction_pct_20260902.tif' for s in ('ssp126', 'ssp245', 'ssp585')],
    ],
    'future_conservation_period_matched_mpei_20260903/tables': ['Supplementary_Table_S5_20260906.csv'],
    'figureS4_mpei_mhw_exposure_20260807/tables': ['regional_mhw_tercile_plot_data_20260807.csv'],
}

def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True, help='The accepted output/01_mainline directory')
    args = parser.parse_args()
    target = ROOT / 'data/source'
    target.mkdir(parents=True, exist_ok=True)
    records = []
    for folder, names in PACKAGES.items():
        for name in names:
            source = args.source_root / folder / name
            destination = target / name
            if destination.exists() and sha256(source) != sha256(destination):
                raise RuntimeError(f'Existing bundled source differs: {name}')
            shutil.copyfile(source, destination)
            record = {'file': f'data/source/{name}', 'project_relative_source': f'output/01_mainline/{folder}/{name}',
                      'sha256': sha256(source), 'bytes': source.stat().st_size,
                      'kind': 'derived_project_result', 'unchanged_copy': True}
            if source.suffix == '.csv':
                with source.open(encoding='utf-8-sig', newline='') as stream:
                    rows = list(csv.DictReader(stream))
                    record.update(rows=len(rows), columns=list(rows[0]) if rows else [])
            if sha256(destination) != record['sha256']:
                raise RuntimeError(f'Copy hash mismatch: {name}')
            records.append(record)
            browser_source = ROOT / 'public/data/source'
            browser_source.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(destination, browser_source / name)
    manifest = {'schema_version': 1, 'curated_date': '2026-10-06',
                'selection': 'Complete small accepted result tables; three accepted Figure 4 GeoTIFFs.',
                'historical_lineage': 'S4: static 54-tracer MPEI 1993–2019, HSI/MHW 1993–2022.',
                'future_lineage': 'Figure 4: period-matched 36-tracer s1–s4 MPEI; conservation repaired 2026-09-06.',
                'files': records}
    (ROOT / 'data/provenance.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
    print(f'Imported and SHA-256 checked {len(records)} real source files.')

if __name__ == '__main__':
    main()
