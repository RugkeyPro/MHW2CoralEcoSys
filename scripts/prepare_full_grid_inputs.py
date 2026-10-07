"""Copy and hash full native projection grids; parent source files remain read-only."""
import argparse
import json
import shutil
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.modules.physiology_maxent_demo.full_grid import SCENARIOS, digest

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--project', type=Path, required=True)
args = parser.parse_args()
source = args.project / '解压整理/01_ecoMarine/02_data/ecoMarine_0907_20260913_02_data_004/filesystem/D/ecoMarine/data/processed/scenarios_shared_cmip6_ensemble'
destination = ROOT / 'analysis/modules/physiology_maxent_demo/data/full_grid'
records = []
for scenario in ['static'] + SCENARIOS:
    source_folder = source / ('static_regional' if scenario == 'static' else f'members/gfdl_esm4_r1i1p1f1_gr/{scenario}_regional')
    variables = ['bathymetry', 'rugosity'] if scenario == 'static' else ['bottomT', 'thetao', 'si', 'so', 'zos', 'omega_arag']
    for variable in variables:
        original = source_folder / f'{variable}.tif'
        target = destination / scenario / original.name
        target.parent.mkdir(parents=True, exist_ok=True)
        original_hash = digest(original)
        if not target.is_file() or digest(target) != original_hash:
            shutil.copyfile(original, target)
        if digest(target) != original_hash:
            raise ValueError('Full-grid input copy mismatch')
        records.append({'file': target.relative_to(destination).as_posix(),
                        'original_source': original.relative_to(args.project).as_posix(),
                        'bytes': target.stat().st_size, 'sha256': original_hash})
        print(f'Full input verified: {scenario}/{variable}', flush=True)
manifest = {'member': 'gfdl_esm4_r1i1p1f1_gr', 'spatial_subsampling': False, 'files': records}
(destination / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
(ROOT / 'analysis/modules/physiology_maxent_demo/full_grid_sources.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
print(f'Prepared {len(records)} verified full source rasters')
