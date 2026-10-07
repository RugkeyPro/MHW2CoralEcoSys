"""Recompute the small demonstration from bundled results, without source rasters."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/source'
SCENARIOS = ('ssp126_2050', 'ssp245_2050', 'ssp585_2050')
REGIONS = ('Global', 'Caribbean', 'GBR', 'Southeast_Asia')
SCENARIO_NAMES = dict(zip(SCENARIOS, ('SSP1-2.6', 'SSP2-4.5', 'SSP5-8.5')))

def read_csv(name):
    with (SOURCE / name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def close(actual, expected, label):
    if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-8):
        raise ValueError(f'{label}: {actual} != {expected}')

def verify_hashes():
    manifest = json.loads((ROOT / 'data/provenance.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        file = ROOT / item['file']
        digest = hashlib.sha256(file.read_bytes()).hexdigest()
        if digest != item['sha256']:
            raise ValueError(f'Source hash mismatch: {item["file"]}')
    return manifest

def assemble():
    provenance = verify_hashes()
    exposure = read_csv('figure4_period_matched_mpei_change_ensemble_20260902.csv')
    quality = read_csv('figure4_quality_reduction_ensemble_20260902.csv')
    members = read_csv('figure4_period_matched_mpei_change_by_member_20260902.csv')
    quality_members = read_csv('figure4_quality_reduction_by_member_20260902.csv')
    conservation = read_csv('Supplementary_Table_S5_20260906.csv')
    historical = read_csv('regional_mhw_tercile_plot_data_20260807.csv')
    expected_keys = {(s, r) for s in SCENARIOS for r in REGIONS}
    def indexed(rows):
        result = {(r['scenario'], r['region']): r for r in rows}
        if len(result) != len(rows) or set(result) != expected_keys:
            raise ValueError('Missing or duplicate scenario-region source rows')
        return result
    eindex, qindex = indexed(exposure), indexed(quality)
    results = []
    checks = 0
    for scenario in SCENARIOS:
        for region in REGIONS:
            e, q = eindex[scenario, region], qindex[scenario, region]
            em = [r for r in members if (r['scenario'], r['region']) == (scenario, region)]
            qm = [r for r in quality_members if (r['scenario'], r['region']) == (scenario, region)]
            if len(em) != 3 or len(qm) != 3 or len({r['member'] for r in em}) != 3:
                raise ValueError('Every scenario-region must contain three distinct climate members')
            metrics = {}
            for short, column, member_rows, ensemble in (
                ('burden', 'dynamic_mpei_burden_change_pct', em, e),
                ('pressure', 'dynamic_mpei_pressure_change_pct', em, e),
                ('quality', 'continuous_quality_reduction_pct', qm, q),
            ):
                values = [float(r[column]) for r in member_rows]
                metrics[short] = {'mean': mean(values), 'min': min(values), 'max': max(values)}
                for suffix, actual in metrics[short].items():
                    close(actual, float(ensemble[f'{column}_{suffix}']), f'{scenario}/{region}/{column}/{suffix}')
                    checks += 1
            member_records = []
            for row in em:
                for kind in ('burden', 'pressure'):
                    past, future = float(row[f'historical_mpei_{kind}_' + ('kg_m3_km2' if kind == 'burden' else 'kg_m3')]), float(row[f'future_mpei_{kind}_' + ('kg_m3_km2' if kind == 'burden' else 'kg_m3')])
                    if past <= 0:
                        raise ValueError('Non-positive historical denominator')
                    close(100 * (future / past - 1), float(row[f'dynamic_mpei_{kind}_change_pct']), f'{kind} change')
                    checks += 1
                close(float(row['future_mpei_burden_kg_m3_km2']) / float(row['future_hsi_weighted_area_km2']), float(row['future_mpei_pressure_kg_m3']), 'burden / weighted habitat = pressure')
                checks += 1
                quality_row = next(x for x in qm if x['member'] == row['member'])
                close(100 * (1 - float(quality_row['future_continuous_quality_km2_hsi']) / float(quality_row['historical_continuous_quality_km2_hsi'])), float(quality_row['continuous_quality_reduction_pct']), 'quality reduction')
                checks += 1
                member_records.append({'member': row['member'], 'burden': float(row['dynamic_mpei_burden_change_pct']), 'pressure': float(row['dynamic_mpei_pressure_change_pct']), 'quality': float(quality_row['continuous_quality_reduction_pct'])})
            results.append({'scenario': scenario, 'region': region, **metrics, 'members': member_records,
                            'historicalPressure': float(e['historical_mpei_pressure_kg_m3_mean']),
                            'futurePressure': float(e['future_mpei_pressure_kg_m3_mean']),
                            'futureHabitat': float(e['future_habitat_area_km2_mean'])})
    cindex = {(r['scenario'], r['region'].replace(' ', '_')): r for r in conservation}
    if len(cindex) != 12:
        raise ValueError('Invalid conservation row count')
    conservation_results = []
    for scenario in SCENARIOS:
        for region in REGIONS:
            row = cindex[SCENARIO_NAMES[scenario], region]
            record = {'scenario': scenario, 'region': region}
            for metric, prefix, lo, hi in (
                ('reallocation', 'future_mpei_vs_climate_area_changed_pct', 'area_changed_lower_95', 'area_changed_upper_95'),
                ('gain', 'future_mpei_vs_climate_mean_mpei_change_pct', 'mpei_change_lower_95', 'mpei_change_upper_95'),
            ):
                value, lower, upper = (float(row[k]) for k in (prefix, lo, hi))
                if not all(math.isfinite(x) for x in (value, lower, upper)) or lower > upper:
                    raise ValueError('Invalid conservation interval')
                record[metric] = {'mean': value, 'min': lower, 'max': upper}
            conservation_results.append(record)
    history_results = []
    for row in historical:
        record = {'region': row['region'].replace(' ', '_'), 'taxon': row['response'], 'group': int(row['group_order']),
                  'exposure': float(row['mean_unit_hsi_weighted_mpei_index']), 'change': float(row['contrast_vs_low_pct']),
                  'min': float(row['ci95_low_pct']), 'max': float(row['ci95_high_pct'])}
        if not all(math.isfinite(record[k]) for k in ('exposure', 'change', 'min', 'max')):
            raise ValueError('Non-finite historical value')
        history_results.append(record)
    if len(history_results) != 60:
        raise ValueError('Expected 4 regions × 5 responses × 3 MHW groups')
    return {'schemaVersion': 1, 'lineageDate': '2026-09-07', 'curatedDate': '2026-10-06',
            'future': results, 'conservation': conservation_results, 'historical': history_results,
            'provenance': provenance, 'checks': {'sourceHashes': len(provenance['files']), 'numericComparisons': checks,
            'ensembleRows': len(results), 'conservationRows': len(conservation_results), 'historicalRows': len(history_results)}}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'outputs')
    parser.add_argument('--publish', action='store_true', help='Refresh the browser dataset')
    args = parser.parse_args()
    data = assemble()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / 'demo_results.json').write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    (args.output / 'validation.json').write_text(json.dumps(data['checks'], indent=2) + '\n', encoding='utf-8')
    if args.publish:
        target = ROOT / 'public/data'
        target.mkdir(parents=True, exist_ok=True)
        (target / 'demo.json').write_text(json.dumps(data, ensure_ascii=False) + '\n', encoding='utf-8')
    print(json.dumps(data['checks'], indent=2))

if __name__ == '__main__':
    main()
