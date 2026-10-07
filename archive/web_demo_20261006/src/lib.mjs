export const SCENARIOS = ['ssp126_2050', 'ssp245_2050', 'ssp585_2050'];
export const REGIONS = ['Global', 'Caribbean', 'GBR', 'Southeast_Asia'];
export const MEMBERS = ['gfdl_esm4_r1i1p1f1_gr', 'cnrm_esm2_1_r1i1p1f2_gr', 'ipsl_cm6a_lr_r1i1p1f1_gr'];

/** A numerical zero may be a small signed roundoff from a bootstrap quantile. */
export function intervalIncludesZero(lower, upper) {
  return lower <= 1e-10 && upper >= -1e-10;
}

/** Filter only; these region records overlap and must never be summed. */
export function findResult(rows, scenario, region) {
  const found = rows.filter(row => row.scenario === scenario && row.region === region);
  if (found.length !== 1) throw new Error(`Expected one source result for ${scenario}/${region}`);
  return found[0];
}

export function formatPct(value, signed = false) {
  if (!Number.isFinite(value)) return '—';
  const clean = Math.abs(value) < 0.05 ? 0 : value;
  return `${signed && clean > 0 ? '+' : ''}${clean.toFixed(1)}%`;
}

export function rowsToCsv(rows) {
  if (!rows.length) return '';
  const columns = Object.keys(rows[0]);
  const escape = value => `"${String(value ?? '').replaceAll('"', '""')}"`;
  return [columns.map(escape).join(','), ...rows.map(row => columns.map(c => escape(row[c])).join(','))].join('\r\n');
}

export function mapColor(value) {
  if (value < -5) return '#397f90';
  if (value < 0) return '#8eb7b1';
  if (value < 20) return '#dbb878';
  if (value < 40) return '#d68c57';
  if (value < 60) return '#c96243';
  return '#8f342f';
}
