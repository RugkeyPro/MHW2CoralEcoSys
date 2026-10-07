import { useMemo, useState } from 'react';
import { geoNaturalEarth1, geoPath } from 'd3-geo';
import { feature } from 'topojson-client';
import world from 'world-atlas/countries-110m.json';
import { mapColor, formatPct } from './lib.mjs';
import type { Metric, MapCell } from './types';

const labels = { Caribbean: 'Caribbean', GBR: 'Great Barrier Reef', Southeast_Asia: 'Southeast Asia' };
const coordinates: Record<string, [number, number]> = { Caribbean: [-73, 18], GBR: [147, -19], Southeast_Asia: [115, 7] };
export function ReefMap({ cells, region, onRegion, english }: { cells: MapCell[]; region: string; onRegion: (r: string) => void; english: boolean }) {
  const [hover, setHover] = useState<MapCell | null>(null);
  const projection = useMemo(() => geoNaturalEarth1().scale(156).translate([465, 240]), []);
  const path = useMemo(() => geoPath(projection), [projection]);
  const land = useMemo(() => feature(world as never, world.objects.countries as never), []);
  const features = (land as unknown as { features: object[] }).features;
  return <div className="map-wrap">
    <svg viewBox="0 0 930 470" className="reef-map" role="img" aria-label={english ? 'Actual habitat quality reduction raster, aggregated to one-degree display cells' : '真实栖息地质量降低栅格，一度显示网格'}>
      <defs><pattern id="map-grid" width="52" height="52" patternUnits="userSpaceOnUse"><path d="M 52 0 L 0 0 0 52" fill="none" stroke="#a3c7c3" strokeWidth=".5" opacity=".25"/></pattern></defs>
      <rect width="930" height="470" fill="#e7f0ed"/><rect width="930" height="470" fill="url(#map-grid)"/>
      <path d={path({ type: 'Sphere' }) ?? ''} fill="none" stroke="#d6e1dc"/>
      {features.map((f, i) => <path key={i} d={path(f as never) ?? ''} fill="#f8f8ee" stroke="#ccd7c9" strokeWidth=".65"/>)}
      {cells.map(cell => { const xy = projection([cell[0], cell[1]]); return xy && <circle key={`${cell[0]}:${cell[1]}`} cx={xy[0]} cy={xy[1]} r="1.8" fill={mapColor(cell[2])} opacity=".84" onMouseEnter={() => setHover(cell)} onMouseLeave={() => setHover(null)}><title>{cell[0]}°, {cell[1]}° · {formatPct(cell[2])} · n={cell[3]}</title></circle>; })}
      {Object.entries(coordinates).map(([name, coord]) => {const xy=projection(coord)!; const selected = region === name; return <g key={name} className="map-region" role="button" tabIndex={0} aria-label={name} onClick={() => onRegion(name)} onKeyDown={e => {if(e.key==='Enter' || e.key===' ') {e.preventDefault();onRegion(name);}}}><circle cx={xy[0]} cy={xy[1]} r={selected ? 12 : 8} fill="white" fillOpacity=".8" stroke={selected ? '#12685c' : '#4f8074'} strokeWidth={selected ? 3 : 1.5}/><circle cx={xy[0]} cy={xy[1]} r="2.6" fill="#12685c"/><text x={xy[0]+13} y={xy[1]-12} fill="#214e45" fontSize="11" fontWeight="600" paintOrder="stroke" stroke="#f5f9f1" strokeWidth="3">{labels[name as keyof typeof labels]}</text></g>;})}
      <text x="24" y="446" fontSize="10" fill="#648178">Natural Earth · EPSG:4326 source · 1° display bins</text>
    </svg>
    <div className="map-chip"><span className="pulse"/>{english ? 'Actual model raster' : '实际模型栅格'}<span>{cells.length.toLocaleString()} {english ? 'display cells' : '显示格点'}</span></div>
    <div className="map-legend"><div><span>{english ? 'Quality reduction (%)' : '质量降低（%）'}</span><span>{english ? 'negative = increase' : '负值表示质量提高'}</span></div><div className="legend-gradient"/><div className="legend-ticks"><span>−45</span><span>0</span><span>20</span><span>40</span><span>60</span><span>85</span></div></div>
    {hover && <div className="map-tooltip">{hover[0]}°, {hover[1]}° · {formatPct(hover[2])} · {hover[3]} {english ? 'source cells' : '原始格点'}</div>}
  </div>;
}

export function IntervalChart({ rows, signed = true, color = '#278071', label }: { rows: { name: string; metric: Metric; active?: boolean }[]; signed?: boolean; color?: string; label: string }) {
  const values = rows.flatMap(row => [row.metric.min, row.metric.max]);
  const low = Math.min(0, ...values), high = Math.max(1, ...values), span = high-low;
  const x = (value: number) => 135 + (value-low)/span * 330;
  const height = 52 + rows.length * 58;
  return <svg viewBox={`0 0 565 ${height}`} className="interval-chart" role="img" aria-label={label}>
    {[0, .25, .5, .75, 1].map(t => <g key={t}><line x1={135+t*330} y1="12" x2={135+t*330} y2={height-35} stroke="#e9eee8"/><text x={135+t*330} y={height-13} textAnchor="middle" fontSize="10" fill="#7e8c83">{Math.round(low+t*span)}%</text></g>)}
    <line x1={x(0)} y1="12" x2={x(0)} y2={height-35} stroke="#a6b8ad" strokeDasharray="3 4"/>
    {rows.map(({name,metric,active},i) => { const y=35+i*58; return <g key={name}><text x="0" y={y+4} fontSize="12" fill={active?'#124b40':'#647568'} fontWeight={active?700:400}>{name}</text><line x1={x(metric.min)} x2={x(metric.max)} y1={y} y2={y} stroke={color} strokeWidth="3" opacity={active?1:.6}/><line x1={x(metric.min)} x2={x(metric.min)} y1={y-5} y2={y+5} stroke={color}/><line x1={x(metric.max)} x2={x(metric.max)} y1={y-5} y2={y+5} stroke={color}/><circle cx={x(metric.mean)} cy={y} r={active?6:4.5} fill={color} stroke="white" strokeWidth="2"/><text x="490" y={y+4} fontSize="13" fill="#294737" fontWeight="600">{formatPct(metric.mean,signed)}</text></g>;})}
  </svg>;
}

export function ScenarioBars({ rows, selected, label }: { rows: { name: string; metric: Metric }[]; selected: number; label: string }) {
  const max = Math.max(...rows.map(r=>r.metric.max), 45);
  const y = (v: number) => 166-(v/max)*134;
  return <svg viewBox="0 0 435 221" className="bar-chart" role="img" aria-label={label}>
    {[0,10,20,30,40].map(v=><g key={v}><line x1="38" x2="425" y1={y(v)} y2={y(v)} stroke="#e8ede6"/><text x="24" y={y(v)+4} textAnchor="end" fontSize="10" fill="#819083">{v}</text></g>)}
    {rows.map((row,i)=>{const x=70+i*125;return <g key={row.name}><rect x={x} y={y(row.metric.mean)} width="61" height={166-y(row.metric.mean)} rx="5" fill={i===selected?'#277d69':'#b9cfbe'}/><line x1={x+30.5} x2={x+30.5} y1={y(row.metric.min)} y2={y(row.metric.max)} stroke="#36594b" strokeWidth="1.5"/><path d={`M${x+24},${y(row.metric.min)}h13 M${x+24},${y(row.metric.max)}h13`} stroke="#36594b"/><text x={x+30.5} y={y(row.metric.max)-11} textAnchor="middle" fontSize="13" fontWeight="600" fill="#284c3b">{formatPct(row.metric.mean)}</text><text x={x+30.5} y="191" textAnchor="middle" fontSize="11" fill="#657769">{row.name}</text></g>;})}
    <text x="3" y="14" fontSize="10" fill="#819083">%</text>
  </svg>;
}
