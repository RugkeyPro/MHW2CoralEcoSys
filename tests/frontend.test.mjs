import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { findResult, formatPct, intervalIncludesZero, rowsToCsv, mapColor } from '../src/lib.mjs';
const data=JSON.parse(readFileSync(new URL('../public/data/demo.json',import.meta.url)));
test('Published numbers match the accepted global SSP585 result',()=>{
  const row=findResult(data.future,'ssp585_2050','Global');
  assert.equal(formatPct(row.quality.mean),'31.6%');
  assert.equal(formatPct(row.pressure.mean,true),'+47.9%');
  assert.equal(row.members.length,3);
});
test('Region changes expose actual GBR pressure, not a global rescaling',()=>{
  assert.equal(formatPct(findResult(data.future,'ssp126_2050','GBR').pressure.mean,true),'+209.6%');
  assert.throws(()=>findResult(data.future,'unknown','Global'));
});
test('Scenario changes are different records',()=>{
  const a=findResult(data.future,'ssp126_2050','Southeast_Asia');
  const b=findResult(data.future,'ssp585_2050','Southeast_Asia');
  assert.notEqual(a.pressure.mean,b.pressure.mean);
  assert.equal(formatPct(a.quality.mean),'30.3%');
});
test('Signed roundoff does not become false exclusion of zero',()=>{
  assert.equal(intervalIncludesZero(-1.8e-14,32.1),true);
  assert.equal(intervalIncludesZero(0.029,2.92),false);
});
test('Exports preserve commas, quotes, zero and missingness',()=>{
  assert.equal(rowsToCsv([{label:'A,"B"',value:0},{label:'C',value:null}]),'"label","value"\r\n"A,""B""","0"\r\n"C",""');
});
test('Negative map values retain their own display class',()=>{
  assert.notEqual(mapColor(-40),mapColor(40));
});
