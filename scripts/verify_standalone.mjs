import { chromium } from '@playwright/test';
import { readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
const url='http://127.0.0.1:8086/';
const manifest=JSON.parse(await readFile('data/provenance.json','utf8'));
for(const source of manifest.files){
  const response=await fetch(`${url}data/source/${source.file.split('/').at(-1)}`);
  if(!response.ok)throw new Error(`Download HTTP ${response.status}: ${source.file}`);
  const sha=createHash('sha256').update(Buffer.from(await response.arrayBuffer())).digest('hex');
  if(sha!==source.sha256)throw new Error(`HTTP source hash mismatch: ${source.file}`);
}
const browser=await chromium.launch({channel:process.env.CI?undefined:'msedge',headless:true});
const page=await browser.newPage({viewport:{width:1365,height:900},reducedMotion:'reduce'});
const externalRequests=[];
page.on('request',r=>{if(!r.url().startsWith(url)&&!r.url().startsWith('data:'))externalRequests.push(r.url());});
await page.goto(url);
await page.getByTestId('metric-pressure').waitFor();
if(await page.getByTestId('metric-pressure').textContent()!=='+47.9%')throw new Error('Standalone metric mismatch');
await page.waitForFunction(()=>document.querySelectorAll('.reef-map circle').length>2000);
if(externalRequests.length)throw new Error(`External runtime requests: ${externalRequests.join(',')}`);
await browser.close();
const report={standaloneBrowser:'passed',httpSourceHashMatches:manifest.files.length,externalRuntimeRequests:0,overviewPressure:'+47.9%'};
await writeFile('outputs/standalone_validation.json',JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify(report));
