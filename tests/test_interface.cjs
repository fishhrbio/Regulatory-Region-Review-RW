// Render/event smoke tests using a minimal DOM adapter; no browser required.
// This checks generated markup and application behavior, not visual layout.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const html=fs.readFileSync('docs/index.html','utf8');
const payload=html.match(/<script id="atlas-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const nodes=new Map(),downloads=[];
function element(id){if(!nodes.has(id))nodes.set(id,{value:'',textContent:'',innerHTML:'',handlers:{},addEventListener(name,fn){this.handlers[name]=fn;},click(){downloads.push(this);}});return nodes.get(id);}
element('atlas-data').textContent=payload;
element('mask').value='primary';element('enhancer').value='strong_enhancer';element('start').value='20000000';element('end').value='20500000';
const blobs=[];
const context={document:{getElementById:element,createElement:()=>element('download-'+downloads.length)},Blob,URL:{createObjectURL:b=>{blobs.push(b);return 'blob:test';},revokeObjectURL(){}},setTimeout:fn=>fn(),console};
vm.createContext(context);vm.runInContext(fs.readFileSync('web/app.js','utf8'),context);
assert.match(element('metrics').innerHTML,/47.4%/);assert.match(element('metrics').innerHTML,/7.3%/);
assert.equal((element('matrix').innerHTML.match(/<button /g)||[]).length,225);
assert.equal(element('locus-error').textContent,'');assert.match(element('locus').innerHTML,/<rect /);
for(const mask of ['primary','retain_repetitive'])for(const enhancer of ['strong_enhancer','broad_enhancer']){
 element('mask').value=mask;element('enhancer').value=enhancer;element('mask').handlers.change();
 assert.match(element('insight').textContent,/less shared/);
 assert.match(element('coverage').innerHTML,/<svg /);
 assert.ok(!element('sharing').innerHTML.includes('NaN'));
}
element('matrix').handlers.click({target:{closest:()=>({dataset:{i:'0',j:'0'}})}});
assert.match(element('matrix-detail').textContent,/265,800 bp/);
element('start').value='10';element('end').value='5';element('show-locus').handlers.click();
assert.match(element('locus-error').textContent,/Use integer/);
assert.equal(element('locus').innerHTML,'');
element('start').value='0';element('end').value='100';element('show-locus').handlers.click();
assert.match(element('locus-summary').textContent,/undefined/);
element('download-csv').handlers.click();element('download-svg').handlers.click();
assert.equal(downloads.length,2);
(async()=>{const csv=await blobs[0].text(),svg=await blobs[1].text();assert.match(csv,/Robert Waite/);assert.match(csv,/retain_repetitive/);assert.match(svg,/M.S. candidate in Biotechnology/);assert.match(svg,/Broad enhancer/);console.log('Interface render/event tests passed: settings, 225 matrix cells, locus validation, empty region, CSV and SVG exports.');})();
