// Independent JS check against Python-derived whole-chromosome results.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const {locusStats,validRegion}=require('../web/app.js');
const html=fs.readFileSync('docs/index.html','utf8');
const data=JSON.parse(html.match(/<script id="atlas-data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
for(const scenario of Object.values(data.scenarios))for(const [name,states] of Object.entries(data.config.classes)){
 const actual=locusStats(data.atoms,0,1e9,scenario.excluded_states,states),expected=scenario.overlap[name];
 assert.equal(actual.eligible,scenario.eligible_bp);
 assert.equal(actual.gmOnly,expected.gm_only_bp);assert.equal(actual.shared,expected.shared_bp);assert.equal(actual.kOnly,expected.k_only_bp);assert.equal(actual.jaccard,expected.jaccard);
}
assert.equal(locusStats([[0,10,1,1]],20,30,[],[1]).jaccard,null);
assert.equal(locusStats([[0,10,4,5]],2,7,[],[4,5]).shared,5);
assert.equal(locusStats([[0,10,4,15]],2,7,[14,15],[4,5]).eligible,0);
assert.equal(validRegion(0,10,10),true);assert.equal(validRegion(0,11,10),false);assert.equal(validRegion(.5,10,20),false);
console.log('Dashboard calculations: six scenarios, clipping, masks, empty unions, coordinate validation passed.');
