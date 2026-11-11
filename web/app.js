/* Copyright (c) 2026 Robert Waite. MIT. */
"use strict";
const PALETTE = ["#c34154","#df7b88","#9766b0","#e49330","#d5a43a","#d5bc57","#dcd488","#5688b0","#468b78","#286b5b","#9bc6a4","#7f8091","#d3d7db","#e9e9eb","#f2f2f4"];
const CREDIT = "Robert Waite | M.S. candidate in Biotechnology | USF";
function locusStats(atoms,start,end,excluded,states){
  const blocked=new Set(excluded), selected=new Set(states);
  let eligible=0,gmOnly=0,shared=0,kOnly=0;
  for(const [a,b,g,k] of atoms){
    const bp=Math.max(0,Math.min(b,end)-Math.max(a,start));
    if(!bp||blocked.has(g)||blocked.has(k))continue;
    eligible+=bp;
    if(selected.has(g)&&selected.has(k))shared+=bp;
    else if(selected.has(g))gmOnly+=bp;
    else if(selected.has(k))kOnly+=bp;
  }
  const union=gmOnly+shared+kOnly;
  return {eligible,gmOnly,shared,kOnly,union,jaccard:union?shared/union:null};
}
function validRegion(start,end,maxEnd){return Number.isSafeInteger(start)&&Number.isSafeInteger(end)&&start>=0&&start<end&&end<=maxEnd;}
function table(headers,rows){return '<table><thead><tr>'+headers.map(h=>`<th scope="col">${h}</th>`).join('')+'</tr></thead><tbody>'+rows.map(r=>'<tr>'+r.map(v=>`<td>${v}</td>`).join('')+'</tr>').join('')+'</tbody></table>';}
if(typeof module!=="undefined")module.exports={locusStats,validRegion};
if(typeof document!=="undefined"){
const DATA=JSON.parse(document.getElementById("atlas-data").textContent);
const el=id=>document.getElementById(id), num=n=>n.toLocaleString("en-US"), pct=n=>n===null?"undefined":(100*n).toFixed(1)+"%";
const maxEnd=Math.max(...Object.values(DATA.qc).map(q=>q.selected_span_end));
function settings(){return {data:DATA.scenarios[el("mask").value],enhancer:el("enhancer").value};}
function render(){
 const {data,enhancer}=settings(),p=data.overlap.active_promoter,e=data.overlap[enhancer],delta=p.jaccard-e.jaccard;
 el("metrics").innerHTML=[['Eligible territory',(data.eligible_bp/1e6).toFixed(2)+' Mb','Common covered bases after the selected mask'],['Promoter Jaccard',pct(p.jaccard),'Shared / union · active promoter state 1'],['Enhancer Jaccard',pct(e.jaccard),'Shared / union · selected enhancer definition'],['Sharing gap',(delta*100).toFixed(1)+' pp','Promoter Jaccard minus enhancer Jaccard']].map(([a,b,c])=>`<div class="metric"><label>${a}</label><div class="value">${b}</div><p>${c}</p></div>`).join('');
 el("insight").textContent=`${enhancer==='strong_enhancer'?'Strong':'Broad'} enhancer annotations are ${delta>0?'less':delta<0?'more':'equally'} shared than active promoters on the selected territory. ${num(data.excluded_bp)} common-covered bp are excluded by this mask. This is a descriptive comparison of two cell lines.`;
 el("sharing").innerHTML=[["Active promoter · state 1",p],[enhancer==='strong_enhancer'?"Strong enhancer · states 4+5":"Broad enhancer · states 4–7",e]].map(([name,r])=>`<div class="barlabel"><strong>${name}</strong><span>J ${r.jaccard.toFixed(3)}</span></div><div class="bar" role="img" aria-label="${name}: ${pct(r.gm_only_bp/r.union_bp)} GM12878 only, ${pct(r.shared_bp/r.union_bp)} shared, ${pct(r.k_only_bp/r.union_bp)} K562 only">`+["gm_only_bp","shared_bp","k_only_bp"].map((k,i)=>`<div class="segment" style="width:${100*r[k]/r.union_bp}%;background:${['#276f83','#373d59','#da9b3a'][i]}" title="${num(r[k])} bp">${r[k]/r.union_bp>.07?pct(r[k]/r.union_bp):''}</div>`).join('')+`</div><div class="note">Union ${num(r.union_bp)} bp · shared ${num(r.shared_bp)} bp<br>Retained in other line: GM12878 ${pct(r.gm_retained)} · K562 ${pct(r.k_retained)}</div>`).join('');
 let coverage='<svg class="chart" viewBox="0 0 620 165" role="img" aria-label="Chromatin state coverage for both cell lines">';
 for(const [idx,key,name] of [[0,'gm_fraction','GM12878'],[1,'k_fraction','K562']]){
   coverage+=`<text x="0" y="${28+idx*65}" fill="#172e39" font-size="14">${name}</text>`;
   let x=100;
   for(const r of data.coverage){const width=r[key]*500;coverage+=`<rect x="${x}" y="${10+idx*65}" width="${width}" height="30" fill="${PALETTE[r.state-1]}"><title>State ${r.state}: ${pct(r[key])}, ${num(r[key==='gm_fraction'?'gm_bp':'k_bp'])} bp</title></rect>`;x+=width;}
 }
 for(const t of [0,25,50,75,100])coverage+=`<text x="${100+t*5}" y="135" text-anchor="middle" fill="#566b73" font-size="12">${t}%</text>`;
 el("coverage").innerHTML=coverage+'</svg>';
 el("state-key").innerHTML=DATA.state_names.map((name,i)=>`<span><i class="swatch" style="background:${PALETTE[i]}"></i>${i+1} ${name}</span>`).join('');
 let matrix='<span></span>'+Array.from({length:15},(_,j)=>`<span class="axis">${j+1}</span>`).join('');
 data.state_matrix_bp.forEach((row,i)=>{const total=row.reduce((a,b)=>a+b,0);matrix+=`<span class="axis row">${i+1}</span>`;row.forEach((bp,j)=>{const f=total?bp/total:null;const c=f===null?'#d3d7db':`rgb(${Math.round(244-213*f)},${Math.round(248-143*f)},${Math.round(248-126*f)})`;matrix+=`<button data-i="${i}" data-j="${j}" style="background:${c}" aria-label="GM12878 state ${i+1} to K562 state ${j+1}: ${num(bp)} bp, ${pct(f)} of row" title="${num(bp)} bp">${f!==null&&f>.1?'<span style="color:'+(f>.5?'white':'#172e39')+'">'+Math.round(f*100)+'</span>':''}</button>`;});});
 el("matrix").innerHTML=matrix;
 el("matrix-detail").textContent='Select a matrix cell to inspect its denominator.';
 const rows=[];
 for(const [name,s] of Object.entries(DATA.scenarios))for(const c of ['strong_enhancer','broad_enhancer'])rows.push([name==='primary'?'Exclude 14/15':'Retain 14/15',c==='strong_enhancer'?'Strong 4+5':'Broad 4–7',pct(s.overlap.active_promoter.jaccard),pct(s.overlap[c].jaccard),(100*(s.overlap.active_promoter.jaccard-s.overlap[c].jaccard)).toFixed(1)]);
 el("sensitivity").innerHTML=table(['Territory','Enhancer','Promoter J','Enhancer J','Gap (pp)'],rows);
 el("audit").innerHTML=table(['Cell line','Whole-file rows','chr22 rows','Covered bp','Outside common bp','Internal gap bp','Duplicates'],Object.entries(DATA.qc).map(([c,q])=>[c,num(q.whole_file_rows),num(q.selected_rows),num(q.selected_bp),num(q.selected_bp_outside_common),num(q.internal_gap_bp),q.duplicate_selected_rows]));
 drawLocus();
}
function drawLocus(){
 const start=Number(el('start').value),end=Number(el('end').value);
 if(!el('start').value.trim()||!el('end').value.trim()||!validRegion(start,end,maxEnd)){el('locus-error').textContent=`Use integer coordinates with 0 ≤ start < end ≤ ${num(maxEnd)}.`;el('locus').innerHTML='';el('locus-summary').textContent='';return;}
 el('locus-error').textContent='';
 const {data,enhancer}=settings(),blocked=new Set(data.excluded_states);
 let svg=`<svg class="chart" viewBox="0 0 1120 185" role="img" aria-label="Exact hg19 chr22 state annotations from ${start} to ${end}"><rect x="110" y="22" width="980" height="105" fill="#fff" stroke="#ccd6d6"/><text x="0" y="56" font-size="15" fill="#172e39">GM12878</text><text x="0" y="111" font-size="15" fill="#172e39">K562</text>`;
 for(const [a,b,g,k] of DATA.atoms){const left=Math.max(a,start),right=Math.min(b,end);if(left>=right||blocked.has(g)||blocked.has(k))continue;const x=110+(left-start)/(end-start)*980,w=(right-left)/(end-start)*980;for(const [s,y] of [[g,30],[k,85]])svg+=`<rect x="${x}" y="${y}" width="${w}" height="34" fill="${PALETTE[s-1]}"><title>chr22:${left}-${right} [0-based), state ${s}: ${DATA.state_names[s-1]}</title></rect>`;}
 for(let i=0;i<=4;i++)svg+=`<text x="${110+i*245}" y="151" text-anchor="middle" font-size="12" fill="#566b73">${num(Math.round(start+i*(end-start)/4))}</text>`;
 el('locus').innerHTML=svg+'</svg>';
 const e=locusStats(DATA.atoms,start,end,data.excluded_states,DATA.config.classes[enhancer]),p=locusStats(DATA.atoms,start,end,data.excluded_states,DATA.config.classes.active_promoter);
 el('locus-summary').textContent=`Region width ${num(end-start)} bp · eligible ${num(e.eligible)} bp · promoter J ${pct(p.jaccard)} (union ${num(p.union)} bp) · enhancer J ${pct(e.jaccard)} (union ${num(e.union)} bp). Undefined means the class union is empty.`;
}
function download(content,type,name){const url=URL.createObjectURL(new Blob([content],{type}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
el('mask').addEventListener('change',render);el('enhancer').addEventListener('change',render);el('show-locus').addEventListener('click',drawLocus);
for(const [id,direction] of [['previous',-1],['next',1]])el(id).addEventListener('click',()=>{const a=Number(el('start').value),b=Number(el('end').value);if(!validRegion(a,b,maxEnd)){drawLocus();return;}const width=b-a,newStart=Math.max(0,Math.min(maxEnd-width,a+direction*width));el('start').value=newStart;el('end').value=newStart+width;drawLocus();});
el('matrix').addEventListener('click',event=>{const b=event.target.closest('button');if(!b)return;const {data}=settings(),i=Number(b.dataset.i),j=Number(b.dataset.j),bp=data.state_matrix_bp[i][j],total=data.state_matrix_bp[i].reduce((a,b)=>a+b,0);el('matrix-detail').textContent=`GM12878 ${i+1} (${DATA.state_names[i]}) → K562 ${j+1} (${DATA.state_names[j]}): ${num(bp)} bp / ${num(total)} row bp = ${total?pct(bp/total):'undefined'}.`;});
el('download-csv').addEventListener('click',()=>{const {data,enhancer}=settings();let csv='author,assembly,chromosome,mask,class,eligible_bp,gm_only_bp,shared_bp,k_only_bp,union_bp,jaccard\n';for(const name of ['active_promoter',enhancer]){const r=data.overlap[name];csv+=[CREDIT,'hg19','chr22',el('mask').value,name,data.eligible_bp,r.gm_only_bp,r.shared_bp,r.k_only_bp,r.union_bp,r.jaccard].join(',')+'\n';}download(csv,'text/csv','waite-chromatin-selected-results.csv');});
el('download-svg').addEventListener('click',()=>{const {data,enhancer}=settings();let svg='<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="360" viewBox="0 0 1000 360"><rect width="1000" height="360" fill="white"/><g font-family="sans-serif" fill="#172e39"><text x="40" y="42" font-size="25">Waite Chromatin Atlas | Regulatory sharing</text><text x="40" y="72" font-size="14">hg19 chr22 | '+num(data.eligible_bp)+' eligible bp | '+el('mask').value+'</text>';
[['active_promoter','Active promoter (1)'],[enhancer,enhancer==='strong_enhancer'?'Strong enhancer (4+5)':'Broad enhancer (4–7)']].forEach(([key,label],i)=>{const r=data.overlap[key],y=120+i*86;svg+=`<text x="40" y="${y-12}" font-size="15">${label} | J = ${r.jaccard.toFixed(3)} | union ${num(r.union_bp)} bp</text>`;let x=40;for(const [j,k] of ['gm_only_bp','shared_bp','k_only_bp'].entries()){const width=920*r[k]/r.union_bp;svg+=`<rect x="${x}" y="${y}" width="${width}" height="40" fill="${['#276f83','#373d59','#da9b3a'][j]}"/>`;if(width>65)svg+=`<text x="${x+width/2}" y="${y+25}" text-anchor="middle" fill="white" font-size="14">${pct(r[k]/r.union_bp)}</text>`;x+=width;}});
svg+='<text x="40" y="288" font-size="13">GM12878 only (teal) · shared (navy) · K562 only (gold). Each bar uses its own class union.</text><text x="40" y="312" font-size="12">Published ENCODE/Broad ChromHMM annotations; descriptive overlap, not functional validation.</text><text x="40" y="340" font-size="13">'+CREDIT+'</text></g></svg>';download(svg,'image/svg+xml','waite-chromatin-sharing.svg');});
render();
}
