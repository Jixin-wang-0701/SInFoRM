import fs from 'node:fs/promises';
import fsSync from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {pathToFileURL} from 'node:url';

const ROOT='D:\\Phd_project\\revision\\model';
const HERE=path.join(ROOT,'outputs','pfam_rebuild_npnt');
const BUILD=path.join(HERE,'figure_build');
const SOURCE=path.join(ROOT,'figures','pptx','rebuilt_pfam_vs_sqtl_two_line_rank_enrichment_pp_gt_0p9_editable.pptx');
const FINAL=path.join(ROOT,'figures','pptx','rebuilt_pfam_six_panel_four_sided_editable.pptx');
const SKILL='C:\\Users\\18553\\.codex\\plugins\\cache\\openai-primary-runtime\\presentations\\26.1022.11436\\skills\\presentations';
const LIB='C:\\Users\\18553\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\node\\node_modules\\@oai\\artifact-tool\\dist\\artifact_tool.mjs';
const PY='C:\\Users\\18553\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe';
process.env.RUNTIME_NODE_MODULES='C:\\Users\\18553\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\node\\node_modules';
const {Presentation,PresentationFile}=await import(pathToFileURL(LIB).href);
const {finalizePresentation}=await import(pathToFileURL(path.join(SKILL,'container_tools','artifact_tool_utils.mjs')).href);
fs.realpath=async p=>fsSync.realpathSync(p);
fs.link=async (source,target)=>fs.copyFile(source,target,fsSync.constants.COPYFILE_EXCL);
function csvRows(s){const lines=s.trim().split(/\r?\n/),h=lines.shift().split(',');return lines.map(l=>Object.fromEntries(l.split(',').map((v,i)=>[h[i],v])));}
const [sqstats,curves,conditional,loco]=await Promise.all([
 'rank_enrichment_model_vs_sqtl_stats.csv','rank_enrichment_model_vs_sqtl_curves.csv','global_top2_figure_stats.csv','loco_five_power_p_sweep.csv'
].map(n=>fs.readFile(path.join(HERE,n),'utf8').then(csvRows)));
const P=Presentation.create({slideSize:{width:1800,height:1900}});
const BLACK='#181818',GRAY='#555555',LIGHT='#D7D7D7',RED='#BA4045',BLUE='#2E729C',GOLD='#B78A24',GREEN='#548064',SOFT='#F6F6F3',FONT='Times New Roman';
let slide;
const slideRefs=[];
function newSlide(){slide=P.slides.add();slide.background.fill='#FFFFFF';slideRefs.push(slide);return slide;}
function rect(x,y,w,h,fill='none',stroke='none',sw=0,name=''){return slide.shapes.add({geometry:'rect',name,position:{left:x,top:y,width:w,height:h},fill,line:{style:'solid',fill:stroke,width:sw}});}
function line(x1,y1,x2,y2,color=BLACK,w=2,name=''){
 if(Math.abs(x1-x2)<.01)return rect(x1-w/2,Math.min(y1,y2),w,Math.abs(y2-y1),color,'none',0,name);
 if(Math.abs(y1-y2)<.01)return rect(Math.min(x1,x2),y1-w/2,Math.abs(x2-x1),w,color,'none',0,name);
 return slide.shapes.add({geometry:'line',name,position:{left:Math.min(x1,x2),top:Math.min(y1,y2),width:Math.abs(x2-x1),height:Math.abs(y2-y1),verticalFlip:(x2-x1)*(y2-y1)<0},fill:'none',line:{style:'solid',fill:color,width:w}});
}
function txt(s,x,y,w,h,size=20,color=BLACK,bold=false,align='left',name=''){
 const sh=slide.shapes.add({geometry:'textbox',name,position:{left:x,top:y,width:w,height:h},fill:'none',line:{style:'solid',fill:'none',width:0}});
 sh.text.style={fontSize:size,color,bold,alignment:align,verticalAlignment:'middle',autoFit:'shrinkText',wrap:'none',insets:0};
 sh.text=[[{run:String(s),style:{typeface:FONT,fontSize:`${size}px`,bold,color}}]];sh.text.insets=0;return sh;
}
function dot(x,y,r,color){slide.shapes.add({geometry:'ellipse',position:{left:x-r,top:y-r,width:r*2,height:r*2},fill:color,line:{style:'solid',fill:color,width:1}});}
function pct(x,d=1){return `${(100*Number(x)).toFixed(d)}%`;}
function pv(x){const v=Number(x);return v<.001?v.toExponential(2):v.toFixed(v<.01?4:3);}
function key(endpoint){return endpoint==='ICGC_COPD_LTRC'?'ICGC COPD × LTRC lung':'FEV1/FVC × GTEx lung';}
function getSq(endpoint,method){return sqstats.find(r=>r.endpoint===endpoint&&+r.posterior_threshold===.9&&r.kind==='ranking'&&r.method===method);}
function getContrast(endpoint){return sqstats.find(r=>r.endpoint===endpoint&&+r.posterior_threshold===.9&&r.kind==='contrast'&&r.method==='primary_full_vs_lead_beta');}
function getCond(endpoint){return conditional.find(r=>r.endpoint===endpoint&&+r.posterior_threshold===.9&&r.panel==='C_D_within_MAGMA_global_top2');}
function thresholdLabel(p){return p===5e-4?'<5×10⁻⁴':p===5e-5?'<5×10⁻⁵':p===5e-6?'<5×10⁻⁶':p===5e-7?'<5×10⁻⁷':'<5×10⁻⁸';}
const ends=['ICGC_COPD_LTRC','FEV1FVC_GTEx'];
const cuts=[5e-4,5e-5,5e-6,5e-7,5e-8];
function getL(model,cut,block){const row=loco.find(r=>r.model===model&&+r.magma_p_cutoff===cut&&r.block===block);if(!row)throw Error(`missing ${model} ${cut} ${block}`);return row;}
function legendLine(x,y,c,label){line(x,y+12,x+48,y+12,c,5);txt(label,x+60,y-5,330,36,21);}

// Slide 1: the two curves and the two conditional MAGMA panels share a slide,
// while their distinct analysis universes are stated directly on the canvas.
newSlide();
txt('Colocalization enrichment: rebuilt Pfam model',75,16,1650,50,36,BLACK,false,'center');
txt('Expression + GO + rebuilt Pfam, posterior > 0.9',150,68,1500,32,23,GRAY,false,'center');
legendLine(500,108,RED,'Functional score');legendLine(975,108,BLUE,'Lead sQTL |β|');
function curvePanel(endpoint,letter,x0,x1){
 const model=getSq(endpoint,'rebuilt_full'),base=getSq(endpoint,'lead_sqtl_abs_beta'),test=getContrast(endpoint);
 txt(`${letter}   ${key(endpoint)}`,x0-28,158,x1-x0+90,38,27,BLACK,true);
 txt(`Common set: ${(+model.common_n).toLocaleString()} junctions, ${model.positive_n} hits`,x0-28,198,x1-x0+90,29,19,GRAY);
 const top=254,bottom=556,w=x1-x0,h=bottom-top;
 for(const tick of [0,.25,.5,.75,1]){const y=bottom-h*tick;line(x0,y,x1,y,tick===0?BLACK:LIGHT,tick===0?2:1);txt(`${Math.round(tick*100)}%`,x0-63,y-15,53,30,16,GRAY,false,'right');}
 for(const tick of [0,.25,.5,.75,1]){const x=x0+w*tick;line(x,bottom,x,bottom+7,BLACK,1.4);txt(`${Math.round(tick*100)}%`,x-29,bottom+10,58,27,16,GRAY,false,'center');}
 line(x0,top,x0,bottom,BLACK,2);
 line(x0,top,x1,top,BLACK,2);line(x1,top,x1,bottom,BLACK,2);
 for(const [method,color] of [['lead_sqtl_abs_beta',BLUE],['rebuilt_full',RED]]){
  const k=curves.filter(r=>r.endpoint===endpoint&&+r.posterior_threshold===.9&&r.method===method);
  for(let i=1;i<k.length;i++)line(x0+w*(+k[i-1].rank_fraction),bottom-h*(+k[i-1].cumulative_recall),x0+w*(+k[i].rank_fraction),bottom-h*(+k[i].cumulative_recall),color,3.2);
 }
 txt('Fraction of ranked junctions',x0+90,602,w-180,27,18,GRAY,false,'center');
 txt(`Curve area: ${(+model.cumulative_area).toFixed(3)} vs ${(+base.cumulative_area).toFixed(3)}`,x0-28,643,w+55,28,20);
 txt(`AUPRC: ${(+model.auprc).toFixed(4)} vs ${(+base.auprc).toFixed(4)}`,x0-28,677,w+55,28,20);
 txt(`Regional P=${pv(test.block_signflip_p_one_sided)} (Holm ${pv(test.holm_p_two_endpoints)})`,x0-28,711,w+55,28,18,GRAY);
}
curvePanel(ends[0],'A',105,800);curvePanel(ends[1],'B',1000,1695);
txt('Cumulative share of colocalization hits',42,229,690,25,17,GRAY);
line(70,768,1730,768,LIGHT,2);
txt('C   Within global MAGMA top 2%: hit rate',58,783,780,43,27,BLACK,true);
txt('D   Conditional odds ratio within MAGMA top 2%',950,783,820,43,27,BLACK,true);
txt('MAGMA and model top 2% each select 1,915 of 95,750 junctions globally',60,826,1630,28,18,GRAY);
// C: direct conditional counts, red is model-selected intersection; blue is other MAGMA top 2%.
rect(180,862,18,18,RED);txt('Model ∩ MAGMA',205,857,300,29,19);
rect(520,862,18,18,BLUE);txt('Other MAGMA top 2%',545,857,330,29,19);
const rateX=260,rateW=565,rateMax=.13;
// Match the four-sided plotting frames in the reference Figures 2 and 4.
for(const [i,e] of ends.entries()){
 const r=getCond(e),y=920+i*120,a=+r.model_positive_n/+r.model_selected_n,b=+r.other_magma_positive_n/+r.other_magma_n;
 txt(key(e),55,y-12,200,30,19);
 rect(rateX,y,rateW*a/rateMax,24,RED);rect(rateX,y+34,rateW*b/rateMax,24,BLUE);
 txt(`${pct(a)} (${r.model_positive_n}/${r.model_selected_n})`,rateX+rateW*a/rateMax+8,y-2,250,28,18);
 txt(`${pct(b)} (${Math.round(+r.other_magma_positive_n)}/${Math.round(+r.other_magma_n)})`,rateX+rateW*b/rateMax+8,y+32,250,28,18);
}
line(rateX,900,890,900,BLACK,2);line(890,900,890,1115,BLACK,2);
line(rateX,1115,890,1115,BLACK,2);line(rateX,900,rateX,1115,BLACK,2);
// D: forest plot on a log OR scale.
const fx0=1290,fx1=1690,fy0=893,fy1=1095,lmin=.5,lmax=30;
const coord=v=>fx0+(Math.log(v)-Math.log(lmin))/(Math.log(lmax)-Math.log(lmin))*(fx1-fx0);
for(const tick of [.5,1,2,5,10,30]){const x=coord(tick);line(x,fy0,x,fy1,tick===1?GRAY:LIGHT,tick===1?2:1);txt(String(tick),x-27,fy1+9,58,26,16,GRAY,false,'center');}
for(const [i,e] of ends.entries()){
 const r=getCond(e),y=i?1055:936,or=+r.event_or,lo=+r.event_or_ci_low,hi=+r.event_or_ci_high;
 line(coord(lo),y,coord(hi),y,RED,3);line(coord(lo),y-9,coord(lo),y+9,RED,2);line(coord(hi),y-9,coord(hi),y+9,RED,2);dot(coord(or),y,8,RED);
 const ty=i?1017:898;
 txt(i?'FEV1/FVC':'ICGC COPD',951,ty,315,23,18,BLACK,true);
 txt(`OR ${or.toFixed(2)} [${lo.toFixed(2)}, ${hi.toFixed(2)}]`,951,ty+23,315,23,17);
 txt(`P(event)=${pv(r.event_p_one_sided)}`,951,ty+46,315,22,16,GRAY);
 txt(`P(region)=${pv(r.region_p_one_sided)}`,951,ty+68,315,22,16,GRAY);
}
line(fx0,fy0,fx1,fy0,BLACK,2);line(fx1,fy0,fx1,fy1,BLACK,2);
line(fx0,fy1,fx1,fy1,BLACK,2);line(fx0,fy0,fx0,fy1,BLACK,2);
txt('Odds ratio (log scale)',1310,1130,320,28,18,GRAY,false,'center');
line(70,1171,1730,1171,LIGHT,2);
txt('A–B rank the sQTL common set; C–D condition on the global MAGMA top 2%. Outcome: posterior > 0.9 in GWAS-significant windows.',64,1184,1672,28,18,GRAY);
txt('The rank curves use the full-data fitted score. P in A–B tests curve area by region; P(event) and P(region) in D test conditional hit enrichment.',64,1217,1672,29,18,GRAY);
slide.speakerNotes.text=[
 'A–B source: outputs/pfam_rebuild_npnt/rank_enrichment_model_vs_sqtl_stats.csv and rank_enrichment_model_vs_sqtl_curves.csv.',
 'C–D source: outputs/pfam_rebuild_npnt/global_top2_figure_stats.csv, panel C_D_within_MAGMA_global_top2, posterior threshold 0.9.',
 'Recreates the requested four panels from the two reference PPTX files as native editable slide elements.',
 'A–B use full-data fitted score and regional sign flip P values; the ranking test is on cumulative curve area rather than AUPRC.',
 'C–D condition on globally selected MAGMA top 2%; event-level and region-level P values differ materially.',
 'E–F source: outputs/pfam_rebuild_npnt/loco_five_power_p_sweep.csv. Existing 22-fold LOCO predictions are relabeled at P < 5×10^-4 through P < 5×10^-8 without refitting.',
 'E shows average precision divided by positive prevalence. F shows top-10% positive sensitivity for component-only scores and their signed-sum combination.',
 'The ≥3/4 chromosome feature coverage rule was chosen after looking at previous results; these are exploratory rather than an independent model confirmation.',
 'Published figure conventions informing the display: Saito and Rehmsmeier, PLOS ONE 2015, doi:10.1371/journal.pone.0118432; PoPS, Nature Genetics 2023, doi:10.1038/s41588-023-01443-6; genomic feature-group ablation example, PMC12185218.',
 'The top four panels use external colocalization endpoints and different denominators from the MAGMA-labeled LOCO panels below.'
].join('\n');

// E–F: final two panels use the same out-of-fold scores and five requested
// MAGMA thresholds. They are component-only evaluations, not drop-one ablations.
line(70,1267,1730,1267,LIGHT,2);
txt('E   LOCO AUPRC relative to random ranking',55,1281,810,41,27,BLACK,true);
txt('F   Feature-only sensitivity among the top 10%',947,1281,820,41,27,BLACK,true);
txt('≥3/4 training-chromosome coverage; exploratory',60,1320,840,27,18,GRAY);
txt('Same LOCO fits, changing MAGMA positive threshold',951,1320,816,27,18,GRAY);

const eX0=144,eX1=805,fX0=1040,fX1=1690,pTop=1400,pBottom=1710;
const xAt=(i,x0,x1)=>x0+i*(x1-x0)/4;
const eMin=1,eMax=1.7,eY=v=>pBottom-(v-eMin)/(eMax-eMin)*(pBottom-pTop);
for(const v of [1,1.2,1.4,1.6]){
 const y=eY(v);line(eX0,y,eX1,y,v===1?BLACK:LIGHT,v===1?2:1);txt(`${v.toFixed(1)}×`,eX0-79,y-14,67,28,17,GRAY,false,'right');
}
line(eX0,pTop,eX0,pBottom,BLACK,2);
line(eX0,pTop,eX1,pTop,BLACK,2);line(eX1,pTop,eX1,pBottom,BLACK,2);
line(eX0,pBottom,eX1,pBottom,BLACK,2);
const eVals=cuts.map(c=>+getL('threequarter_training_chromosomes',c,'combined').auprc_over_prevalence);
for(let i=1;i<cuts.length;i++)line(xAt(i-1,eX0,eX1),eY(eVals[i-1]),xAt(i,eX0,eX1),eY(eVals[i]),RED,4);
for(let i=0;i<cuts.length;i++){
 const row=getL('threequarter_training_chromosomes',cuts[i],'combined'),x=xAt(i,eX0,eX1),y=eY(eVals[i]);
 dot(x,y,7,RED);txt(`+${((eVals[i]-1)*100).toFixed(0)}%`,x-54,y-39,108,26,19,RED,true,'center');
 txt(thresholdLabel(cuts[i]),x-67,pBottom+14,134,27,17,BLACK,false,'center');
 txt(`n=${row.positive_n}`,x-67,pBottom+43,134,24,16,GRAY,false,'center');
}
txt('AUPRC / positive prevalence',182,1370,575,26,19,GRAY,false,'center');

const fMin=.08,fMax=.19,fY=v=>pBottom-(v-fMin)/(fMax-fMin)*(pBottom-pTop);
for(const v of [.08,.10,.12,.14,.16,.18]){
 const y=fY(v);line(fX0,y,fX1,y,v===.10?BLACK:LIGHT,v===.10?2:1);txt(`${Math.round(v*100)}%`,fX0-78,y-14,65,28,17,GRAY,false,'right');
}
line(fX0,pTop,fX0,pBottom,BLACK,2);
line(fX0,pTop,fX1,pTop,BLACK,2);line(fX1,pTop,fX1,pBottom,BLACK,2);
line(fX0,pBottom,fX1,pBottom,BLACK,2);
const featureBlocks=[['expression',GRAY,'Expression'],['GO',BLUE,'GO'],['Pfam',GOLD,'Pfam'],['combined',RED,'Combined']];
for(let i=0;i<featureBlocks.length;i++){
 const [block,color,label]=featureBlocks[i];
 line(960+i*198,1367,995+i*198,1367,color,4);txt(label,1001+i*198,1353,159,29,17,BLACK);
 const vals=cuts.map(c=>+getL('threequarter_training_chromosomes',c,block).sensitivity_at_10_expected_random_ties);
 for(let j=1;j<cuts.length;j++)line(xAt(j-1,fX0,fX1),fY(vals[j-1]),xAt(j,fX0,fX1),fY(vals[j]),color,3.2);
 vals.forEach((v,j)=>dot(xAt(j,fX0,fX1),fY(v),5.5,color));
}
for(let i=0;i<cuts.length;i++)txt(thresholdLabel(cuts[i]),xAt(i,fX0,fX1)-67,pBottom+14,134,27,17,BLACK,false,'center');
txt('Reference line: random ranking = 10%',1040,1784,690,27,17,GRAY);
line(70,1810,1730,1810,LIGHT,2);
txt('E–F: 95,750 junctions; MAGMA P = upper-tail(ZSTAT). Each threshold changes labels only.',68,1818,1670,23,17,GRAY);
txt('Component score = |block contribution|; combined = |signed sum|. No genetic score enters the LOCO ranking.',68,1842,1670,23,17,GRAY);
txt('The ≥3/4 chromosome rule was selected post hoc. A–D use separate colocalization denominators.',68,1866,1670,23,17,GRAY);

await fs.mkdir(BUILD,{recursive:true});
const pngs=[];
for(let i=0;i<1;i++){
 const s=slideRefs[i],png=path.join(ROOT,'figures',`rebuilt_pfam_six_panel_four_sided.png`);
 const preview=await P.export({slide:s,format:'png',scale:1});await fs.writeFile(png,new Uint8Array(await preview.arrayBuffer()));pngs.push(png);
}
const candidate=path.join(BUILD,'candidate_rebuilt_pfam_six_panel.pptx');
await (await PresentationFile.exportPptx(P)).save(candidate);
const refHash=crypto.createHash('sha256').update(await fs.readFile(SOURCE)).digest('hex');
const result=await finalizePresentation({explicitTotalSlideCount:1,workspaceDir:ROOT,candidatePath:candidate,finalPath:FINAL,
 pythonExecutable:PY,integrityValidatorPath:path.join(SKILL,'container_tools','inspect_presentation_package_integrity.py'),
 layoutValidatorPath:path.join(SKILL,'container_tools','inspect_presentation_layout_geometry.py'),layoutArgs:['--validate-heading-fit'],
 requiredNativeTableOwnerSlides:[],fontPolicy:{basis:'reference',families:[FONT],referencePath:SOURCE,referenceSha256:refHash},
 verifyArtifactToolImport:true,receiptPath:path.join(BUILD,'rebuilt_pfam_six_panel_four_sided.validation.json')});
console.log(JSON.stringify({final:FINAL,pngs,result},null,2));
