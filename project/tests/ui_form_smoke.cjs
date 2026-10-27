// Offline form-serialization unit test. Fake DOM and fetch; no browser/network.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict'),path=require('node:path');
const html=fs.readFileSync(path.join(__dirname,'../hexw/ui/index.html'),'utf8');
const elements=new Map();
for(const match of html.matchAll(/<[^>]+\bid="([^"]+)"[^>]*>/g)){
 const value=match[0].match(/\bvalue="([^"]*)"/);
 elements.set(match[1],{value:value?value[1]:'',checked:false,dataset:{},hidden:false,classList:{toggle(){}}});
}
const blocks=['xrd','mossbauer','sem','eds','fmr'].map(method=>({dataset:{method},hidden:false}));
const fields=new Map();
elements.get('forecastForm').elements={namedItem(name){if(!fields.has(name))fields.set(name,{value:''});return fields.get(name)}};
const calls=[];
const responses=new Map();
const context={document:{getElementById:id=>{assert(elements.has(id),'missing element '+id);return elements.get(id)},
 querySelectorAll:selector=>selector==='[data-method]'?blocks:[]},setInterval(){},
 FormData:class{constructor(form){this.form=form}get(name){return this.form.elements.namedItem(name).value}},
 fetch:async(url,options)=>{if(options?.method==='POST'){calls.push({url,data:JSON.parse(options.body)});return{ok:true,json:async()=>({id:'synthetic'})};}
 if(responses.has(url))return responses.get(url);
 return{ok:true,json:async()=>({jobs:[],runs:[],root:'SYNTHETIC'})};}};
vm.runInNewContext(html.match(/<script>([\s\S]*?)<\/script>/)[1],context);
const el=id=>elements.get(id);
async function submit(method,file){el('method').value=method;el('method').onchange();el('analysisFile').value=file;el('analysisOptions').value='{}';
 const old=calls.length;await el('analysisForm').onsubmit({preventDefault(){}});assert.equal(calls.length,old+1,el('message').textContent);
 assert.equal(calls.at(-1).url,'/api/analyze');return calls.at(-1).data;}
(async()=>{
 el('xWave').value='1.54184';let d=await submit('xrd','"C:\\synthetic\\scan.mdi"');assert.equal(d.file,'C:\\synthetic\\scan.mdi');assert.equal(d.options.wavelength_angstrom,1.54184);
 d=await submit('mossbauer','C:\\synthetic\\spectrum.spc');assert.deepEqual(d.options.n_components,[1,2,3]);assert(!d.options.velocity_calibration);
 el('mCal').checked=true;el('mRef').value='independent alpha-Fe';el('mSlope').value='.05';el('mZero').value='-12';d=await submit('mossbauer','spectrum.spc');assert.equal(d.options.velocity_calibration.mm_s_per_channel,.05);
 el('sScale').value='.1';el('sForeground').value='bright';el('sRoi').value='0, 0, 100, 100';d=await submit('sem','image.png');assert.deepEqual(d.options.roi_px,[0,0,100,100]);assert.equal(d.options.um_per_pixel,.1);
 el('eReference').value='Sr';d=await submit('eds','table.csv');assert.equal(d.options.normalize_element,'Sr');
 el('fUnit').value='GHz';el('fSignal').value='derivative';el('fBaseline').value='0';d=await submit('fmr','fmr.csv');assert.equal(d.options.signal_kind,'derivative');assert.equal(d.options.axis_unit,'GHz');
 const before=calls.length;el('coznExample').onclick();assert.equal(calls.length,before,'preset must not launch a calculation');
 await el('forecastForm').onsubmit({preventDefault(){},target:el('forecastForm')});
 d=calls.at(-1).data;assert.equal(calls.at(-1).url,'/api/predict');assert.equal(d.formula,'SrCo1Zn1Fe16O27');assert.deepEqual(d.fields_T,[0,.4,.8]);assert.equal(d.frequencies_GHz.length,401);assert.equal(d.n_draws,48);assert.equal(d.g,2);assert.equal(d.alpha,.02);
 const count=calls.length;fields.get('f1').value='0.01';await el('forecastForm').onsubmit({preventDefault(){},target:el('forecastForm')});assert.equal(calls.length,count,'descending frequency grid must not submit');
 el('coznExample').onclick();el('emOptions').value='{"formula":"different"}';await el('forecastForm').onsubmit({preventDefault(){},target:el('forecastForm')});assert.equal(calls.length,count,'EM options must not silently override the main form');
 el('emOptions').value='[]';await el('forecastForm').onsubmit({preventDefault(){},target:el('forecastForm')});assert.equal(calls.length,count);
 responses.set('/results/predictions/missing/',{ok:false,json:async()=>({error:'Synthetic missing result'})});
 await context.showResult('predictions/missing');assert.equal(el('message').textContent,'Synthetic missing result');
 responses.set('/results/predictions/demo/',{ok:true,json:async()=>({files:['summary.json']})});
 responses.set('/results/predictions/demo/summary.json',{ok:true,json:async()=>({request:{label:'TEST <script>'},status:'complete_unvalidated',field_sweep:[{T_K:300,mu0H_T:0,f_low_median_GHz:7.375,f_low_lo95_GHz:.73,f_low_hi95_GHz:40.9,n_valid:48,n_total:48}]})});
 await context.showResult('predictions/demo');assert(el('resultDetail').innerHTML.includes('TEST &lt;script&gt;'));assert(el('resultDetail').innerHTML.includes('<table>'));assert(!el('resultDetail').innerHTML.includes('TEST <script>'));
 const versions=[{run:'analysis/recent',label:'scan.mdi',method:'xrd',status:'complete',latest_for_source:true,created_utc:'2026-10-27T00:00:00Z'},{run:'analysis/earlier',label:'scan.mdi',method:'xrd',status:'complete',latest_for_source:false,created_utc:'2026-10-26T00:00:00Z'}];
 el('showArchives').checked=false;context.renderAnalyses(versions);assert(el('analyses').innerHTML.includes('analysis/recent'));assert(!el('analyses').innerHTML.includes('analysis/earlier'));
 el('showArchives').checked=true;context.renderAnalyses(versions);assert(el('analyses').innerHTML.includes('analysis/earlier'));
 console.log('6 analysis-form cases plus Co-Zn preset/submission, invalid ranges/JSON, readable summary and result-error handling passed; no network or browser used');
})().catch(error=>{console.error(error);process.exitCode=1});
