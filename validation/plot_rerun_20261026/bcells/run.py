from pathlib import Path
import concurrent.futures, subprocess, os, time, json, hashlib
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
CASE=ROOT/'cases/bcells_and_CLL'
ENV=dict(os.environ,MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/bcells-rerun-mpl',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',EPINFO_ROOT=str(ROOT))
records=[]
def run(relative):
 script=CASE/relative; log=OUT/'logs'/(script.stem+'.log'); start=time.time()
 before={str(p):p.stat().st_mtime_ns for p in CASE.rglob('*') if p.suffix in ('.pdf','.png','.svg')}
 with log.open('w') as handle:
  result=subprocess.run(['python',str(script)],cwd=ROOT,env=ENV,stdout=handle,stderr=subprocess.STDOUT)
 changed=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in CASE.rglob('*') if p.suffix in ('.pdf','.png','.svg') and p.stat().st_mtime_ns!=before.get(str(p))]
 record=dict(script=str(script.relative_to(ROOT)),command=f'EPINFO_ROOT="$PWD" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python {script.relative_to(ROOT)}',returncode=result.returncode,seconds=round(time.time()-start,2),log=str(log.relative_to(ROOT)),outputs=changed)
 print(script.stem, result.returncode, round(time.time()-start,2),flush=True)
 return record
independent=['dmr_repro_cemt/scripts/plot_dmr_volcano.py','dmr_repro_cemt/scripts/plot_igh_tracks_corrected.py','scripts/plot_de_male_ighv_volcano_clear.py','scripts/plot_gsea_eight_panel_v2_18x10cm.py','scripts/plot_gsea_representative_bars.py','scripts/plot_gsea_top20_bars.py','scripts/validate_male_de_in_females.py','scripts/plot_abc_triplet_classes.py']
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for record in pool.map(run,independent):records.append(record)
for path in ['scripts/plot_abc_pooled_combined.py','scripts/plot_abc_pooled_compact_preview.py','global_oinfo_and_triplets/plot_global_oinfo_and_ABC_triplets.py']:
 records.append(run(path))
# A concurrent output snapshot may include another task's file; report generation assigns exact output families.
(OUT/'runs.json').write_text(json.dumps(records,indent=2)+'\n')
