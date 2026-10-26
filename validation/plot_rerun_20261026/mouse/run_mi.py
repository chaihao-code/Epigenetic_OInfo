from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess,os,time,json,sys
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
BASE=ROOT/'cases/mouse_preimp/mutual_information'
os.environ.update(MPLCONFIGDIR='/tmp/epinfo-mouse-mpl',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
records=[]
def run(script,args=(),key=None):
 key=key or Path(script).stem
 cmd=[sys.executable,str(BASE/'scripts'/script),*map(str,args)]
 start=time.time()
 with (OUT/'logs'/(key+'.log')).open('w') as log:
  proc=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=os.environ)
 files=[str(p.relative_to(ROOT)) for p in (BASE/'plots').rglob('*') if p.is_file() and p.stat().st_mtime>=start-1]
 rec=dict(id=key,script=str((BASE/'scripts'/script).relative_to(ROOT)),command=cmd,returncode=proc.returncode,seconds=round(time.time()-start,2),outputs=files,figure_files=[p for p in files if p.endswith(('.png','.pdf'))],log=str((OUT/'logs'/(key+'.log')).relative_to(ROOT)),input_type='intermediate MI CSVs; plot summaries only')
 records.append(rec);(OUT/'mi_results.json').write_text(json.dumps(records,indent=2)+'\n'); print(key,proc.returncode,len(rec['figure_files']),flush=True)
with ThreadPoolExecutor(max_workers=3) as pool:
 list(pool.map(lambda n:run('plot_by_block_all_anchors.py',('--n-bins',n),f'mi_anchors_n{n}'),[3,5,10]))
for n in [3,5,10]:run('plot_mi_ratio_hist.py',('--n-bins',n),f'mi_hist_n{n}')
for s in ['plot_joint_vs_pairwise_ratio_heatmap_v3.py','plot_composite_trend_hist_n3_n5_n10.py','plot_composite_trend_hist_n3_n5_n10_transposed.py']:run(s)
