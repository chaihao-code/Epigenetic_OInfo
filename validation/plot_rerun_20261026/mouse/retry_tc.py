from pathlib import Path
import subprocess,sys,os,time,json
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
os.environ.update(MPLCONFIGDIR='/tmp/epinfo-mouse-mpl',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
records=json.loads((OUT/'development_results.json').read_text())
for key in ['plot_total_oinfo_trajectories_18cm']:
 r=next(r for r in records if r['id']==key);start=time.time()
 with (OUT/'logs'/(key+'.log')).open('w') as log:p=subprocess.run(r['command'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=os.environ)
 r['initial_returncode']=r['returncode'];r['returncode']=p.returncode;r['seconds']=round(time.time()-start,2)
 (OUT/'development_results.json').write_text(json.dumps(records,indent=2)+'\n');print(key,p.returncode,flush=True)
