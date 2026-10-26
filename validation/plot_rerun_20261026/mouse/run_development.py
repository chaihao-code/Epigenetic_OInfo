from pathlib import Path
import subprocess,os,time,json,sys
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parent; BASE=ROOT/'cases/mouse_development'
os.environ.update(MPLCONFIGDIR='/tmp/epinfo-mouse-mpl',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
records=[]
tasks=[('script/plot_triplet_counts_corrected_by_stage.py',[]),('script/plot_triplet_A_B_C_counts_trajectories.py',[]),('script/plot_special_triplet_f1.py',[]),('script/plot_triplet_ABC_selected_regions.py',[]),('script/plot_synergistic_triplet_mark_trajectories.py',['--plot-only']),('script/plot_triplet_ABC_and_F1_panels.py',[]),('script/plot_triplet_A_C_D_counts_trajectories.py',[]),('plot_tc_dtc/plot_combined_hoi_trajectories_2x4.py',['--output-dir',str(OUT/'figures/tc_dtc_original')]),('plot_tc_dtc/plot_tc_dtc_stacked_18cm.py',[]),('plot_tc_dtc/plot_total_oinfo_trajectories_18cm.py',[]),('redundancy_core/resize_redundancy_core.py',[])]
for script,args in tasks:
 key=Path(script).stem;cmd=[sys.executable,str(BASE/script),*args];start=time.time()
 with (OUT/'logs'/(key+'.log')).open('w') as log:proc=subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=os.environ)
 files=[str(p.relative_to(ROOT)) for p in BASE.rglob('*') if p.is_file() and p.stat().st_mtime>=start-1 and p.suffix in ['.pdf','.png','.tsv','.json','.csv']]
 rec=dict(id=key,script=str((BASE/script).relative_to(ROOT)),command=cmd,returncode=proc.returncode,seconds=round(time.time()-start,2),outputs=files,figure_files=[p for p in files if p.endswith(('.png','.pdf'))],log=str((OUT/'logs'/(key+'.log')).relative_to(ROOT)),input_type='previous result CSV/TSV summaries',note='PDF vector reconstruction from existing figure, not intermediate tables' if 'resize_' in script else '')
 records.append(rec);(OUT/'development_results.json').write_text(json.dumps(records,indent=2)+'\n');print(key,proc.returncode,len(rec['figure_files']),flush=True)
