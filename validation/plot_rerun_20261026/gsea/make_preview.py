from pathlib import Path
import json
from PIL import Image,ImageDraw
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
r=json.loads((OUT/'report.json').read_text());figs=[];examples=[]
for e in r['entries']:
 these=[dict(term=f['term'],kind=g['kind'],region=Path(g['input_folder']).parent.name,path=next(o['path'] for o in f['outputs'] if o['path'].endswith('.png'))) for g in e['groups'] for f in g['figures']]
 figs+=these
 if these:
  examples += [these[round(i*(len(these)-1)/11)] for i in range(12)]
def sheet(records,target,w,h,cols):
 rows=(len(records)+cols-1)//cols;canvas=Image.new('RGB',(w*cols,h*rows),'#eeeeee');draw=ImageDraw.Draw(canvas)
 for i,f in enumerate(records):
  im=Image.open(ROOT/f['path']).convert('RGB');im.thumbnail((w-8,h-38));x=(i%cols)*w+(w-im.width)//2;y=(i//cols)*h+4;canvas.paste(im,(x,y));label=f'{i+1}. {f["term"]}'
  draw.text(((i%cols)*w+4,(i//cols)*h+h-30),label[:int(w/6)],fill='black')
  draw.text(((i%cols)*w+4,(i//cols)*h+h-17),f['kind']+' '+f['region'][:int(w/6)-16],fill='black')
 canvas.save(OUT/target,quality=90)
sheet(figs,'contact_sheet_all.jpg',160,150,20);sheet(examples,'contact_sheet_examples.jpg',500,410,4)
(OUT/'preview_index.json').write_text(json.dumps(dict(all=figs,examples=examples),ensure_ascii=False,indent=2)+'\n')
p=OUT/'REPORT.md';s=p.read_text();s+='\n[全部572幅缩略图](contact_sheet_all.jpg) · [24幅代表图预览](contact_sheet_examples.jpg) · [预览索引](preview_index.json)\n';p.write_text(s)
print('Previewed',len(figs),'figures')
