#!/usr/bin/env python3
"""Refresh comparative MD dashboards using completed SCF/MD records only."""
import argparse,csv,hashlib,json,os,re
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-round24-25-dashboard')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from io import BytesIO
ROOT=Path(__file__).resolve().parent
F=r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[EeDd][-+]?\d+)?'
def scalar(pattern,text):
    x=re.findall(pattern,text);return float(x[-1].replace('D','E')) if x else float('nan')
def parse(path,n,template):
    t=path.read_text(errors='replace');cell=np.array([[float(v) for v in line.split()] for line in re.search(r'(?m)^CELL_PARAMETERS[^\n]*\n((?:[^\n]+\n){3})',template)[1].splitlines()])
    nat=int(re.search(r'number of atoms/cell\s*=\s*(\d+)',t)[1])
    alat=scalar(r'lattice parameter \(alat\)\s*=\s*('+F+')',t)*0.529177210903
    tau=np.array(re.findall(r'tau\(\s*\d+\s*\)\s*=\s*\(\s*('+F+')\s+('+F+')\s+('+F+')',t)[:nat],float)
    initial=tau*alat@np.linalg.inv(cell);current=initial.copy()
    bounds=list(re.finditer(r'Entering Dynamics:\s*iteration\s*=\s*(\d+)',t))
    rows=[];positions=[];spins=[];prev=0
    for i,m in enumerate(bounds):
        before=t[prev:m.start()];after=t[m.end():bounds[i+1].start() if i+1<len(bounds) else len(t)];prev=m.end()
        if 'convergence has been achieved' not in before:raise ValueError('Dynamics without convergence')
        vm=re.findall(r'total magnetization\s*=\s*('+F+')\s+('+F+')\s+('+F+')',before)
        spin=np.array(re.findall(r'(?m)^\s*magnetization\s*:\s*('+F+')\s+('+F+')\s+('+F+')',before)[-nat:],float)
        assert spin.shape==(nat,3)
        v=np.array(vm[-1],float);absolute=scalar(r'absolute magnetization\s*=\s*('+F+')',before)
        force=scalar(r'Total force\s*=\s*('+F+')',before);corr=scalar(r'Total SCF correction\s*=\s*('+F+')',before)
        post=after.split('Self-consistent Calculation')[0]
        row={'step':int(m[1]),'temperature_K':scalar(r'\n\s*temperature\s*=\s*('+F+')',post),'pressure_GPa':scalar(r'P=\s*('+F+')',before)/10,'scf_iterations':scalar(r'convergence has been achieved in\s*(\d+)',before),'net_moment_muB_atom':np.linalg.norm(v)/nat,'absolute_moment_muB_atom':absolute/nat,'local_moment_muB':np.linalg.norm(spin,axis=1).mean(),'force_correction_percent':100*corr/force,'energy_Ry':scalar(r'!\s*total energy\s*=\s*('+F+')',before),'Fermi_eV':scalar(r'the Fermi energy is\s*('+F+')',before)}
        positions.append(current.copy());spins.append(spin);rows.append(row)
        prop=re.search(r'(?m)^ATOMIC_POSITIONS \(crystal\)\n((?:Fe\d+[^\n]*\n){'+str(nat)+'})',post)
        if prop:current=np.array([[float(v) for v in line.split()[1:4]] for line in prop[1].splitlines()])
        else:raise ValueError('Incomplete propagated coordinates')
    return {'round':n,'natoms':nat,'cell':cell,'positions':np.array(positions),'spins':np.array(spins),'rows':rows,'source':str(path),'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'job_done':'JOB DONE' in t}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in [24,25]:p.add_argument(f'--round{n}-output',type=Path)
    p.add_argument('--no-gif',action='store_true');a=p.parse_args()
    dest=ROOT/'dashboard';dest.mkdir(exist_ok=True);data=[]
    for n in [24,25]:
        path=getattr(a,f'round{n}_output')
        if path is None:
            folders=[f for f in ROOT.parent.glob(f'round{n}_*') if f.is_dir()]
            if len(folders)!=1:raise ValueError(f'Pass --round{n}-output')
            path=folders[0]/'md_400steps.run001.out'
        d=parse(path,n,(ROOT/f'round{n}/md_template.in').read_text());data.append(d)
        with (dest/f'round{n}_steps.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=d['rows'][0]);w.writeheader();w.writerows(d['rows'])
        np.savez_compressed(dest/f'round{n}_frames.npz',cell_A=d['cell'],force_aligned_positions_crystal=d['positions'],local_moments_muB=d['spins'],step=[r['step'] for r in d['rows']])
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(15,11),facecolor='#f4f7fb')
    grid=fig.add_gridspec(4,2,left=.08,right=.96,bottom=.085,top=.85,hspace=.7,wspace=.28)
    structures=[fig.add_subplot(grid[0,j],projection='3d') for j in range(2)]
    specs=[('temperature_K','Temperature after propagation','K'),('pressure_GPa','Electronic stress pressure','GPa'),('scf_iterations','SCF iterations','Iterations'),('net_moment_muB_atom','Net magnetization per atom','μB / Fe'),('force_correction_percent','SCF correction / total force','%'),('local_moment_muB','Mean local moment','μB / Fe')]
    axes=[fig.add_subplot(grid[1+i//2,i%2]) for i in range(6)];cursors=[]
    colors=['#2563eb','#dc751e'];maxstep=max(len(d['rows']) for d in data)
    for ax,(key,title,ylabel) in zip(axes,specs):
        for d,color in zip(data,colors):ax.plot([r['step'] for r in d['rows']],[r[key] for r in d['rows']],'.-',color=color,ms=3,label=f'Round{d["round"]} ({d["natoms"]} Fe)')
        ax.set(title=title,xlabel='MD step',ylabel=ylabel,xlim=(.5,maxstep+.5));ax.grid(alpha=.2)
        cursors.append(ax.axvline(1,color='#64748b',alpha=.6,lw=1))
    axes[0].axhline(4000,color='#64748b',ls='--',lw=1);axes[0].legend(fontsize=8)
    fig.text(.06,.955,'ROUND24 + ROUND25  |  MD progress dashboard',fontsize=22,weight='bold',color='#12304d')
    fig.text(.06,.922,'4000 K target · nraise 20 · λ 0.1 Ry · constrained noncollinear Fe · no SOC',fontsize=11)
    status=fig.text(.06,.887,'',fontsize=11,color='#12304d')
    fig.text(.06,.028,'Completed records only; local snapshots do not prove live job status. Positions/spins correspond to force evaluation; T is post-propagation.\nDifferent initial states and trajectory lengths. Electronic bands require the separate checkpoint calculation; no band results are fabricated.',fontsize=9,color='#526278')
    images=[]
    for frame in range(maxstep):
        for ax,d in zip(structures,data):
            idx=min(frame,len(d['rows'])-1);ax.clear();pos=(d['positions'][idx]%1)@d['cell'];spin=d['spins'][idx];u=spin/np.linalg.norm(spin,axis=1)[:,None]
            ax.scatter(*pos.T,c=u[:,2],cmap='coolwarm',vmin=-1,vmax=1,s=22)
            ax.quiver(*pos.T,*u.T,length=.5,color='#334155',linewidth=.5)
            ax.set(xlim=(0,d['cell'][0,0]),ylim=(0,d['cell'][1,1]),zlim=(0,d['cell'][2,2]))
            ax.set_box_aspect(np.diag(d['cell']).copy());ax.view_init(22,-55);ax.tick_params(labelsize=6)
            ax.set_title(f'Round{d["round"]}: {d["natoms"]} Fe · step {d["rows"][idx]["step"]}'+(' (last available)' if frame>=len(d['rows']) else ''),fontsize=10)
        status.set_text('Available: '+ '  |  '.join(f'Round{d["round"]}: {len(d["rows"])} / 400 steps; '+('JOB DONE' if d['job_done'] else 'partial output') for d in data))
        for cursor in cursors:cursor.set_xdata([frame+1]*2)
        if not a.no_gif:
            buffer=BytesIO();fig.savefig(buffer,format='png',dpi=100);buffer.seek(0);images.append(Image.open(buffer).convert('RGB').quantize(colors=256))
        if frame==maxstep-1:fig.savefig(dest/'round24_round25_dashboard.png',dpi=170)
    if images:images[0].save(dest/'round24_round25_dashboard.gif',save_all=True,append_images=images[1:],duration=250,loop=0,disposal=2)
    plt.close(fig)
    summaries=[{k:d[k] for k in ['round','natoms','source','source_sha256','job_done']}|{'completed_steps':len(d['rows'])} for d in data]
    (dest/'provenance.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (dest/'index.html').write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>Round24 / Round25 dashboard</title><style>body{font:16px system-ui;background:#f4f7fb;color:#12304d;max-width:1500px;margin:30px auto;padding:15px}img{width:100%;border-radius:12px}button,a{padding:10px;cursor:pointer}section{background:white;padding:20px;border-radius:12px;margin:20px 0}</style><h1>Round24 + Round25</h1><p>MD snapshots and electronic-band workflow</p><section><button onclick="document.getElementById('md').src='round24_round25_dashboard.gif'">Play MD animation</button><button onclick="document.getElementById('md').src='round24_round25_dashboard.png'">Final summary</button><img id="md" src="round24_round25_dashboard.png"><p><a href="round24_steps.csv">Round24 measurements</a><a href="round25_steps.csv">Round25 measurements</a><a href="provenance.json">Data provenance</a></p></section><section><h2>Electronic bands</h2><p id="pending">Awaiting Jakar calculations. Run plot_bands.py after both band jobs to populate this panel. Curves will be folded supercell bands, not primitive BCC bands.</p><img src="electronic_bands.png" alt="Electronic band results not computed yet" onload="document.getElementById('pending').textContent='Frozen-snapshot spinor bands, referenced to each snapshot SCF Fermi level.'" onerror="this.style.display='none'"></section></html>''')
    print('Dashboard saved:',dest)
if __name__=='__main__':main()
