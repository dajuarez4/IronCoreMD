#!/usr/bin/env python3
"""Prepare fixed-snapshot SCF + supercell bands from an inactive QE MD checkpoint."""
import argparse,fcntl,hashlib,json,re,shutil,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
BOHR=0.529177210903

def xml_root(path):
    r=ET.parse(path).getroot()
    for el in r.iter():el.tag=el.tag.split('}')[-1]
    return r

def system_from(text,nbnd):
    block=re.search(r'(?ms)^&SYSTEM\n.*?^/',text)[0]
    return block[:-1]+f'   nbnd={nbnd}\n/'

def render_inputs(template,prefix,cell,positions,names,nbnd,intervals):
    species=re.search(r'(?ms)^ATOMIC_SPECIES\n(.*?)^ATOMIC_POSITIONS',template)[1]
    cards='ATOMIC_SPECIES\n'+species+'ATOMIC_POSITIONS angstrom\n'
    cards+=''.join(n+' '+' '.join(f'{v:.12f}' for v in xyz)+'\n' for n,xyz in zip(names,positions))
    cards+='CELL_PARAMETERS angstrom\n'+''.join(' '.join(f'{v:.12f}' for v in row)+'\n' for row in cell)
    control=f"&CONTROL\n calculation='scf'\n restart_mode='from_scratch'\n prefix='{prefix}'\n outdir='./qe_tmp'\n pseudo_dir='../pseudo'\n disk_io='low'\n verbosity='high'\n max_seconds=160000\n/\n"
    syst=system_from(template,nbnd)+'\n'
    elec=re.search(r'(?ms)^&ELECTRONS\n.*?^/',template)[0]+'\n'
    scf=control+syst+elec+cards+'K_POINTS automatic\n1 1 1 0 0 0\n'
    nodes=np.array([[0,0,0],[.5,0,0],[.5,.5,0],[0,.5,0],[0,0,0],[0,0,.5]])
    ks=np.concatenate([np.linspace(a,b,intervals,endpoint=False) for a,b in zip(nodes[:-1],nodes[1:])]+[nodes[-1:]])
    bands=control.replace("calculation='scf'","calculation='bands'")+syst
    bands+="&ELECTRONS\n startingpot='file'\n startingwfc='atomic+random'\n diagonalization='cg'\n diago_thr_init=1.0d-6\n diago_cg_maxiter=200\n diago_full_acc=.true.\n/\n"
    bands+=cards+f'K_POINTS crystal\n{len(ks)}\n'+''.join(' '.join(f'{v:.10f}' for v in k)+' 1.0\n' for k in ks)
    post=f"&BANDS\n prefix='{prefix}'\n outdir='./qe_tmp'\n filband='bands.dat'\n lsym=.false.\n no_overlap=.true.\n/\n"
    reciprocal=2*np.pi*np.linalg.inv(cell).T
    dist=np.r_[0,np.cumsum(np.linalg.norm(np.diff(ks@reciprocal,axis=0),axis=1))]
    path={'labels':['Γ','X','S','Y','Γ','Z'],'tick_indices':[i*intervals for i in range(6)],'kpoints_crystal':ks.tolist(),'distance_inv_A':dist.tolist(),'description':'Geometric supercell reciprocal path; no primitive-cell unfolding or symmetry claim'}
    return scf,bands,post,path

def signature(root):
    return {str(p.relative_to(root)):(p.stat().st_size,p.stat().st_mtime_ns) for p in root.rglob('*') if p.is_file()}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--round',type=int,choices=[24,25],required=True)
    p.add_argument('--source-case',type=Path,required=True,help='Inactive MD case directory containing qe_tmp')
    p.add_argument('--points-per-segment',type=int,default=2)
    args=p.parse_args()
    if args.points_per_segment<1:raise ValueError('points-per-segment must be positive')
    source=args.source_case.resolve();case=ROOT/f'round{args.round}'
    config=json.loads((case/'config.json').read_text());prefix=config['prefix']
    source_tmp=source/'qe_tmp';save=source_tmp/(prefix+'.save');schema=save/'data-file-schema.xml'
    if not schema.is_file():raise FileNotFoundError(schema)
    if not (save/'charge-density.dat').is_file():raise FileNotFoundError('Missing saved charge-density.dat')
    if not (list(save.glob('wfc*.dat'))+list(source_tmp.glob(prefix+'.wfc*'))):raise FileNotFoundError('Missing wavefunction files; supply full qe_tmp')
    lockpath=source/'.workflow.lock'
    lock=lockpath.open('r') if lockpath.exists() else None
    if lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('MD is active. Wait for a clean stop before copying its checkpoint.')
    else:raise RuntimeError('Expected workflow lock file not found: provide the original round24/25 case directory after its workflow has stopped')
    try:
        if (case/'qe_tmp').exists() or (case/'checkpoint.json').exists():raise RuntimeError('Checkpoint already imported; refusing overwrite')
        before=signature(source_tmp)
        root=xml_root(schema);out=root.find('output')
        if out is None:raise ValueError('No output state in checkpoint XML')
        state=out.find('atomic_structure');nat=int(state.attrib['nat'])
        if nat!=config['natoms']:raise ValueError('Wrong atom count')
        cell=np.array([[float(v) for v in state.find('cell/'+a).text.split()] for a in ['a1','a2','a3']])*BOHR
        atoms=state.find('atomic_positions').findall('atom');names=[a.attrib['name'] for a in atoms]
        positions=np.array([[float(v) for v in a.text.split()] for a in atoms])*BOHR
        if not np.allclose(cell,np.diag(config['cell_A']),atol=1e-5):raise ValueError('Unexpected checkpoint cell')
        if names!=config['atom_labels']:raise ValueError('Checkpoint species order differs from source template')
        bs=out.find('band_structure');nbnd=int(bs.findtext('nbnd'))
        if abs(float(bs.findtext('nelec'))-16*nat)>1e-6:raise ValueError('Unexpected electron count')
        if bs.findtext('noncolin','false').strip().lower()!='true':raise ValueError('Expected spinor checkpoint')
        for key,expected in [('ecutwfc',35.5),('ecutrho',248.0)]:
            if not np.isclose(float(out.findtext('basis_set/'+key)),expected):raise ValueError('Checkpoint cutoff differs from template')
        pseudoname='Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'
        stored=list(save.glob('*.UPF'))
        if not stored:raise ValueError('Missing UPF in save directory; cannot verify checkpoint pseudopotential')
        expected=hashlib.sha256((ROOT/'pseudo'/pseudoname).read_bytes()).hexdigest()
        if not any(hashlib.sha256(f.read_bytes()).hexdigest()==expected for f in stored):raise ValueError('Pseudopotential mismatch')
        template=(case/'md_template.in').read_text()
        scf,bands,post,path=render_inputs(template,prefix,cell,positions,names,nbnd,args.points_per_segment)
        shutil.copytree(source_tmp,case/'qe_tmp')
        if before!=signature(source_tmp):raise RuntimeError('Source changed during copy; do not run this imported checkpoint')
        for exitfile in (case/'qe_tmp').glob('*.EXIT'):exitfile.unlink()
        for filename,text in [('01_snapshot_scf.in',scf),('02_bands.in',bands),('03_bands_post.in',post)]: (case/filename).write_text(text)
        (case/'kpath.json').write_text(json.dumps(path,indent=2)+'\n')
        provenance={'source_case':str(source),'source_xml_sha256':hashlib.sha256(schema.read_bytes()).hexdigest(),'natoms':nat,'nbnd':nbnd,'nk':len(path['kpoints_crystal']),'source_scf_converged':out.findtext('convergence_info/scf_conv/convergence_achieved'),'snapshot_source':'output/atomic_structure from copied XML','method':'reconverge fixed snapshot from copied density and Gamma wavefunctions; new wavefunctions along band path','source_file_count':len(before)}
        (case/'checkpoint.json').write_text(json.dumps(provenance,indent=2)+'\n')
        print(f'Prepared {case}: {nbnd} spinor bands, {len(path["kpoints_crystal"])} k points')
    finally:
        if lock:lock.close()
if __name__=='__main__':main()
