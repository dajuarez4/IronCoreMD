#!/usr/bin/env python3
from pathlib import Path
import importlib.util,re,json,csv,shutil
import numpy as np
ROOT=Path(__file__).resolve().parent
R23=ROOT.parent/'round23_32atom_4x2x2_round22settings_nraise20_md100_jakar'
R19=ROOT.parent/'round19_48atom_6x2x2_qe75_ntyp128_lambda0p100_md1000_jakar'
CASE='conv1e-4_beta1e-2_48atom'
PREFIX='Fe_bcc48_r24_conv1e4_beta1e2'
def main():
    if list((ROOT/'cases').rglob('*.out')) or list((ROOT/'cases').rglob('*.save')):
        raise RuntimeError('Refusing to overwrite calculation results')
    spec=importlib.util.spec_from_file_location('r19',R19/'generate_round19_jakar.py')
    h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    src=R23/'cases/conv1e-4_beta1e-2_32atom'
    lines=(src/'md_100steps.in').read_text().splitlines()
    pr=h.input_block(lines,'ATOMIC_POSITIONS crystal',32)
    vr=h.input_block(lines,'ATOMIC_VELOCITIES { a.u }',32)
    pos=np.array([[float(x) for x in row.split()[1:4]] for row in pr])
    vel=np.array([[float(x) for x in row.split()[1:4]] for row in vr])
    first=pos.copy();first[:,0]*=2/3
    third=pos[:16].copy();third[:,0]=(third[:,0]*2+2)/3
    positions=np.concatenate([first,third])
    velocities,vs=h.rescale_velocities(np.concatenate([vel,np.roll(vel[:16],9,axis=0)[:,[1,2,0]]]))
    vs['source']='Round23 initial 32-atom velocities plus row/axis-permuted first 16; COM removed and full 48-atom set rescaled to 4000 K'
    extra=np.array(h.LABELS[32:],dtype=object);np.random.default_rng(240048).shuffle(extra)
    labels=[r.split()[0] for r in pr]+extra.tolist()
    for name in ['pseudo','logs',f'cases/{CASE}']:(ROOT/name).mkdir(parents=True,exist_ok=True)
    shutil.copy2(R23/'pseudo'/h.PSEUDO,ROOT/'pseudo'/h.PSEUDO)
    for f in sorted(src.glob('*.in')):
        text=f.read_text().split('ATOMIC_SPECIES')[0]
        text=text.replace('nat=32','nat=48').replace('ntyp=32','ntyp=48').replace('Fe_bcc32_r23_conv1e4_beta1e2',PREFIX).replace('nstep=100','nstep=400')
        additions=''
        for i,(theta,phi) in enumerate(h.ANGLES[32:],33):
            additions+=f'   starting_magnetization({i})=0.12500000\n   angle1({i})={theta:.10f}\n   angle2({i})={phi:.10f}\n'
        match=re.search(r'(?ms)^&SYSTEM\n.*?^/',text)
        text=text[:match.end()-1]+additions+text[match.end()-1:]
        md=f.name.startswith('md_')
        text+='\n'.join(h.cards(positions,labels,velocities if md else None))
        (ROOT/'cases'/CASE/f.name.replace('100steps','400steps')).write_text(text)
    (ROOT/'cases.txt').write_text(CASE+'\n')
    vectors=np.array([h.direction(*a) for a in h.ANGLES])
    ms={'natoms':48,'unique_directions':48,'target_moment_per_atom_muB':2.0,'target_vector_sum_muB':(2*vectors.sum(axis=0)).tolist(),'position_labels':labels,'extra_label_seed':240048}
    for name,data in [('magnetic_seed_summary.json',ms),('velocity_summary.json',vs),('structure_summary.json',{'source':'Round23 original 32-atom configuration extended by its first 16-atom tile along x','cell_A':[14.34,4.78,4.78],'a_A':2.39,'supercell':'6x2x2','natoms':48})]:
        (ROOT/name).write_text(json.dumps(data,indent=2)+'\n')
    meta=json.loads((src/'case_metadata.json').read_text());meta.update(case=CASE,source_round=23,supercell='6x2x2',natoms=48,ntyp=48,md_steps=400)
    (ROOT/'cases'/CASE/'case_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    with (ROOT/'case_manifest.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=meta);w.writeheader();w.writerow(meta)
    print('Generated 48 atoms, 400 MD steps; target spin sum:',ms['target_vector_sum_muB'])
if __name__=='__main__':main()
