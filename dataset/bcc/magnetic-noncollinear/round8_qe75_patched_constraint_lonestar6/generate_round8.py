#!/usr/bin/env python3
"""Generate patched-QE Round-8 local-moment validation SCFs."""
import csv, json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
PSEUDO='Fe.pbe-spn-kjpaw_psl.1.0.0.UPF'
LABELS=[f'Fe{i:02d}' for i in range(1,9)]
CASES=(
 ('lambda0p010',(0.0025,0.005,0.010)),
 ('lambda0p020',(0.0050,0.010,0.020)),
 ('lambda0p050',(0.0050,0.020,0.050)),
 ('lambda0p100',(0.0050,0.020,0.100)),
)
STAGES=(('stage1_4e-3',4e-3,500),('stage2_1e-3',1e-3,500),('stage3_3e-4',3e-4,600))
ANGLES=((72.8157590495,72.7795960677),(73.5303944933,92.3546343919),(33.6579186653,331.8844103957),(95.5478199306,0.9655109597),(107.1842409505,252.7795960677),(106.4696055067,272.3546343919),(146.3420813347,151.8844103957),(84.4521800694,180.9655109597))
POSITIONS='''ATOMIC_POSITIONS crystal
Fe01 0.012027047575 0.073702611029 0.968028113246
Fe02 0.284393727779 0.288436353207 0.473334848881
Fe03 0.984382264316 0.580798864365 0.031577926129
Fe04 0.243689253926 0.798996508121 0.501007080078
Fe05 0.466723054647 0.925974182785 0.995024635457
Fe06 0.772856652737 0.200137659907 0.519648075104
Fe07 0.513922214508 0.431229859591 0.018644046038
Fe08 0.722005784512 0.700724005699 0.492735266685'''.splitlines()
CELL='''CELL_PARAMETERS angstrom
4.780000000000 0.000000000000 0.000000000000
0.000000000000 4.780000000000 0.000000000000
0.000000000000 0.000000000000 2.390000000000'''.splitlines()

def make_input(case,path,stage_i,pos,cell):
 stem,thr,maxstep=STAGES[stage_i]; lam=path[stage_i]
 x=["&CONTROL","   calculation='scf'","   restart_mode='from_scratch'",f"   prefix='Fe_bcc8_r8_{case}'","   pseudo_dir='../../pseudo'","   outdir='./qe_tmp'","   disk_io='low'","   verbosity='high'","   tstress=.true.","   tprnfor=.true.","/","&SYSTEM","   ibrav=0","   nat=8","   ntyp=8","   ecutwfc=71.0","   ecutrho=496.0","   occupations='smearing'","   smearing='fd'","   degauss=0.025334","   nosym=.true.","   noncolin=.true.","   lspinorb=.false.","   constrained_magnetization='atomic'",f"   lambda={lam:.8f}","   report=1"]
 for i,(theta,phi) in enumerate(ANGLES,1): x += [f"   starting_magnetization({i})=0.12500000",f"   angle1({i})={theta:.10f}",f"   angle2({i})={phi:.10f}"]
 x += ["/","&ELECTRONS",f"   startingpot='{'atomic' if stage_i==0 else 'file'}'",f"   startingwfc='{'atomic+random' if stage_i==0 else 'file'}'",f"   conv_thr={thr:.1e}",f"   electron_maxstep={maxstep}","   scf_must_converge=.true.","   mixing_mode='local-TF'","   mixing_beta=0.010d0","   mixing_ndim=20","   diagonalization='cg'","   diago_thr_init=1.0d-3","   diago_cg_maxiter=200","/","ATOMIC_SPECIES"]
 x += [f'{label} 55.845 {PSEUDO}' for label in LABELS] + pos + cell + ['K_POINTS automatic','1 1 1 0 0 0','']
 return '\n'.join(x)

def main():
 for d in ('cases','logs','pseudo'): (ROOT/d).mkdir(parents=True,exist_ok=True)
 if not (ROOT/'pseudo'/PSEUDO).is_file(): raise FileNotFoundError(f'Missing bundled pseudopotential: {PSEUDO}')
 pos,cell=POSITIONS,CELL; rows=[]
 for idx,(label,path) in enumerate(CASES):
  case=f'a2p39_m2_fixed_{label}'; folder=ROOT/'cases'/case; folder.mkdir(exist_ok=True)
  if list(folder.glob('*.out')) or (folder/'qe_tmp').exists(): raise RuntimeError(f'Results exist in {folder}')
  for s,(stem,_,_) in enumerate(STAGES): (folder/f'{stem}.in').write_text(make_input(case,path,s,pos,cell))
  row={'serial_index':idx,'case':case,'lattice_a_A':2.39,'intended_target_Bohr':2.0,'starting_magnetization':0.125,'lambda_stage1':path[0],'lambda_stage2':path[1],'lambda_stage3':path[2],'qe_requirement':'patched qe-7.5 mcons=zv*starting_magnetization'}
  rows.append(row); (folder/'case_metadata.json').write_text(json.dumps(row,indent=2)+'\n')
 with (ROOT/'case_manifest.csv').open('w',newline='') as h: w=csv.DictWriter(h,fieldnames=list(rows[0]),lineterminator='\n'); w.writeheader(); w.writerows(rows)
 (ROOT/'cases.txt').write_text('\n'.join(r['case'] for r in rows)+'\n'); print(f'Generated {len(rows)} Round-8 cases')
if __name__=='__main__': main()
