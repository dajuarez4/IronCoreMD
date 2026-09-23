#!/usr/bin/env python3
"""Local structural and mocked-checkpoint tests; does not execute Quantum ESPRESSO."""
from pathlib import Path
import fcntl,json,re,shutil,subprocess,sys,tempfile,xml.etree.ElementTree as E
import numpy as np
from ase.io import read
from prepare_checkpoint import BOHR,xml_root
from plot_bands import read_bands
ROOT=Path(__file__).resolve().parent

def main():
 for n in [24,25]:
  case=ROOT/f'round{n}';config=json.loads((case/'config.json').read_text());nat=config['natoms']
  for name in ['01_snapshot_scf.in.example','02_bands.in.example']:
   atoms=read(case/name,format='espresso-in');assert len(atoms)==nat
   assert np.allclose(atoms.cell,np.diag(config['cell_A']))
   text=(case/name).read_text();template=(case/'md_template.in').read_text()
   s=re.search(r'(?ms)^&SYSTEM\n.*?^/',text)[0];s=re.sub(r'^.*nbnd=.*\n','',s,flags=re.M)
   assert s==re.search(r'(?ms)^&SYSTEM\n.*?^/',template)[0]
  path=json.loads((case/'kpath.example.json').read_text());assert len(path['kpoints_crystal'])==11
  assert path['tick_indices']==[0,2,4,6,8,10]
  subprocess.run(['bash','-n',str(case/'run_jakar.sbatch')],check=True)
 # Standalone mock copies of both rounds test importer and file isolation.
 with tempfile.TemporaryDirectory(prefix='bands_import_test_') as temp:
  root=Path(temp)/'package';root.mkdir()
  shutil.copy2(ROOT/'prepare_checkpoint.py',root/'prepare_checkpoint.py')
  shutil.copytree(ROOT/'pseudo',root/'pseudo')
  for n in [24,25]:
   case=root/f'round{n}';case.mkdir()
   for name in ['config.json','md_template.in']:shutil.copy2(ROOT/f'round{n}'/name,case/name)
   c=json.loads((case/'config.json').read_text());nat=c['natoms'];source=Path(temp)/f'source{n}';source.mkdir()
   (source/'.workflow.lock').touch();save=source/'qe_tmp'/(c['prefix']+'.save');save.mkdir(parents=True)
   xml=E.Element('espresso');out=E.SubElement(xml,'output');state=E.SubElement(out,'atomic_structure',nat=str(nat))
   cell=E.SubElement(state,'cell')
   for name,row in zip(['a1','a2','a3'],np.diag(c['cell_A'])/BOHR):E.SubElement(cell,name).text=' '.join(map(str,row))
   positions=E.SubElement(state,'atomic_positions')
   for i,label in enumerate(c['atom_labels']):E.SubElement(positions,'atom',name=label,index=str(i+1)).text=f'{i*.01} 0 0'
   bands=E.SubElement(out,'band_structure')
   for name,value in [('nbnd',922 if n==24 else 614),('nelec',16*nat),('noncolin','true')]:E.SubElement(bands,name).text=str(value)
   basis=E.SubElement(out,'basis_set')
   for name,val in [('ecutwfc',35.5),('ecutrho',248.0)]:E.SubElement(basis,name).text=str(val)
   E.ElementTree(xml).write(save/'data-file-schema.xml');(save/'charge-density.dat').write_bytes(b'MOCK density');(save/'wfc1.dat').write_bytes(b'MOCK wavefunctions')
   pseudo='Fe.pbe-spn-kjpaw_psl.1.0.0.UPF';shutil.copy2(ROOT/'pseudo'/pseudo,save/pseudo)
   cmd=[sys.executable,str(root/'prepare_checkpoint.py'),'--round',str(n),'--source-case',str(source)]
   with (source/'.workflow.lock').open('r') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    result=subprocess.run(cmd,capture_output=True,text=True);assert result.returncode!=0 and 'MD is active' in result.stderr
   subprocess.run(cmd,check=True)
   assert (case/'qe_tmp'/save.name/'wfc1.dat').read_bytes()==b'MOCK wavefunctions'
   (case/'qe_tmp'/save.name/'wfc1.dat').write_bytes(b'CHANGED copy')
   assert (save/'wfc1.dat').read_bytes()==b'MOCK wavefunctions'
   text=(case/'01_snapshot_scf.in').read_text();assert "startingwfc='file'" in text
   assert "startingwfc='atomic+random'" in (case/'02_bands.in').read_text()
   assert read(case/'01_snapshot_scf.in',format='espresso-in').positions[1,0]==float(f'{.01*BOHR:.12f}')
   assert subprocess.run(cmd,capture_output=True).returncode!=0
  test=Path(temp)/'bands.dat';test.write_text('&plot nbnd=2, nks=2 /\n0 0 0\n1 2\n.5 0 0\n3 4\n')
  k,e=read_bands(test);assert np.array_equal(e,[[1,2],[3,4]])
 print('PASS: example input parsing, parameter retention, path count, shell syntax, active-lock rejection, mock import, independent copies, overwrite refusal, band parser.')
if __name__=='__main__':main()
