#!/usr/bin/env python3
"""Incrementally update BCC TDEP curves, writing all calculations outside the repo."""
from pathlib import Path
import argparse,csv,fcntl,hashlib,json,os,re,shutil,subprocess,sys,tempfile
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/tmp/matplotlib')
from npz_to_tdep import write_tdep_folder
from run_harmonic_tdep import ensure_forceconstant_link,run_logged,refresh_single_temperature_plots,refresh_comparison_plot
from tdep_common import find_tdep_root,discover_npz_files,read_free_energy,read_u0_second_order,classify_free_energy
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
INPUTS=['infile.ucposcar','infile.ssposcar','infile.positions','infile.forces','infile.stat','infile.meta','infile.qpoints_dispersion']
OUTPUTS=['outfile.forceconstant','outfile.U0','outfile.free_energy','outfile.dispersion_relations','outfile.dispersion_relations.gnuplot','outfile.phonon_dos']
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--temperature',type=int,default=4000)
 p.add_argument('--source-dir',type=Path,default=REPO/'dataset/bcc/non-mag')
 p.add_argument('--output-dir',type=Path,default=REPO.parent/'dataset/bcc/non-mag')
 p.add_argument('--tdep-root',type=Path)
 p.add_argument('--skip',type=int,default=0)
 p.add_argument('--every',type=int,default=1)
 p.add_argument('--rc2',type=float,default=5.)
 p.add_argument('--dos-qgrid',type=int,nargs=3,default=[32,32,32])
 p.add_argument('--force',action='store_true',help='Recompute all selected points')
 args=p.parse_args();source=args.source_dir.resolve();dest=args.output_dir.resolve()
 if dest==REPO or REPO in dest.parents:raise SystemExit('Output must be outside IronCoreMD; no TDEP folders will be created in the repo.')
 if args.skip<0 or args.every<1:raise SystemExit('Require skip>=0 and every>=1')
 dest.mkdir(parents=True,exist_ok=True);temp=str(args.temperature)
 lock=(dest/f'.update_{temp}K.lock').open('w')
 fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 # Keep canonical NPZ copies outside the repo so all plotting utilities share one dataset.
 files=discover_npz_files(source,temp,'bcc')
 if not files:raise SystemExit('No matching NPZ inputs')
 for f in files:
  target=dest/f.name
  if f!=target and (not target.exists() or digest(f)!=digest(target)):shutil.copy2(f,target)
 binaries=find_tdep_root(args.tdep_root)
 extract=binaries/'extract_forceconstants/extract_forceconstants';disp=binaries/'phonon_dispersion_relations/phonon_dispersion_relations'
 settings={'temperature':args.temperature,'skip':args.skip,'every':args.every,'rc2':args.rc2,'qgrid':args.dos_qgrid,'extract_sha256':digest(extract),'dispersion_sha256':digest(disp),'converter_sha256':digest(HERE/'npz_to_tdep.py'),'phases_sha256':digest(HERE/'tdep_phases.py')}
 records=[];failures=[]
 for f in files:
  npz=dest/f.name;folder=dest/('tdep_'+f.stem);stamp=folder/'update_fingerprint.json';pending=folder/'.update_in_progress'
  key={'npz_sha256':digest(npz),'settings':settings}
  complete=all((folder/n).is_file() for n in OUTPUTS)
  prior=json.loads(stamp.read_text()) if stamp.exists() else {}
  action='cached'
  try:
   match=complete and not pending.exists() and prior.get('key')==key and prior.get('output_hashes')=={n:digest(folder/n) for n in OUTPUTS}
   if not match or args.force:
    # Validate legacy inputs byte-for-byte before adopting existing calculations.
    with tempfile.TemporaryDirectory(prefix='.tdep_check_',dir=dest) as td:
     generated=Path(td)
     write_tdep_folder(npz_path=npz,outdir=generated,phase='bcc',supercell=(4,4,4),temperature_override=args.temperature,skip=args.skip,every=args.every,max_frames=0,keep_invalid=False)
     old_options=(args.skip==0 and args.every==1 and args.rc2==5. and args.dos_qgrid==[32,32,32])
     legacy=(complete and not pending.exists() and not stamp.exists() and old_options and all((folder/n).exists() and digest(folder/n)==digest(generated/n) for n in INPUTS))
     if legacy and not args.force:
      action='adopted_existing';print('[reuse verified inputs]',f.name,flush=True)
     else:
      action='computed';folder.mkdir(exist_ok=True)
      pending.write_text('Recomputation has not completed. Do not adopt previous output files.\n')
      # Invalidate the cache first; unsuccessful calculations must not look current.
      if stamp.exists():stamp.unlink()
      for item in generated.iterdir():shutil.copy2(item,folder/item.name)
      print('[TDEP]',f.name,flush=True)
      run_logged([str(extract),'-rc2',str(args.rc2)],folder,folder/'extract_forceconstants.log')
      ensure_forceconstant_link(folder)
      run_logged([str(disp)],folder,folder/'phonon_dispersion_relations.log')
      run_logged([str(disp),'--dos','--qpoint_grid',*map(str,args.dos_qgrid),'--temperature',temp],folder,folder/f'free_energy_{temp}K.log')
    stamp.write_text(json.dumps({'key':key,'output_hashes':{n:digest(folder/n) for n in OUTPUTS},'origin':action},indent=2)+'\n')
    if pending.exists():pending.unlink()
   values=read_free_energy(folder/'outfile.free_energy');u0=read_u0_second_order(folder/'outfile.U0');status=classify_free_energy(*values,u0)
   records.append({'npz':str(npz),'folder':str(folder),'action':action,'free_energy_status':status,'F_vib_eV_atom':values[1],'U0_eV_atom':u0,'F_total_eV_atom':values[1]+u0})
   if action!='cached' or not (folder/'phonon_dispersion_and_dos.png').exists():
    run_logged([sys.executable,str(HERE/'plot_single_tdep_phonons.py'),str(folder),'--title',f'BCC Fe non-magnetic: {f.stem}'],folder,folder/'plot_phonons.log')
  except Exception as e:
   failures.append({'npz':str(npz),'error':str(e)});print('[FAILED]',f.name,e,flush=True)
 report=dest/f'update_{temp}K_report.json';report.write_text(json.dumps({'settings':settings,'source_dir':str(source),'output_dir':str(dest),'points':records,'failures':failures},indent=2)+'\n')
 if failures:raise SystemExit(f'{len(failures)} points failed; inspect {report}. Existing combined plots retained.')
 refresh_single_temperature_plots(dest,'bcc',temp)
 refresh_comparison_plot(dest,'bcc',[])
 print('Updated:',report,flush=True)
if __name__=='__main__':main()
