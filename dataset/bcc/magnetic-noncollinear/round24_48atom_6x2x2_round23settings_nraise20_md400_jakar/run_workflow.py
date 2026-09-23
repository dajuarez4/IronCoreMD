#!/usr/bin/env python3
"""Run/continue the staged SCF and 400-step MD within one allocation."""
import argparse,fcntl,re,subprocess,time
from pathlib import Path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seconds',type=int,default=165600,help='QE time budget for this allocation (default 46 hours)')
    p.add_argument('command',nargs=argparse.REMAINDER)
    args=p.parse_args();command=args.command
    if command[:1]==['--']:command=command[1:]
    if not command:raise SystemExit('Provide launcher and pw.x after --')
    root=Path(__file__).resolve().parent
    case=root/'cases'/(root/'cases.txt').read_text().strip()
    with (case/'.workflow.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        deadline=time.monotonic()+args.seconds
        for stem in ['stage1_4e-3','stage2_1e-3','stage3_1e-4','md_400steps']:
            old=sorted(case.glob(stem+'.run*.out'))
            histories=[f.read_text(errors='replace') for f in old]
            # A hard-killed run is not a verified restart checkpoint.
            if any('JOB DONE' not in t or 'Error in routine' in t or 'convergence NOT achieved' in t for t in histories):
                raise RuntimeError(f'Inspect unsuccessful {stem} output before restarting; automatic continuation refused')
            md=stem=='md_400steps'
            completed=sum(t.count('Entering Dynamics:') for t in histories) if md else 0
            if md and completed>=400:
                print(f'MD complete: {completed} steps',flush=True);return
            if not md and histories and 'convergence has been achieved' in histories[-1]:continue
            budget=int(deadline-time.monotonic())
            if budget<300:
                print('Time budget exhausted; submit the same sbatch file again.',flush=True);return
            text=(case/(stem+'.in')).read_text()
            if old:
                if not (case/'qe_tmp/Fe_bcc48_r24_conv1e4_beta1e2.save').is_dir():
                    raise RuntimeError('Missing checkpoint directory')
                text=text.replace("restart_mode='from_scratch'","restart_mode='restart'")
            if md:text=re.sub(r'nstep=400',f'nstep={400-completed}',text)
            text=text.replace('&CONTROL','&CONTROL\n   max_seconds='+str(budget),1)
            attempt=len(old)+1
            inp=case/f'{stem}.run{attempt:03d}.in';out=case/f'{stem}.run{attempt:03d}.out'
            if inp.exists() or out.exists():raise RuntimeError('Existing attempt files; inspect before retrying')
            inp.write_text(text)
            print(f'Running {inp.name}, budget {budget}s; previous MD steps {completed}',flush=True)
            with out.open('x') as stream:
                result=subprocess.run(command+['-in',inp.name],cwd=case,stdout=stream)
            output=out.read_text(errors='replace')
            if result.returncode or 'JOB DONE' not in output or 'Error in routine' in output or 'convergence NOT achieved' in output:
                raise RuntimeError(f'{out.name} did not finish cleanly; inspect output')
            finished=(completed+output.count('Entering Dynamics:')>=400) if md else 'convergence has been achieved' in output
            if not finished:
                print('Clean time-limit stop. Submit the same sbatch file to continue.',flush=True);return
        print('All preparation stages and 400 MD steps complete.',flush=True)
if __name__=='__main__':main()
