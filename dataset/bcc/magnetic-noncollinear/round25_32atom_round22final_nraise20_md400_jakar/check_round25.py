#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parent
case=root/'cases'/(root/'cases.txt').read_text().strip()
steps=0
for p in sorted(case.glob('*.run*.out')):
    t=p.read_text(errors='replace');n=t.count('Entering Dynamics:');steps+=n
    state='CLEAN STOP' if 'JOB DONE' in t else 'RUNNING OR INTERRUPTED'
    if 'convergence NOT achieved' in t or 'Error in routine' in t:state='ERROR'
    print(p.name,n,'MD steps',state)
print(f'Total MD steps recorded: {steps}/400')
