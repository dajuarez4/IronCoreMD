#!/usr/bin/env python3
"""Plot QE bands.dat for both snapshots, using their SCF Fermi energies."""
import os,json,re,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/mpl-round24-25-dashboard')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent

def read_bands(path):
    text=path.read_text();nb=int(re.search(r'nbnd\s*=\s*(\d+)',text)[1]);nk=int(re.search(r'nks\s*=\s*(\d+)',text)[1])
    vals=np.fromstring(text.split('/',1)[1].replace('D','E'),sep=' ')
    if vals.size!=nk*(nb+3):raise ValueError(f'Incomplete bands file {path}')
    values=vals.reshape(nk,nb+3)
    return values[:,:3],values[:,3:]

def main():
    fig,axes=plt.subplots(1,2,figsize=(14,6),layout='constrained');ready=0
    for n,ax in zip([24,25],axes):
        case=ROOT/f'round{n}';ax.set_title(f'Round{n} — electronic bands')
        if not (case/'post.complete').exists():
            ax.text(.5,.5,'Awaiting completed band calculation',ha='center',transform=ax.transAxes);continue
        coords,energies=read_bands(case/'bands.dat');path=json.loads((case/'kpath.json').read_text())
        xml=ET.parse(case/'snapshot_scf.xml').getroot()
        for e in xml.iter():e.tag=e.tag.split('}')[-1]
        ef=float(xml.findtext('output/band_structure/fermi_energy'))*27.211386245988
        x=np.array(path['distance_inv_A']);assert len(x)==len(energies)
        assert np.isfinite(energies).all()
        relative=energies-ef
        visible=(relative.min(axis=0)<5)&(relative.max(axis=0)>-8)
        ax.plot(x,relative[:,visible],lw=.55,color='#2563eb' if n==24 else '#dc751e',alpha=.5)
        ticks=np.array(path['tick_indices']);ax.set_xticks(x[ticks],path['labels']);ax.set(xlim=(x[0],x[-1]),ylim=(-8,5),ylabel='E − SCF EF (eV)',xlabel='Supercell path (distance in Å⁻¹)')
        for tick in x[ticks]:ax.axvline(tick,color='gray',lw=.4)
        ax.axhline(0,color='black',ls='--',lw=.7)
        ax.set_title(f'Round{n}: {energies.shape[1]} spinor bands / {len(x)} k points')
        np.savez_compressed(case/'electronic_bands.npz',distance_inv_A=x,energies_eV=energies,energies_minus_EF_eV=relative,SCF_fermi_eV=ef,kpoints_crystal=path['kpoints_crystal'])
        ready+=1
    if not ready:raise SystemExit('No completed band results yet; no plot generated.')
    fig.suptitle('Frozen-snapshot, energy-ordered supercell bands\nNo unfolding; each panel uses its own SCF Fermi energy',fontsize=12)
    (ROOT/'dashboard').mkdir(exist_ok=True);fig.savefig(ROOT/'dashboard/electronic_bands.png',dpi=180)
    print('Saved dashboard/electronic_bands.png')
if __name__=='__main__':main()
