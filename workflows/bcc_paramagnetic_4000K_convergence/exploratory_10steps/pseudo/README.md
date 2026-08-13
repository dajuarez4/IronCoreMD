# Pseudopotential staging directory

The large pseudopotential is not duplicated in this workflow. Before copying
the workflow to Jakar, stage the repository's existing Fe PAW file here:

```bash
cp ../../../../dataset/hcp/hcp_mag/pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF .
```

The QE input expects:

```text
pseudo/Fe.pbe-spn-kjpaw_psl.1.0.0.UPF
```
