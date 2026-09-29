# HCP Fe data

This directory retains the existing AIMD archives and source data. TDEP calculations and tooling are stored outside the repository:

`/Users/dajuarez4/Documents/Fe/dataset/hcp`

- Corrected fits: `tdep_corrected/`
- Workflow scripts: `codes/tdep_workflow/`
- Regression tests: `tests/test_hcp_reference.py`

All 21 existing AIMD archives were refitted with consistent HCP TDEP reference cells and explicit reciprocal paths. The raw simulations are unchanged. Current summary figures are in the repository's `assets/` directory; see [the correction report](TDEP_REFERENCE_CORRECTION.md) for results, limitations, and reproduction commands.
