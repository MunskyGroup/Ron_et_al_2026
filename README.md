# Ron et al. (2026) Analysis Codes

This repository contains the code, model files, parameter tables, and saved outputs used for the analyses and figures in the Ron et al. 2026 manuscript.

## Manuscript And SI Drafts

Current draft files are included in this repository:

- `Ron_et_al_2026_MAIN_DRAFT_7_1_2026.pdf`
- `Ron_et_al_2026_SI_DRAFT_7_1_2026.pdf`

These local draft references will later be replaced with a bioRxiv link.

## Install SSIT

This project relies on the SSIT package. See the official SSIT README for full setup details:

- https://github.com/MunskyGroup/SSIT/

Simple clone step:

```bash
git clone https://github.com/MunskyGroup/SSIT.git
```

After cloning, follow the remaining installation/setup instructions in the SSIT repository README. Specifically, navigate to the SSIT directory and run the command:
'''MATLAB
install(0,0)
'''

## Reproduce Analyses And Figures

Use `Main.m` as the main entry point to rerun analyses and regenerate manuscript figures.

From MATLAB, open this repository and run:

```matlab
Main
```

`Main.m` is intended to reproduce the analyses and figures reported in the manuscript.