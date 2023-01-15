# Network_comparison
This repository serves as a research archive for the mini-project "Comparison of Gaussian graphical models (GGM) and Directed Cyclic Graph (DCG) Models as Causal Discovery Tools", 
which is an extension of my Master's thesis [(mid-report PDF here)](https://github.com/KyuriP/Markup_manuscript/blob/main/main.pdf). 
This repository contains all necessary files to replicate the [original simulation study](https://rpubs.com/KyuriP/992071) 
and [the follow-up analyses](https://rpubs.com/KyuriP/992072).

## Contents
The file structure is shown as below.

```bash
├── LICENSE
├── Network_comparison.Rproj
├── code
│   ├── R
│   │   ├── CCD_fnc.R
│   │   ├── data_generating_fnc.R
│   │   ├── dsep_fnc.R
│   │   ├── equivset_fnc.R
│   │   ├── eval_metric_fnc.R
│   │   ├── plot_fnc.R
│   │   ├── searchAM_KP_fnc.R
│   │   └── variation_fnc.R
│   ├── empircal_example.R
│   ├── simulation.R
│   └── variation_DCG.R
├── data
│   ├── McNally.csv
│   ├── equiv4p.RData
│   ├── equiv4p_high.RData
│   ├── equiv5p.RData
│   ├── equiv5p_high.RData
│   ├── equiv6p.RData
│   └── equiv6p_high.RData
└── manuscript
    ├── apa.csl
    ├── img
    │   ├── CCDoverallsummary.png
    │   ├── CCDsummary.png
    │   ├── degreevariation.png
    │   ├── densityvariation.png
    │   └── truemodels.png
    ├── follow-up.html
    ├── follow-up.qmd
    ├── original_simulation.html
    ├── original_simulation.qmd
    ├── references.bib
    ├── style.css
    └── style2.css

5 directories, 35 files
```

### CODE
- `code` folder contains all the code used for the original simulation study and follow-up analyses.
- Three R scripts contains the code for simulating models in the original study and creating figures presented in the follow-up report.
  - `simulation.R`: code to simulate all six different models.
  - `variation_DCG.R`: code to create the density plots (Figure 2) and degree centrality plots (Figure 3) in the follow-up report.
  - `empircal_example.R`: code to create the PAG (Figure 4) and GGM (Figure 5) on empirical data in the follow-up report.
- `R` folder contains all supporting functions required to run the R scripts mentioned above.

### DATA
- `data` folder contains six equivalence classes of directed cyclic graphs (DCG) from each of the simulated models saved as `.Rdata`.
- Also, an example empirical data `McNally.csv` can be found in this folder.

### MANUSCRIPT
- `original_simulation.html` is the manuscript for the main analyses.
- `original_simulation.qmd` renders the main manuscript.
- `follow-up.html` is the manuscript for the follow-up analyses.
- `follow-up.qmd` renders the follow-up manuscript.
- other supporting files/folders:
  - `img`: img folder contains all necessary image files to render the manuscripts.
  - `references.bib`: BibTex for references cited in the manuscripts.
  - `apa.csl`: csl file to format the references as per APA guidelines.
  - `.css` files: for the minor styling and layout of the manuscripts.
  
### RECOMMENDATION
1) Check the manuscripts. (`original_simulation.html` and `follow-up.html`)
2) To reproduce the original simulation study, either check the code presented in `original_simulation.html` (click `See code here`) or check the [quarto document](https://github.com/KyuriP/Network_comparison/blob/main/manuscript/original_simulation.qmd) where all code is self-contained.
3) To reproduce the follow-up analyses, run the following R scripts: `variation_DCG.R` and `empircal_example.R` for the first sub-analysis and the second sub-analysis, respectively. 
4) For any further information, see the detailed description attached to each R script.

---
For any help with the files in this archive, please contact Kyuri Park (k.park@uu.nl). 
