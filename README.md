# Biochar_ASI
An open-source decision-support framework implementing a two-stage evaluation process—combining regulatory chemical safety screening and multi-criteria ranking process —to identify fit-for-purpose biochars for targeted agricultural systems.

📋 Overview

This repository contains the complete analytical R script and supplementary dataset associated with the research paper evaluating 36 biochars derived from green wastes (GW) and composting residues (CR) produced across varying pyrolysis conditions (summarized by the Heat-Treatment Severity Index, $HTSI = T \times t$).
The computational pipeline calculates the Aggregated Suitability Index (ASI) by synthesizing multi-dimensional physical, chemical, and environmental safety traits:
  pH ($w = 0.290$) 
  EC ($w = 0.290$)
  Biochar Stability Index (BSI) ($w = 0.145$)
  Water-Holding Capacity (WHC) ($w = 0.097$)
  Cation Exchange Capacity (CEC) ($w = 0.072$)
  Nutrient Content Index (NUTC) ($w = 0.058$)
  Heavy Metal Load Index (HMC) ($w = 0.048$)

🛠️ Framework Workflow

Stage 1: Regulatory Safety Screening: Filters biochars against strict trace-element thresholds (EU REACH / EBC / IBI standards).
Stage 2: Agronomic Trait Scoring: Maps physical and chemical traits into normalized scores $[0, 1]$ using piecewise agronomic threshold functions (pH, EC, BSI, WHC, CEC) and min-max relative normalization (NUTC, HMC).
Stage 3: ASI Aggregation: Multiplies normalized scores by criteria weights ($CR < 0.1$, $\sum w_i = 1.000$) to yield the composite Aggregated Suitability Index (ASI) for target agricultural contexts (e.g., acidic soils subject to drought risk).

📂 Repository Structure
```
├── data/
│   ├── raw_biochar_dataset.csv       # Raw physicochemical and trace-element properties
│   └── processed_scores.csv          # Calculated normalized scores and ASI ranks
├── biochar_ASI.R                     # Single master R script containing all analyses
├── outputs/
│   ├── figures/                      # Exported vector PDFs (e.g., Figure 1 flowchart)
│   └── tables/                       # Ranked output tables (e.g., Table S3 top biochars)
├── renv.lock                         # Package dependency lockfile for precise reproducibility
├── LICENSE                           # MIT License
└── README.md                         # Repository documentation

```
🚀 Getting Started & Execution
Prerequisites
R (version 4.4.3 or higher recommended)
RStudio
Reproduce the Environment
Restore the exact package versions used in the publication via `renv`:
```
# Open the project in RStudio and run:
if (!require("renv")) install.packages("renv")
renv::restore()

```
Run the Analysis
All data cleaning, safety screening, multivariate profiling, AHP weighting, ASI calculation, and vector figure exports are compiled into a single master R script.
Execute `biochar_analysis_master.R` inside RStudio, or from your system terminal:
```
Rscript biochar_ASI.R

```
⚙️ Adapting for Alternative Scenarios
While the default script is parameterized for acidic soils with moderate-to-high drought risk, the code is fully adaptable for other pedoclimatic contexts (e.g., alkaline soils in Mediterranean agricultural systems).

To customize the framework:
  Open `biochar_analysis_master.R`.
  Locate the Weights Vector block at the top of the script.
  Modify the priority vector ($w_i$) or piecewise scoring threshold functions to align with your target soil conditions.
  Re-run the script to update the composite ASI rankings.

📖 Citation
If you use this codebase or multi-criteria decision framework in your research, please cite:
> **Giannopoulos, G. et al.,** (2026). *Towards Fit-for-Purpose Biochar: A Multi‑Criteria Framework to Select Appropriate Biochar for Agricultural Applications* Journal Name. \[DOI Link\]

📄 License
Distributed under the MIT License. See `LICENSE` for more information.
```
eof

```
