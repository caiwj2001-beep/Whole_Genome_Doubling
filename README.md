# Whole-Genome Doubling and Cancer Evolution

Analysis code for the pan-cancer WGD study using the MSK-IMPACT 50K clinical sequencing cohort (53,654 tumors, 66 cancer types).



## Repository Contents

| File | Description |
|---|---|
| `01_data_loading_and_qc.R` | Data loading, quality control, pan-cancer WGD frequency |
| `02_genomic_correlates.R` | Driver gene associations, TP53/RB1, mutational signatures |
| `03_survival_analysis.R` | Kaplan-Meier, Cox regression, cancer-type-specific HRs |
| `04_figure_generation.R` | All 8 publication figures (PDF/PNG/TIFF) |

## Data Source

MSK-IMPACT 50K Clinical Sequencing Cohort (Bandlamudi et al., *Cancer Cell* 2026)
- cBioPortal: https://www.cbioportal.org/study/summary?id=msk_impact_50k_2026
- PMID: 41895280


## Author

Wenjie Cai (蔡文杰), Department of Radiation Oncology, First Hospital of Quanzhou Affiliated to Fujian Medical University
