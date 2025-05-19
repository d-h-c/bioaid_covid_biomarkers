
This is code for the analysis performed by Heather R. Jackson for the paper:

# Differentiation of COVID-19 from other emergency infectious disease presentations using whole blood transcriptomics then rapid qPCR: a case-control and observational cohort study

Ho Kwong Li, Heather R. Jackson, Luca Miglietta, Dominic Habgood-Coote, Ewurabena Mills, Ravi Mehta, Ali Hamady, Anna Haber, Maisarah Amran, Robert Hammond, Dominique Arancon, Graham Cooke, Mahdad Noursadeghi, Peter J.M. Openshaw, Jesus Rodriguez-Manzano, Myrsini Kaforou, Shiranee Sriskandan

Correspondence: Jesus Rodriguez-Manzano j.rodriguez-manzano@imperial.ac.uk; Myrsini Kaforou m.kaforou@imperial.ac.uk; Shiranee Sriskandan s.sriskandan@imperial.ac.uk

## Code:

- Signature_Normalisation_Discovery.R -- the discovery
- RT_PCR_Results.R -- RT-PCR validation
- McClain_validation.R -- validation in Mcclain RNAseq dataset

## data:

### discovery:

- discovery_cohort_phenotypes.csv - Phenotypes in the discovery cohort
- betas_discovery.csv  -- model coefficients in the discovery RNA-Seq data
- genes_remove.csv -- non-coding genes removed from the discovery

### Validation:

- CT_data.csv - Cycle threshold values for the RT-PCR validation cohort
- RT-PCR_validation_phenotypes.csv - Phenotypes in the RT-PCR validation cohort
- crpwcc_validation.csv - C-reactive protein and white cell count in the validation cohort
- samples_no_gapdh.csv -- samples removed from the RT-PCR validation due to the failure of GAPDH reaction



