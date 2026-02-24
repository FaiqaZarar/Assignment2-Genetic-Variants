# Genetic Variant Analysis of Selected Genetic Disorders and Rare Diseases

<div align="center">

[![NUST](https://img.shields.io/badge/NUST-SINES-003087?style=flat-square&logoColor=white)](https://www.nust.edu.pk)
[![Course](https://img.shields.io/badge/Course-Special%20Topics%20in%20Bioinformatics-0057A8?style=flat-square)](https://www.nust.edu.pk)
[![Assembly](https://img.shields.io/badge/Reference%20Genome-GRCh38%2Fhg38-2E8B57?style=flat-square)](https://genome.ucsc.edu)
[![ClinVar](https://img.shields.io/badge/Database-ClinVar-E87722?style=flat-square)](https://www.ncbi.nlm.nih.gov/clinvar/)
[![License](https://img.shields.io/badge/License-Academic-lightgrey?style=flat-square)](#)

</div>

---

## Project Information

| Field | Details |
|:------|:--------|
| **Institution** | National University of Sciences and Technology (NUST) |
| **Department** | School of Interdisciplinary Engineering & Sciences (SINES) |
| **Course** | Special Topics in Bioinformatics |
| **Assignment** | #2 — Genetic Disorders & Rare Diseases Variant Analysis |
| **Authors** | Faiqa Zarar Noor (471543), Pukhraj Tahir (467407) |
| **Submission Date** | 1st March 2026 |

---

## Table of Contents

- [Objective](#objective)
- [Diseases and Variants Selected](#diseases-and-variants-selected)
- [Methodology](#methodology)
- [Disease 1: Cystic Fibrosis](#disease-1-cystic-fibrosis)
- [Disease 2: Sickle Cell Anemia](#disease-2-sickle-cell-anemia)
- [Disease 3: Huntingtons Disease](#disease-3-huntingtons-disease)
- [VCF File and Annotation](#vcf-file-and-annotation)
- [Repository Structure](#repository-structure)
- [Tools and Databases](#tools-and-databases)
- [References](#references)

---

## Objective

The objective of this assignment is to identify and characterize pathogenic genetic variants associated with three selected diseases using established clinical bioinformatics databases and tools. For each disease, the analysis covers variant identification from ClinVar, phenotypic characterization from OMIM, in silico pathogenicity scoring from UCSC Genome Browser (AlphaMissense and REVEL), variant classification using ACMG/AMP 2015 guidelines, and VCF file creation with Ensembl VEP annotation.

---

## Diseases and Variants Selected

| # | Disease | Category | Gene | Variant | ClinVar ID | ACMG Classification |
|:--|:--------|:---------|:-----|:--------|:-----------|:--------------------|
| 1 | Cystic Fibrosis | Genetic Disorder | *CFTR* | NM_000492.4:c.1521_1523delCTT (p.Phe508del) | VCV000007107 | Pathogenic |
| 2 | Sickle Cell Anemia | Genetic Disorder | *HBB* | NM_000518.5:c.20A>T (p.Glu7Val) | VCV000015280 | Pathogenic |
| 3 | Huntington's Disease | Rare Disease | *HTT* | NM_002111.8:c.53CAG[>35] | VCV000525889 | Pathogenic |

---

## Methodology

The following workflow was applied consistently across all three diseases:

```
1. ClinVar Search
   └── Searched disease name → filtered by Pathogenic significance
   └── Recorded: HGVS notation, dbSNP ID, ClinVar ID, chromosome position,
       mutation type, review status, and clinical significance

2. Literature Review (Explanation Field)
   └── Read ClinVar submissions, linked publications, and study observations
   └── Summarized molecular mechanism and clinical evidence

3. OMIM Phenotype Collection
   └── Retrieved phenotype, inheritance pattern, and clinical features
   └── Recorded OMIM accession number

4. UCSC Genome Browser — Pathogenicity Scoring
   └── Navigated to GRCh38 chromosome position for each variant
   └── Enabled AlphaMissense track (AI-based pathogenicity)
   └── Enabled REVEL Scores track (ensemble pathogenicity)
   └── Captured screenshots for Excel sheet

5. ACMG/AMP Classification
   └── Applied Richards et al. 2015 criteria
   └── Recorded supporting evidence codes (PVS1, PS, PM, PP)

6. VCF File Creation
   └── Compiled all 3 variants into standard VCFv4.2 format
   └── Annotated using Ensembl Variant Effect Predictor (VEP)
```

---

## Disease 1: Cystic Fibrosis

### Variant Details

| Field | Value |
|:------|:------|
| **Gene** | *CFTR* (Cystic Fibrosis Transmembrane Conductance Regulator) |
| **HGVS Notation** | NM_000492.4(CFTR):c.1521_1523delCTT |
| **Protein Change** | p.Phe508del |
| **Mutation Type** | 3 bp in-frame deletion |
| **Chromosome Position** | chr7:117,548,628 (GRCh38) |
| **dbSNP ID** | rs113993960 |
| **ClinVar ID** | VCV000007107 |
| **Review Status** | Practice guideline (highest confidence) |
| **Clinical Significance** | Pathogenic |

### Molecular Mechanism

The F508del variant results in deletion of phenylalanine at position 508 of the CFTR protein. The misfolded protein is retained in the endoplasmic reticulum and targeted for proteasomal degradation, preventing it from reaching the apical membrane of epithelial cells. This causes near-complete loss of chloride and bicarbonate channel activity, leading to thickened secretions in the lungs, pancreas, liver, and intestine. F508del accounts for approximately 70% of CF alleles globally and is present in at least one allele in ~90% of CF patients.

### OMIM Phenotype — #219700

Autosomal recessive multisystem disorder characterized by chronic obstructive pulmonary disease, recurrent pulmonary infections (*Pseudomonas aeruginosa*, *Staphylococcus aureus*), bronchiectasis, exocrine pancreatic insufficiency, elevated sweat chloride concentration (>60 mmol/L), and male infertility due to congenital bilateral absence of the vas deferens (CBAVD). Meconium ileus presents in ~15% of affected neonates. With modern CFTR modulator therapy (elexacaftor/tezacaftor/ivacaftor), median predicted survival now exceeds 50 years.

### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
|:----------|:---------|
| **PVS1** | Loss-of-function variant in *CFTR* where LOF is an established disease mechanism |
| **PS3** | Functional studies demonstrate absent chloride channel activity (Dalemans et al., 1991) |
| **PP5** | >500 independent ClinVar submissions, all classified as Pathogenic |
| **PM3** | Detected in trans with other pathogenic *CFTR* variants in affected individuals |

---

## Disease 2: Sickle Cell Anemia

### Variant Details

| Field | Value |
|:------|:------|
| **Gene** | *HBB* (Hemoglobin Subunit Beta) |
| **HGVS Notation** | NM_000518.5(HBB):c.20A>T |
| **Protein Change** | p.Glu7Val |
| **Mutation Type** | Single nucleotide variant — missense |
| **Chromosome Position** | chr11:5,227,002 (GRCh38) |
| **dbSNP ID** | rs334 |
| **ClinVar ID** | VCV000015280 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |

### Molecular Mechanism

The c.20A>T transversion substitutes valine for glutamic acid at position 7 of the beta-globin subunit, producing Hemoglobin S (HbS). Unlike glutamic acid, valine is hydrophobic and, under deoxygenated conditions, promotes intermolecular polymerization of HbS tetramers into rigid rod-like fibers. These fibers distort erythrocytes into the characteristic sickle morphology, reducing deformability, increasing adhesion to the vascular endothelium, and causing microvascular occlusion. Heterozygous carriers (HbAS) are largely asymptomatic and demonstrate resistance to severe *Plasmodium falciparum* malaria, explaining the high allele frequency in malaria-endemic regions.

### OMIM Phenotype — #603903

Autosomal recessive hemoglobinopathy presenting with chronic hemolytic anemia (hemoglobin ~6–9 g/dL), episodic vaso-occlusive pain crises, acute chest syndrome, ischemic stroke (occurring in ~11% of children without prophylaxis), splenic sequestration, progressive renal impairment, and pulmonary hypertension. Management includes hydroxyurea, chronic transfusion therapy, and hematopoietic stem cell transplantation. Recently approved gene therapies — Casgevy (CRISPR/Cas9-based) and Lyfgenia (lentiviral vector) — offer potential curative options.

### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
|:----------|:---------|
| **PS1** | Identical amino acid change (Glu→Val) previously established as pathogenic |
| **PS3** | Functional studies confirm HbS polymerization and erythrocyte sickling under deoxygenation |
| **PS4** | Variant prevalence significantly elevated in SCA cohorts compared to general population |
| **PM1** | Variant located in the critical oxygen-binding domain of the beta-globin protein |
| **PP5** | Classified Pathogenic by multiple independent expert clinical sources |

---

## Disease 3: Huntington's Disease

### Variant Details

| Field | Value |
|:------|:------|
| **Gene** | *HTT* (Huntingtin) |
| **HGVS Notation** | NM_002111.8(HTT):c.53CAG[>35] |
| **Variant Type** | CAG trinucleotide repeat expansion — exon 1 |
| **Chromosome Position** | chr4:3,074,877 (GRCh38) |
| **dbSNP ID** | rs193922927 |
| **ClinVar ID** | VCV000525889 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |
| **Normal Range** | 10–35 CAG repeats |
| **Reduced Penetrance** | 36–39 CAG repeats |
| **Full Penetrance** | ≥40 CAG repeats |

### Molecular Mechanism

Huntington's disease arises from the expansion of a CAG repeat tract in exon 1 of the *HTT* gene, encoding a polyglutamine (polyQ) stretch in the huntingtin protein. When the polyQ tract exceeds 35 residues, the protein misfolds, oligomerizes, and forms cytotoxic intranuclear inclusions. Mutant huntingtin disrupts transcriptional regulation, impairs mitochondrial function, and interferes with vesicular trafficking and synaptic signaling, leading to preferential degeneration of medium spiny neurons in the striatum and cortical neurons. Age of onset is inversely correlated with repeat length. Paternal transmission is associated with intergenerational repeat instability and anticipation.

### OMIM Phenotype — #143100

Autosomal dominant progressive neurodegenerative disorder with a clinical triad of: (1) motor dysfunction — chorea, dystonia, bradykinesia, and oculomotor abnormalities; (2) cognitive decline — executive dysfunction, memory impairment, and dementia; (3) neuropsychiatric symptoms — depression, irritability, obsessive-compulsive behavior, and psychosis. Typical age of onset is 30–50 years. Disease duration averages 15–20 years after symptom onset. No disease-modifying therapy is currently approved; management is symptomatic (tetrabenazine or deutetrabenazine for chorea).

### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
|:----------|:---------|
| **PVS1** | Repeat expansion in *HTT* where this mechanism is the established cause of HD |
| **PS4** | Repeat expansion found exclusively in HD-affected individuals (≥40 repeats) |
| **PS3** | Functional studies confirm toxic polyQ gain-of-function and neurodegeneration |
| **PP1** | Robust cosegregation with HD across multiple multigenerational pedigrees |
| **PP5** | Consistent Pathogenic classification across all reputable clinical laboratories |

---

## VCF File and Annotation

### VCF File (variants.vcf)

```
##fileformat=VCFv4.2
##fileDate=20260223
##source=Assignment2_GeneticVariantAnalysis
##reference=GRCh38/hg38
##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC gene symbol">
##INFO=<ID=DISEASE,Number=1,Type=String,Description="Associated disease">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance">
##INFO=<ID=CLNID,Number=1,Type=String,Description="ClinVar accession">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM  POS        ID           REF  ALT        QUAL  FILTER  INFO
7       117548628  rs113993960  CTT  .          .     PASS    GENE=CFTR;DISEASE=Cystic_Fibrosis;CLNSIG=Pathogenic;CLNID=VCV000007107
11      5227002    rs334        A    T          .     PASS    GENE=HBB;DISEASE=Sickle_Cell_Anemia;CLNSIG=Pathogenic;CLNID=VCV000015280
4       3074877    rs193922927  CAG  CAGCAGCAG  .     PASS    GENE=HTT;DISEASE=Huntingtons_Disease;CLNSIG=Pathogenic;CLNID=VCV000525889
```

### Ensembl VEP Annotation Summary

The VCF file was submitted to Ensembl Variant Effect Predictor (GRCh38) at https://asia.ensembl.org/Tools/VEP.

| Metric | Result |
|:-------|:-------|
| Variants processed | 2 |
| Variants filtered out | 0 |
| Novel / existing variants | 1 (50%) / 1 (50%) |
| Overlapped genes | 7 |
| Overlapped transcripts | 25 |
| Overlapped regulatory features | 0 |

**Key consequences identified:**

| Variant | Gene | Primary Consequence | Biotype | Exon |
|:--------|:-----|:--------------------|:--------|:-----|
| rs193922927 | HTT | inframe_insertion | protein_coding | 1/67 |
| rs193922927 | HTT-AS | upstream_gene_variant | lncRNA | — |
| rs334 | HBB | synonymous_variant | protein_coding | 1/3 |

---

## Repository Structure

```
Assignment2-Genetic-Variants/
├── Assignment2.xlsx          # Main data sheet: variant details, OMIM phenotypes,
│                             # UCSC AlphaMissense & REVEL screenshots, ACMG criteria
├── variants.vcf              # VCFv4.2 file containing all 3 pathogenic variants
├── Assignment2_Report.docx   # Full written report with figures and analysis
└── README.md                 # This file
```

---

## Tools and Databases

| Tool / Database | Version / Build | Purpose |
|:----------------|:----------------|:--------|
| NCBI ClinVar | Accessed February 2026 | Variant identification, clinical significance, review status |
| OMIM | Accessed February 2026 | Phenotype, inheritance, clinical features |
| UCSC Genome Browser | GRCh38/hg38 | AlphaMissense and REVEL pathogenicity track visualization |
| AlphaMissense | Google DeepMind (2023) | AI-based proteome-wide missense pathogenicity prediction |
| REVEL | v1.3 | Ensemble-based missense variant pathogenicity scoring |
| Ensembl VEP | GRCh38.p14 | VCF annotation, consequence and transcript overlap prediction |
| ACMG/AMP Guidelines | Richards et al. 2015 | Variant pathogenicity classification framework |

---

## References

1. Richards S, et al. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine*, 17(5), 405–423. https://doi.org/10.1038/gim.2015.30

2. Riordan JR, et al. (1989). Identification of the cystic fibrosis gene: cloning and characterization of complementary DNA. *Science*, 245(4922), 1066–1073. https://doi.org/10.1126/science.2475911

3. Cutting GR (2015). Cystic fibrosis genetics: from molecular understanding to clinical application. *Nature Reviews Genetics*, 16(1), 45–56. https://doi.org/10.1038/nrg3849

4. Ingram VM (1956). A specific chemical difference between the globins of normal human and sickle-cell anaemia haemoglobin. *Nature*, 178(4537), 792–794. https://doi.org/10.1038/178792a0

5. Pauling L, et al. (1949). Sickle cell anemia, a molecular disease. *Science*, 110(2865), 543–548. https://doi.org/10.1126/science.110.2865.543

6. The Huntington's Disease Collaborative Research Group (1993). A novel gene containing a trinucleotide repeat that is expanded and unstable on Huntington's disease chromosomes. *Cell*, 72(6), 971–983. https://doi.org/10.1016/0092-8674(93)90585-E

7. Bates GP, et al. (2015). Huntington disease. *Nature Reviews Disease Primers*, 1, 15005. https://doi.org/10.1038/nrdp.2015.5

8. Cheng J, et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*, 381(6660), eadg7492. https://doi.org/10.1126/science.adg7492

---

*National University of Sciences and Technology — School of Interdisciplinary Engineering & Sciences — 2026*
