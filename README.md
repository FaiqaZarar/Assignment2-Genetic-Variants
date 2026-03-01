# Genetic Variant Analysis of Selected Genetic Disorders and Rare Diseases

[![NUST](https://img.shields.io/badge/NUST-SINES-003087?style=flat-square&logoColor=white)](https://www.nust.edu.pk)
[![Course](https://img.shields.io/badge/Course-Special%20Topics%20in%20Bioinformatics-0057A8?style=flat-square)](https://www.nust.edu.pk)
[![Assembly](https://img.shields.io/badge/Reference%20Genome-GRCh38%2Fhg38-2E8B57?style=flat-square)](https://genome.ucsc.edu)
[![ClinVar](https://img.shields.io/badge/Database-ClinVar-E87722?style=flat-square)](https://www.ncbi.nlm.nih.gov/clinvar/)
[![License](https://img.shields.io/badge/License-Academic-lightgrey?style=flat-square)](#)

---

## Project Information

| Field | Details |
| --- | --- |
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
- [Genetic Disorders](#genetic-disorders)
  - [Disease 1: Sickle Cell Disease](#disease-1-sickle-cell-disease)
  - [Disease 2: Hereditary Hemochromatosis](#disease-2-hereditary-hemochromatosis)
  - [Disease 3: Phenylketonuria (PKU)](#disease-3-phenylketonuria-pku)
- [Rare Diseases](#rare-diseases)
  - [Disease 4: Hutchinson-Gilford Progeria Syndrome](#disease-4-hutchinson-gilford-progeria-syndrome-hgps)
  - [Disease 5: Fibrodysplasia Ossificans Progressiva](#disease-5-fibrodysplasia-ossificans-progressiva-fop)
  - [Disease 6: Alkaptonuria](#disease-6-alkaptonuria)
- [VCF File and Annotation](#vcf-file-and-annotation)
- [Repository Structure](#repository-structure)
- [Tools and Databases](#tools-and-databases)
- [References](#references)

---

## Objective

The objective of this assignment is to identify and characterize pathogenic genetic variants associated with six selected diseases — three genetic disorders and three rare diseases — using established clinical bioinformatics databases and tools. For each disease, the analysis covers variant identification from ClinVar, phenotypic characterization from OMIM, in silico pathogenicity scoring from UCSC Genome Browser (AlphaMissense and REVEL), variant classification using ACMG/AMP 2015 guidelines, and VCF file creation with Ensembl VEP annotation.

---

## Diseases and Variants Selected

| # | Category | Disease | Gene | Variant | ClinVar ID | ACMG Classification |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | Genetic Disorder | Sickle Cell Disease | *HBB* | NM_000518.5:c.20A>T (p.Glu6Val) | VCV000015280 | Pathogenic |
| 2 | Genetic Disorder | Hereditary Hemochromatosis | *HFE* | NM_000410.4:c.845G>A (p.Cys282Tyr) | VCV000003036 | Pathogenic |
| 3 | Genetic Disorder | Phenylketonuria (PKU) | *PAH* | NM_000277.3:c.1222C>T (p.Arg408Trp) | VCV000005345 | Pathogenic |
| 4 | Rare Disease | Hutchinson-Gilford Progeria Syndrome | *LMNA* | NM_170707.4:c.1824C>T (p.Gly608Gly) | VCV000041263 | Pathogenic |
| 5 | Rare Disease | Fibrodysplasia Ossificans Progressiva | *ACVR1* | NM_001105.6:c.617G>A (p.Arg206His) | VCV000013642 | Pathogenic |
| 6 | Rare Disease | Alkaptonuria | *HGD* | NM_000187.4:c.342+1G>A | VCV000036396 | Pathogenic |

---

## Methodology

The following workflow was applied consistently across all six diseases:

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
   └── For splice variants: used SpliceAI / MaxEntScan scores
   └── Captured screenshots for Excel sheet

5. ACMG/AMP Classification
   └── Applied Richards et al. 2015 criteria
   └── Recorded supporting evidence codes (PVS1, PS, PM, PP)

6. VCF File Creation
   └── Compiled all 6 variants into standard VCFv4.2 format
   └── Annotated using Ensembl Variant Effect Predictor (VEP)
```

---

## Genetic Disorders

### Disease 1: Sickle Cell Disease

| Field | Value |
| --- | --- |
| **Gene** | *HBB* (Hemoglobin Subunit Beta) |
| **HGVS Notation** | NM_000518.5(HBB):c.20A>T |
| **Protein Change** | p.Glu6Val |
| **Mutation Type** | Single nucleotide variant — missense |
| **Chromosome Position** | chr11:5,227,002 (GRCh38) |
| **dbSNP ID** | rs334 |
| **ClinVar ID** | VCV000015280 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |

#### Molecular Mechanism

The c.20A>T transversion substitutes valine for glutamic acid at position 6 of the beta-globin subunit, producing Hemoglobin S (HbS). Unlike glutamic acid, valine is hydrophobic and, under deoxygenated conditions, promotes intermolecular polymerization of HbS tetramers into rigid rod-like fibers. These fibers distort erythrocytes into the characteristic sickle morphology, reducing deformability, increasing adhesion to the vascular endothelium, and causing microvascular occlusion. Heterozygous carriers (HbAS) demonstrate resistance to severe *Plasmodium falciparum* malaria, explaining the high allele frequency in malaria-endemic regions (carrier rate up to 40% in sub-Saharan Africa; ~300,000 affected births annually).

#### OMIM Phenotype — #603903

Autosomal recessive hemoglobinopathy with chronic hemolytic anemia, recurrent vaso-occlusive pain crises, acute chest syndrome, ischemic stroke, splenic sequestration, progressive renal impairment, priapism, leg ulcers, and pulmonary hypertension. Recently approved gene therapies — Casgevy (CRISPR/Cas9-based) and Lyfgenia (lentiviral vector) — offer potential curative options.

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PS3** | Well-established functional studies (70+ years) confirming HbS polymerization and sickling |
| **PS4** | High prevalence in affected cohorts vs. controls |
| **PM1** | Critical surface domain of beta-globin protein |
| **PM2** | Rare in homozygous controls |
| **PP3** | Computational pathogenicity evidence (AlphaMissense 0.874, REVEL 0.924) |
| **PP5** | Classified Pathogenic by multiple reputable sources |

---

### Disease 2: Hereditary Hemochromatosis

| Field | Value |
| --- | --- |
| **Gene** | *HFE* (Homeostatic Iron Regulator) |
| **HGVS Notation** | NM_000410.4(HFE):c.845G>A |
| **Protein Change** | p.Cys282Tyr |
| **Mutation Type** | Single nucleotide variant — missense |
| **Chromosome Position** | chr6:26,093,141 (GRCh38) |
| **dbSNP ID** | rs1800562 |
| **ClinVar ID** | VCV000003036 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |

#### Molecular Mechanism

The C282Y variant disrupts a critical disulfide bond in the HFE protein, preventing it from folding correctly and binding β2-microglobulin. This abolishes HFE's regulation of transferrin receptor 1-mediated iron uptake, resulting in uncontrolled intestinal iron absorption. Found in 83% of hereditary hemochromatosis patients vs. 0.01% of controls (Feder et al., 1996). Carrier frequency is ~10% in Northern European populations. Penetrance is incomplete (~10% for clinically significant disease), with environmental factors influencing expression.

#### OMIM Phenotype — #235200

Autosomal recessive iron overload disorder causing cirrhosis, hepatocellular carcinoma, diabetes mellitus ("bronze diabetes"), cardiomyopathy, arthropathy, hypogonadism, and skin hyperpigmentation. Symptoms typically begin age 40–60 in males and later in females. Treatment with therapeutic phlebotomy is highly effective when initiated before end-organ damage.

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PS3** | Functional studies demonstrate abolished β2-microglobulin binding and iron dysregulation |
| **PS4** | Variant enriched in 83% of hemochromatosis cases vs. 0.01% controls |
| **PM1** | Critical cysteine residue required for disulfide bond formation |
| **PM2** | Low homozygous frequency in general population |
| **PP3** | Computational predictions (AlphaMissense 0.783, REVEL 0.856) |
| **PP5** | Consistent Pathogenic classification across multiple sources |

---

### Disease 3: Phenylketonuria (PKU)

| Field | Value |
| --- | --- |
| **Gene** | *PAH* (Phenylalanine Hydroxylase) |
| **HGVS Notation** | NM_000277.3(PAH):c.1222C>T |
| **Protein Change** | p.Arg408Trp |
| **Mutation Type** | Single nucleotide variant — missense |
| **Chromosome Position** | chr12:102,836,891 (GRCh38) |
| **dbSNP ID** | rs5030858 |
| **ClinVar ID** | VCV000005345 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |

#### Molecular Mechanism

The R408W variant causes protein misfolding and accelerated proteasomal degradation of phenylalanine hydroxylase (PAH), reducing residual enzyme activity to 3–8% of normal. This prevents conversion of phenylalanine to tyrosine, leading to toxic phenylalanine accumulation in blood and brain (levels >1200 μmol/L in classical PKU; normal <120 μmol/L). R408W accounts for ~15–20% of PKU chromosomes in Northern Europe and produces the most severe classical PKU phenotype. Incidence is approximately 1:10,000 births worldwide.

#### OMIM Phenotype — #261600

Autosomal recessive inborn error of metabolism causing severe intellectual disability, seizures, behavioral problems, eczema, hypopigmentation, and a musty body odor if untreated. Brain MRI shows white matter abnormalities. With early newborn screening and lifelong low-phenylalanine diet (and/or sapropterin/pegvaliase therapy), outcomes are excellent.

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PS3** | Well-established functional studies showing protein misfolding and 3–8% residual activity |
| **PS4** | High prevalence in affected alleles (15–20% in Northern Europe) |
| **PM1** | Critical catalytic domain of PAH enzyme |
| **PM2** | Rare in general population databases |
| **PM3** | Detected in trans with other pathogenic *PAH* variants |
| **PP3** | Computational predictions (AlphaMissense 0.912, REVEL 0.947) |
| **PP5** | Literature and database consensus |

---

## Rare Diseases

### Disease 4: Hutchinson-Gilford Progeria Syndrome (HGPS)

| Field | Value |
| --- | --- |
| **Gene** | *LMNA* (Lamin A/C) |
| **HGVS Notation** | NM_170707.4(LMNA):c.1824C>T |
| **Protein Change** | p.Gly608Gly (synonymous — cryptic splice site) |
| **Mutation Type** | Splice region variant |
| **Chromosome Position** | chr1:156,105,681 (GRCh38) |
| **dbSNP ID** | rs121912706 |
| **ClinVar ID** | VCV000041263 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |
| **Incidence** | <1 in 4,000,000 |

#### Molecular Mechanism

Despite being synonymous at the amino acid level, this c.1824C>T variant activates a cryptic splice donor site in exon 11 of *LMNA*, producing a truncated prelamin A protein (progerin) with a 50-amino-acid internal deletion. Progerin retains a farnesyl group normally removed in wild-type processing, causing it to permanently anchor to the nuclear envelope. Accumulation of progerin disrupts nuclear architecture, impairs DNA repair, and accelerates cellular senescence. Found in ~90% of classical progeria cases (Eriksson et al., 2003); de novo in most patients. Median survival is 14.6 years, predominantly due to premature cardiovascular disease.

#### OMIM Phenotype — #176670

Autosomal dominant (de novo) premature aging syndrome: growth retardation from infancy, alopecia, lipodystrophy, scleroderma-like skin changes, prominent scalp veins, craniofacial disproportion, loss of subcutaneous fat, and rapidly progressive atherosclerosis leading to myocardial infarction or stroke. Intellect is preserved.

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PVS1** | Activates cryptic splice site — null/loss-of-function mechanism established |
| **PS3** | Functional studies confirm progerin production and nuclear abnormalities |
| **PS4** | Found in ~90% of all classical HGPS cases |
| **PM2** | Extremely rare in population databases |
| **PP1** | Cosegregation with disease in families |
| **PP5** | Classified Pathogenic by all reputable clinical laboratories |

---

### Disease 5: Fibrodysplasia Ossificans Progressiva (FOP)

| Field | Value |
| --- | --- |
| **Gene** | *ACVR1* (Activin A Receptor Type 1) |
| **HGVS Notation** | NM_001105.6(ACVR1):c.617G>A |
| **Protein Change** | p.Arg206His |
| **Mutation Type** | Single nucleotide variant — missense |
| **Chromosome Position** | chr2:157,779,160 (GRCh38) |
| **dbSNP ID** | rs121913490 |
| **ClinVar ID** | VCV000013642 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |
| **Incidence** | 1 in 2,000,000 |

#### Molecular Mechanism

The R206His variant causes gain-of-function constitutive activation of the ACVR1 (ALK2) BMP type I receptor kinase, located in the glycine-serine (GS) activation domain. This drives inappropriate chondrogenesis and heterotopic osteogenesis in soft tissues — particularly following injury or inflammation. The receptor also becomes aberrantly responsive to Activin A (normally an inhibitor), further amplifying ossification signals. Found in >95% of all classical FOP cases worldwide (Shore et al., 2006). De novo in ~95% of cases; incidence is 1 in 2 million with no ethnic predilection.

#### OMIM Phenotype — #135100

Autosomal dominant (mostly de novo) disorder with progressive heterotopic ossification of skeletal muscle, tendons, ligaments, and fascia. A diagnostic hallmark is congenital malformation of the great toes. Painful soft tissue flare-ups triggered by trauma, infection, or intramuscular injections lead to irreversible bone formation. Most patients lose ambulation by age 30. No approved disease-modifying therapy exists (though Palovarotene received FDA approval in 2023 for flare management).

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PS3** | Functional studies confirm constitutive BMP signaling activation |
| **PS4** | Present in >95% of classical FOP cases |
| **PM1** | GS activation domain — critical for receptor regulation |
| **PM2** | Absent in population control databases |
| **PM5** | Different amino acid changes at same position also pathogenic |
| **PP3** | Computational predictions (AlphaMissense 0.891, REVEL 0.940) |
| **PP5** | Consistent Pathogenic classification across all reputable sources |

---

### Disease 6: Alkaptonuria

| Field | Value |
| --- | --- |
| **Gene** | *HGD* (Homogentisate 1,2-Dioxygenase) |
| **HGVS Notation** | NM_000187.4(HGD):c.342+1G>A |
| **Protein Change** | Splice donor disruption — exon 6 |
| **Mutation Type** | Splice donor variant |
| **Chromosome Position** | chr3:120,312,459 (GRCh38) |
| **dbSNP ID** | rs121907954 |
| **ClinVar ID** | VCV000036396 |
| **Review Status** | Reviewed by expert panel |
| **Clinical Significance** | Pathogenic |
| **Prevalence** | 1 in 250,000–1,000,000 |

#### Molecular Mechanism

The c.342+1G>A variant disrupts the invariant GT splice donor dinucleotide at the intron 6 boundary of *HGD*, completely abolishing normal splicing of exon 6. This results in intron retention or exon skipping, producing a non-functional truncated homogentisate 1,2-dioxygenase enzyme. Loss of HGD activity prevents catabolism of homogentisic acid (HGA) in the tyrosine degradation pathway, leading to HGA accumulation and excretion (4–8 g/day in urine). This variant accounts for ~12% of mutant alleles in European alkaptonuria populations (Zatkova et al., 2000).

#### OMIM Phenotype — #203500

Autosomal recessive inborn error of tyrosine metabolism characterized by homogentisic aciduria (urine darkens on standing — first sign in infancy), ochronosis (bluish-black pigmentation of connective tissue, visible in sclerae and ear cartilage by age 30), and ochronotic arthropathy (severe degenerative arthritis of the spine and large joints by age 40–50). Cardiac valve involvement and renal/prostate calculi may also occur. Nitisinone (NTBC) reduces HGA excretion and is used off-label.

#### ACMG/AMP Classification — Pathogenic

| Criterion | Evidence |
| --- | --- |
| **PVS1** | Canonical splice site ±1/2 variant — loss of function established mechanism |
| **PS3** | Complete loss of HGD enzyme activity demonstrated in functional studies |
| **PM2** | Absent in population control databases |
| **PM3** | Detected in trans with other pathogenic *HGD* variants in affected individuals |
| **PP4** | Phenotype highly specific to alkaptonuria (dark urine, ochronosis) |
| **PP5** | Classified Pathogenic by reputable clinical sources |

---

## VCF File and Annotation

### VCF File (variants.vcf)

```
##fileformat=VCFv4.2
##fileDate=20260301
##source=Assignment2_GeneticVariantAnalysis
##reference=GRCh38/hg38
##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC gene symbol">
##INFO=<ID=DISEASE,Number=1,Type=String,Description="Associated disease">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance">
##INFO=<ID=CLNID,Number=1,Type=String,Description="ClinVar accession">
##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Variant consequence type">
##INFO=<ID=PROTEIN,Number=1,Type=String,Description="Protein change (HGVS)">
##INFO=<ID=OMIM,Number=1,Type=String,Description="OMIM phenotype accession">
##FILTER=<ID=PASS,Description="All filters passed">
#CHROM  POS        ID           REF  ALT  QUAL  FILTER  INFO
11      5227002    rs334        A    T    .     PASS    GENE=HBB;DISEASE=Sickle_Cell_Disease;CLNSIG=Pathogenic;CLNID=VCV000015280;CONSEQUENCE=missense_variant;PROTEIN=p.Glu6Val;OMIM=603903
6       26093141   rs1800562    G    A    .     PASS    GENE=HFE;DISEASE=Hereditary_Hemochromatosis;CLNSIG=Pathogenic;CLNID=VCV000003036;CONSEQUENCE=missense_variant;PROTEIN=p.Cys282Tyr;OMIM=235200
12      102836891  rs5030858    C    T    .     PASS    GENE=PAH;DISEASE=Phenylketonuria;CLNSIG=Pathogenic;CLNID=VCV000005345;CONSEQUENCE=missense_variant;PROTEIN=p.Arg408Trp;OMIM=261600
1       156105681  rs121912706  C    T    .     PASS    GENE=LMNA;DISEASE=Hutchinson-Gilford_Progeria_Syndrome;CLNSIG=Pathogenic;CLNID=VCV000041263;CONSEQUENCE=splice_region_variant;PROTEIN=p.Gly608Gly;OMIM=176670
2       157779160  rs121913490  G    A    .     PASS    GENE=ACVR1;DISEASE=Fibrodysplasia_Ossificans_Progressiva;CLNSIG=Pathogenic;CLNID=VCV000013642;CONSEQUENCE=missense_variant;PROTEIN=p.Arg206His;OMIM=135100
3       120312459  rs121907954  G    A    .     PASS    GENE=HGD;DISEASE=Alkaptonuria;CLNSIG=Pathogenic;CLNID=VCV000036396;CONSEQUENCE=splice_donor_variant;PROTEIN=c.342+1G>A;OMIM=203500
```

### Ensembl VEP Annotation

Submit `variants.vcf` to Ensembl VEP (GRCh38) at https://asia.ensembl.org/Tools/VEP using the following settings:
- Species: *Homo sapiens*
- Assembly: GRCh38/hg38
- Transcript database: Ensembl/GENCODE
- Enable: SIFT, PolyPhen, regulatory features, population frequencies (gnomAD)

---

## Repository Structure

```
Assignment2-Genetic-Variants/
├── Assignment2.xlsx     # Main data sheet: variant details, OMIM phenotypes,
│                        # UCSC AlphaMissense & REVEL screenshots, ACMG criteria
├── variants.vcf         # VCFv4.2 file containing all 6 pathogenic variants
└── README.md            # This file
```

---

## Tools and Databases

| Tool / Database | Version / Build | Purpose |
| --- | --- | --- |
| NCBI ClinVar | Accessed March 2026 | Variant identification, clinical significance, review status |
| OMIM | Accessed March 2026 | Phenotype, inheritance, clinical features |
| UCSC Genome Browser | GRCh38/hg38 | AlphaMissense and REVEL pathogenicity track visualization |
| AlphaMissense | Google DeepMind (2023) | AI-based proteome-wide missense pathogenicity prediction |
| REVEL | v1.3 | Ensemble-based missense variant pathogenicity scoring |
| SpliceAI | v1.3 | Deep learning splice variant pathogenicity prediction |
| MaxEntScan | — | Splice site strength scoring |
| Ensembl VEP | GRCh38.p14 | VCF annotation, consequence and transcript overlap prediction |
| ACMG/AMP Guidelines | Richards et al. 2015 | Variant pathogenicity classification framework |

---

## Steps to Reproduce

1. Clone this repository:
   ```bash
   git clone https://github.com/FaiqaZarar/Assignment2-Genetic-Variants.git
   cd Assignment2-Genetic-Variants
   ```

2. Open `Assignment2.xlsx` to view all variant details, phenotype data, pathogenicity scores, and ACMG classifications.

3. To annotate the VCF file using Ensembl VEP:
   - Go to https://asia.ensembl.org/Tools/VEP
   - Upload `variants.vcf`
   - Select GRCh38/hg38 assembly
   - Run annotation and download results

4. For UCSC Genome Browser screenshots:
   - Go to https://genome.ucsc.edu
   - Select GRCh38/hg38 assembly
   - Navigate to each chromosomal position listed in the VCF
   - Enable AlphaMissense and REVEL tracks under "Phenotype and Literature"
   - Capture screenshots and insert into `Assignment2.xlsx`

---

## References

1. Richards S, et al. (2015). Standards and guidelines for the interpretation of sequence variants. *Genetics in Medicine*, 17(5), 405–423. https://doi.org/10.1038/gim.2015.30
2. Pauling L, et al. (1949). Sickle cell anemia, a molecular disease. *Science*, 110(2865), 543–548. https://doi.org/10.1126/science.110.2865.543
3. Ingram VM (1956). A specific chemical difference between the globins of normal human and sickle-cell anaemia haemoglobin. *Nature*, 178(4537), 792–794. https://doi.org/10.1038/178792a0
4. Piel FB, et al. (2017). Global epidemiology of sickle haemoglobin in neonates. *PLOS Medicine*. https://doi.org/10.1371/journal.pmed.1001484
5. Feder JN, et al. (1996). A novel MHC class I–like gene is mutated in patients with hereditary haemochromatosis. *Nature Genetics*, 13(4), 399–408. https://doi.org/10.1038/ng0896-399
6. Allen KJ, et al. (2008). Iron-overload–related disease in HFE hereditary hemochromatosis. *NEJM*, 358(3), 221–230. https://doi.org/10.1056/NEJMoa073286
7. DiLella AG, et al. (1986). Molecular structure and polymorphic map of the human phenylalanine hydroxylase gene. *Biochemistry*, 25(4), 743–749. https://doi.org/10.1021/bi00352a001
8. Guldberg P, et al. (1998). A European multicenter study of phenylalanine hydroxylase deficiency. *European Journal of Pediatrics*, 157, S27–S33. https://doi.org/10.1007/PL00014320
9. Eriksson M, et al. (2003). Recurrent de novo point mutations in lamin A cause Hutchinson-Gilford progeria syndrome. *Nature*, 423(6937), 293–298. https://doi.org/10.1038/nature01629
10. De Sandre-Giovannoli A, et al. (2003). Lamin A truncation in Hutchinson-Gilford progeria. *Science*, 300(5628), 2055. https://doi.org/10.1126/science.1084125
11. Shore EM, et al. (2006). A recurrent mutation in the BMP type I receptor ACVR1 causes inherited and sporadic fibrodysplasia ossificans progressiva. *Nature Genetics*, 38(5), 525–527. https://doi.org/10.1038/ng1783
12. Kaplan FS, et al. (2009). Classic and atypical fibrodysplasia ossificans progressiva (FOP) phenotypes are caused by mutations in the bone morphogenetic protein (BMP) type I receptor ACVR1. *Human Mutation*, 30(3), 379–390. https://doi.org/10.1002/humu.20868
13. Fernández-Cañón JM, et al. (1996). The molecular basis of alkaptonuria. *Nature Genetics*, 14(1), 19–24. https://doi.org/10.1038/ng0996-19
14. Zatkova A, et al. (2000). Identification of 11 novel alkaptonuria mutations in Central European patients. *Human Mutation*, 15(5), 394–402. https://doi.org/10.1002/(SICI)1098-1004(200005)15:5<394::AID-HUMU2>3.0.CO;2-9
15. Cheng J, et al. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. *Science*, 381(6660), eadg7492. https://doi.org/10.1126/science.adg7492

---

*National University of Sciences and Technology — School of Interdisciplinary Engineering & Sciences — 2026*
