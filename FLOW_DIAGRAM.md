# KOIOS-VRS Pipeline - Flow Diagram

## Main Pipeline Flow

```mermaid
flowchart TD
    Start([Start Pipeline]) --> Input[Input Files]

    Input --> VCF[VCF File<br/>hg19 or hg38]
    Input --> Clinical[Clinical CSV<br/>diagnosis, sex, age]

    VCF --> PlatformDetect{Platform<br/>Detection}

    PlatformDetect -->|Ion Torrent| IonParams[Ion Torrent Parameters<br/>AF≥1%, DP≥100, FAO≥10]
    PlatformDetect -->|Illumina| IllParams[Illumina Parameters<br/>AF≥5%, DP≥50]
    PlatformDetect -->|Others| DefaultParams[Default Parameters]

    IonParams --> Preprocess
    IllParams --> Preprocess
    DefaultParams --> Preprocess

    Preprocess[VCF Preprocessing] --> LiftCheck{Genome<br/>hg19?}

    LiftCheck -->|Yes| Liftover[Liftover to hg38<br/>AnnotationHub]
    LiftCheck -->|No| QualFilter
    Liftover --> QualFilter[Quality Filters]

    QualFilter --> FilterSteps[• Remove undefined genotypes<br/>• Filter by AF<br/>• Filter by DP<br/>• Filter by FAO]

    FilterSteps --> MultiAllelic{Multi-allelic<br/>Variants?}

    MultiAllelic -->|Yes| SelectBest[Select allele with<br/>maximum FAO/AD]
    MultiAllelic -->|No| SeparateTypes
    SelectBest --> SeparateTypes

    SeparateTypes[Type Separation] --> SNP[SNP/INDEL<br/>→ Main pipeline]
    SeparateTypes --> CNV[CNV<br/>→ Separate file]

    SNP --> Annotation[KOIOS Annotation]

    Annotation --> AnnotSteps[• Load OMOP concepts<br/>• Generate HGVSg notation<br/>• ClinGen/ClinVar classification<br/>• Extract protein HGVS]

    AnnotSteps --> PopFilter[Population Filter<br/>gnomAD AF > 1%]

    PopFilter --> VRS[VRS Integration]

    VRS --> VRSSteps[• Query VRS Registry API<br/>• Fetch VRS Allele IDs<br/>• Retry logic<br/>• Add Sequence IDs]

    Clinical --> Phenopacket
    VRSSteps --> Phenopacket[Phenopacket Generation]

    Phenopacket --> Outputs[Output Generation]

    Outputs --> JSON[Phenopackets JSON<br/>GA4GH v2]
    Outputs --> VUS[VUS CSV<br/>Unknown Variants]
    Outputs --> Known[Known CSV<br/>Significant Variants]
    Outputs --> GeneSummary[Gene Summary CSV<br/>Counts per gene]
    Outputs --> AnnotVCF[Annotated VCF<br/>VRS + OMOP]

    JSON --> End([End Pipeline])
    VUS --> End
    Known --> End
    GeneSummary --> End
    AnnotVCF --> End
    CNV --> End
```

## Platform Detection Detail

```mermaid
flowchart TD
    Start[VCF File] --> ReadHeader[Read VCF Header<br/>first 200 lines]

    ReadHeader --> Parse[Parse Metadata]

    Parse --> CheckTVC{Contains<br/>TVC?}
    Parse --> CheckIon{Contains<br/>Ion Torrent?}
    Parse --> CheckOncomine{Contains<br/>Oncomine?}

    CheckTVC -->|Yes| IonTorrent[Platform:<br/>Ion Torrent]
    CheckIon -->|Yes| IonTorrent
    CheckOncomine -->|Yes| IonTorrent

    Parse --> CheckGATK{Contains<br/>GATK?}
    Parse --> CheckStrelka{Contains<br/>Strelka?}
    Parse --> CheckMuTect{Contains<br/>MuTect?}

    CheckGATK -->|Yes| Illumina[Platform:<br/>Illumina]
    CheckStrelka -->|Yes| Illumina
    CheckMuTect -->|Yes| Illumina

    CheckTVC -->|No| Unknown
    CheckIon -->|No| Unknown
    CheckOncomine -->|No| Unknown
    CheckGATK -->|No| Unknown
    CheckStrelka -->|No| Unknown
    CheckMuTect -->|No| Unknown[Platform:<br/>Unknown]

    IonTorrent --> SetParams1[Set Parameters<br/>AF≥1%, DP≥100<br/>FAO≥10]
    Illumina --> SetParams2[Set Parameters<br/>AF≥5%, DP≥50<br/>FAO=NA]
    Unknown --> SetParams3[Default Parameters<br/>AF≥5%, DP≥50]

    SetParams1 --> Return[Return Parameters]
    SetParams2 --> Return
    SetParams3 --> Return
```

## VRS API Query Detail

```mermaid
flowchart TD
    Start[Annotated Variant] --> PrepareQuery[Prepare Query<br/>HGVSg notation]

    PrepareQuery --> Query[HTTP GET Request<br/>reg.genome.network/allele]

    Query --> Response{HTTP<br/>Status}

    Response -->|200 OK| Parse[Parse JSON Response]
    Response -->|404| Retry{Retry<br/>Count < 3?}
    Response -->|500| Retry
    Response -->|Timeout| Retry

    Retry -->|Yes| Wait[Wait<br/>exponential backoff]
    Wait --> Query

    Retry -->|No| Failed[Mark as<br/>VRS lookup failed]

    Parse --> Extract[Extract VRS Fields]

    Extract --> Fields[• VRS Allele ID<br/>• VRS Sequence ID<br/>• Location info]

    Fields --> AddToVariant[Add to Variant]
    Failed --> AddToVariant

    AddToVariant --> CheckMore{More<br/>variants?}

    CheckMore -->|Yes| PrepareQuery
    CheckMore -->|No| Complete[VRS Integration<br/>Complete]
```

## Phenopacket Generation Detail

```mermaid
flowchart TD
    Start[Annotated Data] --> CreatePacket[Create Phenopacket Object]

    CreatePacket --> AddID[Add UUID<br/>unique identifier]

    AddID --> AddSubject[Add Subject]
    AddSubject --> SubjectFields[• sample_id<br/>• sex PATO ID<br/>• age of onset ISO]

    SubjectFields --> AddDisease[Add Disease]
    AddDisease --> DiseaseFields[• diagnosis label<br/>• NCIT term ID<br/>• onset age]

    DiseaseFields --> AddInterpretations[Add Interpretations]

    AddInterpretations --> LoopVariants{For each<br/>variant}

    LoopVariants --> CreateInterp[Create Interpretation]

    CreateInterp --> InterpFields[• Variant HGVS<br/>• VRS Allele ID<br/>• Gene symbol<br/>• ClinVar class<br/>• OMOP concepts]

    InterpFields --> AddToList[Add to List]

    AddToList --> LoopVariants

    LoopVariants -->|Complete| AddMetadata[Add Metadata]

    AddMetadata --> MetadataFields[• Created timestamp<br/>• Software version<br/>• Resources used]

    MetadataFields --> Serialize[Serialize to JSON]

    Serialize --> ValidateSchema{Phenopackets v2<br/>Schema<br/>valid?}

    ValidateSchema -->|Yes| WriteFile[Write JSON File]
    ValidateSchema -->|No| Error[Validation Error]

    WriteFile --> Complete[Phenopacket Generated]
```

## Data Flow: Input → Output

```mermaid
flowchart LR
    subgraph Input["INPUT FILES"]
        VCF[VCF File<br/>Genomic variants]
        CSV[Clinical CSV<br/>Patient metadata]
    end

    subgraph Processing["PROCESSING STAGES"]
        Stage1[1. Platform Detection<br/>& Parameter Setup]
        Stage2[2. Preprocessing<br/>& Quality Filtering]
        Stage3[3. KOIOS Annotation<br/>OMOP/ClinVar]
        Stage4[4. VRS Integration<br/>API Queries]
        Stage5[5. Phenopacket<br/>Generation]
    end

    subgraph Output["OUTPUT FILES"]
        JSON[Phenopackets JSON<br/>GA4GH Standard]
        VUS[VUS CSV<br/>Unknown Variants]
        KNOWN[Known CSV<br/>Significant Variants]
        GENE[GeneSummary CSV<br/>Gene-level Stats]
        AVCF[Annotated VCF<br/>Enhanced VCF]
        CNV[CNV VCF<br/>Copy Number]
    end

    VCF --> Stage1
    CSV --> Stage5
    Stage1 --> Stage2
    Stage2 --> Stage3
    Stage3 --> Stage4
    Stage4 --> Stage5

    Stage5 --> JSON
    Stage5 --> VUS
    Stage5 --> KNOWN
    Stage5 --> GENE
    Stage5 --> AVCF
    Stage2 --> CNV
```

## Modular Architecture

```mermaid
flowchart TD
    subgraph Main["R/koios_pipeline.R"]
        Pipeline[koios_vrs_pipeline]
        Preprocess[preprocess_vcf]
        Extract[extract_protein_hgvs]
        Fetch[fetch_vrs_fields]
        Generate[generate_phenopacket]
    end

    subgraph Platform["R/platform_detection.R"]
        Detect[detect_sequencing_platform]
        GetParams[get_platform_parameters]
        PrintSummary[print_platform_summary]
    end

    subgraph Helper["create_clinical_csv.R"]
        CreateCSV[create_clinical_csv]
        Interactive[interactive_mode]
        CommandLine[command_line_mode]
    end

    subgraph External["External Dependencies"]
        VcfR[vcfR<br/>VCF parsing]
        KOIOS[KOIOS<br/>OMOP annotation]
        RTL[rtracklayer<br/>Liftover]
        HTTR[httr<br/>VRS API]
        JSON[jsonlite<br/>JSON I/O]
    end

    Pipeline --> Preprocess
    Pipeline --> Detect
    Pipeline --> Extract
    Pipeline --> Fetch
    Pipeline --> Generate

    Detect --> GetParams
    Detect --> PrintSummary

    CreateCSV --> Interactive
    CreateCSV --> CommandLine

    Pipeline --> VcfR
    Pipeline --> KOIOS
    Preprocess --> RTL
    Fetch --> HTTR
    Generate --> JSON
```

## Standards and Ontologies Used

```mermaid
flowchart LR
    subgraph Standards["Genomic Standards"]
        VRS[VRS<br/>Variant Representation<br/>Specification]
        HGVS[HGVS<br/>Variant<br/>Nomenclature]
        VCF[VCF<br/>Variant Call<br/>Format]
        Pheno[Phenopackets v2<br/>GA4GH]
    end

    subgraph Ontologies["Ontologies"]
        OMOP[OMOP<br/>Medical Concepts]
        NCIT[NCI Thesaurus<br/>Cancer Diagnosis]
        PATO[PATO<br/>Phenotype/Sex]
        ClinVar[ClinVar/ClinGen<br/>Clinical Significance]
    end

    subgraph Databases["Population Databases"]
        gnomAD[gnomAD<br/>Allele Frequencies<br/>Populations]
    end

    Pipeline[KOIOS-VRS<br/>Pipeline] --> Standards
    Pipeline --> Ontologies
    Pipeline --> Databases
```

---

## Technical Notes

### Main Technologies
- **Language**: R (≥ 4.0)
- **Type**: Bioconductor package
- **API**: VRS Registry (reg.genome.network)
- **Standards**: GA4GH Phenopackets v2, VRS, HGVS

### Supported Platforms
- Ion Torrent (Oncomine panels, TVC caller)
- Illumina (GATK, Strelka, MuTect2, VarScan)
- Other sequencers (default parameters)

### Genome Formats
- **Input**: hg19 or hg38
- **Output**: hg38 (automatic liftover if needed)

### Use Cases
- Cancer genomics (melanoma, lung, colorectal)
- Clinical variant reporting
- OMOP CDM integration
- GA4GH data interoperability

### Output Files

For each pipeline run, the following files are generated:
1. **`*_Phenopackets.json`** - GA4GH Phenopackets v2 format with clinical and genomic data
2. **`*_VUS.csv`** - Variants of Unknown Significance (concept_id = 1028197)
3. **`*_Known.csv`** - Clinically significant variants (known pathogenic/benign)
4. **`*_GeneSummary.csv`** - Gene-level mutation counts and summary statistics
5. **`*_annotated.vcf`** - Final VCF with VRS and OMOP annotations in INFO field
6. **`*_preprocessed.vcf`** - Quality-filtered intermediate VCF (after QC, before annotation)
7. **`*_CNV.vcf`** - Copy number variants (separated from main pipeline)
8. **`*_hg38.vcf`** - Lifted-over VCF (generated only if input was hg19)

---

*Flow diagram for KOIOS-VRS Pipeline v1.0*
*DOI: 10.5281/zenodo.17476991*
