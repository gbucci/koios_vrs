# KOIOS-VRS Pipeline - Diagramma di Flusso

## Panoramica del Flusso Principale

```mermaid
flowchart TD
    Start([Inizio Pipeline]) --> Input[Input Files]

    Input --> VCF[File VCF<br/>hg19 o hg38]
    Input --> Clinical[Dati Clinici CSV<br/>diagnosi, sesso, età]

    VCF --> PlatformDetect{Rilevamento<br/>Piattaforma}

    PlatformDetect -->|Ion Torrent| IonParams[Parametri Ion Torrent<br/>AF≥1%, DP≥100, FAO≥10]
    PlatformDetect -->|Illumina| IllParams[Parametri Illumina<br/>AF≥5%, DP≥50]
    PlatformDetect -->|Altri| DefaultParams[Parametri Default]

    IonParams --> Preprocess
    IllParams --> Preprocess
    DefaultParams --> Preprocess

    Preprocess[Preprocessing VCF] --> LiftCheck{Genoma<br/>hg19?}

    LiftCheck -->|Sì| Liftover[Liftover a hg38<br/>AnnotationHub]
    LiftCheck -->|No| QualFilter
    Liftover --> QualFilter[Filtri di Qualità]

    QualFilter --> FilterSteps[• Rimuovi genotipi indefiniti<br/>• Filtra per AF<br/>• Filtra per DP<br/>• Filtra per FAO]

    FilterSteps --> MultiAllelic{Varianti<br/>Multi-alleliche?}

    MultiAllelic -->|Sì| SelectBest[Seleziona allele<br/>con FAO/AD massimo]
    MultiAllelic -->|No| SeparateTypes
    SelectBest --> SeparateTypes

    SeparateTypes[Separazione Tipi] --> SNP[SNP/INDEL<br/>→ Pipeline principale]
    SeparateTypes --> CNV[CNV<br/>→ File separato]

    SNP --> Annotation[Annotazione KOIOS]

    Annotation --> AnnotSteps[• Carica concetti OMOP<br/>• Genera notazione HGVSg<br/>• Classifica ClinGen/ClinVar<br/>• Estrai HGVS proteico]

    AnnotSteps --> PopFilter[Filtro Popolazioni<br/>gnomAD AF > 1%]

    PopFilter --> VRS[Integrazione VRS]

    VRS --> VRSSteps[• Query VRS Registry API<br/>• Recupera VRS Allele IDs<br/>• Logica di retry<br/>• Aggiungi Sequence IDs]

    Clinical --> Phenopacket
    VRSSteps --> Phenopacket[Generazione Phenopacket]

    Phenopacket --> Outputs[Generazione Output]

    Outputs --> JSON[Phenopackets JSON<br/>GA4GH v2]
    Outputs --> VUS[VUS CSV<br/>Varianti Unknown]
    Outputs --> Note[Note CSV<br/>Varianti Significative]
    Outputs --> GeneSummary[Gene Summary CSV<br/>Conteggi per gene]
    Outputs --> AnnotVCF[VCF Annotato<br/>VRS + OMOP]

    JSON --> End([Fine Pipeline])
    VUS --> End
    Note --> End
    GeneSummary --> End
    AnnotVCF --> End
    CNV --> End
```

## Dettaglio Rilevamento Piattaforma

```mermaid
flowchart TD
    Start[VCF File] --> ReadHeader[Leggi Header VCF<br/>prime 200 righe]

    ReadHeader --> Parse[Parse Metadati]

    Parse --> CheckTVC{Contiene<br/>TVC?}
    Parse --> CheckIon{Contiene<br/>Ion Torrent?}
    Parse --> CheckOncomine{Contiene<br/>Oncomine?}

    CheckTVC -->|Sì| IonTorrent[Piattaforma:<br/>Ion Torrent]
    CheckIon -->|Sì| IonTorrent
    CheckOncomine -->|Sì| IonTorrent

    Parse --> CheckGATK{Contiene<br/>GATK?}
    Parse --> CheckStrelka{Contiene<br/>Strelka?}
    Parse --> CheckMuTect{Contiene<br/>MuTect?}

    CheckGATK -->|Sì| Illumina[Piattaforma:<br/>Illumina]
    CheckStrelka -->|Sì| Illumina
    CheckMuTect -->|Sì| Illumina

    CheckTVC -->|No| Unknown
    CheckIon -->|No| Unknown
    CheckOncomine -->|No| Unknown
    CheckGATK -->|No| Unknown
    CheckStrelka -->|No| Unknown
    CheckMuTect -->|No| Unknown[Piattaforma:<br/>Sconosciuta]

    IonTorrent --> SetParams1[Imposta Parametri<br/>AF≥1%, DP≥100<br/>FAO≥10]
    Illumina --> SetParams2[Imposta Parametri<br/>AF≥5%, DP≥50<br/>FAO=NA]
    Unknown --> SetParams3[Parametri Default<br/>AF≥5%, DP≥50]

    SetParams1 --> Return[Ritorna Parametri]
    SetParams2 --> Return
    SetParams3 --> Return
```

## Dettaglio Query VRS API

```mermaid
flowchart TD
    Start[Variante Annotata] --> PrepareQuery[Prepara Query<br/>HGVSg notation]

    PrepareQuery --> Query[HTTP GET Request<br/>reg.genome.network/allele]

    Query --> Response{Status<br/>HTTP}

    Response -->|200 OK| Parse[Parse JSON Response]
    Response -->|404| Retry{Retry<br/>Count < 3?}
    Response -->|500| Retry
    Response -->|Timeout| Retry

    Retry -->|Sì| Wait[Attendi<br/>backoff esponenziale]
    Wait --> Query

    Retry -->|No| Failed[Segna come<br/>VRS lookup failed]

    Parse --> Extract[Estrai Campi VRS]

    Extract --> Fields[• VRS Allele ID<br/>• VRS Sequence ID<br/>• Location info]

    Fields --> AddToVariant[Aggiungi a Variante]
    Failed --> AddToVariant

    AddToVariant --> CheckMore{Altre<br/>varianti?}

    CheckMore -->|Sì| PrepareQuery
    CheckMore -->|No| Complete[VRS Integration<br/>Completa]
```

## Dettaglio Generazione Phenopacket

```mermaid
flowchart TD
    Start[Dati Annotati] --> CreatePacket[Crea Phenopacket Object]

    CreatePacket --> AddID[Aggiungi UUID<br/>identificatore unico]

    AddID --> AddSubject[Aggiungi Subject]
    AddSubject --> SubjectFields[• sample_id<br/>• sex PATO ID<br/>• age of onset ISO]

    SubjectFields --> AddDisease[Aggiungi Disease]
    AddDisease --> DiseaseFields[• diagnosis label<br/>• NCIT term ID<br/>• onset age]

    DiseaseFields --> AddInterpretations[Aggiungi Interpretations]

    AddInterpretations --> LoopVariants{Per ogni<br/>variante}

    LoopVariants --> CreateInterp[Crea Interpretation]

    CreateInterp --> InterpFields[• Variant HGVS<br/>• VRS Allele ID<br/>• Gene symbol<br/>• ClinVar class<br/>• OMOP concepts]

    InterpFields --> AddToList[Aggiungi a Lista]

    AddToList --> LoopVariants

    LoopVariants -->|Completato| AddMetadata[Aggiungi Metadata]

    AddMetadata --> MetadataFields[• Created timestamp<br/>• Software version<br/>• Resources used]

    MetadataFields --> Serialize[Serializza a JSON]

    Serialize --> ValidateSchema{Schema<br/>Phenopackets v2<br/>valido?}

    ValidateSchema -->|Sì| WriteFile[Scrivi File JSON]
    ValidateSchema -->|No| Error[Errore Validazione]

    WriteFile --> Complete[Phenopacket Generato]
```

## Flusso Dati: Input → Output

```mermaid
flowchart LR
    subgraph Input["INPUT FILES"]
        VCF[VCF File<br/>Varianti genomiche]
        CSV[Clinical CSV<br/>Metadati paziente]
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
        NOTE[Note CSV<br/>Significant Variants]
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
    Stage5 --> NOTE
    Stage5 --> GENE
    Stage5 --> AVCF
    Stage2 --> CNV
```

## Architettura Modulare

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

    subgraph External["Dipendenze Esterne"]
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

## Standard e Ontologie Utilizzate

```mermaid
flowchart LR
    subgraph Standards["Standard Genomici"]
        VRS[VRS<br/>Variant Representation<br/>Specification]
        HGVS[HGVS<br/>Nomenclatura<br/>Varianti]
        VCF[VCF<br/>Variant Call<br/>Format]
        Pheno[Phenopackets v2<br/>GA4GH]
    end

    subgraph Ontologies["Ontologie"]
        OMOP[OMOP<br/>Medical Concepts]
        NCIT[NCI Thesaurus<br/>Diagnosi Cancer]
        PATO[PATO<br/>Fenotipo/Sesso]
        ClinVar[ClinVar/ClinGen<br/>Significato Clinico]
    end

    subgraph Databases["Database Popolazioni"]
        gnomAD[gnomAD<br/>Frequenze Alleliche<br/>Popolazioni]
    end

    Pipeline[KOIOS-VRS<br/>Pipeline] --> Standards
    Pipeline --> Ontologies
    Pipeline --> Databases
```

---

## Note Tecniche

### Tecnologie Principali
- **Linguaggio**: R (≥ 4.0)
- **Tipo**: Pacchetto Bioconductor
- **API**: VRS Registry (reg.genome.network)
- **Standard**: GA4GH Phenopackets v2, VRS, HGVS

### Piattaforme Supportate
- Ion Torrent (Oncomine panels, TVC caller)
- Illumina (GATK, Strelka, MuTect2, VarScan)
- Altri sequenziatori (parametri default)

### Formati Genoma
- **Input**: hg19 o hg38
- **Output**: hg38 (liftover automatico se necessario)

### Casi d'Uso
- Genomica del cancro (melanoma, polmone, colon-retto)
- Reporting clinico varianti
- Integrazione OMOP CDM
- Interoperabilità dati GA4GH

---

*Diagramma generato per KOIOS-VRS Pipeline v1.0*
*DOI: 10.5281/zenodo.17476991*
