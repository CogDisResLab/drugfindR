---
title: "`drugfindR`: Transcriptomic signature analysis and drug repurposing using iLINCS in R"
tags:
  - LINCS
  - Bioinformatics
  - Drug Repurposing
authors:
  - name: Ali Sajid Imami
    affiliation: 1
    orcid: 0000-0003-3684-3539
    email: Ali.Sajid.Imami@gmail.com
    equal-contrib: true
    corresponding: yes
  - name: Smita Sahay
    affiliation: 1
    orcid: 0009-0003-4377-8963
    email: Smita.Sahay@rockets.utoledo.edu
    equal-contrib: true
  - name: Justin Fortune Creeden
    affiliation: 3
    orcid: 0000-0003-3123-8401
    email: Justin.Creeden@rockets.utoledo.edu
  - name: Robert Erne McCullumsmith
    affiliation: "1, 2"
    orcid: 0000-0003-3123-8401
    email: Robert.McCullumsmith@utoledo.edu

affiliations:
  - name: "Department of Neurosciences and Psychiatry, University of Toledo, Toledo, OH, USA"
    index: 1
    ror: 01pbdzh19
  - name: "Promedica Neurosciences Institute, Toledo, OH, USA"
    index: 2
    ror: 01k800776
  - name: "Department of Psychiatry, University of California San Diego, La Jolla, Diego, CA, USA"
    index: 3
    ror: 0168r3w48

date: April 15, 2025
bibliography: paper.bib
---

# Summary

`drugfindR` is an R-based computational framework designed to abstract and automate the process of retrieving and comparing transcriptomic signatures using the Library of Integrated Network-Based Cellular Signatures (LINCS) database [@RN479] via the iLINCS platform [@RN6411]. LINCS is a National Institutes of Health (NIH) backed project with the goal of creating a comprehensive database of cell responses to different perturbations. While the biological importance of the L1000 dataset is well-established, programmatic access that maintains statistical rigor and batch-processing capabilities has historically been a significant barrier for researchers.

`drugfindR` provides a structured environment for comparing Differential Gene Expression (DGE) signatures against thousands of chemical and biological perturbations available in the iLINCS platform. By formalizing the workflow of raw signature preparation to consensus ranking, the package enables systematic investigation of candidate drugs for repurposing and functional genomics discovery while maintaining reproducible environments.

# Statement of Need

Developing a new drug from scratch is extraordinarily time-consuming and expensive, often requiring over a decade and billions of dollars, with a high risk of failure due to toxicity or lack of efficacy [@RN6409; @RN6408]. Drug repurposing, finding new therapeutic uses for existing, approved drugs—offers a powerful and cost-effective alternative. By leveraging known safety profiles, this strategy can dramatically accelerate the path to new treatments, as demonstrated by successful cases like allopurinol (gout) and gemcitabine (cancer) [@RN6410; @RN6412].

A leading computational method for repurposing is transcriptomic signature reversion (TSR). This approach posits that a drug capable of reversing the gene expression signature of a disease may have therapeutic potential. The Library of Integrated Network-Based Cellular Signatures (LINCS) serves as a major repository for this methodology, providing a vast public database of gene expression profiles generated in response to thousands of chemical and genetic perturbations.

For researchers, the primary gateway to this resource is the iLINCS web portal. While suitable for analyzing individual signatures, iLINCS has critical inadequacies for rigorous research: its manual, web-based interface makes batch analysis and large-scale studies impractical, and it inherently lacks reproducibility, creating a significant barrier to the systematic use of LINCS data.

`drugfindR` is designed to fill this gap. It provides a programmable abstraction layer, treating transcriptomic signatures as first-class objects within R and standardizing access to the iLINCS REST API. This enables researchers to build complex, automated, and fully reproducible discovery pipelines, transforming the LINCS resource from a manually queried database into a core component of scalable computational biology workflows.

# State of the Field

While other tools like `signatureSearch` [@RN9956] offer broad connectivity mapping, they often require local hosting of massive datasets or complex pre-processing. `drugfindR` uniquely bridges the gap by providing a high-level R interface to the iLINCS API, allowing for real-time cloud-based queries without the "black box" nature of web GUIs. Unlike the raw API examples provided by the iLINCS consortium, `drugfindR` enforces data integrity via S4 classes, ensuring that biological signatures are statistically comparable across different libraries (KD, OE, CP).

# Software Design and Architecture

The architectural philosophy of `drugfindR` is built upon a dual mandate: to provide **extensible tools** for computational biologists constructing complex pipelines, while ensuring an **accessible interface** for domain scientists seeking actionable biological insights. This balance was achieved by adhering to established R package development best practices and Bioconductor submission guidelines [@RN9955], ensuring robustness, maintainability, and interoperability within the Bioconductor ecosystem.

## Foundational Design Principles

The package is constructed around a functional, pipe-friendly (`|>`) interface. This design allows outputs from one core function to seamlessly serve as inputs to the next, enabling clear, readable, and modular workflow construction. The five core modular functions form the essential units of this pipeline, each adhering to the **single responsibility principle**. Each function encapsulates a discrete, critical step in the signature analysis process, from data ingestion to result synthesis. This high degree of composability allows advanced users to customize, intercept, or extend any stage of the analysis.

A key architectural decision was the strict **separation of concerns** between functions that perform local data manipulation and those that call external resources. Specifically, functions that interact with the iLINCS REST API (e.g., `getSignature()`, `getConcordants()`) are isolated from core analytical functions. This design prevents unintended server load, simplifies error handling for network issues, and enhances the overall stability and reliability of the package during batch processing.

## Implementation via S4 Object-Oriented System

`drugfindR` leverages R's **S4 object-oriented system** to ensure rigorous data structure validation and formal method dispatch. This implementation provides a stable and predictable framework for handling the complex data objects central to transcriptomic analysis. A significant advantage of this approach is native interoperability with common Bioconductor classes and widespread data container classes. The package is designed to work fluidly with `base::data.frame`, `tibble::tibble`, and the Bioconductor-standard `S4Vectors::DataFrame`, allowing researchers to integrate `drugfindR` into diverse existing workflows with minimal data reformatting overhead.

## Core API and Wrapper Functions

To serve both programmable and rapid-prototyping use cases, `drugfindR` offers two distinct layers of interaction: a granular core API and simplified wrapper functions.

**Core Modular Functions:** These five functions provide complete control over the signature matching pipeline.

| Function | Type | Description |
| :--- | :--- | :--- |
| `getSignature()` | Core | Retrieves a standardized L1000 gene expression signature from the iLINCS database using a unique identifier. |
| `prepareSignature()` | Core | Processes a user-supplied differential gene expression (DGE) list into a curated, analysis-ready L1000 signature object. |
| `filterSignature()` | Core | Applies user-defined thresholds (e.g., p-value, fold-change) to focus subsequent analysis on the most statistically robust genes in a signature. |
| `getConcordants()` | Core | Queries the iLINCS API to retrieve a list of all LINCS perturbation signatures that are statistically concordant (positively or inversely correlated) with the input signature. |
| `consensusConcordants()` | Core | Integrates results from multiple related queries to generate a consensus-ranked list of top candidate drugs or genes, aggregating evidence for robust prioritization. |

**Higher-Level Wrapper Functions:** These functions encapsulate the complete core pipeline with sensible, validated defaults, enabling users to perform sophisticated analyses with minimal code.

| Function | Type | Description |
| :--- | :--- | :--- |
| `investigateTarget()` | Wrapper | Executes a full analysis on all signatures associated with a given drug or gene target (e.g., "vorinostat"), identifying other perturbations with similar mechanistic effects. |
| `investigateSignature()` | Wrapper | The primary entry point for user data. Accepts a raw DGE signature, automatically processes it through the standard pipeline (`prepare` -> `filter` -> `query` -> `rank`), and returns a list of candidate repurposing drugs. |

This two-tiered architecture ensures that `drugfindR` is both a library for building custom bioinformatics pipelines and an efficient toolkit for direct therapeutic hypothesis generation.

# Research Impact and Significance

`drugfindR` emerged as the formalization of methods developed for research projects investigating therapies for COVID-19 [@RN3234; @RN654; @RN614]. The identification of the anti-inflammatory properties of fluoxetine in COVID-19 was particularly significant since this finding was later validated in in-vitro studies [@RN9954], representing a major shift in our understanding of both antidepressants and depression.

Since then, evolving versions of `drugfindR` have supported investigations and research projects across fields as diverse as cancer research and neuropsychiatry [@RN9507; @RN9403; @RN9505; @RN9504; @RN9952].

# Open Source Practice

`drugfindR` is built for open-source from the beginning. The codebase is licensed under the GNU General Public License version 3 (GPLv3-or-later), allowing freedom to obtain, use and extend the codebase. In addition, the project also follows the package guidelines published by the Bioconductor projects for inclusion into their repository.
This includes running continuous integration checks on every push and open development using git and GitHub.

# AI Usage Disclosure

No generative AI tools were used in the development of the software code. Generative AI (Model: Gemini) was used for copy-editing and structural refinement of the manuscript text.

# Acknowledgements

We acknowledge Dr. Sinead O'Donovan for her early conceptual contributions. ASI expresses gratitude to Saadia Farooq, Romana Imami, and Samara Imami for their patience and support during the development of this software.

# Funding

This work was supported by NIH grants NIMH R01 MH107487 and NIMH R01 MH121102.

# References
