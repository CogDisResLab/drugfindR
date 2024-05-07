## code to prepare `metadata_tables` dataset goes here
suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(tidyverse)
})

renameColumns <- function(inputNames) {
    if ("Perturbagen" %in% inputNames) {
        colNames <- c(
            "SourceSignature", "Source", "SourceConcentration",
            "SourceCellLine", "SourceTime"
        )
    } else {
        colNames <- c(
            "SourceSignature", "Source",
            "SourceCellLine", "SourceTime"
        )
    }

    colNames
}

oeMetadata <- read_tsv("raw/LINCS-OE-Metadata.xls.gz") %>%
    dplyr::select(-CGSID, -is_exemplar) %>%
    rename_with(renameColumns)

kdMetadata <- read_tsv("raw/LINCS-KD-Metadata.xls.gz") %>%
    dplyr::select(-CGSID, -is_exemplar) %>%
    rename_with(renameColumns)

cpMetadata <- read_tsv("raw/LINCS-Perturbagen-Metadata.xls.gz") %>%
    dplyr::select(-GeneTargets, -is_exemplar) %>%
    rename_with(renameColumns)

l1000List <- read_tsv("raw/l1000_genes.tsv.gz")

l1000 <- l1000List %>%
    pull(HGNC) %>%
    unique() %>%
    AnnotationDbi::select(org.Hs.eg.db,
        keys = .,
        columns = c("ENTREZID", "ALIAS", "SYMBOL"),
        keytype = "SYMBOL"
    ) %>%
    inner_join(l1000List, by = c(SYMBOL = "HGNC")) %>%
    dplyr::select(ENTREZID, L1000, SYMBOL, ALIAS)

usethis::use_data(oeMetadata, kdMetadata, cpMetadata, l1000,
    internal = TRUE, overwrite = TRUE, compress = "bzip2", version = 2L
)
