# Table 2: Results ----
# Check dependencies
check_packages(
  cran_packages = c("ggplot2"),
  bioc_packages = c("UniProt.ws")
)

# Get data
table2_results_1 <- subset(
  merge(
    tableS2_xlsx$Model1,
    tableS2_xlsx$Model2,
    by = c("name1", "panel", "uniprot_id", "olink_id", "LOD", "below_LOD_freq", "data_type")
  ),
  p.x < 0.05/122
)

# Get ENSG from uniprot
table2_up <- UniProt.ws::UniProt.ws(taxId = 9606)
table2_keys <- table2_results_1$uniprot_id
table2_columns <- c("ENSEMBL")
table2_kt <- "UNIPROTKB"
table2_upres <- UniProt.ws::select(table2_up, table2_keys, table2_columns, table2_kt)

table2_results_2 <- merge(
  table2_upres,
  table2_results_1,
  by.x = "UNIPROTKB",
  by.y = "uniprot_id"
)

# Load human protein atlas 20.0
proteinatlas <- read.delim(
  "~/projekt_data/2019-11-19_SCD-FIA3_data/extdata/proteinatlas.tsv"
)

# Merge tables
table2 <- merge(
  proteinatlas,
  table2_results_2,
  by.x = "Ensembl",
  by.y = "ENSEMBL"
)

# Extract predicted location
Predicted.location <- sapply(
  table2$Protein.class, 
  function(x) {
    strsplit(x, ", ")
  }
)

Predicted.location2 <- lapply(
  Predicted.location, 
  function(x) {
  grep("^Predicted", x, value = TRUE)
  }
)

Predicted.location3 <- lapply(
  Predicted.location2, 
  function(x) {
  paste(x, collapse = ", ")
  }
)

Predicted.location4 <- Reduce(
  rbind,
  Predicted.location3
)

Predicted.location5 <- sapply(
  Predicted.location4,
  function(x) {
    y <- gsub("Predicted ", "", x)
    z <- gsub(" proteins", "", y)
    return(z)
  }
)
  
table2$Predicted.location <- Predicted.location5

# Finalize table
table2.out.vars <- c(
  "Gene", 
  "Gene.synonym",
  "Gene.description",
  #"Biological.process", 
  #"Molecular.function",
  "Secretome.location"
)
table2.fin  <- table2[table2.out.vars]

# Format tissue specificity
table2.fin$RNA.tissue.specificity <- mapply(
  function(x, y) {
    out <- paste0(x, " ", y)
    return(out)
  },
  table2$RNA.tissue.specificity,
  table2$RNA.tissue.specific.NX
)

# Format model statistics
table2.fin$OR.1 <- mapply(
  function(x, y, z) {
    out <- paste0(x, " [", y, " - ", z, "]")
    return(out)
  },
  table2$OR.x,
  table2$ci_low.x,
  table2$ci_high.x
)
table2.fin$p.1 <- table2$p.x
table2.fin$OR.2 <- mapply(
  function(x, y, z) {
    out <- paste0(x, " [", y, " - ", z, "]")
    return(out)
  },
  table2$OR.y,
  table2$ci_low.y,
  table2$ci_high.y
)
table2.fin$p.2 <- table2$p.y

# Reorder table
table2.fin <- table2.fin[order(table2.fin$p.2), ]
