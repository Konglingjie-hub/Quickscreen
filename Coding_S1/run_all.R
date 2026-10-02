# Run the retained FoG example by default. --all runs models on user-supplied CSVs.
if (.Platform$OS.type == "windows") suppressWarnings(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent=TRUE))
a <- grep("^--file=", commandArgs(FALSE), value=TRUE)
root <- if(length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), winslash="/")) else getwd()
setwd(root)
needed <- c("metafor","forestplot")
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Install required packages: ",paste(missing,collapse=", "))
dir.create("outputs",showWarnings=FALSE)
source("scripts/SMD_metaanlysis.R", local=new.env())
if("--all" %in% commandArgs(TRUE)) {
  required <- c("motor_overall.csv","motor_medication.csv","motor_high_risk_excluded.csv","frequency_extracted.csv")
  absent <- required[!file.exists(file.path("data",required))]
  if(length(absent)) stop("Provide adapted input datasets in data/: ",paste(absent,collapse=", "))
  source("scripts/export_embedded_fog_data.R")
  source("scripts/analyze.R")
}
writeLines(capture.output(sessionInfo()),"outputs/sessionInfo.txt")
cat("Example analysis completed. Adapt the data and models before scientific use.\n")
