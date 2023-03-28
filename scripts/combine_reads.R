combine_reads <- function(samplesheet_path) {
  samples <-
    read.delim(samplesheet_path) |>
    dplyr::mutate(combine_reads_command =
                    paste0("zcat /data/GridION/GridIONOutput/",
                           run_name, "/", flowcell_name, "/*/fastq_pass/",
                           # barcode_name, "/*.fastq.gz > /analyses/Diagnostics/enterovirus_typing/ONT_amplicon/reads_raw/",
                           barcode_name, "/*.fastq.gz > reads_raw/",
                           sample_name, "_fastq_pass.fastq"))
  
  purrr::walk(samples$combine_reads_command, system)
}

combine_reads(snakemake@input[["samplesheet_path"]])