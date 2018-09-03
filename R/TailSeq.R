# Build and Reload Package: 'Cmd + Shift + B' Check Package:
# 'Cmd + Shift + E' Test Package: 'Cmd + Shift + T'

#' @title TailSeq Calculating Trunncation, Extension, and Tail of RNA Molecules
#'
#' @description Current packages that take care of sequence-based RNA or DNA sequence although very effective in many of applications, namely Biostrings; in our experience, they cannot capture the specific nuisances 3?\200\231 sequence analysis up to the details we had in our mind. To address this issue, we developed an R package called TailSeq. Please refer to the package manual for details.
#'
#' @param file_path, primerSeq, adaptorSeq, distPrimerFromEoG, tailNucleotide, tailThreshold
#'
#' @return dat_normal
#'
#' @export
TailSeq <- function(file_path, primerSeq, adaptorSeq, distPrimerFromEoG, 
    tailNucleotide, tailThreshold) {
    dat <- read.table(file_path, header = F, fill = T)
    colnames(dat) <- "OriginalSequence"
    dat$OriginalSequence <- as.character(dat$OriginalSequence)
    primerSeq <- primerSeq
    adaptorSeq <- adaptorSeq
    
    pbapply::pboptions(type = "timer", char = "=")
    
    cat("Calculate Primer Incidence and Locations")
    
    dat$NumberOfPrimer <- pbapply::pbsapply(dat$OriginalSequence, 
        function(x) dim(as.data.frame(stringr::str_locate_all(x, 
            primerSeq)))[1])
    
    cat("Calculate Adaptor Incidence and Locations")
    
    dat$NumberOfAdaptor <- pbapply::pbsapply(dat$OriginalSequence, 
        function(x) dim(as.data.frame(stringr::str_locate_all(x, 
            adaptorSeq)))[1])
    
    dat$NormalReads <- (dat$NumberOfPrimer == 1) == TRUE & (dat$NumberOfAdaptor == 
        1) == TRUE
    
    pdf(paste0(zonator::file_path_sans_ext(file_path), ".pdf", 
        collapse = ""), onefile = T)
    pie(table(dat$NumberOfAdaptor), main = "Number of Adaptor per reads")
    pie(table(dat$NumberOfPrimer), main = "Number of Primer per read")
    pie(table(dat$NormalReads), labels = 2, main = "Number of Good Reads")
    dev.off()
    
    cat(c("Number of Adaptor per Read: ", table(dat$NumberOfAdaptor), 
        "\n"))
    cat(c("Number of Primer per Read: ", table(dat$NumberOfPrimer), 
        "\n"))
    cat(c("Number of Normal Reads: ", table(dat$NormalReads), 
        "\n"))
    dat_normal <- dat[dat$NormalReads, ]
    cat("\n")
    cat(paste0(c("Number of Reads Excluded Due to Multiple Adaptor or Primer", 
        dim(dat)[1] - dim(dat_normal)[1]), collapse = " : "))
    dat_normal$end_of_primer <- pbapply::pbsapply(dat_normal$OriginalSequence, 
        function(x) stringr::str_locate_all(x, primerSeq)[[1]][2])
    dat_normal$start_of_adaptor <- pbapply::pbsapply(dat_normal$OriginalSequence, 
        function(x) stringr::str_locate_all(x, paste0(adaptorSeq, 
            ".*", collapse = ""))[[1]][1])
    flank <- function(string, eOFp, sOFa) {
        stringr::str_sub(string = string, start = eOFp, end = sOFa)
    }
    
    dat_normal$flank <- base::mapply(flank, dat_normal$OriginalSequence, 
        dat_normal$end_of_primer + 1, dat_normal$start_of_adaptor - 
            1)
    tailPattern <- paste0(c(rep(tailNucleotide, tailThreshold), 
        ".*"), collapse = "")
    dat_normal$woTail <- gsub(tailPattern, "", dat_normal$flank)
    
    dat_normal$sizeOfTail <- nchar(dat_normal$flank) - nchar(dat_normal$woTail)
    
    dat_normal$sizeofExt <- nchar(dat_normal$woTail) - distPrimerFromEoG
    
    write.csv(dat_normal, paste0(file_path_sans_ext(file_path), 
        "processed_TailSeqv2.csv", collapse = ""))
    write.csv(table(dat_normal$sizeofExt), paste0(file_path_sans_ext(file_path), 
        "SizeOfExt_processed_TailSeqv2.csv", collapse = ""))
    write.csv(table(dat_normal$sizeOfTail), paste0(file_path_sans_ext(file_path), 
        "SizeOfTail_processed_TailSeqv2.csv", collapse = ""))
    return(dat_normal)
}

