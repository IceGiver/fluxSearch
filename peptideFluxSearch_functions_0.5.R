require(fuzzyjoin)
require(data.table)
require(dplyr)

#====================================================================================
# Build the isotopomer reference files (13 carbon)
# Protien/peptide files only have positive mode
refLib_pep <- function(refInput) {
        refPeaklist <- fread(refInput)
        # extract all the protein/peptide info for the reference file
        # Re-organize the refPeaklist
        mz2 <- refPeaklist$`Observed m/z`
        rt2 <- (refPeaklist$`Retention Time (sec)`)/60
        Charge <- refPeaklist$`Spectrum charge`
        Mass <- refPeaklist$`Actual peptide mass (AMU)`
        Seqec <- refPeaklist$`Peptide sequence`
        Modification <- refPeaklist$`Variable modifications identified by spectrum`
        Proteins <- refPeaklist$`Protein accession numbers`
        Protein_names <- refPeaklist$`Protein name`
        
        refPeaklist <- data.frame(mz2, rt2, Charge, Mass, Seqec, Modification,
                                  Proteins, Protein_names)
        # Make sure the type of all the numeric cols
        refPeaklist$mz2 <- as.numeric(refPeaklist$mz2)
        refPeaklist$rt2 <- as.numeric(refPeaklist$rt2)
        refPeaklist$Charge <- as.numeric(refPeaklist$Charge)
        refPeaklist$Mass <- as.numeric(refPeaklist$Mass)
        
        return(refPeaklist)
}

# Create 13 reference files for 13C protein/peptide flux search
#=====================================================================================
# All the features are positive
ref_13CPep <- function(refInput, refOutput) {
        C13 <- 1.0033548378
        ref_pep <- refLib_pep(refInput)
        ref_12C <- data.frame(ref_pep, Annotation = "12C")
        ref_13C1 <- data.frame(ref_pep, Annotation = "13C_01")
        ref_13C1$mz2 <- ref_pep$mz2 + (C13/ref_pep$Charge)
        ref_13C2 <- data.frame(ref_pep, Annotation = "13C_02")
        ref_13C2$mz2 <- ref_pep$mz2 + (2*C13/ref_pep$Charge)
        ref_13C3 <- data.frame(ref_pep, Annotation = "13C_03")
        ref_13C3$mz2 <- ref_pep$mz2 + (3*C13/ref_pep$Charge)
        ref_13C4 <- data.frame(ref_pep, Annotation = "13C_04")
        ref_13C4$mz2 <- ref_pep$mz2 + (4*C13/ref_pep$Charge)
        ref_13C5 <- data.frame(ref_pep, Annotation = "13C_05")
        ref_13C5$mz2 <- ref_pep$mz2 + (5*C13/ref_pep$Charge)
        ref_13C6 <- data.frame(ref_pep, Annotation = "13C_06")
        ref_13C6$mz2 <- ref_pep$mz2 + (6*C13/ref_pep$Charge)
        ref_13C7 <- data.frame(ref_pep, Annotation = "13C_07")
        ref_13C7$mz2 <- ref_pep$mz2 + (7*C13/ref_pep$Charge)
        ref_13C8 <- data.frame(ref_pep, Annotation = "13C_08")
        ref_13C8$mz2 <- ref_pep$mz2 + (8*C13/ref_pep$Charge)
        ref_13C9 <- data.frame(ref_pep, Annotation = "13C_09")
        ref_13C9$mz2 <- ref_pep$mz2 + (9*C13/ref_pep$Charge)
        ref_13C10 <- data.frame(ref_pep, Annotation = "13C_10")
        ref_13C10$mz2 <- ref_pep$mz2 + (10*C13/ref_pep$Charge)
        ref_13C11 <- data.frame(ref_pep, Annotation = "13C_11")
        ref_13C11$mz2 <- ref_pep$mz2 + (11*C13/ref_pep$Charge)
        ref_13C12 <- data.frame(ref_pep, Annotation = "13C_12")
        ref_13C12$mz2 <- ref_pep$mz2 + (12*C13/ref_pep$Charge)
        ref_13C13 <- data.frame(ref_pep, Annotation = "13C_13")
        ref_13C13$mz2 <- ref_pep$mz2 + (13*C13/ref_pep$Charge)
        
        ref_pep <- rbind(ref_12C, ref_13C1, ref_13C2, ref_13C3, ref_13C4, ref_13C5,
                         ref_13C6, ref_13C7, ref_13C8, ref_13C9, ref_13C10, ref_13C11,
                         ref_13C12, ref_13C13)
        return(ref_pep)
        
}


#======================================================================================
# The fuzzy join merge parameter settings

mmf <- function(x, y) {
        # The differnce between data vs. reference
        delta_mz <- abs(x[,1] - y[,1])
        mz_dist <- 1000000*abs(x[,1] - y[,1])/abs(x[,1])
        rt_dist <- abs(x[,2] - y[,2])
        out <- data_frame(merge = rt_dist <= 2 & mz_dist <= 10,
                          Score = 1 - sqrt((mz_dist*0.1)^2 + (rt_dist)^2),
                          delmz_ppm = mz_dist)
        return(out)
}

#Conditional matching the mz and RT from peaklist to the personalized references

msMatch <- function(mzFile, refFile) {
        resultFile <- fuzzy_join(
                mzFile, 
                refFile, 
                multi_by = c("mz1" = "mz2", "rt1" = "rt2"), 
                multi_match_fun = mmf,
                mode = "full")
        #Debug using        
        if(length(resultFile$score) != 0) {
                res_order <- resultFile[order(-resultFile$score,
                                              resultFile$delmz_ppm),]
                resultReducedFile <- na.omit(res_order)
                uniResult <- resultReducedFile[!duplicated(resultReducedFile[, c(8, 9, 10, 11)]),]
                return(uniResult)
        } else {
                score <- NA
                delmz_ppm <- NA
                uniResult_na <- data.frame(resultFile, score, delmz_ppm)
                return(uniResult_na)
        }
}

#======================================================================================
# Process the feature searching and wrap the result files

resultWrap <- function(inputFile) {
        inputPeaklist <- fread(inputFile)
        mz1 <- inputPeaklist$"m/z"
        rt1 <- (inputPeaklist$"RT")/60
        Intensity <- inputPeaklist$"max_int"
        mzList <- data.frame(mz1, rt1, Intensity)
        
        #==========================================
        # create all the positive mode results
        result_12C <- msMatch(mzList, ref_12C)
        result_13C1 <- msMatch(mzList, ref_13C1)
        result_13C2 <- msMatch(mzList, ref_13C2)
        result_13C3 <- msMatch(mzList, ref_13C3)
        result_13C4 <- msMatch(mzList, ref_13C4)
        result_13C5 <- msMatch(mzList, ref_13C5)
        result_13C6 <- msMatch(mzList, ref_13C6)
        result_13C7 <- msMatch(mzList, ref_13C7)
        result_13C8 <- msMatch(mzList, ref_13C8)
        result_13C9 <- msMatch(mzList, ref_13C9)
        result_13C10 <- msMatch(mzList, ref_13C10)
        result_13C11 <- msMatch(mzList, ref_13C11)
        result_13C12 <- msMatch(mzList, ref_13C12)
        result_13C13 <- msMatch(mzList, ref_13C13)
        
        
        result_all <- rbind(result_12C, result_13C1, result_13C2,
                            result_13C3, result_13C4, result_13C5,
                            result_13C6, result_13C7, result_13C8,
                            result_13C9, result_13C10, result_13C11,
                            result_13C12, result_13C13)
        return(result_all)
}

#=====================================================================================
# Wrap all the files together and scoring

FLux_result <- function(input_pep, refer_pep, score = 0.5) {
        flux_pep <- resultWrap(input_pep, refer_pep)
        output_pep <- Sgrade(flux_pep)
        # Select the featues with good grades
        output_pep <- subset(output_pep, output_pep$Score >= score)
        
        return(output_pep)
        
}

# Grading function

Sgrade <- function(inputResult) {
        inputResult$Grades <- cut(inputResult$Score, 
                                  breaks = c(-2, 0, 0.5, 0.7, 0.8, 0.9, 1),
                                  labels = c("F", "D", "C", "B", "A", "A+"))
        return(inputResult)
}












