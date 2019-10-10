require(fuzzyjoin)
require(data.table)
require(dplyr)





#====================================================================
# Build the isotopomer reference files (11 carbon)
# Positive mode
refLibPos_ele <- function (refInput) {
        refPeaklist <- fread(refInput)
        refPeaklist <- subset(refPeaklist, refPeaklist$'ID Score' >= 0.7)
        # extract the positive peaks for reference positive files
        refPeaklist_pos <- subset(refPeaklist, refPeaklist$Charge == "1")
        # re-organize the refPeaklist_pos
        mz2 <- refPeaklist_pos$'m/z'
        rt2 <- refPeaklist_pos$'Aligned RT'
        Metabolite <- refPeaklist_pos$'Metabolite'
        Accession <- refPeaklist_pos$'Accession'
        Theoretical_mz <- refPeaklist_pos$'Theoretical m/z'
        Adduct <- refPeaklist_pos$'Adduct'
        Charge <- refPeaklist_pos$'Charge'
        
        ref_pos <- data.frame(mz2, rt2, Metabolite, Accession, 
                              Theoretical_mz, Adduct, Charge)
        ref_pos$mz2 <- as.numeric(ref_pos$mz2)
        ref_pos$rt2 <- as.numeric(ref_pos$rt2)
        return(ref_pos)
}

# Negative mode
refLibNeg_ele <- function (refInput) {
        refPeaklist <- fread(refInput)
        refPeaklist <- subset(refPeaklist, refPeaklist$'ID Score' >= 0.7)
        # extract the positive peaks for reference positive files
        refPeaklist_neg <- subset(refPeaklist, refPeaklist$Charge == "-1")
        # re-organize the refPeaklist_neg
        mz2 <- refPeaklist_neg$'m/z'
        rt2 <- refPeaklist_neg$'Aligned RT'
        Metabolite <- refPeaklist_neg$'Metabolite'
        Accession <- refPeaklist_neg$'Accession'
        Theoretical_mz <- refPeaklist_neg$'Theoretical m/z'
        Adduct <- refPeaklist_neg$'Adduct'
        Charge <- refPeaklist_neg$'Charge'
        
        ref_neg <- data.frame(mz2, rt2, Metabolite, Accession, 
                              Theoretical_mz, Adduct, Charge)
        ref_neg$mz2 <- as.numeric(ref_neg$mz2)
        ref_neg$rt2 <- as.numeric(ref_neg$rt2)
        return(ref_neg)
}

# Create 11 reference files for 13C metabolite flux search
#===============================================================================
# negative references
ref_13C_neg <- function(refInput, refOutput) {
        C13 <- 1.0033548378
        ref_neg <- refLibNeg_ele(refInput)
        ref_nC12 <- data.frame(ref_neg, Annotation = "12C")
        ref_n1C13 <- data.frame(ref_neg, Annotation = "13C_01")
        ref_n1C13$mz2 <- ref_neg$mz2 + C13
        ref_n2C13 <- data.frame(ref_neg, Annotation = "13C_02")
        ref_n2C13$mz2 <- ref_neg$mz2 + C13*2
        ref_n3C13 <- data.frame(ref_neg, Annotation = "13C_03")
        ref_n3C13$mz2 <- ref_neg$mz2 + C13*3
        ref_n4C13 <- data.frame(ref_neg, Annotation = "13C_04")
        ref_n4C13$mz2 <- ref_neg$mz2 + C13*4
        ref_n5C13 <- data.frame(ref_neg, Annotation = "13C_05")
        ref_n5C13$mz2 <- ref_neg$mz2 + C13*5
        ref_n6C13 <- data.frame(ref_neg, Annotation = "13C_06")
        ref_n6C13$mz2 <- ref_neg$mz2 + C13*6
        ref_n7C13 <- data.frame(ref_neg, Annotation = "13C_07")
        ref_n7C13$mz2 <- ref_neg$mz2 + C13*7
        ref_n8C13 <- data.frame(ref_neg, Annotation = "13C_08")
        ref_n8C13$mz2 <- ref_neg$mz2 + C13*8
        ref_n9C13 <- data.frame(ref_neg, Annotation = "13C_09")
        ref_n9C13$mz2 <- ref_neg$mz2 + C13*9
        ref_n10C13 <- data.frame(ref_neg, Annotation = "13C_10")
        ref_n10C13$mz2 <- ref_neg$mz2 + C13*10
        ref_n11C13 <- data.frame(ref_neg, Annotation = "13C_11")
        ref_n11C13$mz2 <- ref_neg$mz2 + C13*11

        ref_neg_13C <- rbind(ref_nC12, ref_n1C13, ref_n2C13, ref_n3C13, ref_n4C13, ref_n5C13,
                             ref_n6C13, ref_n7C13, ref_n8C13, ref_n9C13, ref_n10C13, ref_n11C13)
        return(ref_neg_13C)
}

# Positive references
ref_13C_pos <- function(refInput, refOutput) {
        C13 <- 1.0033548378
        ref_pos <- refLibPos_ele(refInput)
        ref_pC12 <- data.frame(ref_pos, Annotation = "12C")
        ref_p1C13 <- data.frame(ref_pos, Annotation = "13C_01")
        ref_p1C13$mz2 <- ref_pos$mz2 + C13
        ref_p2C13 <- data.frame(ref_pos, Annotation = "13C_02")
        ref_p2C13$mz2 <- ref_pos$mz2 + C13*2
        ref_p3C13 <- data.frame(ref_pos, Annotation = "13C_03")
        ref_p3C13$mz2 <- ref_pos$mz2 + C13*3
        ref_p4C13 <- data.frame(ref_pos, Annotation = "13C_04")
        ref_p4C13$mz2 <- ref_pos$mz2 + C13*4
        ref_p5C13 <- data.frame(ref_pos, Annotation = "13C_05")
        ref_p5C13$mz2 <- ref_pos$mz2 + C13*5
        ref_p6C13 <- data.frame(ref_pos, Annotation = "13C_06")
        ref_p6C13$mz2 <- ref_pos$mz2 + C13*6
        ref_p7C13 <- data.frame(ref_pos, Annotation = "13C_07")
        ref_p7C13$mz2 <- ref_pos$mz2 + C13*7
        ref_p8C13 <- data.frame(ref_pos, Annotation = "13C_08")
        ref_p8C13$mz2 <- ref_pos$mz2 + C13*8
        ref_p9C13 <- data.frame(ref_pos, Annotation = "13C_09")
        ref_p9C13$mz2 <- ref_pos$mz2 + C13*9
        ref_p10C13 <- data.frame(ref_pos, Annotation = "13C_10")
        ref_p10C13$mz2 <- ref_pos$mz2 + C13*10
        ref_p11C13 <- data.frame(ref_pos, Annotation = "13C_11")
        ref_p11C13$mz2 <- ref_pos$mz2 + C13*11
        
        ref_pos_13C <- rbind(ref_pC12, ref_p1C13, ref_p2C13, ref_p3C13, ref_p4C13, ref_p5C13,
                             ref_p6C13, ref_p7C13, ref_p8C13, ref_p9C13, ref_p10C13, ref_p11C13)
        return(ref_pos_13C)
}

#==========================================================================================
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

# Conditional matching the mz and RT from peaklist to the personalized references

msMatch <- function(mzFile, refFile) {
        resultFile <- fuzzy_join(
                mzFile,
                refFile,
                multi_by = c("mz1" = "mz2", "rt1" = "rt2"),
                multi_match_fun = mmf,
                mode = "full"
        )
        
        if(length(resultFile$Score) != 0) {
                res_order <- resultFile[order(-resultFile$Score,
                                              resultFile$delmz_ppm), ]
                resultReducedFile <- na.omit(res_order)
                uniResult <- resultReducedFile[!duplicated(resultReducedFile[, c(7, 9, 10)]),]
                return(uniResult)
        } else {
                Score <- NA
                delmz_ppm <- NA
                uniResult_na <- data.frame(resultFile, Score, delmz_ppm)
                return(uniResult_na)
        }
        
}


#===========================================================================================
# Process the feature searching and wrap the result files
# Negative files mz matching
resultWrap_neg <- function(inputFile_neg, refInput) {
        inputPeaklist <- fread(inputFile_neg)
        reference <- ref_13C_neg(refInput)
        
        mz1 <- as.numeric(inputPeaklist$"m/z")
        rt1 <- (as.numeric(inputPeaklist$"RT")/60)
        Intensity <- inputPeaklist$"max_int"
        mzList <- data.frame(mz1, rt1, Intensity)
        
        # fuzzy join all the negative files
        
        result_nC12 <- msMatch(mzList, reference[reference$Annotation == "12C",])
        result_n1C13 <- msMatch(mzList, reference[reference$Annotation == "13C_01", ])
        result_n2C13 <- msMatch(mzList, reference[reference$Annotation == "13C_02", ])
        result_n3C13 <- msMatch(mzList, reference[reference$Annotation == "13C_03", ])
        result_n4C13 <- msMatch(mzList, reference[reference$Annotation == "13C_04", ])
        result_n5C13 <- msMatch(mzList, reference[reference$Annotation == "13C_05", ])
        result_n6C13 <- msMatch(mzList, reference[reference$Annotation == "13C_06", ])
        result_n7C13 <- msMatch(mzList, reference[reference$Annotation == "13C_07", ])
        result_n8C13 <- msMatch(mzList, reference[reference$Annotation == "13C_08", ])
        result_n9C13 <- msMatch(mzList, reference[reference$Annotation == "13C_09", ])
        result_n10C13 <- msMatch(mzList, reference[reference$Annotation == "13C_10", ])
        result_n11C13 <- msMatch(mzList, reference[reference$Annotation == "13C_11", ])
        
        #Wrap all the results as one output dataframe
        result_all_neg <- rbind(result_nC12, result_n1C13, result_n2C13, result_n3C13, result_n4C13,
                                result_n5C13, result_n6C13, result_n7C13, result_n8C13, result_n9C13,
                                result_n10C13, result_n11C13)
        return(result_all_neg)
}

# Positive files mz matching
resultWrap_pos <- function(inputFile_pos, refInput) {
        inputPeaklist <- fread(inputFile_pos)
        reference <- ref_13C_pos(refInput)
        
        mz1 <- as.numeric(inputPeaklist$"m/z")
        rt1 <- (as.numeric(inputPeaklist$"RT")/60)
        Intensity <- inputPeaklist$"max_int"
        mzList <- data.frame(mz1, rt1, Intensity)
        
        # fuzzy join all the negative files
        result_pC12 <- msMatch(mzList, reference[reference$Annotation == "12C", ])
        result_p1C13 <- msMatch(mzList, reference[reference$Annotation == "13C_01", ])
        result_p2C13 <- msMatch(mzList, reference[reference$Annotation == "13C_02", ])
        result_p3C13 <- msMatch(mzList, reference[reference$Annotation == "13C_03", ])
        result_p4C13 <- msMatch(mzList, reference[reference$Annotation == "13C_04", ])
        result_p5C13 <- msMatch(mzList, reference[reference$Annotation == "13C_05", ])
        result_p6C13 <- msMatch(mzList, reference[reference$Annotation == "13C_06", ])
        result_p7C13 <- msMatch(mzList, reference[reference$Annotation == "13C_07", ])
        result_p8C13 <- msMatch(mzList, reference[reference$Annotation == "13C_08", ])
        result_p9C13 <- msMatch(mzList, reference[reference$Annotation == "13C_09", ])
        result_p10C13 <- msMatch(mzList, reference[reference$Annotation == "13C_10", ])
        result_p11C13 <- msMatch(mzList, reference[reference$Annotation == "13C_11", ])

        #Wrap all the results as one output dataframe
        result_all_pos <- rbind(result_pC12, result_p1C13, result_p2C13, result_p3C13, result_p4C13,
                                result_p5C13, result_p6C13, result_p7C13, result_p8C13, result_p9C13,
                                result_p10C13, result_p11C13)
        return(result_all_pos)
}



#===================================================================================================
# Wrap all positive and negative files together

Flux_result <- function(input_negative, input_positive, referInput, score = 0.8) {
        flux_neg <- resultWrap_neg(input_negative, referInput)
        flux_pos <- resultWrap_pos(input_positive, referInput)
        output_all <- rbind(flux_neg, flux_pos)
        output_all <- Sgrade(output_all)
        # select feature with good grades
        output_all <- subset(output_all, output_all$Score >= score)
        
        return(output_all)
}


# Grading function
Sgrade <- function(inputResult) {
        inputResult$Grades <- cut(inputResult$Score,
                                  breaks = c(-2, 0, 0.5, 0.7, 0.8, 0.9, 1),
                                  labels = c("F", "D", "C", "B", "A", "A+"))
        return(inputResult)
}























