require(fuzzyjoin)
require(data.table)
require(dplyr)





#====================================================================
# Build the isotopomer reference files (40 carbon)
# Positive mode
refLibPos_ls <- function (refInput) {
        refPeaklist <- fread(refInput)
        # extract the positive peaks for reference positive files
        refPeaklist_pos <- subset(refPeaklist, refPeaklist$Pol == "P")
        # re-organize the refPeaklist_pos
        mz2 <- refPeaklist_pos$ObsMz
        rt2 <- refPeaklist_pos$Rt
        Lipid <- refPeaklist_pos$LipidIon
        fattyAcid <- refPeaklist_pos$FattyAcid
        lipidClass <- refPeaklist_pos$Class
        lipidForm <- refPeaklist_pos$Formula
        Theoretical_mz <- refPeaklist_pos$CalcMz
        Adduct <- refPeaklist_pos$Ion
        Charge <- refPeaklist_pos$Pol
        z <- refPeaklist_pos$z
        Grade_ls <- refPeaklist_pos$Grade
        
        ref_pos <- data.frame(mz2, rt2, Lipid, fattyAcid, lipidClass, lipidForm,
                              Theoretical_mz, Adduct, Charge, z, Grade_ls)
        ref_pos$mz2 <- as.numeric(ref_pos$mz2)
        ref_pos$rt2 <- as.numeric(ref_pos$rt2)
        ref_pos$z <- as.numeric(ref_pos$z)
        return(ref_pos)
}

# Negative mode
refLibNeg_ls <- function (refInput) {
        refPeaklist <- fread(refInput)
        # extract the negative peaks for reference negative files
        refPeaklist_neg <- subset(refPeaklist, refPeaklist$Pol == "N")
        # re-organize the refPeaklist_neg
        mz2 <- refPeaklist_neg$ObsMz
        rt2 <- refPeaklist_neg$Rt
        Lipid <- refPeaklist_neg$LipidIon
        fattyAcid <- refPeaklist_neg$FattyAcid
        lipidClass <- refPeaklist_neg$Class
        lipidForm <- refPeaklist_neg$Formula
        Theoretical_mz <- refPeaklist_neg$CalcMz
        Adduct <- refPeaklist_neg$Ion
        Charge <- refPeaklist_neg$Pol
        z <- refPeaklist_neg$z
        Grade_ls <- refPeaklist_neg$Grade
        
        ref_neg <- data.frame(mz2, rt2, Lipid, fattyAcid, lipidClass, lipidForm,
                              Theoretical_mz, Adduct, Charge, z, Grade_ls)
        ref_neg$mz2 <- as.numeric(ref_neg$mz2)
        ref_neg$rt2 <- as.numeric(ref_neg$rt2)
        ref_neg$z <- as.numeric(ref_neg$z)
        return(ref_neg)
}

# Create 40 reference files for 13C lipid flux search
#===============================================================================
# negative references
ref_13C_neg <- function(refInput, refOutput) {
       C13 <- 1.0033548378
       ref_neg <- refLibNeg_ls(refInput)
       ref_nC12 <- data.frame(ref_neg, Annotation = "12C")
       ref_n1C13 <- data.frame(ref_neg, Annotation = "13C_01")
       ref_n1C13$mz2 <- ref_neg$mz2 + (C13/ref_neg$z)
       ref_n2C13 <- data.frame(ref_neg, Annotation = "13C_02")
       ref_n2C13$mz2 <- ref_neg$mz2 + (2*(C13/ref_neg$z))
       ref_n3C13 <- data.frame(ref_neg, Annotation = "13C_03")
       ref_n3C13$mz2 <- ref_neg$mz2 + (3*(C13/ref_neg$z))
       ref_n4C13 <- data.frame(ref_neg, Annotation = "13C_04")
       ref_n4C13$mz2 <- ref_neg$mz2 + (4*(C13/ref_neg$z))
       ref_n5C13 <- data.frame(ref_neg, Annotation = "13C_05")
       ref_n5C13$mz2 <- ref_neg$mz2 + (5*(C13/ref_neg$z))
       ref_n6C13 <- data.frame(ref_neg, Annotation = "13C_06")
       ref_n6C13$mz2 <- ref_neg$mz2 + (6*(C13/ref_neg$z))
       ref_n7C13 <- data.frame(ref_neg, Annotation = "13C_07")
       ref_n7C13$mz2 <- ref_neg$mz2 + (7*(C13/ref_neg$z))
       ref_n8C13 <- data.frame(ref_neg, Annotation = "13C_08")
       ref_n8C13$mz2 <- ref_neg$mz2 + (8*(C13/ref_neg$z))
       ref_n9C13 <- data.frame(ref_neg, Annotation = "13C_09")
       ref_n9C13$mz2 <- ref_neg$mz2 + (9*(C13/ref_neg$z))
       ref_n10C13 <- data.frame(ref_neg, Annotation = "13C_10")
       ref_n10C13$mz2 <- ref_neg$mz2 + (10*(C13/ref_neg$z))
       ref_n11C13 <- data.frame(ref_neg, Annotation = "13C_11")
       ref_n11C13$mz2 <- ref_neg$mz2 + (11*(C13/ref_neg$z))
       ref_n12C13 <- data.frame(ref_neg, Annotation = "13C_12")
       ref_n12C13$mz2 <- ref_neg$mz2 + (12*(C13/ref_neg$z))
       ref_n13C13 <- data.frame(ref_neg, Annotation = "13C_13")
       ref_n13C13$mz2 <- ref_neg$mz2 + (13*(C13/ref_neg$z))
       ref_n14C13 <- data.frame(ref_neg, Annotation = "13C_14")
       ref_n14C13$mz2 <- ref_neg$mz2 + (14*(C13/ref_neg$z))
       ref_n15C13 <- data.frame(ref_neg, Annotation = "13C_15")
       ref_n15C13$mz2 <- ref_neg$mz2 + (15*(C13/ref_neg$z))
       ref_n16C13 <- data.frame(ref_neg, Annotation = "13C_16")
       ref_n16C13$mz2 <- ref_neg$mz2 + (16*(C13/ref_neg$z))
       ref_n17C13 <- data.frame(ref_neg, Annotation = "13C_17")
       ref_n17C13$mz2 <- ref_neg$mz2 + (17*(C13/ref_neg$z))
       ref_n18C13 <- data.frame(ref_neg, Annotation = "13C_18")
       ref_n18C13$mz2 <- ref_neg$mz2 + (18*(C13/ref_neg$z))
       ref_n19C13 <- data.frame(ref_neg, Annotation = "13C_19")
       ref_n19C13$mz2 <- ref_neg$mz2 + (19*(C13/ref_neg$z))
       ref_n20C13 <- data.frame(ref_neg, Annotation = "13C_20")
       ref_n20C13$mz2 <- ref_neg$mz2 + (20*(C13/ref_neg$z))
       ref_n21C13 <- data.frame(ref_neg, Annotation = "13C_21")
       ref_n21C13$mz2 <- ref_neg$mz2 + (21*(C13/ref_neg$z))
       ref_n22C13 <- data.frame(ref_neg, Annotation = "13C_22")
       ref_n22C13$mz2 <- ref_neg$mz2 + (22*(C13/ref_neg$z))
       ref_n23C13 <- data.frame(ref_neg, Annotation = "13C_23")
       ref_n23C13$mz2 <- ref_neg$mz2 + (23*(C13/ref_neg$z))
       ref_n24C13 <- data.frame(ref_neg, Annotation = "13C_24")
       ref_n24C13$mz2 <- ref_neg$mz2 + (24*(C13/ref_neg$z))
       ref_n25C13 <- data.frame(ref_neg, Annotation = "13C_25")
       ref_n25C13$mz2 <- ref_neg$mz2 + (25*(C13/ref_neg$z))
       ref_n26C13 <- data.frame(ref_neg, Annotation = "13C_26")
       ref_n26C13$mz2 <- ref_neg$mz2 + (26*(C13/ref_neg$z))
       ref_n27C13 <- data.frame(ref_neg, Annotation = "13C_27")
       ref_n27C13$mz2 <- ref_neg$mz2 + (27*(C13/ref_neg$z))
       ref_n28C13 <- data.frame(ref_neg, Annotation = "13C_28")
       ref_n28C13$mz2 <- ref_neg$mz2 + (28*(C13/ref_neg$z))
       ref_n29C13 <- data.frame(ref_neg, Annotation = "13C_29")
       ref_n29C13$mz2 <- ref_neg$mz2 + (29*(C13/ref_neg$z))
       ref_n30C13 <- data.frame(ref_neg, Annotation = "13C_30")
       ref_n30C13$mz2 <- ref_neg$mz2 + (30*(C13/ref_neg$z))
       ref_n31C13 <- data.frame(ref_neg, Annotation = "13C_31")
       ref_n31C13$mz2 <- ref_neg$mz2 + (31*(C13/ref_neg$z))
       ref_n32C13 <- data.frame(ref_neg, Annotation = "13C_32")
       ref_n32C13$mz2 <- ref_neg$mz2 + (32*(C13/ref_neg$z))
       ref_n33C13 <- data.frame(ref_neg, Annotation = "13C_33")
       ref_n33C13$mz2 <- ref_neg$mz2 + (33*(C13/ref_neg$z))
       ref_n34C13 <- data.frame(ref_neg, Annotation = "13C_34")
       ref_n34C13$mz2 <- ref_neg$mz2 + (34*(C13/ref_neg$z))
       ref_n35C13 <- data.frame(ref_neg, Annotation = "13C_35")
       ref_n35C13$mz2 <- ref_neg$mz2 + (35*(C13/ref_neg$z))
       ref_n36C13 <- data.frame(ref_neg, Annotation = "13C_36")
       ref_n36C13$mz2 <- ref_neg$mz2 + (36*(C13/ref_neg$z))
       ref_n37C13 <- data.frame(ref_neg, Annotation = "13C_37")
       ref_n37C13$mz2 <- ref_neg$mz2 + (37*(C13/ref_neg$z))
       ref_n38C13 <- data.frame(ref_neg, Annotation = "13C_38")
       ref_n38C13$mz2 <- ref_neg$mz2 + (38*(C13/ref_neg$z))
       ref_n39C13 <- data.frame(ref_neg, Annotation = "13C_39")
       ref_n39C13$mz2 <- ref_neg$mz2 + (39*(C13/ref_neg$z))
       ref_n40C13 <- data.frame(ref_neg, Annotation = "13C_40")
       ref_n40C13$mz2 <- ref_neg$mz2 + (40*(C13/ref_neg$z))
       
       ref_neg_13C <- rbind(ref_nC12, ref_n1C13, ref_n2C13, ref_n3C13, ref_n4C13, ref_n5C13,
                            ref_n6C13, ref_n7C13, ref_n8C13, ref_n9C13, ref_n10C13, ref_n11C13,
                            ref_n12C13, ref_n13C13, ref_n14C13, ref_n15C13, ref_n16C13,
                            ref_n17C13, ref_n18C13, ref_n19C13, ref_n20C13, ref_n21C13,
                            ref_n22C13, ref_n23C13, ref_n24C13, ref_n25C13, ref_n26C13,
                            ref_n27C13, ref_n28C13, ref_n29C13, ref_n30C13, ref_n31C13,
                            ref_n32C13, ref_n33C13, ref_n34C13, ref_n35C13, ref_n36C13,
                            ref_n37C13, ref_n38C13, ref_n39C13, ref_n40C13)
       return(ref_neg_13C)
}

# Positive references
ref_13C_pos <- function(refInput, refOutput) {
       C13 <- 1.0033548378
       ref_pos <- refLibPos_ls(refInput)
       ref_pC12 <- data.frame(ref_pos, Annotation = "12C")
       ref_p1C13 <- data.frame(ref_pos, Annotation = "13C_01")
       ref_p1C13$mz2 <- ref_pos$mz2 + (C13/ref_pos$z)
       ref_p2C13 <- data.frame(ref_pos, Annotation = "13C_02")
       ref_p2C13$mz2 <- ref_pos$mz2 + (2*(C13/ref_pos$z))
       ref_p3C13 <- data.frame(ref_pos, Annotation = "13C_03")
       ref_p3C13$mz2 <- ref_pos$mz2 + (3*(C13/ref_pos$z))
       ref_p4C13 <- data.frame(ref_pos, Annotation = "13C_04")
       ref_p4C13$mz2 <- ref_pos$mz2 + (4*(C13/ref_pos$z))
       ref_p5C13 <- data.frame(ref_pos, Annotation = "13C_05")
       ref_p5C13$mz2 <- ref_pos$mz2 + (5*(C13/ref_pos$z))
       ref_p6C13 <- data.frame(ref_pos, Annotation = "13C_06")
       ref_p6C13$mz2 <- ref_pos$mz2 + (6*(C13/ref_pos$z))
       ref_p7C13 <- data.frame(ref_pos, Annotation = "13C_07")
       ref_p7C13$mz2 <- ref_pos$mz2 + (7*(C13/ref_pos$z))
       ref_p8C13 <- data.frame(ref_pos, Annotation = "13C_08")
       ref_p8C13$mz2 <- ref_pos$mz2 + (8*(C13/ref_pos$z))
       ref_p9C13 <- data.frame(ref_pos, Annotation = "13C_09")
       ref_p9C13$mz2 <- ref_pos$mz2 + (9*(C13/ref_pos$z))
       ref_p10C13 <- data.frame(ref_pos, Annotation = "13C_10")
       ref_p10C13$mz2 <- ref_pos$mz2 + (10*(C13/ref_pos$z))
       ref_p11C13 <- data.frame(ref_pos, Annotation = "13C_11")
       ref_p11C13$mz2 <- ref_pos$mz2 + (11*(C13/ref_pos$z))
       ref_p12C13 <- data.frame(ref_pos, Annotation = "13C_12")
       ref_p12C13$mz2 <- ref_pos$mz2 + (12*(C13/ref_pos$z))
       ref_p13C13 <- data.frame(ref_pos, Annotation = "13C_13")
       ref_p13C13$mz2 <- ref_pos$mz2 + (13*(C13/ref_pos$z))
       ref_p14C13 <- data.frame(ref_pos, Annotation = "13C_14")
       ref_p14C13$mz2 <- ref_pos$mz2 + (14*(C13/ref_pos$z))
       ref_p15C13 <- data.frame(ref_pos, Annotation = "13C_15")
       ref_p15C13$mz2 <- ref_pos$mz2 + (15*(C13/ref_pos$z))
       ref_p16C13 <- data.frame(ref_pos, Annotation = "13C_16")
       ref_p16C13$mz2 <- ref_pos$mz2 + (16*(C13/ref_pos$z))
       ref_p17C13 <- data.frame(ref_pos, Annotation = "13C_17")
       ref_p17C13$mz2 <- ref_pos$mz2 + (17*(C13/ref_pos$z))
       ref_p18C13 <- data.frame(ref_pos, Annotation = "13C_18")
       ref_p18C13$mz2 <- ref_pos$mz2 + (18*(C13/ref_pos$z))
       ref_p19C13 <- data.frame(ref_pos, Annotation = "13C_19")
       ref_p19C13$mz2 <- ref_pos$mz2 + (19*(C13/ref_pos$z))
       ref_p20C13 <- data.frame(ref_pos, Annotation = "13C_20")
       ref_p20C13$mz2 <- ref_pos$mz2 + (20*(C13/ref_pos$z))
       ref_p21C13 <- data.frame(ref_pos, Annotation = "13C_21")
       ref_p21C13$mz2 <- ref_pos$mz2 + (21*(C13/ref_pos$z))
       ref_p22C13 <- data.frame(ref_pos, Annotation = "13C_22")
       ref_p22C13$mz2 <- ref_pos$mz2 + (22*(C13/ref_pos$z))
       ref_p23C13 <- data.frame(ref_pos, Annotation = "13C_23")
       ref_p23C13$mz2 <- ref_pos$mz2 + (23*(C13/ref_pos$z))
       ref_p24C13 <- data.frame(ref_pos, Annotation = "13C_24")
       ref_p24C13$mz2 <- ref_pos$mz2 + (24*(C13/ref_pos$z))
       ref_p25C13 <- data.frame(ref_pos, Annotation = "13C_25")
       ref_p25C13$mz2 <- ref_pos$mz2 + (25*(C13/ref_pos$z))
       ref_p26C13 <- data.frame(ref_pos, Annotation = "13C_26")
       ref_p26C13$mz2 <- ref_pos$mz2 + (26*(C13/ref_pos$z))
       ref_p27C13 <- data.frame(ref_pos, Annotation = "13C_27")
       ref_p27C13$mz2 <- ref_pos$mz2 + (27*(C13/ref_pos$z))
       ref_p28C13 <- data.frame(ref_pos, Annotation = "13C_28")
       ref_p28C13$mz2 <- ref_pos$mz2 + (28*(C13/ref_pos$z))
       ref_p29C13 <- data.frame(ref_pos, Annotation = "13C_29")
       ref_p29C13$mz2 <- ref_pos$mz2 + (29*(C13/ref_pos$z))
       ref_p30C13 <- data.frame(ref_pos, Annotation = "13C_30")
       ref_p30C13$mz2 <- ref_pos$mz2 + (30*(C13/ref_pos$z))
       ref_p31C13 <- data.frame(ref_pos, Annotation = "13C_31")
       ref_p31C13$mz2 <- ref_pos$mz2 + (31*(C13/ref_pos$z))
       ref_p32C13 <- data.frame(ref_pos, Annotation = "13C_32")
       ref_p32C13$mz2 <- ref_pos$mz2 + (32*(C13/ref_pos$z))
       ref_p33C13 <- data.frame(ref_pos, Annotation = "13C_33")
       ref_p33C13$mz2 <- ref_pos$mz2 + (33*(C13/ref_pos$z))
       ref_p34C13 <- data.frame(ref_pos, Annotation = "13C_34")
       ref_p34C13$mz2 <- ref_pos$mz2 + (34*(C13/ref_pos$z))
       ref_p35C13 <- data.frame(ref_pos, Annotation = "13C_35")
       ref_p35C13$mz2 <- ref_pos$mz2 + (35*(C13/ref_pos$z))
       ref_p36C13 <- data.frame(ref_pos, Annotation = "13C_36")
       ref_p36C13$mz2 <- ref_pos$mz2 + (36*(C13/ref_pos$z))
       ref_p37C13 <- data.frame(ref_pos, Annotation = "13C_37")
       ref_p37C13$mz2 <- ref_pos$mz2 + (37*(C13/ref_pos$z))
       ref_p38C13 <- data.frame(ref_pos, Annotation = "13C_38")
       ref_p38C13$mz2 <- ref_pos$mz2 + (38*(C13/ref_pos$z))
       ref_p39C13 <- data.frame(ref_pos, Annotation = "13C_39")
       ref_p39C13$mz2 <- ref_pos$mz2 + (39*(C13/ref_pos$z))
       ref_p40C13 <- data.frame(ref_pos, Annotation = "13C_40")
       ref_p40C13$mz2 <- ref_pos$mz2 + (40*(C13/ref_pos$z))
       
       ref_pos_13C <- rbind(ref_pC12, ref_p1C13, ref_p2C13, ref_p3C13, ref_p4C13, ref_p5C13,
                            ref_p6C13, ref_p7C13, ref_p8C13, ref_p9C13, ref_p10C13, ref_p11C13,
                            ref_p12C13, ref_p13C13, ref_p14C13, ref_p15C13, ref_p16C13,
                            ref_p17C13, ref_p18C13, ref_p19C13, ref_p20C13, ref_p21C13,
                            ref_p22C13, ref_p23C13, ref_p24C13, ref_p25C13, ref_p26C13,
                            ref_p27C13, ref_p28C13, ref_p29C13, ref_p30C13, ref_p31C13,
                            ref_p32C13, ref_p33C13, ref_p34C13, ref_p35C13, ref_p36C13,
                            ref_p37C13, ref_p38C13, ref_p39C13, ref_p40C13)
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
        result_n12C13 <- msMatch(mzList, reference[reference$Annotation == "13C_12", ])
        result_n13C13 <- msMatch(mzList, reference[reference$Annotation == "13C_13", ])
        result_n14C13 <- msMatch(mzList, reference[reference$Annotation == "13C_14", ])
        result_n15C13 <- msMatch(mzList, reference[reference$Annotation == "13C_15", ])
        result_n16C13 <- msMatch(mzList, reference[reference$Annotation == "13C_16", ])
        result_n17C13 <- msMatch(mzList, reference[reference$Annotation == "13C_17", ])
        result_n18C13 <- msMatch(mzList, reference[reference$Annotation == "13C_18", ])
        result_n19C13 <- msMatch(mzList, reference[reference$Annotation == "13C_19", ])
        result_n20C13 <- msMatch(mzList, reference[reference$Annotation == "13C_20", ])
        result_n21C13 <- msMatch(mzList, reference[reference$Annotation == "13C_21", ])
        result_n22C13 <- msMatch(mzList, reference[reference$Annotation == "13C_22", ])
        result_n23C13 <- msMatch(mzList, reference[reference$Annotation == "13C_23", ])
        result_n24C13 <- msMatch(mzList, reference[reference$Annotation == "13C_24", ])
        result_n25C13 <- msMatch(mzList, reference[reference$Annotation == "13C_25", ])
        result_n26C13 <- msMatch(mzList, reference[reference$Annotation == "13C_26", ])
        result_n27C13 <- msMatch(mzList, reference[reference$Annotation == "13C_27", ])
        result_n28C13 <- msMatch(mzList, reference[reference$Annotation == "13C_28", ])
        result_n29C13 <- msMatch(mzList, reference[reference$Annotation == "13C_29", ])
        result_n30C13 <- msMatch(mzList, reference[reference$Annotation == "13C_30", ])
        result_n31C13 <- msMatch(mzList, reference[reference$Annotation == "13C_31", ])
        result_n32C13 <- msMatch(mzList, reference[reference$Annotation == "13C_32", ])
        result_n33C13 <- msMatch(mzList, reference[reference$Annotation == "13C_33", ])
        result_n34C13 <- msMatch(mzList, reference[reference$Annotation == "13C_34", ])
        result_n35C13 <- msMatch(mzList, reference[reference$Annotation == "13C_35", ])
        result_n36C13 <- msMatch(mzList, reference[reference$Annotation == "13C_36", ])
        result_n37C13 <- msMatch(mzList, reference[reference$Annotation == "13C_37", ])
        result_n38C13 <- msMatch(mzList, reference[reference$Annotation == "13C_38", ])
        result_n39C13 <- msMatch(mzList, reference[reference$Annotation == "13C_39", ])
        result_n40C13 <- msMatch(mzList, reference[reference$Annotation == "13C_40", ])
        
        #Wrap all the results as one output dataframe
        result_all_neg <- rbind(result_nC12, result_n1C13, result_n2C13, result_n3C13, result_n4C13,
                                result_n5C13, result_n6C13, result_n7C13, result_n8C13, result_n9C13,
                                result_n10C13, result_n11C13, result_n12C13, result_n13C13,
                                result_n14C13, result_n15C13, result_n16C13, result_n17C13,
                                result_n18C13, result_n19C13, result_n20C13, result_n21C13,
                                result_n22C13, result_n23C13, result_n24C13, result_n25C13,
                                result_n26C13, result_n27C13, result_n28C13, result_n29C13,
                                result_n30C13, result_n31C13, result_n32C13, result_n33C13,
                                result_n34C13, result_n35C13, result_n36C13, result_n37C13,
                                result_n38C13, result_n39C13, result_n40C13)
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
        result_p12C13 <- msMatch(mzList, reference[reference$Annotation == "13C_12", ])
        result_p13C13 <- msMatch(mzList, reference[reference$Annotation == "13C_13", ])
        result_p14C13 <- msMatch(mzList, reference[reference$Annotation == "13C_14", ])
        result_p15C13 <- msMatch(mzList, reference[reference$Annotation == "13C_15", ])
        result_p16C13 <- msMatch(mzList, reference[reference$Annotation == "13C_16", ])
        result_p17C13 <- msMatch(mzList, reference[reference$Annotation == "13C_17", ])
        result_p18C13 <- msMatch(mzList, reference[reference$Annotation == "13C_18", ])
        result_p19C13 <- msMatch(mzList, reference[reference$Annotation == "13C_19", ])
        result_p20C13 <- msMatch(mzList, reference[reference$Annotation == "13C_20", ])
        result_p21C13 <- msMatch(mzList, reference[reference$Annotation == "13C_21", ])
        result_p22C13 <- msMatch(mzList, reference[reference$Annotation == "13C_22", ])
        result_p23C13 <- msMatch(mzList, reference[reference$Annotation == "13C_23", ])
        result_p24C13 <- msMatch(mzList, reference[reference$Annotation == "13C_24", ])
        result_p25C13 <- msMatch(mzList, reference[reference$Annotation == "13C_25", ])
        result_p26C13 <- msMatch(mzList, reference[reference$Annotation == "13C_26", ])
        result_p27C13 <- msMatch(mzList, reference[reference$Annotation == "13C_27", ])
        result_p28C13 <- msMatch(mzList, reference[reference$Annotation == "13C_28", ])
        result_p29C13 <- msMatch(mzList, reference[reference$Annotation == "13C_29", ])
        result_p30C13 <- msMatch(mzList, reference[reference$Annotation == "13C_30", ])
        result_p31C13 <- msMatch(mzList, reference[reference$Annotation == "13C_31", ])
        result_p32C13 <- msMatch(mzList, reference[reference$Annotation == "13C_32", ])
        result_p33C13 <- msMatch(mzList, reference[reference$Annotation == "13C_33", ])
        result_p34C13 <- msMatch(mzList, reference[reference$Annotation == "13C_34", ])
        result_p35C13 <- msMatch(mzList, reference[reference$Annotation == "13C_35", ])
        result_p36C13 <- msMatch(mzList, reference[reference$Annotation == "13C_36", ])
        result_p37C13 <- msMatch(mzList, reference[reference$Annotation == "13C_37", ])
        result_p38C13 <- msMatch(mzList, reference[reference$Annotation == "13C_38", ])
        result_p39C13 <- msMatch(mzList, reference[reference$Annotation == "13C_39", ])
        result_p40C13 <- msMatch(mzList, reference[reference$Annotation == "13C_40", ])
        
        #Wrap all the results as one output dataframe
        result_all_pos <- rbind(result_pC12, result_p1C13, result_p2C13, result_p3C13, result_p4C13,
                                result_p5C13, result_p6C13, result_p7C13, result_p8C13, result_p9C13,
                                result_p10C13, result_p11C13, result_p12C13, result_p13C13,
                                result_p14C13, result_p15C13, result_p16C13, result_p17C13,
                                result_p18C13, result_p19C13, result_p20C13, result_p21C13,
                                result_p22C13, result_p23C13, result_p24C13, result_p25C13,
                                result_p26C13, result_p27C13, result_p28C13, result_p29C13,
                                result_p30C13, result_p31C13, result_p32C13, result_p33C13,
                                result_p34C13, result_p35C13, result_p36C13, result_p37C13,
                                result_p38C13, result_p39C13, result_p40C13)
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























