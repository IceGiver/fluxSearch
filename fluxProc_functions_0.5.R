require(fuzzyjoin)
require(data.table)
require(dplyr)



#====================================================================
# Isotopomer natural abundance correction

IsoCorr <- function(fluxResultFile) {
        fluxResult <- fread(fluxResultFile)
        # Missing value debug using
        fluxResult <- fluxResult[(fluxResult$'lipidForm' != 'NA'),]
        
        # natC13 <- 0.01109
        # Extract the total number of each lipid as well as the 13C labeling number
        
        fluxResult$'C_num' <- as.numeric(sub('C([0-9]+).*', '\\1', fluxResult$'lipidForm'))
        fluxResult$'lab_num' <- NA
        fluxResult$'lab_num' <- as.numeric(fluxResult$'lab_num')
        fluxResult$'lab_num'[fluxResult$Annotation == "12C"] <- 0
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_01"] <- 1
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_02"] <- 2
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_03"] <- 3
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_04"] <- 4
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_05"] <- 5
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_06"] <- 6
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_07"] <- 7
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_08"] <- 8
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_09"] <- 9
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_10"] <- 10
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_11"] <- 11
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_12"] <- 12
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_13"] <- 13
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_14"] <- 14
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_15"] <- 15
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_16"] <- 16
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_17"] <- 17
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_18"] <- 18
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_19"] <- 19
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_20"] <- 20
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_21"] <- 21
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_22"] <- 22
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_23"] <- 23
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_24"] <- 24
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_25"] <- 25
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_26"] <- 26
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_27"] <- 27
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_28"] <- 28
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_29"] <- 29
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_30"] <- 30
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_31"] <- 31
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_32"] <- 32
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_33"] <- 33
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_34"] <- 34
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_35"] <- 35
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_36"] <- 36
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_37"] <- 37
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_38"] <- 38
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_39"] <- 39
        fluxResult$'lab_num'[fluxResult$Annotation == "13C_40"] <- 40
        
        # Extract a list with only 12C unlabeled features and rename the 'Intensity' col
        fluxRes12C <- fluxResult[(fluxResult$Annotation == "12C"),]
        
        names(fluxRes12C)[names(fluxRes12C) == "Intensity"] <- "Int_12C"
        fluxRes12C_sub <- subset(fluxRes12C, select = c("Int_12C", "Lipid",
                                                        "Theoretical_mz", "Charge"))
        
        # Extract a list without any 12C unlabeled features and rename the 'Intensity' col
        fluxRes13C <- fluxResult
        names(fluxRes13C)[names(fluxRes13C) == "Intensity"] <- "Int_13C"
        
        # Merge the labeled and unlabeled features in a parallel wide list
        fluxRes12C13C <- inner_join(fluxRes13C, fluxRes12C_sub)
        
        # Check all the vector type
        fluxRes12C13C$Int_13C <- as.numeric(fluxRes12C13C$Int_13C)
        fluxRes12C13C$Int_12C <- as.numeric(fluxRes12C13C$Int_12C)
        
        # Using the nature abundance correction function to calculate the isotopomer nature abundances
        fluxRes12C13C$Int_13C_corr <- NA
        fluxRes12C13C$Int_13C_corr <- 
                IsoNat(fluxRes12C13C$Int_13C,
                       fluxRes12C13C$Int_12C,
                       fluxRes12C13C$C_num,
                       fluxRes12C13C$lab_num)
        # Remove the corrected intensities with negative values
        fluxRes12C13C <- fluxRes12C13C[fluxRes12C13C$Int_13C_corr > 0, ]

        return(fluxRes12C13C)
}


#================================================================================
# The function using to calculate the 13C nature isotopic abundances
# 'Moseley BMC Bioinformatics 2010, 11:139'

IsoNat <- function(Int_13C, Int_12C, C_num, lab_num) {

        # The natural abundant of 13C ~ 1.109%
        natC13 <- 0.01109

        # Binomial distribution model of 13C isotopomer distributions
        
        Int_13C_corr <-
                Int_13C - Int_12C*(choose(C_num, lab_num))*(natC13^lab_num)*((1-natC13)^(C_num - lab_num))

}

#================================================================================
# The function using to calculate the MDV (mass distribution vactor) of the lipid
IsoLipMDV <- function (resultFile) {
        # Process the isotopic pattern correction first
        df <- IsoCorr(resultFile)
        # Remove all the NA rows
        df_comp <- df[complete.cases(df),]
        # Remove all the replicates
        df_uni <- df_comp[!duplicated(df_comp[order(-df_comp$Annotation),]),]
        # Calculate the mass distribution vactor (MDV)
        df_MDV <- df_uni[, MVD := Int_13C_corr/sum(Int_13C_corr), by = Lipid]
        
        return(df_MDV)
        
} 


#================================================================================
# The function using to calculate the FC (fractional contribution) of the lipid
IsoLipFC <- function (resultFile) {
        # Process the isotopic pattern correction first
        df <- IsoCorr(resultFile)
        # Remove all the NA rows
        df_comp <- df[complete.cases(df),]
        # Remove all the replicates
        df_uni <- df_comp[!duplicated(df_comp[order(-df_comp$Annotation),]),]
        # Calculate the mass distribution vactor (MDV)
        df_MDV <- df_uni[, MVD := Int_13C_corr/sum(Int_13C_corr), by = Lipid]
        # Calculate the fractional contribution (frac/FC)
        df_frac <- df_uni[, Frac := sum(lab_num*MVD)/C_num, by = Lipid]
        
        return(df_frac)
}













