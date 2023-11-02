library("MSGFplus")
library("MSnID")
library(xlsx)

#BiocManager::install("Rcpp", version = '1.0')
#remove.packages("Rcpp")

#BiocManager::install("MSnID", force=TRUE)

#rcpp <- "http://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.6.tar.gz"
#install.packages(rcpp, repos=NULL, type="source")
#devtools::install_github("RcppCore/Rcpp")
#install.packages("Rcpp", repos="https://RcppCore.github.io/drat")


setwd('../data/spectra/roots/MS_MS/MGF/')


par <- msgfPar()

db(par) <- '../data/IPG_db/IPG_rice_chunk0.fasta'

tolerance(par) <- '0.2 Da'      # Set parent ion tolerance
chargeRange(par) <- c(1, 6)     # Set the range of charge states to look after
lengthRange(par) <- c(4, 30)    # Set the range of peptide length to look after
instrument(par) <- 'TOF'  # Set the instrument used for acquisition
enzyme(par) <- 'Trypsin'        # Set the enzyme used for digestion
fragmentation(par) <- 0         # Set the fragmentation method
protocol(par) <- 0              # Set the protocol type
isotopeError(par) <- c(0,0)     # Set the isotope error
matches(par) <- 2               # Set the number of matches to report per scan
ntt(par) <- 1                   # Set number of tolerable termini
tda(par) <- TRUE

mods(par)[[1]] <- msgfParModification(name = 'Carbamidomethyl',
                                      composition='C2H3N1O1',
                                      residues = 'C', 
                                      type = 'opt', 
                                      position = 'any')
mods(par)[[2]] <- msgfParModification(name = 'Oxidation', 
                                      mass=15.9949,
                                      residues = 'M', 
                                      type = 'opt', 
                                      position = 'any')
nMod(par) <- 2                  # Set max number of modifications per peptide



#msnid <- MSnID(".")
#msnid <- read_mzIDs(msnid,'C:/Users/HP/Desktop/PPB/2020_oxy_proteome/MS/rice/samples/7_1.mzid')
#msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
#msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")

#msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
#msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
#filtObj <- MSnIDFilter(msnid)
#filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
#filtObj$msmsScore <- list(comparison=">", threshold=8.0)
#evaluate_filter(msnid, filtObj, level="PSM")
#evaluate_filter(msnid, filtObj, level="peptide")
#evaluate_filter(msnid, filtObj, level="accession")

#filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
#                                method="Grid", level="peptide",
#                                n.iter=500)

#filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01,
#                              method="Nelder-Mead", level="peptide",
#                              n.iter=500)

#evaluate_filter(msnid, filtObj, level="peptide")
#evaluate_filter(msnid, filtObj.grid, level="peptide")

#evaluate_filter(msnid, filtObj.nm, level="accession")

#msnid <- apply_filter(msnid, filtObj.nm)
#msnid <- apply_filter(msnid, "isDecoy == FALSE")
#msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
#psm.df <- psms(msnid)


dir('C:/Users/HP/Desktop/PPB/2021_proteome/rice_MS_MS/iter_2')
setwd('C:/Users/HP/Desktop/PPB/2021_proteome/rice_MS_MS/iter_2')


for (file in dir()){
  runMSGF(par, file)
}

#library("MSnID")

df_tot=data.frame()


for (file in dir()){
  if (grepl( "mzid", file, fixed = TRUE)){
    msnid <- MSnID(".")
    msnid <- read_mzIDs(msnid,file)
    msnid <- assess_termini(msnid, validCleavagePattern="[KR]\\.[^P]")
    msnid <- assess_missed_cleavages(msnid, missedCleavagePattern="[KR](?=[^P$])")
    
    msnid$msmsScore <- -log10(msnid$`MS-GF:SpecEValue`)
    msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
    filtObj <- MSnIDFilter(msnid)
    filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=5.0)
    filtObj$msmsScore <- list(comparison=">", threshold=8.0)
    evaluate_filter(msnid, filtObj, level="PSM")
    evaluate_filter(msnid, filtObj, level="peptide")
    evaluate_filter(msnid, filtObj, level="accession")
    
    filtObj.grid <- optimize_filter(filtObj, msnid, fdr.max=0.01,
                                    method="Grid", level="peptide",
                                    n.iter=500)
    
    filtObj.nm <- optimize_filter(filtObj.grid, msnid, fdr.max=0.01,
                                  method="Nelder-Mead", level="peptide",
                                  n.iter=500)
    
    msnid <- apply_filter(msnid, filtObj.nm)
    msnid <- apply_filter(msnid, "isDecoy == FALSE")
    msnid <- apply_filter(msnid, "!grepl('Contaminant',accession)")
    psm.df <- psms(msnid)
    df_tot <- rbind(df_tot,psm.df[1,])
  }
}



df_tot$splot=sapply(strsplit(df_tot[,2], "_"), head, 1)

df_not_na <- df_tot[!is.na(df_tot$spectrumID),]
df_exp <- df_not_na[,c(36,22,24,34,31)]
df_exp$splot <- as.numeric(df_exp$splot )
write.xlsx(df_exp, 'MS_MS_MSGF_annotation_iter_2.xlsx', sheetName = "Sheet1", 
           col.names = TRUE, row.names = FALSE, append = FALSE)
