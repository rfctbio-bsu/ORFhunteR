Vectorising_RNAs_V = function(x) {

RNAs <- x

stats <- data.table(index = (1:length(RNAs)), ID = "",
                    #monograms
                    freqA = 0, freqT = 0, freqG = 0, freqC = 0,
                    #bigrams
                    freqAA = 0, freqAT = 0, freqAG = 0, freqAC = 0,
                    freqTA = 0, freqTT = 0, freqTG = 0, freqTC = 0,
                    freqGA = 0, freqGT = 0, freqGG = 0, freqGC = 0,
                    freqCA = 0, freqCT = 0, freqCG = 0, freqCC = 0,
                    #trigrams
                    freqAAA = 0, freqAAT = 0, freqAAG = 0, freqAAC = 0,
                    freqATA = 0, freqATT = 0, freqATG = 0, freqATC = 0,
                    freqAGA = 0, freqAGT = 0, freqAGG = 0, freqAGC = 0,
                    freqACA = 0, freqACT = 0, freqACG = 0, freqACC = 0,
                    freqTAA = 0, freqTAT = 0, freqTAG = 0, freqTAC = 0,
                    freqTTA = 0, freqTTT = 0, freqTTG = 0, freqTTC = 0,
                    freqTGA = 0, freqTGT = 0, freqTGG = 0, freqTGC = 0,
                    freqTCA = 0, freqTCT = 0, freqTCG = 0, freqTCC = 0,
                    freqGAA = 0, freqGAT = 0, freqGAG = 0, freqGAC = 0,
                    freqGTA = 0, freqGTT = 0, freqGTG = 0, freqGTC = 0,
                    freqGGA = 0, freqGGT = 0, freqGGG = 0, freqGGC = 0,
                    freqGCA = 0, freqGCT = 0, freqGCG = 0, freqGCC = 0,
                    freqCAA = 0, freqCAT = 0, freqCAG = 0, freqCAC = 0,
                    freqCTA = 0, freqCTT = 0, freqCTG = 0, freqCTC = 0,
                    freqCGA = 0, freqCGT = 0, freqCGG = 0, freqCGC = 0,
                    freqCCA = 0, freqCCT = 0, freqCCG = 0, freqCCC = 0,
                    
                    #segmental probabilities of correlation factors
                    corrFactAT = 0, corrFactAG = 0, corrFactAC = 0, corrFactTG = 0, corrFactTC = 0, corrFactGC = 0,
                    #string lengthes
                    length = 0, logLength = 0,
                    
                    #Bao metrics
                    HRR = 0, HRY = 0, HYR = 0, HYY = 0,
                    HMM = 0, HMK = 0, HKM = 0, HKK = 0,
                    HWW = 0, HWS = 0, HSW = 0, HSS = 0)

#stats$ID = names(RNAs)
stats$ID = 1:length(RNAs)

monograms = c("A","T","G","C");


#stats[,3] <- lapply(RNAs,getPatternCount_V,monograms[1])
# frequencies of single nucleotides (bigrams)
for (i in 1 : 4) stats[,i+2] <- sapply(RNAs,getPatternCount_V,monograms[i])/
  sapply(RNAs,str_length)

# frequencies of nucleotid pairs (bigrams)
bigrams = c("AA","AT","AG","AC","TA","TT","TG","TC","GA","GT","GG","GC","CA",
            "CT","CG","CC");
for (i in 1 : 16) stats[,i+6] <- sapply(RNAs,getPatternCount_V,bigrams[i])/
  (sapply(RNAs,str_length)-1)

#frequencies of nucleotide triplets (trigrams)
trigrams = c("AAA","AAT","AAG","AAC","ATA","ATT","ATG","ATC","AGA","AGT","AGG",
             "AGC","ACA","ACT","ACG","ACC",
             "TAA","TAT","TAG","TAC","TTA","TTT","TTG","TTC","TGA","TGT","TGG",
             "TGC","TCA","TCT","TCG","TCC",
             "GAA","GAT","GAG","GAC","GTA","GTT","GTG","GTC","GGA","GGT","GGG",
             "GGC","GCA","GCT","GCG","GCC",
             "CAA","CAT","CAG","CAC","CTA","CTT","CTG","CTC","CGA","CGT","CGG",
             "CGC","CCA","CCT","CCG","CCC")
for (i in 1 : 64) stats[,i+22] <- sapply(RNAs,getPatternCount_V,trigrams[i])/
  (sapply(RNAs,str_length)-2)

#correlation factors
corrFactors = c("AT","AG","AC","TG","TC","GC")
for (i in 1 : 6) stats[,i+86] <- sapply(RNAs,getCorrelationFactors,corrFactors[i],10)

#lengthes
stats$length <- sapply(RNAs,str_length)
stats$logLength <- log10(sapply(RNAs,str_length))

#Bao metrics
testRowsNum = length(RNAs)
bao <- sapply(RNAs,getBaoMetrics,testRowsNum,testRowsNum)
bao <- data.frame(t(bao))

rm(RNAs)

# add Bao features to existing features set
stats$HRR <- as.numeric(bao$HRR)
stats$HRY <- as.numeric(bao$HRY)
stats$HYR <- as.numeric(bao$HYR)
stats$HYY <- as.numeric(bao$HYY)
stats$HMM <- as.numeric(bao$HMM)
stats$HMK <- as.numeric(bao$HMK)
stats$HKM <- as.numeric(bao$HKM)
stats$HKK <- as.numeric(bao$HKK)
stats$HWW <- as.numeric(bao$HWW)
stats$HWS <- as.numeric(bao$HWS)
stats$HSW <- as.numeric(bao$HSW)
stats$HSS <- as.numeric(bao$HSS)

#  row.names(stats) <- stats$ID
#  stats$index <- NULL
#  stats$ID <- NULL

RNAs_vctd <- data.frame(stats)
row.names(RNAs_vctd) <- RNAs_vctd$ID
RNAs_vctd$index <- NULL
RNAs_vctd$ID <- NULL

#t2 <- Sys.time( )
#print(t2-t1)

#write.table(RNAs_vctd,'lncRNAs_vctd.txt')
return(RNAs_vctd)
}
