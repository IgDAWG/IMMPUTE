SCOREPGM <- function(fileNAME,AncLevel,DataDir,GlobalDir,Mask=F,Bar="Down",Copy=FALSE,Res=2) {
  
  require(gplots)
  
  # Script is limited by locus: A,B,C,DRB1,DQB1
  
  # Argument Definitions and acceptable values:
  # fileName - string, file name of text data file, UTF-8, tab delimited.
  # AncLevel - string, desired ancestry level to analyze: Pop, Geo, AL1, AL2.
  # DataDir - string, path for text data file.
  # GlobalDir - string, path for accessory files.
  # Mask - logical, should subjects with untrained alleles be masked? Default = F.
  # Bar - string, position of color bar on plots: Up or Down. Default="Down".
  # Copy - logical, should same genotype level probabilities be used for each allele?
  # Res - Numeric, Resolution Desired: 1 or 2. Default=2.
  #
  # 1) This script reads in accessory files that must be contained within a defined common Global directory: GlobalDir
  #     - scoringSource_Functions.R
  #    Empty versions of following text files provided on GitHub
  #     - Known_HLA_Genotypes.txt
  #       header: Known.id,Locus,Allele 1,Allele 2,Allele 1 2D,Allele 2 2D,Allele 1 4D,Allele 2 4D,Allele1G,Allele2G,Allele1P,Allele2P
  #     - HLA_G_Groups.txt & HLA_P_Groups.txt (Database Version: 2013-04-17 at time of script construction)
  #       header: Known.id,Locus,Allele
  #     - Known_Demographics.txt
  #       header: Known.id,Population,Geography,AL1,AL2
  #     - Known_Homozygotes_1F.txt (1-Field) and Known_Homozygotes_2F.txt (2-Field)
  #       header: Known.id,Locus,Allele
  #     - Known_MASK.txt
  #       header: sampleid,A,B,C,DRB1,DQB1
  #       Mask is a file of 1's and/or NA's whether or not you want a specific subject to be masked for that locus.
  #       1 = No Masking, NA = Masking
  # 2) Header definitions:
  #       known.id - Sample ID for known genotypes
  #       Locus - HLA Locus
  #       Allele - Full field allele name
  #       Allele 1/Allele 2 - Full field name for 1st and 2nd allele of a given locus genotype
  #       Allele 1 2D/Allele 2 2D - 2 digit (1-Field) allele name
  #       Allele 1 4D/Allele 2 4D - 4 digit (2-Field) allele name
  #       Allele1G / Allele2G - 1-Field and 2-Field G group designations respectively
  #       Allele1P / Allele2P - 1-Field and 2-Field P group designations respectively
  #       Population, Geography, AL1 (Ancestry Level 1), AL2 (Ancestry Level 2) - user defined ancestry information
  # 3) Required package: gplots

#################################################################################### Data and Variable Definition ####
# Global Variables ----------------------------------------------------------
 
  # Reading in common functions and data
  # Builds: GTYPE, DEMO, Loci, KnownUNAll, GGROUPS 
  
  MainDir <- getwd()
  
  setwd(GlobalDir)
  source('scoringSOURCE_Functions.R')
  
  # Read in full Known genotype table with PG conversions
  GTYPE <- read.table(file = "Known_HLA_Genotypes.txt",sep="\t", header=T, stringsAsFactors=F)
  GTYPE <- GTYPE[which(GTYPE$Allele.1!="****"),]
    
  # Read in P/G group extended conversion tables
  # Database Version: 2013-04-17
  # origin: http://hla.alleles.org/wmda/hla_nom_g.txt
  # author: WHO, Steven G. E. Marsh (steven.marsh@ucl.ac.uk)
  if(Res==2) {
    GGROUPS <- read.table(file = "HLA_G_Groups.txt",sep="\t", header=T,fill=T)
      GGROUPS <- GGROUPS[grep("G",GGROUPS[,5],fixed=T),]
    PGROUPS <- read.table(file = "HLA_P_Groups.txt",sep="\t", header=T,fill=T)
      PGROUPS <- PGROUPS[grep("P",PGROUPS[,5],fixed=T),] }
    
  # Demographic File corresponding to Known
  DEMO <- read.table("Known_Demographics.txt",sep="\t", header=T,colClasses="character")
  
  # Loci List
  Loci <- c("A","B","C","DRB1","DQB1")
  
  # List of Unique Known ID/Locus combinations without full training set overlap
  if(Res==2) { MASK <- read.table("Known_MASK.txt",sep="\t", header=T) }

  # Protected objects list for clean-up
  SAFE <- c("SAFE","Loci","AncLevel","DEMO","DataDir","GlobalDir", "GTYPE","IMP","IMPsub","SCORES","MATCH",
            "fileNAME","is.between","percent","FixField","colorBAR","GGROUPS","PGROUPS","KnownUNAll","Mask",
            "MASK","Bar","CheckMatch","CheckMatchv2","Copy","CopyProb","Zipper","HLA.Mirror","Res","MainDir")


# DATA READ IN ----------------------------------------------------------

  # Imputation Analysis and Scoring
  # Builds: IMP, IMPsub, DEMO
  
  setwd(DataDir)
  
  # Read in imputed data ("IMP") and keep only allele call columns ("IMPsub")
  colHeaders <- c("Known.id","A.1","A.2","A.Prob.1","A.Prob.2","B.1","B.2","B.Prob.1","B.Prob.2","C.1","C.2","C.Prob.1","C.Prob.2","DRB1.1","DRB1.2","DRB1.Prob.1","DRB1.Prob.2","DQB1.1","DQB1.2","DQB1.Prob.1","DQB1.Prob.2")
  IMP <- read.table(file = fileNAME, sep="\t", header=T)
  if(Copy) {
    IMP <- CopyProb(IMP)
    colnames(IMP) <- colHeaders
  } else { colnames(IMP) <- colHeaders }
  IMPsub <- IMP[,c(1:3,6:7,10:11,14:15,18:19)]

  # Ensure Mask Sorts the same as IMPsub
  if(Res==2) {
    MASK <- MASK[match(IMPsub[,'Known.id'],MASK[,'Known.id']),]
    rownames(MASK) <- NULL }

  #Fix 1 field data to ensure HLA alleles contain leading zeros
  if(Res==1) {
    for(i in 2:ncol(IMPsub)) {
      IMPsub[,i] <- sapply(IMPsub[,i],FixField)
    }; rm(i) }

  # Build DEMO file specific for the indicated run
  DEMO <- DEMO[match(IMPsub$Known.id,DEMO$Known.id),]
  
  # Define Output Directories Based on Resolution
  if(Res==1) { OutDir <- paste(DataDir,"/Score_1F",sep="") }
  if (Res==2) {
    if(Mask) {
      OutDir <- paste(DataDir,"/Score_2F_FullMask",sep="")
    } else {
      OutDir <- paste(DataDir,"/Score_2F_NoMask",sep="")
    }     
  }
  dir.create(OutDir)
  setwd(OutDir)
  
  #Write Variables to File for tracking
  if(Res==1) { Vars <- c(fileNAME,AncLevel,Bar,Copy,Res) }
  if(Res==2) { Vars <- c(fileNAME,AncLevel,Mask,Bar,Copy,Res) }
  write.table(Vars, "variables.txt", sep=",",col.names=F,row.names=F,quote=F)
  
  rm(list=ls()[!(ls() %in% SAFE)])


#################################################################################### Begin Analysis Script ####
# 1. BUILD MATCH AND SCORE MATRICES ----------------------------------------------------------

  # Current scores for A,B,C,DRB1. Must updated k range to reflect inclusion of DQB1 when available
  # Requires: Loci
  # Builds: MATCH, SCORES
  # Writes: Match_out.txt, Scores_out.txt
  
  MATCH <- data.frame(SampleID=character(),
                      A1=character(),
                      A2=character(), 
                      B1=character(),
                      B2=character(),
                      C1=character(),
                      C2=character(),
                      DRB1.1=character(),
                      DRB1.2=character(),
                      DQB1.1=character(),
                      DQB1.2=character(),
                      stringsAsFactors=FALSE)
  
  # Loop through Sample IDs (SID)
  for (i in 1:nrow(IMPsub)) {
    
    SID <- as.numeric(IMPsub[i,1])
    if(Res==1) { Tab <- GTYPE[which(GTYPE$Known.id == SID),c('Known.id','Locus','Allele.1.2D','Allele.2.2D')] }
    if(Res==2) {
      Tab <- GTYPE[which(GTYPE$Known.id == SID),c('Known.id','Locus','Allele.1.4D','Allele.2.4D')]; rownames(Tab) <- NULL
      TabPG <- GTYPE[which(GTYPE$Known.id == SID),c('Known.id','Locus','Allele1G','Allele2G','Allele1P','Allele2P')]; rownames(TabPG) <- NULL
    }
    
    # Loop through Loci
    for (k in 1:length(Loci)){
      
      #Locus
      Locus <- Loci[k]
      
      First <- k*2
      Second <- First + 1
      
      MATCH[i,1] <- SID
      
      if(Res==1) { MaxLen <- as.numeric(length(which(Tab$Locus == Locus))) }
      if(Res==2) { MaxLen <- as.numeric(max(length(which(Tab$Locus==Locus)),length(which(TabPG$Locus==Locus)))) }
      
      if(MaxLen==0){
        
        MATCH[i,First] <- NA
        MATCH[i,Second] <- NA
        
      } else {
        
        if(Res==1){
          Calls <- unique(Tab[which(Tab$Locus==Locus),c('Allele.1.2D','Allele.2.2D')]); rownames(Calls) <- NULL
          
          # All Calls
          Calls.grp <- unique(as.character(unlist(sapply(as.character(unlist(Calls)),strsplit,split="/"))))
          
          # All Zipped Calls
          Calls.fix <- Zipper(Calls)
          
          #IMP Calls
          IMPsub.Calls <-IMPsub[i,First:Second]
          IMPsub.Calls.grp1 <- as.character(unlist(sapply(as.character(unlist(IMPsub.Calls[,1])),strsplit,split="/")))
          IMPsub.Calls.grp2 <- as.character(unlist(sapply(as.character(unlist(IMPsub.Calls[,2])),strsplit,split="/")))
          IMPsub.Calls.fix <- Zipper(IMPsub.Calls)
          
          #Add in reverse ordered calls to control for sorting
          IMPsub.Calls.fix <- unique(rbind(IMPsub.Calls.fix,
                                           HLA.Mirror(IMPsub.Calls.fix)))
        }
        
        
        if(Res==2){
          SubTab <- unique(Tab[which(Tab$Locus==Locus),]); rownames(SubTab) <- NULL
          SubTabPG <- unique(TabPG[which(TabPG$Locus==Locus),]); rownames(SubTabPG) <- NULL
          
          Calls <- list(C=SubTab[,c('Allele.1.4D','Allele.2.4D')],
                        G=SubTabPG[,c('Allele1G','Allele2G')],
                        P=SubTabPG[,c('Allele1P','Allele2P')])
          
          colnames(Calls$C) <- c("Allele.1","Allele.2")
          colnames(Calls$G) <- c("Allele.1","Allele.2")
          colnames(Calls$P) <- c("Allele.1","Allele.2")
          Calls <- unique(do.call(rbind,Calls)); rownames(Calls) <- NULL
          
          # All Calls
          Calls.grp <- unique(as.character(unlist(sapply(as.character(unlist(Calls)),strsplit,split="/"))))
          Calls.grp <- unique(c(Calls.grp,sapply(Calls.grp,gsub,pattern='[[:alpha:]]',replacement="")))
          
          # All Zipped Calls
          Calls.fix <- Zipper(Calls)
          Calls.fix.noalpha <- matrix(sapply(Calls.fix,gsub,pattern='[[:alpha:]]',replacement=""),ncol=1)
          colnames(Calls.fix.noalpha) <- colnames(Calls.fix)
          Calls.fix <- unique(rbind(Calls.fix,Calls.fix.noalpha))
          
          #IMP Calls
          IMPsub.Calls <-IMPsub[i,First:Second]
          
          #IMP G Conversion
            GGROUPS.sub <- GGROUPS[which(GGROUPS$Locus==Locus),]
            #A1
            A1.g <-  unique(as.character(GGROUPS.sub[grep(IMPsub.Calls[,1],GGROUPS.sub[,'Allele.2F'],fixed=T),5]))[1]
            if(is.na(A1.g)) { A1.g <- as.character(IMPsub.Calls[,1]) }
          
            #A2
            A2.g <-  unique(as.character(GGROUPS.sub[grep(IMPsub.Calls[,2],GGROUPS.sub[,'Allele.2F'],fixed=T),5]))[1]
            if(is.na(A2.g)) { A2.g <- as.character(IMPsub.Calls[,2]) }
          
            #Call
            tmp.g <- unique(cbind(A1.g,A2.g))
            colnames(tmp.g) <- colnames(IMPsub.Calls)
          
          #IMP P Conversion
            PGROUPS.sub <- PGROUPS[which(PGROUPS$Locus==Locus),]
            #A1
            A1.p <-  unique(as.character(PGROUPS.sub[grep(IMPsub.Calls[,1],PGROUPS.sub[,4],fixed=T),5]))[1]
            if(is.na(A1.p)) { A1.p <- as.character(IMPsub.Calls[,1]) }
  
            #A2
            A2.p <- unique(as.character(PGROUPS.sub[grep(IMPsub.Calls[,2],PGROUPS.sub[,4],fixed=T),5]))[1]
            if(is.na(A2.p)) { A2.p <- as.character(IMPsub.Calls[,2]) }
          
            #Call
            tmp.p <- unique(cbind(A1.p,A2.p))
            colnames(tmp.p) <- colnames(IMPsub.Calls)
          
          #IMP Calls Combined
          IMPsub.Calls.pg <- unique(rbind(IMPsub.Calls,tmp.g,tmp.p))
          IMPsub.Calls.grp1 <- as.character(unlist(sapply(as.character(unlist(IMPsub.Calls.pg[,1])),strsplit,split="/")))
          IMPsub.Calls.grp2 <- as.character(unlist(sapply(as.character(unlist(IMPsub.Calls.pg[,2])),strsplit,split="/")))
          IMPsub.Calls.fix <- Zipper(IMPsub.Calls.pg)
          
          #Add in reverse ordered calls to control for sorting
          IMPsub.Calls.fix <- unique(rbind(IMPsub.Calls.fix,HLA.Mirror(IMPsub.Calls.fix)))
        }
        
        #Final Scoring
        
          # Conservative Matching for genotype exclusivity
          if(sum(IMPsub.Calls.fix %in% Calls.fix)>0) {
            
            MATCH[i,First] <- 1
            MATCH[i,Second] <- 1
            
          } else {
          
          # Partial Matching for single match
          if(CheckMatchv2(IMPsub.Calls.grp1,Calls.grp)==1) {
            MATCH[i,First] <- 1
            MATCH[i,Second] <- 0
          } else if (CheckMatchv2(IMPsub.Calls.grp2,Calls.grp)==1) {
            MATCH[i,First] <- 0
            MATCH[i,Second] <- 1
          } else {
            MATCH[i,First] <- 0
            MATCH[i,Second] <- 0
          }
        
        }

      }
      
    }; rm(k)
    
  }; rm(i)
  
  SCORES <- list(A = rowSums(data.matrix(MATCH[,2:3])),
                 B = rowSums(data.matrix(MATCH[,4:5])),
                 C = rowSums(data.matrix(MATCH[,6:7])),
                 DRB1 = rowSums(data.matrix(MATCH[,8:9])),
                 DQB1 = rowSums(data.matrix(MATCH[,10:11])))

  SCORES <- do.call(cbind,SCORES)
  SCORES <- cbind(IMP$Known.id,SCORES)
  
  if (Mask) {
    MATCH <- cbind(MATCH[,1],
                   as.numeric(MATCH[,2])*MASK[,2],
                   as.numeric(MATCH[,3])*MASK[,2],
                   as.numeric(MATCH[,4])*MASK[,3],
                   as.numeric(MATCH[,5])*MASK[,3],
                   as.numeric(MATCH[,6])*MASK[,4],
                   as.numeric(MATCH[,7])*MASK[,4],
                   as.numeric(MATCH[,8])*MASK[,5],
                   as.numeric(MATCH[,9])*MASK[,5],
                   as.numeric(MATCH[,10])*MASK[,6],
                   as.numeric(MATCH[,11])*MASK[,6])
    MATCH <- matrix(as.numeric(MATCH),ncol=ncol(MATCH))
    colnames(MATCH) <- c("SampleID","A1","A2","B1","B2","C1","C2","DRB1.1","DRB1.2","DQB1.1","DQB1.2")
    SCORES <- cbind(SCORES[,1], SCORES[,2:6]*MASK[,2:6])
  }
  
  write.table(MATCH,file="MATCH_out.txt",sep="\t",quote=F,col.names=T,row.names=F)
  write.table(SCORES,file="SCORES_out.txt",sep="\t",quote=F,col.names=T,row.names=F)
  
  rm(list=ls()[!(ls() %in% SAFE)])


# 2. SCORE BY COUNT ----------------------------------------------------------
  
  # Requires: Loci, SCORES
  # Writes: Scores_byCount_output.txt
  
  Output <- data.frame(Locus=character(),
                       noIND=character(),
                       Match.Overall=character(),
                       Match.0=character(),
                       Match.1=character(),
                       Match.2=character(),
                       stringsAsFactors=FALSE)
  Outputp <- Output
  
  for(i in 1:length(Loci)){
    
    Locus <- Loci[i]
    noIND <- as.numeric(sum(SCORES[,Locus]>=0, na.rm=TRUE))  #Subjects with available genotype
    noAll <- as.numeric(noIND*2) #Total Alleles from Subjects
    
    #Correct overall allele imputations
    List.All <- sum(na.omit(SCORES[,i+1]))
    List.All.pc <-  percent(List.All/(noAll),2)
    
    #Correct Across Individual Calls
    #No Correct calls
    List.0MATCH <- length(which(SCORES[,i+1]==0))
    List.0MATCH.pc <- percent(List.0MATCH/(noIND),2)
    
    #At least 1 correct calls
    List.1MATCH <- length(which(SCORES[,i+1]==1))
    List.1MATCH.pc <- percent(List.1MATCH/(noIND),2)
    
    #At least 2 correct calls
    List.2MATCH <- length(which(SCORES[,i+1]==2))
    List.2MATCH.pc <- percent(List.2MATCH/(noIND),2)
    
    Output[i,] <- c(Locus,
                    noIND,
                    List.All,
                    List.0MATCH,
                    List.1MATCH,
                    List.2MATCH)
    
    Outputp[i,] <- c(Locus,
                     noIND,
                     List.All.pc,
                     List.0MATCH.pc,
                     List.1MATCH.pc,
                     List.2MATCH.pc)
    
  }
  
  SEP <- rbind(Output,rep("-",ncol(Output)),Outputp)
  write.table(SEP,file="SCORES_byCount_output.txt",sep="\t",quote=F,row.names=F,col.names=T,append=F)
  
  # Subjects with correct genotypes at all loci
  FullGTYPE <- rbind(rep("-",2),
                     c("Correctly Imputed Subjects:",
                       length(which(rowSums(SCORES[,Loci])==length(Loci)*2))),
                     c("Correctly Imputed Subjects (%):",
                       percent(length(which(rowSums(SCORES[,Loci])==length(Loci)*2))/nrow(na.omit(SCORES[,Loci])),2)),
                     c("Subjects with full genotype:",
                       nrow(na.omit(SCORES[,Loci]))))
  write.table(FullGTYPE,file="SCORES_byCount_output.txt",sep="\t",quote=F,row.names=F,col.names=F,append=T)
  
  # Subjects with correct genotypes at all loci
  Loci.sub <- Loci[1:4]
  FullGTYPE.4L <- rbind(rep("-",2),
                        c("Correctly Imputed Subjects (-DQB1):",
                          length(which(rowSums(SCORES[,Loci.sub])==length(Loci.sub)*2))),
                        c("Correctly Imputed Subjects (%) (-DQB1):",
                          percent(length(which(rowSums(SCORES[,Loci.sub])==length(Loci.sub)*2))/nrow(na.omit(SCORES[,Loci.sub])),2)),
                        c("Subjects with full genotype (-DQB1):",
                          nrow(na.omit(SCORES[,Loci.sub]))))
  write.table(FullGTYPE.4L,file="SCORES_byCount_output.txt",sep="\t",quote=F,row.names=F,col.names=F,append=T)

  
  rm(list=ls()[!(ls() %in% SAFE)])
  


# 3. SCORE BY PROBABILITIES ----------------------------------------------------------

  # Requires: Loci, IMP, SCORES
  # Writes: Scores_byProb_output.txt
  
  Output.byProb <- data.frame(Locus=character(),
                              Bin=character(),
                              BinCount=character(),
                              Match.Overall=character(),
                              Match.0=character(),
                              Match.1=character(),
                              Match.2=character(),
                              stringsAsFactors=FALSE)
  
  for(i in 1:length(Loci)) {
    
    Locus <- Loci[i]
    
    #IMP Table extracting allele specific Prob columns
    Prob1 <- as.numeric(which(colnames(IMP)==paste(Locus,".Prob.1",sep=""))) #Allele 1 Probs
    Prob2 <- as.numeric(which(colnames(IMP)==paste(Locus,".Prob.2",sep=""))) #Allele 2 Probs
    
    IMPprob <- IMP[,Prob1]*IMP[,Prob2]
    
    Bin1 <- which(IMPprob<=0.25)
    Bin2 <- which(is.between(IMPprob,0.25,0.50))
    Bin3 <- which(is.between(IMPprob,0.50,0.75))
    Bin4 <- which(is.between(IMPprob,0.75,1))
    
    SCORESsub <- SCORES[,i+1]
    
    #Possible Imputations
    List.Poss.byProb <- list(length(na.omit(SCORESsub[Bin1])),
                             length(na.omit(SCORESsub[Bin2])),
                             length(na.omit(SCORESsub[Bin3])),
                             length(na.omit(SCORESsub[Bin4])))
    
    #Correct overall allele imputations
    List.All.byProb <- list(sum(na.omit(SCORESsub[Bin1])),
                            sum(na.omit(SCORESsub[Bin2])),
                            sum(na.omit(SCORESsub[Bin3])),
                            sum(na.omit(SCORESsub[Bin4])))
    
    #Per Individual Allele 0 Match
    List.0MATCH.byProb <- list(length(which(SCORESsub[Bin1]==0)),
                               length(which(SCORESsub[Bin2]==0)),
                               length(which(SCORESsub[Bin3]==0)),
                               length(which(SCORESsub[Bin4]==0)))
    
    #Per Individual Allele 1 Matches
    List.1MATCH.byProb <- list(length(which(SCORESsub[Bin1]==1)),
                               length(which(SCORESsub[Bin2]==1)),
                               length(which(SCORESsub[Bin3]==1)),
                               length(which(SCORESsub[Bin4]==1)))
    
    #Per Individual Allele 2 Matches
    List.2MATCH.byProb <- list(length(which(SCORESsub[Bin1]==2)),
                               length(which(SCORESsub[Bin2]==2)),
                               length(which(SCORESsub[Bin3]==2)),
                               length(which(SCORESsub[Bin4]==2)))
    
    Tab <- cbind(rep(Locus,4),
                 c(0.25,0.50,0.75,1),
                 List.Poss.byProb,
                 List.All.byProb,
                 List.0MATCH.byProb,
                 List.1MATCH.byProb,
                 List.2MATCH.byProb)
    colnames(Tab) <- colnames(Output.byProb)
    Output.byProb <- as.matrix(rbind(Output.byProb,Tab))
    
  }
  
  
  write.table(as.matrix(Output.byProb),file="SCORES_byProb_output.txt",sep="\t",quote=F,row.names=F,col.names=T,append=F)
  
  rm(list=ls()[!(ls() %in% SAFE)])  


# 4. SCORE BY DEMOGRAPHY ----------------------------------------------------------

  # Requires: DEMO, SCORES
  # Writes: SCORES_byDEMO_output.txt
  
  switch(AncLevel,
         Pop = assign("Col",2),
         Geo = assign("Col",3),
         AL1 = assign("Col",4),
         AL2 = assign("Col",5))
  
  AncUn <- sort(unique(DEMO[,Col]))
  
  Output.Demo <- data.frame(Ancestry=character(),
                            Count=character(),
                            Locus=character(),
                            Match.Overall=character(),
                            Match.0=character(),
                            Match.1=character(),
                            Match.2=character(),
                            stringsAsFactors=FALSE)
  
  for(k in 1:length(AncUn)) {
    
    Grp <- AncUn[k]
    Range <- which(DEMO[,Col]==Grp)
    
    #Correct overall allele imputations
    List.All.Demo <- list(sum(na.omit(SCORES[Range,2])),
                          sum(na.omit(SCORES[Range,3])),
                          sum(na.omit(SCORES[Range,4])),
                          sum(na.omit(SCORES[Range,5])),
                          sum(na.omit(SCORES[Range,6])))
    
    #Per Individual Allele 0 Match
    List.0MATCH.Demo <- list(length(which(SCORES[Range,2]==0)),
                             length(which(SCORES[Range,3]==0)),
                             length(which(SCORES[Range,4]==0)),
                             length(which(SCORES[Range,5]==0)),
                             length(which(SCORES[Range,6]==0)))
    
    #Per Individual Allele 1 Matches
    List.1MATCH.Demo <- list(length(which(SCORES[Range,2]==1)),
                             length(which(SCORES[Range,3]==1)),
                             length(which(SCORES[Range,4]==1)),
                             length(which(SCORES[Range,5]==1)),
                             length(which(SCORES[Range,6]==1)))
    
    #Per Individual Allele 2 Matches
    List.2MATCH.Demo <- list(length(which(SCORES[Range,2]==2)),
                             length(which(SCORES[Range,3]==2)),
                             length(which(SCORES[Range,4]==2)),
                             length(which(SCORES[Range,5]==2)),
                             length(which(SCORES[Range,6]==2)))
    
    Tab <- as.data.frame(cbind(Grp,
                               length(Range),
                               c("A","B","C","DRB1","DBQ1"),
                               List.All.Demo,
                               List.0MATCH.Demo,
                               List.1MATCH.Demo,
                               List.2MATCH.Demo))
    colnames(Tab) <- colnames(Output.Demo)
    Output.Demo <- rbind(Output.Demo,Tab)
    
  }; rm(k)
  
  nameDEMO <- paste("SCORES_byDEMO_", AncLevel, "_output.txt",sep="")
  write.table(as.matrix(Output.Demo),file=nameDEMO,sep="\t",quote=F,row.names=F,col.names=T,append=F)
  
  rm(list=ls()[!(ls() %in% SAFE)])


# 5. PERFORMANCE METHOD: HARD THRESHOLDING ----------------------------------------------------------
  
  # Graph Performance of Imputation by incremented cut points M1
  # cut points from list of unique probabilities
  # Allows for individual probabilities for each allele at a given locus
  # Requires: IMP, Loci, MATCH
  # Writes: 4D_*NAME*_Pcurve_A.txt, *NAME*_Pcurve_B.txt, *NAME*_Pcurve_C.txt, *NAME*_Pcurve_DRB1.txt, *NAME*_Pcurve_DQB1.txt
  
  fileSTEM <- gsub(".txt","",fileNAME)
  
  for(k in 1:length(Loci)) {
    
    Locus <- Loci[k]
    Tab <- data.frame(Cut=character(),
                      Alleles=character(),
                      AllelesIMP=character(),
                      Score.All=character(), 
                      Score.0=character(),
                      Score.1=character(),
                      Score.2=character(),
                      stringsAsFactors=FALSE)
    
    #IMPUTED Table, extracting allele specific Columns
    CN.IMP <- colnames(IMP)
    
    Grp.IMP <- grep(paste(Locus,".",sep=""),CN.IMP,fixed=T)
    IMPsub.A <- IMP[,c(1,Grp.IMP)]
    
    #Match Tables, extracting allele specific columns
    CN.MATCH <- colnames(MATCH)
    Grp.MATCH <- grep(Locus,CN.MATCH,fixed=T)[1:2]
    MATCHsub <- data.matrix(MATCH[,Grp.MATCH])
    
    #Adjusting for unknown genotypes
    noIMP.start <- nrow(na.omit(MATCHsub))
    IMPsub.A[setdiff(1:nrow(MATCH),which(MATCHsub[,1]>=0)),4] <- NA
    IMPsub.A[setdiff(1:nrow(MATCH),which(MATCHsub[,2]>=0)),5] <- NA
        
    #Cutpoint vector
    Probs <- c(0,sort(unique(union(IMPsub.A[,4],IMPsub.A[,5])),decreasing=F))
    
    for(j in 1:length(Probs)) {
      Cut <- Probs[j]
      RangeCut1 <- which(IMPsub.A[,4]<=Cut) #Allele 1 Calls Cut
      RangeCut2 <- which(IMPsub.A[,5]<=Cut) #Allele 2 Calls Cut
      
      if(length(RangeCut1)>0) { MATCHsub[RangeCut1,1] <- NA }
      if(length(RangeCut2)>0) { MATCHsub[RangeCut2,2] <- NA }
      
      # Number of imputed alleles
      noImp <- (noIMP.start - length(RangeCut1)) + (noIMP.start - length(RangeCut2))
      
      # Correct overall allele imputations
      Match.Overall <- as.numeric(sum(rowSums(MATCHsub,na.rm=TRUE))) / noImp
      
      # No Correct Calls
      Match.0 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==0))) / length(rowSums(MATCHsub,na.rm=T))
      
      # At least 1 correct calls
      Match.1 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==1))) / length(rowSums(MATCHsub,na.rm=T))
      
      # At least 2 correct calls
      Match.2 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==2))) / length(rowSums(MATCHsub,na.rm=T))
      
      if(noImp==0){
        Tab[j,] <- c(Cut,                 
                     rep(0,6))        
      } else {
        Tab[j,] <- c(Cut,
                     noImp,
                     (100*noImp) / (noIMP.start*2),
                     Match.Overall,
                     Match.0,
                     Match.1,
                     Match.2) }
    }; rm(j)
    
    
    CurveName=paste(fileSTEM,"_PcurveM1alt_",Locus,".txt",sep="")
    write.table(Tab,file=CurveName,sep="\t",quote=F,col.names=T,row.names=F)
  }; rm(k)
  
  Names <- list(paste(fileSTEM,"_PcurveM1alt_","A.txt",sep=""),
                paste(fileSTEM,"_PcurveM1alt_","B.txt",sep=""),
                paste(fileSTEM,"_PcurveM1alt_","C.txt",sep=""),
                paste(fileSTEM,"_PcurveM1alt_","DRB1.txt",sep=""),
                paste(fileSTEM,"_PcurveM1alt_","DQB1.txt",sep=""))
  
  ATab <- read.table(file = Names[[1]], sep="\t", header=T)
  BTab <- read.table(file = Names[[2]], sep="\t", header=T)
  CTab <- read.table(file = Names[[3]], sep="\t", header=T)
  DRTab <- read.table(file = Names[[4]], sep="\t", header=T)
  DQTab <- read.table(file = Names[[5]], sep="\t", header=T)
  
  #Plot Cut Tables. Assumes reverse x-Axis order
  palette(rev(rich.colors(32)))
  
  ImageName=paste(fileSTEM,"_Performance_HARD.pdf",sep="")
  cairo_pdf(filename = ImageName, width = 9, height = 5, onefile=TRUE)
  
  Ymax <- round(max(ATab$Score.All,BTab$Score.All,CTab$Score.All,DRTab$Score.All,DQTab$Score.All),1)+0.1
  Ymin <- round(min(min(ATab$Score.All[which(ATab$Score.All!=0)]),
                    min(BTab$Score.All[which(BTab$Score.All!=0)]),
                    min(CTab$Score.All[which(CTab$Score.All!=0)]),
                    min(DRTab$Score.All[which(DRTab$Score.All!=0)]),
                    min(DQTab$Score.All[which(DQTab$Score.All!=0)])),1)-0.1
  if(Ymax > 1) { Ymax <- 1 }
  if(Ymin < 0) { Ymin <- 0 }
  
  END <- max(which(ATab$AllelesIMP!=0))
  x <- ATab$AllelesIMP[1:END]
  y <- ATab$Score.All[1:END]
  Colmap=1 + 31*ATab$Cut
  plot(y ~ x,cex=0.8,ylim=c(Ymin,Ymax),xlim=c(100,0),xlab="Call Rate",ylab="Performance",main="Imputation Performance by Locus\nHard Threshold",col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"A",font=2,adj=0)
  
  END <- max(which(BTab$AllelesIMP!=0))
  x <- BTab$AllelesIMP[1:END]
  y <- BTab$Score.All[1:END]
  Colmap=1 + 31*BTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0.5,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"B",font=2,adj=0)
  
  END <- max(which(CTab$AllelesIMP!=0))
  x <- CTab$AllelesIMP[1:END]
  y <- CTab$Score.All[1:END]
  Colmap=1 + 31*CTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0.5,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"C",font=2,adj=0)
  
  END <- max(which(DRTab$AllelesIMP!=0))
  x <- DRTab$AllelesIMP[1:END]
  y <- DRTab$Score.All[1:END]
  Colmap=1 + 31*DRTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0.5,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DRB1",font=2,adj=0)

  END <- max(which(DQTab$AllelesIMP!=0))
  x <- DQTab$AllelesIMP[1:END]
  y <- DQTab$Score.All[1:END]
  Colmap=1 + 31*DQTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0.5,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DQB1",font=2,adj=0)
  
  #Labels = color bar properties: min value, max value, description
  Labels <- list(0,1,"Probability Cut Threshold")
  switch(Bar,
         Up = colorBAR(40,60,Ymax-0.01,Ymax-Ymin,Labels),
         Down = colorBAR(40,60,Ymin+0.01,Ymax-Ymin,Labels))
  
  dev.off()
  
  rm(list=ls()[!(ls() %in% SAFE)])


# 6. PERFORMANCE Method: LINEAR BINNING  ----------------------------------------------------------
  
  # Graph Performance of Imputation by incremented cut points M3
  # Allows for individual probabilities for each allele at a given locus
  # Requires: IMP, Loci, MATCH
  # Writes: 2D_*NAME*_Pcurve_A.txt, 2D_*NAME*_Pcurve_B.txt, 2D_*NAME*_Pcurve_C.txt, 2D_*NAME*_Pcurve_DRB1.txt
  
  fileSTEM <- gsub(".txt","",fileNAME)
  
  Increments <- 1000
  
  for(k in 1:length(Loci)){
    
    Locus <- Loci[k]
    Tab <- data.frame(Cut=character(),
                      Alleles=character(),
                      AllelesIMP=character(),
                      Score.All=character(), 
                      Score.0=character(),
                      Score.1=character(),
                      Score.2=character(),
                      stringsAsFactors=FALSE)
    
    #IMPUTED Table, extracting allele specific Columns
    CN.IMP <- colnames(IMP)
    Grp.IMP <- grep(paste(Locus,".",sep=""),CN.IMP,fixed=T)
    IMPsub.A <- IMP[,c(1,Grp.IMP)]
    
    #Match Tables, extracting allele specific columns
    CN.MATCH <- colnames(MATCH)
    Grp.MATCH <- grep(Locus,CN.MATCH,fixed=T)[1:2]
    MATCHsub <- data.matrix(MATCH[,Grp.MATCH])
    
    #Adjusting for unknown genotypes
    noIMP.start <- nrow(na.omit(MATCHsub))
    IMPsub.A[setdiff(1:nrow(MATCH),which(MATCHsub[,1]>=0)),4] <- NA
    IMPsub.A[setdiff(1:nrow(MATCH),which(MATCHsub[,2]>=0)),5] <- NA
    
    #Cutpoint vector
    CutPoints <- c(0,(1:Increments)/Increments)
    
    for(j in 1:length(CutPoints)){
      Cut <- CutPoints[j]
      RangeCut1 <- which(IMPsub.A[,4]<=Cut) #Allele 1 Calls Cut
      RangeCut2 <- which(IMPsub.A[,5]<=Cut) #Allele 2 Calls Cut
      
      if(length(RangeCut1)>0) { MATCHsub[RangeCut1,1] <- NA }
      if(length(RangeCut2)>0) { MATCHsub[RangeCut2,2] <- NA }
      
      # Number of imputed alleles
      noImp <- (noIMP.start - length(RangeCut1)) + (noIMP.start - length(RangeCut2)) 
      
      # Correct overall allele imputations
      Match.Overall <- as.numeric(sum(rowSums(MATCHsub,na.rm=TRUE))) / noImp
      
      # No Correct Calls
      Match.0 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==0))) / length(rowSums(MATCHsub,na.rm=T))
      
      # At least 1 correct calls
      Match.1 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==1))) / length(rowSums(MATCHsub,na.rm=T))
      
      # At least 2 correct calls
      Match.2 <- as.numeric(length(which(rowSums(MATCHsub,na.rm=T)==2))) / length(rowSums(MATCHsub,na.rm=T))
      
      if(noImp==0){
        Tab[j,] <- c(Cut,                 
                     rep(0,6))        
      } else {
        Tab[j,] <- c(Cut,
                     noImp,
                     (100*noImp) / (noIMP.start*2),
                     Match.Overall,
                     Match.0,
                     Match.1,
                     Match.2) }
    }; rm(j)
    
    
    CurveName=paste(fileSTEM,"_PcurveM3_",Locus,".txt",sep="")
    write.table(Tab,file=CurveName,sep="\t",quote=F,col.names=T,row.names=F)
  }; rm(k)
  
  # FLAG: Change to incorporate DQB1
  Names <- list(paste(fileSTEM,"_PcurveM3_","A.txt",sep=""),
                paste(fileSTEM,"_PcurveM3_","B.txt",sep=""),
                paste(fileSTEM,"_PcurveM3_","C.txt",sep=""),
                paste(fileSTEM,"_PcurveM3_","DRB1.txt",sep=""),
                paste(fileSTEM,"_PcurveM3_","DQB1.txt",sep=""))
  
  # FLAG: Change to incorporate DQB1
  ATab <- read.table(file = Names[[1]], sep="\t", header=T)
  BTab <- read.table(file = Names[[2]], sep="\t", header=T)
  CTab <- read.table(file = Names[[3]], sep="\t", header=T)
  DRTab <- read.table(file = Names[[4]], sep="\t", header=T)
  DQTab <- read.table(file = Names[[5]], sep="\t", header=T)
  
  #Plot Cut Tables. Assumes reverse x-Axis order
  palette(rev(rich.colors(32)))
  
  ImageName=paste(fileSTEM,"_Performance_LinBin.pdf",sep="")
  cairo_pdf(filename = ImageName, width = 9, height = 5, onefile=TRUE)
  
  Ymax <- round(max(ATab$Score.All,BTab$Score.All,CTab$Score.All,DRTab$Score.All,DQTab$Score.All),1)+0.1
  Ymin <- round(min(min(ATab$Score.All[which(ATab$Score.All!=0)]),
                    min(BTab$Score.All[which(BTab$Score.All!=0)]),
                    min(CTab$Score.All[which(CTab$Score.All!=0)]),
                    min(DRTab$Score.All[which(DRTab$Score.All!=0)]),
                    min(DQTab$Score.All[which(DQTab$Score.All!=0)])),1)-0.1
  if(Ymax > 1) { Ymax <- 1 }
  if(Ymin < 0) { Ymin <- 0 }
  
  END <- max(which(ATab$AllelesIMP!=0))
  x <- ATab$AllelesIMP[1:END]
  y <- ATab$Score.All[1:END]
  Colmap=1 + 31*ATab$Cut
  plot(y ~ x,cex=0.8,ylim=c(Ymin,Ymax),xlim=c(100,0),xlab="Call Rate",ylab="Performance",main="Imputation Performance by Locus\nLinear Binning",col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"A",font=2,adj=0)
  
  END <- max(which(BTab$AllelesIMP!=0))
  x <- BTab$AllelesIMP[1:END]
  y <- BTab$Score.All[1:END]
  Colmap=1 + 31*BTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"B",font=2,adj=0)
  
  END <- max(which(CTab$AllelesIMP!=0))
  x <- CTab$AllelesIMP[1:END]
  y <- CTab$Score.All[1:END]
  Colmap=1 + 31*CTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"C",font=2,adj=0)
  
  END <- max(which(DRTab$AllelesIMP!=0))
  x <- DRTab$AllelesIMP[1:END]
  y <- DRTab$Score.All[1:END]
  Colmap=1 + 31*DRTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DRB1",font=2,adj=0)
  
  END <- max(which(DQTab$AllelesIMP!=0))
  x <- DQTab$AllelesIMP[1:END]
  y <- DQTab$Score.All[1:END]
  Colmap=1 + 31*DQTab$Cut
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(100,0),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DQB1",font=2,adj=0)
  
  #Labels = color bar properties: min value, max value, description
  Labels <- list(0,1,"Probability Cut Threshold")
  switch(Bar,
         Up = colorBAR(40,60,Ymax-0.01,Ymax-Ymin,Labels),
         Down = colorBAR(40,60,Ymin+0.01,Ymax-Ymin,Labels))
  
  dev.off()
  
  #Plot Cut Tables. Alternative x-axis for graphing data v2.
  palette(rev(rich.colors(32)))
  
  ImageName=paste(fileSTEM,"_Performance_LinBinv2.pdf",sep="")
  cairo_pdf(filename = ImageName, width = 9, height = 5, onefile=TRUE, pointsize=10)
  
  END <- max(which(ATab$AllelesIMP!=0))
  x <- ATab$Cut[1:END]
  y <- ATab$Score.All[1:END]
  Colmap=1 + 31*(ATab$AllelesIMP[1:END])/100
  plot(y ~ x,cex=0.8,ylim=c(Ymin,Ymax),xlim=c(0,1),xlab="Cut Point",ylab="Performance",main="Imputation Performance by Locus\nLinear Binning",col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"A",font=2,adj=0)
  
  END <- max(which(BTab$AllelesIMP!=0))
  x <- BTab$Cut[1:END]
  y <- BTab$Score.All[1:END]
  Colmap=1 + 31*(BTab$AllelesIMP[1:END])/100
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(0,100),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"B",font=2,adj=0)
  
  END <- max(which(CTab$AllelesIMP!=0))
  x <- CTab$Cut[1:END]
  y <- CTab$Score.All[1:END]
  Colmap=1 + 31*(CTab$AllelesIMP[1:END])/100
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(0,100),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"C",font=2,adj=0)
  
  END <- max(which(DRTab$AllelesIMP!=0))
  x <- DRTab$Cut[1:END]
  y <- DRTab$Score.All[1:END]
  Colmap=1 + 31*(DRTab$AllelesIMP[1:END])/100
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(0,100),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DRB1",font=2,adj=0)
  
  END <- max(which(DQTab$AllelesIMP!=0))
  x <- DQTab$Cut[1:END]
  y <- DQTab$Score.All[1:END]
  Colmap=1 + 31*(DQTab$AllelesIMP[1:END])/100
  points(y ~ x,cex=0.8,ylim=c(0,1),xlim=c(0,100),col=Colmap,pch=20)
  text(x[1]-0.025,y[1],"DQB1",font=2,adj=0)
  
  #Labels = color bar properties: min value, max value, description
  Labels <- list(0,1,"Call Rate")
  switch(Bar,
         Up = colorBAR(0.4,0.6,Ymax-0.01,Ymax-Ymin,Labels),
         Down = colorBAR(0.4,0.6,Ymin+0.01,Ymax-Ymin,Labels))
  
  dev.off()
  
  rm(list=ls()[!(ls() %in% SAFE)])
  


# 7. SCORE BY HLA Genotype Calls ----------------------------------------------------------

  # Score outputs of imputed HLA calls by Allele
  # Requires: IMP, Loci, MATCH
  # Writes: SCORES_byAllele_output.txt
  
  Output.byAllele <- data.frame(Locus=character(),
                                Allele=character(),
                                ImpCount=character(),
                                PossCorrect=character(),
                                CorrectCount=character(),
                                CorrectPercent=character())
  
  for(k in 1:length(Loci)) {
    
    Locus <- Loci[k]
    
    # Match Tables, extracting allele specific columns
    CN.MATCH <- colnames(MATCH)
    Grp.MATCH <- grep(Locus,CN.MATCH,fixed=T)[1:2]
    MATCHsub <- data.matrix(MATCH[,Grp.MATCH])
    
    # IMPUTED Table, extracting allele specific Columns
    CN.IMP <- colnames(IMP)
    Grp.IMP <- which(substr(CN.IMP,1,nchar(Locus)+1)==paste(Locus,".",sep=""))
    IMPsub.A <- IMP[,Grp.IMP]
    
    # Generating unique Alleles list
    AllUn <- sort(unique(intersect(IMPsub.A[,1],IMPsub.A[,2])))
    
    #Adjusting for unknown genotypes
    #Error .. adjusts for NA's across all defined columns, not just locus specific columns
    noIND <- as.numeric(nrow(na.omit(MATCHsub)))
    noAll <- noIND * 2    
    
    for(j in 1:length(AllUn)) {
      
      Call <- as.character(AllUn[j])
      
      Grp1 <- which(IMPsub.A[,1]==Call)
      Grp2 <- which(IMPsub.A[,2]==Call)
              
      ImpCount <- as.numeric(length(Grp1) + length(Grp2))
      PossCount <- sum(as.numeric(length(na.omit(MATCHsub[Grp1,1]))),
                       as.numeric(length(na.omit(MATCHsub[Grp2,2])))) 
      CorrectCall <- sum(as.numeric(length(which(MATCHsub[Grp1,1]==1))),
                         as.numeric(length(which(MATCHsub[Grp2,2]==1))))
      CorrectCallp <- percent(CorrectCall/PossCount,2)
            
      Tab <- cbind(Locus,
                   Call,
                   ImpCount,
                   PossCount,
                   CorrectCall,
                   CorrectCallp)
      colnames(Tab) <- colnames(Output.byAllele)
      Output.byAllele <- as.matrix(rbind(Output.byAllele,Tab))
      
    }; rm(j)
    
  }; rm(k)
  
  Output.byAllele <- Output.byAllele[order(Output.byAllele[,1],-as.numeric(Output.byAllele[,3])), ]
  write.table(Output.byAllele,file="SCORES_byImpAllele_output.txt",sep="\t",quote=F,col.names=T,row.names=F)
  
  rm(list=ls()[!(ls() %in% SAFE)])




#  
###################################################################################### End Analysis Script ####

  setwd(MainDir)
  
}




