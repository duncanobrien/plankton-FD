###########################################################################################
### TidyFD function ###
###########################################################################################
#@speciesdat = matrix of species abundances, names across columns, time across rows
#@traitdat = matrix of species traits, names across rows, traits across columns
#@traittype = fuzzy/ordinal
#@trophic.lvl = phyto(plankton), zoo(plankton) or fish. For classification
#@method = which functional diversity metric to be calculated. Multiple can be provided
#@correction = to correct for non-Euclidean distances during construction of trait space
#@ndim = maximum number of dimensions for trait space

tidyFD <- function(speciesdat, traitdat, trophic.lvl = c("phyto","zoo","fish"), traittype = c("fuzzy","ordinal"),
                   method = c("MNND","Rao","FDis","FRic","FEve","FDiv","SRic","SDiv","RARDis","RARScr"),
                   correction = c("sqrt", "cailliez", "lingoes", "none"), ndim = 10){
  require(tibble)
  require(dplyr)
  require(FD)
  require(picante)
  require(ade4)
  require(vegan)
  require(gawdis)
  require(funrar)
  set.seed(1234)
  
source("Code/melodic_fn.R")
  
  #select trait data type. Fuzzy used in this work
  if(traittype == "ordinal"){
    #ensure matched species timeseries and trait data
    species.data <- speciesdat[,colSums(speciesdat != 0) > 0] #drop any full zero timeseries
    traits.data <- traitdat[rowSums(is.na(traitdat)) != (ncol(traitdat)-1),] #drop any full NA species trait vectors                                                                                                      
    
    a <- colnames(species.data)
    b <- traits.data$Species
    c <- intersect(a,b) #match species present in both datasets
    
    species.data<- species.data[, c]
    traits.data <- traits.data[traits.data$Species %in% c, ]
    
    if(length(which(trophic.lvl == "phyto"))==1){
      traits.data <- traits.data %>%
        mutate_at(c("Cell length", "Cell SA","Cell biovol",
                    "Flagellated", "Colonial","Filamentous", "Heterotrophic","Mixotrophic",
                    "Nitrogen"), as.numeric) %>% 
        mutate_at(c("Geometrical shape","Morphoclass"),as.factor) %>%
        column_to_rownames(var = "Species")
      
      traits.data.in <- as.matrix(gowdis(traits.data))
    }
    
    
    if(length(which(trophic.lvl == "zoo"))==1){
      
      traits.data <- traits.data %>% 
        mutate_at(c("Trophic group", "Feeding mode","Reproductive mechanism"), as.factor) %>% 
        mutate_at(c("Mean length μm"),as.numeric) %>%
        #traits.data$`Feeding mode`<- factor(traits.data$`Feeding mode`, levels = c("Cruise","Current","Active ambush"))
        ##traits.data$`Trophic group`<- factor(traits.data$`Trophic group`, levels = c("Mi","Ma","P"), ordered = T)
        #traits.data$`Mean length μm`<- as.numeric(traits.data$`Mean length μm`)
        #traits.data$`Reproductive mechanism`<- factor(traits.data$`Reproductive mechanism`, levels = c("sac","cyc parth","noncyc parth"))
        column_to_rownames( var = "Species")
      
      traits.data.in <- as.matrix(gowdis(traits.data))  
    }
    
    if(length(which(trophic.lvl == "fish"))==1){
      
      traits.data.in <- gowdis(traits.data)
    }
  }
  
  if(traittype=="fuzzy"){
    species.data <- speciesdat[,colSums(speciesdat != 0) > 0] #drop any full zero timeseries
    traits.data <- traitdat[rowSums(is.na(traitdat)) != (ncol(traitdat)-1),] #drop any full NA species trait vectors                                                                                                      
    
    a <- colnames(species.data)
    b <- traits.data$Name
    int <- intersect(a,b) #match species present in both datasets
    
    species.data <- species.data[, int]
    traits.data <- traits.data[traits.data$Name %in% int, ]
    
    if(length(which(trophic.lvl == "phyto"))==1){
      #fuz.trait <-prep.fuzzy.var(df = traits.data[,-c(1:3)], #df = data for  trait (s)
      #col.blocks = c(4,4,4,4,2,3,3), #col.blocks = num categories per trait
      #row.w = rep(1, nrow(traits.data[]))) #row.w = weight of each category
      
      #fuz.lgth <- prep.fuzzy.var(traits.data[,4:7], #df = data for one trait i.e. cell length
      #col.blocks = 4, #col.blocks = num categories for trait
      #row.w = rep(1, nrow(traits.data[,4:7]))) #row.w = weight of each species
      
      #fuz.sa <- prep.fuzzy.var(traits.data[,8:11], col.blocks = 4,row.w = rep(1, nrow(traits.data[,8:11])))
      
      #fuz.vol <- prep.fuzzy.var(traits.data[,12:15],col.blocks = 4,row.w = rep(1, nrow(traits.data[,12:15])))
      
      #fuz.c.ratio <- prep.fuzzy.var(traits.data[,16:19],col.blocks = 4,row.w = rep(1, nrow(traits.data[,16:19])))
      
      #fuz.col <- prep.fuzzy.var(traits.data[,20:21],col.blocks = 2,row.w = rep(1, nrow(traits.data[,20:21])))
      
      #fuz.mob <- prep.fuzzy.var(traits.data[,22:24],col.blocks = 3,row.w = rep(1, nrow(traits.data[,22:24])))
      
      #fuz.troph <- prep.fuzzy.var(traits.data[,25:27],col.blocks = 3,row.w = rep(1, nrow(traits.data[,25:27])))
      
      #fuz.n2 <- prep.fuzzy.var(traits.data[,28:29],col.blocks = 2,row.w = rep(1, nrow(traits.data[,28:29])))
      
      #fuz.sil <- prep.fuzzy.var(traits.data[,30:31],col.blocks = 2,row.w = rep(1, nrow(traits.data[,30:31])))
      
      #ktab <- ktab.list.df(list(fuz.lgth,fuz.sa,fuz.vol,fuz.c.ratio,fuz.col,fuz.mob,fuz.troph,fuz.n2,fuz.sil))
      #ktab<- ktab.list.df(list(fuz.trait))
      #traits.data.in <- dist.ktab(ktab,type = c(rep("F",times=9)),scann = F)
    
      #each group represents similar traits e.g. size is represented by length/SA/vol
      #Can also be used for fuzzy data where a group are the categories for a single trait
      traits.data.in<-gawdis::gawdis(traits.data[,4:20], w.type="equal", 
                                     groups =c(1,1,1,1,2,2,2,2,2,3,4,4,5,6,7,8,9), fuzzy=c(1,2,4)) #estimate gawdis dissimilarity for fuzzy data
      attr(traits.data.in,"Labels") <- int #label rows and columns with species names    
    }
    
    if(length(which(trophic.lvl == "zoo"))==1){
      traits.data.in<-gawdis::gawdis(traits.data[,5:15], w.type="equal", 
                                     groups =c(1,1,1,1,2,2,2,2,2,3,4), fuzzy=c(1,2))
      attr(traits.data.in,"Labels") <- int    
    }
  }
  
  cor<-match.arg(correction, choices=c("sqrt", "cailliez", "lingoes", "none"), several.ok=F)
  func.tmp <- match.fun(cor) #convert string to function from FD package
  meth <- match.arg(method , choices= c("MNND","Rao","FDis","FDiv","FRic","FEve","SRic","SDiv","RARDis","RARScr"), several.ok=T)
  
  if(length(which(meth == "MNND"))==1){
    MNNDout <- picante::mntd(samp = as.matrix(species.data), dis = as.matrix(func.tmp(traits.data.in)), abundance.weighted = T)
    
  }
  
  if(length(which(meth == "Rao"))==1){
    RAOout <- melodic(samp = as.matrix(species.data), dis = as.matrix(func.tmp(traits.data.in)),type="abundance")$abundance$rao
    
  }
  
  if(length(which(grepl("^[RAR]", meth)==1))){
    rel<- funrar::make_relative(as.matrix(species.data))
    funrar <- funrar::funrar(as.matrix(rel),as.matrix(traits.data.in),rel_abund = T) 
    #can't apply correction due to incomaptibility with funrar's calc
    
    if(length(which(meth == "RARDis"))==1){
      RARDISout <- rowMeans(funrar$Di,na.rm=T)
    }
    if(length(which(meth == "RARScr"))==1){
      funrar$Si[is.na(funrar$Si)] <- 1 #if Si is close to 1, the species has relatively low abundances. Therefore 0 abundance would correspond to Si = 1 
      RARSCRout <- rowMeans(funrar$Si,na.rm=T)
    }
  }
  
  if(length(which(grepl("^[F]", meth)==1))){ #if meth starts with the letter 'F'
    tmp.FD <- FD::dbFD(x = traits.data.in,a = as.matrix(species.data), corr = paste(cor), 
                       calc.CWM = F, calc.FRic = T,calc.FDiv = T,m=ndim) #for speed, create single dbFD object
    
    if(length(which(meth == "FDis"))==1){
      FDISout <- tmp.FD$FDis
    }
    if(length(which(meth == "FRic"))==1){
      FRICout <- tmp.FD$FRic
    }
    if(length(which(meth == "FEve"))==1){
      FEVEout <-tmp.FD$FEve
    }
    if(length(which(meth == "FDiv"))==1){
      FDIVout <-tmp.FD$FDiv
    }
  }
  
  if(length(which(meth == "SRic"))==1){
    SRICout <- vegan::specnumber(x = as.matrix(species.data),MARGIN = 1)
  }
  
  if(length(which(meth == "SDiv"))==1){
    SDIVout <- vegan::diversity(x = as.matrix(species.data),index = "shannon",MARGIN = 1,)
  }
  
  out <- data.frame("MNND" = if(length(which(meth=="MNND"))==1){"MNND"=MNNDout}else{NA},
                    "Rao" = if(length(which(meth=="Rao"))==1){"Rao"=RAOout}else{NA},
                    "FRic" = if(length(which(meth=="FRic"))==1){"FRic"=FRICout}else{NA},
                    "FEve" = if(length(which(meth=="FEve"))==1){"FEve"=FEVEout}else{NA},
                    "FDiv" = if(length(which(meth=="FDiv"))==1){"FDiv"=FDIVout}else{NA},
                    "FDis" = if(length(which(meth=="FDis"))==1){"FDis"=FDISout}else{NA},
                    "SRic" = if(length(which(meth=="SRic"))==1){"SRic"=SRICout}else{NA},
                    "SDiv" = if(length(which(meth=="SDiv"))==1){"SDiv"=SDIVout}else{NA},
                    "RARDis" = if(length(which(meth=="RARDis"))==1){"RARDis"=RARDISout}else{NA},
                    "RARScr" = if(length(which(meth=="RARScr"))==1){"RARScr"=RARSCRout}else{NA}
  )
  return(out[,meth])
  
}
