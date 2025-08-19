
#load libraries
my_lib <- "path/to/my/lib"

BiocManager::install("DECIPHER", lib = my_lib)

install.packages("devtools", lib = my_lib)
install_github("janhoo/worms", lib = my_lib)
install_github("pmartinezarbizu/dada2pp", lib = my_lib)

library(DECIPHER, lib.loc = my_lib)
library(devtools, lib.loc = my_lib)
library(dada2pp, lib.loc = my_lib)
library(worms, lib.loc = my_lib)
library(dplyr)
library(tidyverse)
library(data.table)
library(vegan)
library(ggplot2)
library(reshape2)

#load ASV Table
meio<-read.csv('taxonomy_all_worms_curated_final.csv')
nrow(meio)
#[1] 41711

#load reads
read<-read.csv('ASV_table_metazoa_no_rfy.csv')
colnames(read)[1] <- "asv"
#[1] 18925


#add taxonomy to the counts filtering by taxonomy
ASV_meio<- left_join(read, meio, by = c("asv"))

nrow(ASV_meio)
#[1]18925


#######################################################################
#aggregate by NJ 
#load sequences
seqs<-readDNAStringSet('Seqs_nochim_Lukas_Mangan18_eDNA_bulk.fasta')
str(seqs)
#get ASV numbers
ASVno<-data.frame('row'=grep('^*.|',names(seqs)),'asv'=as.character(sapply(strsplit(names(seqs),'|',1),"[[",1)))

#get ASVs of meio
ASVmeio<-ASV_meio[,1]

seqmeio<-seqs[which(ASVno$asv%in%ASVmeio)]

#align the sequences
aligmeio<-AlignSeqs(seqmeio)

#view the alignment
#BrowseSeqs(aligmacro,highlight=0)

#MLTree
d<-DistanceMatrix(aligmeio,correction="Jukes-Cantor",verbose=FALSE,processors=48)

#clustering by 97% and 98% similarity
c<- Clusterize(seqmeio, cutoff=.03, invertCenters = T, processors=1)
#c2<-IdClusters(d,method="NJ",cutoff=.02,showPlot=TRUE,myXStringSet=aligmeio,verbose=FALSE,processors=48)
length(unique(c[,1]))
# 3448


c$asv <- rownames(c)
c$OTU_ID <- paste0("OTU_", as.numeric(as.factor(c[,1])))
head(c)
#order by ASV numbers and aggregate by 97%
#ASV_meio<-ASV_meio[order(ASV_meio$asv),]
ASV_meio <- left_join(ASV_meio, c[, c("asv", "OTU_ID")], by = "asv")

aggNJ_97<-aggregateASV(x=ASV_meio[,2:38],
                       by=ASV_meio$OTU_ID,
                       other_str=ASV_meio[,c(39:49)])

colnames(aggNJ_97)[12]<-"OTU_ID"


nrow(aggNJ_97)
#[1] 3448



table(as.character(aggNJ_97$kingdom))
#Metazoa 
#3448 


#keep only metazoan

aggNJ_97<-aggNJ_97[aggNJ_97$kingdom%in%c("Metazoa"),]

nrow(aggNJ_97)
#[1] 3448


table(as.character(aggNJ_97$Phylum))
#                           Annelida                       Arthropoda                      Brachiopoda 
#1                              736                              747                               14 
#Bryozoa     c("Annelida", "Brachiopoda")        c("Annelida", "Mollusca")       c("Annelida", "Phoronida") 
#31                                4                                7                                1 
#c("Annelida", "Platyhelminthes")       c("Bryozoa", "Arthropoda")        c("Cnidaria", "Porifera")        c("Mollusca", "Annelida") 
#5                                1                                1                                2 
#c("Nematoda", "Tardigrada") c("Platyhelminthes", "Annelida")      c("Tardigrada", "Nematoda")                     Chaetognatha 
#1                                3                                1                                2 
#Chordata                         Cnidaria                    Echinodermata                       Entoprocta 
#24                               48                                6                                1 
#Gastrotricha                     Hemichordata                       Loricifera                         Mollusca 
#34                               17                               15                               77 
#Nematoda                         Nemertea                     Orthonectida                  Platyhelminthes 
#1198                               52                                4                              256 
#Porifera                         Rotifera                       Tardigrada                  Xenacoelomorpha 
#38                                1                               15                               96 
 

#number of ambiguous phylum
#58 

#keep only unique phylums
aggNJ_97<-aggNJ_97[aggNJ_97$Phylum%in%c("Annelida","Arthropoda","Bivalvia","Brachiopoda","Bryozoa","Chaetognatha","Chordata","Cnidaria","Echinodermata",
"Entoprocta","Gastrotricha","Hemichordata","Loricifera","Mollusca","Nematoda","Nemertea","Platyhelminthes","Porifera","Tardigrada","Xenacoelomorpha"),]

nrow(aggNJ_97)
#[1] 3407



table(as.character(aggNJ_97$Phylum))

#Annelida      Arthropoda     Brachiopoda         Bryozoa    Chaetognatha        Chordata        Cnidaria   Echinodermata      Entoprocta 
#736             747              14              31               2              24              48               6               1 
#Gastrotricha    Hemichordata      Loricifera        Mollusca        Nematoda        Nemertea Platyhelminthes        Porifera      Tardigrada 
#34              17              15              77            1198              52             256              38              15 
#Xenacoelomorpha 
#96 

 


table(as.character(aggNJ_97$Class))


#Appendicularia                   Ascidiacea                   Asteroidea                     Bivalvia 
#10                            7                            4                            3                           59 
#Branchiopoda  c("Chromadorea", "Enoplea")         c("Chromadorea", NA)   c("Copepoda", "Crustacea") c("Enopla", "Hoplonemertea") 
#1                            1                           13                            1                            1 
#c("Ostracoda", "") c("Thaliacea", "Thaliaceae")         c(NA, "Chromadorea")                     Calcarea                  Chromadorea 
#1                            1                            6                           22                          704 
#Clitellata                     Copepoda                 Demospongiae                 Eleutherozoa                       Enopla 
#7                          638                           16                            1                            1 
#Enoplea                Enteropneusta                   Gastropoda                 Gymnolaemata             Heterotardigrada 
#391                           17                           18                           23                            8 
#Hexacorallia                Holothuroidea                Hoplonemertea                     Hydrozoa                    Lingulata 
#8                            1                           44                           35                           10 
#Malacostraca                 Octocorallia                  Ophiuroidea                    Ostracoda               Palaeonemertea 
#20                            1                            1                           77                            2 
#Pilidiophora                   Polychaeta                Rhabditophora               Rhynchonellata                  Sagittoidea 
#4                          721                            2                            4                            2 
#Scyphozoa                 Stenolaemata                Tantulocarida                    Teleostei                    Thaliacea 
#4                            8                            1                            1                           10 
#Thaliaceae                  Thecostraca                    Trematoda 
#1                            3                            2 
 


#export tables 
fwrite(aggNJ_97, file ="aggML_no_rfy.csv", row.names = TRUE)


#assign class to the groups with NAs included. note that the groups that had no defined class according to WORMs (e.g. some Platyhelminthes) 'unassigned' has been assigned to them 
#read the tables after cleaning in excel

otu_table_3 <- read.csv("aggNJ_97.csv", row.names = 1)



table(as.character(aggNJ_98$Class))

  Appendicularia       Ascidiacea       Asteroidea         Bivalvia     Branchiopoda         Calcarea 
               6                4                2               95                1               25 
     Chromadorea       Clitellata         Copepoda         Craniata     Demospongiae     Eleutherozoa 
            1221                5             1405                3               16                1 
          Enopla          Enoplea    Enteropneusta       Gastropoda     Gymnolaemata Heterotardigrada 
               3              739               14               16               29               12 
    Hexacorallia    Holothuroidea    Hoplonemertea         Hydrozoa        Lingulata     Malacostraca 
              10                1              105               28               10               21 
    Octocorallia      Ophiuroidea        Ostracoda         Ostreida   Palaeonemertea     Pilidiophora 
               1                1              112                1                1                3 
      Polychaeta    Rhabditophora   Rhynchonellata      Sagittoidea        Scyphozoa     Stenolaemata 
            1805                2                8                7                3               15 
   Tantulocarida        Teleostei        Thaliacea       Thaliaceae      Thecostraca        Trematoda 
               1                1                6                2                3                2 
      unassigned 
             639 

table(as.character(aggNJ_97$Class))

  Appendicularia       Ascidiacea       Asteroidea         Bivalvia     Branchiopoda         Calcarea 
               5                4                2               69                1               17 
     Chromadorea       Clitellata         Copepoda         Craniata     Demospongiae           Enopla 
             858                5              751                3                8                1 
         Enoplea         Enoplea|    Enteropneusta       Gastropoda     Gymnolaemata Heterotardigrada 
             487                1               12               13               18                7 
    Hexacorallia    Holothuroidea    Hoplonemertea         Hydrozoa        Lingulata     Malacostraca 
               8                1               60               22                9               15 
    Octocorallia      Ophiuroidea        Ostracoda   Palaeonemertea     Pilidiophora       Polychaeta 
               1                1               72                1                2             1083 
   Rhabditophora   Rhynchonellata      Sagittoidea        Scyphozoa     Stenolaemata    Tantulocarida 
               2                8                5                2                9                1 
       Teleostei        Thaliacea       Thaliaceae      Thecostraca        Trematoda       unassigned 
               1                6                1                3                2              459 


