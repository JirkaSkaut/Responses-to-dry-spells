#############################################################################################
#           Responses of stem and leaf biomass of temperate conifers to drought spells 
#############################################################################################

#                 Jiří Mašek a*, Isabel Dorado Liñán b, Václav Treml a

#   a Department of Physical Geography and Geoecology, Faculty of Science, Charles University, Albertov 6, 128 43 Prague, Czech Republic
#   b Dpto. de Sistemas y Recursos Naturales, Universidad Politécnica de Madrid, Madrid, Spain.
#   * Corresponding author: jiri.masek@natur.cuni.cz (Jiří Mašek)


# necessary packages
library(dplR)
library(pointRes)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(cowplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(scales)
library(tidytext)

# set working directory
setwd("C:/Users/jirka/Desktop/Pointers/")

##############################################################
#         1. Calculation of TRI chronologies
##############################################################

List_loc<- read.table("Data/List_loc.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8")

Chron<- data.frame(matrix(ncol = 40, nrow = 245)); colnames(Chron)<- List_loc$Code; rownames(Chron)<- c(1778:2022)

for (i in c(1:40)){
  
  serie <- read.rwl(paste("C:/Users/jirka/Desktop/Pointers/Data/", List_loc[i, "Code"], ".rwl", sep = ""))
  
  detrend_serie <- detrend(serie, method = "Spline", nyrs = 30)
  chronology <- chron(detrend_serie, biweight = T)
  
  for (b in rownames(chronology)) {
    
    Chron[rownames(Chron[b,]==rownames(chronology[b,])), colnames(Chron) == List_loc[i, "Code"]] <- chronology[b,1]
    
  }
}

TRI<- Chron[rownames(Chron)<2018 & rownames(Chron)>1984,]

write.table(TRI, "Vystupy/TRI.txt", row.names = T, col.names = T)

PISY_TRI<- TRI[,1:20]
PCAB_TRI<- TRI[,21:40]

##############################################################
#         2. Detrending of NDVI data
##############################################################

NDVI_orig <- read.table("Data/NDVI_orig.txt", check.names=FALSE, row.names = 1, header=TRUE, sep="\t", na.strings="NA", dec=",")

NDVI_det<- data.frame(matrix(ncol = ncol(NDVI_orig), nrow = nrow(NDVI_orig))); colnames(NDVI_det)<- colnames(NDVI_orig); rownames(NDVI_det)<- rownames(NDVI_orig)

for (a in c(1:ncol(NDVI_orig))) {
  
  model<- lm(NDVI_orig[,a] ~as.numeric(rownames(NDVI_orig)))
  
  NDVI_det[,a] <- model$residuals+mean(NDVI_orig[,a])
}

NDVI<- NDVI_det[rownames(NDVI_det)<2018,]

write.table(NDVI, "Vystupy/NDVI.txt", row.names = T, col.names = T)

PISY_NDVI<- NDVI[,1:20]
PCAB_NDVI<- NDVI[,21:40]

##############################################################
#         3. Preparation of climate data from ERA5
##############################################################

# Loading of climadata
List_clima<-read.table("Data/List_clima.txt", check.names=FALSE, dec=",", header = T)

JJA_mean<- data.frame(matrix(ncol = nrow(List_clima), nrow = 33)); colnames(JJA_mean)<- List_clima$Code; rownames(JJA_mean)<- c(1985:2017)

for (i in c(1:nrow(List_clima))){
  
  Clima_var <- read.table(paste("C:/Users/jirka/Desktop/Pointers/Data/", List_clima[i, "Code"], ".txt", sep = ""), dec=",")
  
  # averaging of climadata
  JJA_mean[,i]<- rowMeans(Clima_var[6:8])
  
}

write.table(JJA_mean, "Vystupy/JJA_mean.txt", row.names = T, col.names = T)

# deviding for species
PISY_JJA<- JJA_mean %>% select(contains("PISY"))
colnames(PISY_JJA)<- c("SPEI", "T", "P", "SM", "SR")
PCAB_JJA<- JJA_mean %>% select(contains("PCAB"))
colnames(PCAB_JJA)<- c("SPEI", "T", "P", "SM", "SR")

##############################################################
#    4. Correlation of TRI and NDVI with climatic variables
##############################################################

## for PISY TRI
PISY_TRI_COR<- data.frame(LOK=NA, CLIMA=NA, COR=NA, SIG=NA)
counter<- 1

for (b in c(1:ncol(PISY_TRI))) {
  
  for (c in c(1:ncol(PISY_JJA))) {
    
    Test<- cor.test(PISY_TRI[,b], PISY_JJA[,c])
    
    PISY_TRI_COR[counter, "LOK"]<- colnames(PISY_TRI)[b]
    PISY_TRI_COR[counter, "CLIMA"]<- colnames(PISY_JJA)[c]
    PISY_TRI_COR[counter, "COR"]<- Test$estimate
    PISY_TRI_COR[counter, "SIG"]<- Test$p.value
    
    counter<- counter+1
  }
}

PISY_TRI_COR$SPECIES<- "PISY"
PISY_TRI_COR$VAR<- "TRI"

## for PCAB TRI
PCAB_TRI_COR<- data.frame(LOK=NA, CLIMA=NA, COR=NA, SIG=NA)
counter<- 1

for (b in c(1:ncol(PCAB_TRI))) {
  
  for (c in c(1:ncol(PCAB_JJA))) {
    
    Test<- cor.test(PCAB_TRI[,b], PCAB_JJA[,c])
    
    PCAB_TRI_COR[counter, "LOK"]<- colnames(PCAB_TRI)[b]
    PCAB_TRI_COR[counter, "CLIMA"]<- colnames(PCAB_JJA)[c]
    PCAB_TRI_COR[counter, "COR"]<- Test$estimate
    PCAB_TRI_COR[counter, "SIG"]<- Test$p.value
    
    counter<- counter+1
  }
}

PCAB_TRI_COR$SPECIES<- "PCAB"
PCAB_TRI_COR$VAR<- "TRI"

## for PISY NDVI
PISY_NDVI_COR<- data.frame(LOK=NA, CLIMA=NA, COR=NA, SIG=NA)
counter<- 1

for (b in c(1:ncol(PISY_NDVI))) {
  
  for (c in c(1:ncol(PISY_JJA))) {
    
    Test<- cor.test(PISY_NDVI[,b], PISY_JJA[,c])
    
    PISY_NDVI_COR[counter, "LOK"]<- colnames(PISY_NDVI)[b]
    PISY_NDVI_COR[counter, "CLIMA"]<- colnames(PISY_JJA)[c]
    PISY_NDVI_COR[counter, "COR"]<- Test$estimate
    PISY_NDVI_COR[counter, "SIG"]<- Test$p.value
    
    counter<- counter+1
  }
}

PISY_NDVI_COR$SPECIES<- "PISY"
PISY_NDVI_COR$VAR<- "NDVI"

## for PCAB NDVI
PCAB_NDVI_COR<- data.frame(LOK=NA, CLIMA=NA, COR=NA, SIG=NA)
counter<- 1

for (b in c(1:ncol(PCAB_NDVI))) {
  
  for (c in c(1:ncol(PCAB_JJA))) {
    
    Test<- cor.test(PCAB_NDVI[,b], PCAB_JJA[,c])
    
    PCAB_NDVI_COR[counter, "LOK"]<- colnames(PCAB_NDVI)[b]
    PCAB_NDVI_COR[counter, "CLIMA"]<- colnames(PCAB_JJA)[c]
    PCAB_NDVI_COR[counter, "COR"]<- Test$estimate
    PCAB_NDVI_COR[counter, "SIG"]<- Test$p.value
    
    counter<- counter+1
  }
}

PCAB_NDVI_COR$SPECIES<- "PCAB"
PCAB_NDVI_COR$VAR<- "NDVI"

# combining correlations together
COR<- rbind(PISY_TRI_COR, PISY_NDVI_COR, PCAB_TRI_COR, PCAB_NDVI_COR)
COR$TOPO<- substring(COR$LOK, 3,4)

# aggregating for topographic categories
COR_TOPO<-aggregate(COR$COR, list(COR$SPECIES, COR$VAR, COR$CLIMA, COR$TOPO), mean)
SD<-aggregate(COR$COR, list(COR$SPECIES, COR$VAR, COR$CLIMA, COR$TOPO), sd)
COR_TOPO$SD<- SD$x
colnames(COR_TOPO)<- c("SPECIES", "VAR", "CLIM", "TOPO", "COR", "SD")

##########################################################
#               Fig. S1 (Clima correlations)
##########################################################

order<- c("PL", "SS", "NS", "VA")

ggplot(COR_TOPO, aes(x=factor(TOPO,level=order), y=COR, fill=CLIM))+geom_bar(stat='identity', position=position_dodge())+
  facet_grid2(VAR~SPECIES, axes = "all")+
  geom_errorbar(aes(ymin=COR-SD, ymax=COR+SD), width=.2, position=position_dodge(.9))+
  geom_hline(yintercept=0, color = "black")+
  geom_hline(yintercept=0.346, color = "black",linetype="dashed")+
  geom_hline(yintercept=-0.346, color = "black",linetype="dashed")+
  labs(x= "", y= "Mean correlation")+
  coord_cartesian(ylim= c(-0.6,0.6))+
  
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legend
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  scale_fill_manual(values = c("#3242d1", "#5badf0","#f2e713","#f5a925", "#f0261a"))

ggsave("Grafy/Fig. S2 (Clima correlations).tiff", height = 150, width = 250, units = "mm", dpi = 300)

##############################################################
#    5. Krukal Walis of correlations between site categories
##############################################################

COR$Code<- paste(COR$SPECIES, COR$VAR, COR$CLIMA, sep = "_")
COR <- as.data.frame(lapply(COR, unlist))

Levels1<- as.data.frame(unique(COR$Code)); colnames(Levels1)<- "Code"

Kruskal_Wallis1<- data.frame(matrix(nrow = nrow(Levels1), ncol = 9)); rownames(Kruskal_Wallis1)<- Levels1$Code; colnames(Kruskal_Wallis1)<- c("PL_NS", "SS_NS", "VA_NS", "PL_PL", "SS_PL", "VA_PL", "PL_SS", "SS_SS", "VA_SS")

for (q in c(1:nrow(Levels1))) {
  
  Set1<- COR[COR$Code == Levels1[q, "Code"],]
  
  PostHoc<- pairwise.wilcox.test(Set1$COR,Set1$TOPO)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  Kruskal_Wallis1[q, c(2:10)]<- PH_P_val
  
}

Kruskal_Wallis1<- Kruskal_Wallis1[, -c(4,7,8)]

write.table(Kruskal_Wallis1, "Vystupy/1. KW_topo_cor.txt", row.names = T, col.names = T, sep = "\t")

##############################################################
#    6. Selection of nonconsecutive dry events
##############################################################

JJA_mean$row.names <- as.numeric(rownames(JJA_mean))
JJA_mean<- JJA_mean[, c(11, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]
JJA_mean[,c(3,4,8,9)]<- NULL
JJA_mean$PISY_SR<- JJA_mean$PISY_SR*(-1)
JJA_mean$PCAB_SR<- JJA_mean$PCAB_SR*(-1)

CPY_selection <- function(input.table, cycles = 4, column = 2){
  
  # Vytvorim tabulku pro ulozeni nalezenych pointeru
  result <- data.frame(YEAR = NA, Serie = colnames(input.table)[column], Order = c(1:cycles))
  
  for (i in c(1:cycles)){ # Proved pro zadany pocet cyklu = pocet nejvyznamnejsich pointeru, ktere je potreba najit  
    
    # Ktery rok je v datasetu tim nejvyznamnejsim pointerem?
    pointer <- input.table[min(input.table[, column]) == input.table[, column], 1]
    
    # Nejvyznamnejsi pointer odeberu z datasetu spolu s rokem, ktery po nem nasleduje nebo mu predchazi (pozor, pokud se jedna o posledni nebo prvni rok chronologie)
    if (max(input.table[,1]) > pointer & min(input.table[,1]) < pointer)  {input.table <- input.table[!(input.table[,1] %in% c(pointer - 1, pointer, pointer + 1)),]}
    else { if (max(input.table[,1]) == pointer) {input.table <- input.table[!(input.table[,1] %in% c(pointer - 1, pointer)),]}
      else { if (min(input.table[,1]) == pointer) {input.table <- input.table[!(input.table[,1] %in% c(pointer, pointer + 1)),]}}}
    
    # Ulozim hodnotu roku pro nalezeny pointer
    result[i, "YEAR"] <- pointer
    
    # Pred spustenim dalsiho cyklu radeji zaznam o poslednim pointeru vymazu (ale pripadne by se automaticky prepsal)
    rm(pointer)
  }
  
  return(result)
}


CPY<- data.frame(YEAR=NA, Serie=NA, Order=NA)

counter<- 1

JJA_mean_cpy<- JJA_mean[-c(32,33),]

for (j in c(2:7)) {
  
  Pointers<- CPY_selection(JJA_mean_cpy, 4, j)
  
  
  for (a in c(1:nrow(Pointers))) {
    
    CPY[counter, "YEAR"]<- Pointers[a, "YEAR"]
    CPY[counter, "Serie"]<- Pointers[a, "Serie"]
    CPY[counter, "Order"]<- Pointers[a, "Order"]
    
    counter<- counter+1
    
  }
}


CPY$SPECIES<- lapply(strsplit(as.character(CPY$Serie), "\\_"), "[",1)
CPY$CLIM<- lapply(strsplit(as.character(CPY$Serie), "\\_"), "[",2)
CPY <- as.data.frame(lapply(CPY, unlist))

write.table(CPY, "Vystupy/2. CPY.txt", row.names = F, col.names = T, sep = "\t")


### preparation of both variables data with SR CPY for ggplot
TRI_gg<- TRI
TRI_gg$YEAR<- rownames(TRI_gg)
TRI_gg<- gather(TRI_gg, "LOK", "VAL", 1:40)
TRI_gg[1:660,"SPECIES"]<- "PISY"
TRI_gg[661:1320,"SPECIES"]<- "PCAB"

NDVI_gg<- NDVI
NDVI_gg$YEAR<- rownames(NDVI_gg)
NDVI_gg<- gather(NDVI_gg, "LOK", "VAL", 1:40)
NDVI_gg[1:660,"SPECIES"]<- "PISY"
NDVI_gg[661:1320,"SPECIES"]<- "PCAB"

TRI_NDVI_gg<- rbind(TRI_gg, NDVI_gg)
TRI_NDVI_gg[1:1320, "VAR"]<- "TRI"
TRI_NDVI_gg[1321:2640, "VAR"]<- "NDVI"
TRI_NDVI_gg$TOPO<- substring(TRI_NDVI_gg$LOK,3,4)

##########################################################
#             Fig. 2 (Variables pointers)
##########################################################

SR<- JJA_mean[,c(1,4,7)]
SR<- gather(SR, "SPECIES", "SR", 2:3)
SR[1:33, 2]<- "PISY"
SR[34:66, 2]<- "PCAB"
SR$SR<- SR$SR*(-1)
SR<- rbind(SR, SR)
SR[1:66, "VAR"]<- "TRI"
SR[67:132, "VAR"]<- "NDVI"

min<- min(SR[SR$VAR=="TRI", "SR"])
max<- max(SR[SR$VAR=="TRI", "SR"])
new_min = 0.4
new_max = 1.5

SR[SR$VAR=="TRI", "SR_1"]<- ((SR[SR$VAR=="TRI", "SR"]- min) / (max-min) ) * (new_max - new_min) + new_min
SR[SR$VAR=="TRI", "SR_1"]<- ((SR[SR$VAR=="TRI", "SR"]- min) / (max-min) ) * (new_max - new_min) + new_min

min<- min(SR[SR$VAR=="NDVI", "SR"])
max<- max(SR[SR$VAR=="NDVI", "SR"])
new_min = 0.4
new_max = 0.9

SR[SR$VAR=="NDVI", "SR_1"]<- ((SR[SR$VAR=="NDVI", "SR"]- min) / (max-min) ) * (new_max - new_min) + new_min
SR[SR$VAR=="NDVI", "SR_1"]<- ((SR[SR$VAR=="NDVI", "SR"]- min) / (max-min) ) * (new_max - new_min) + new_min

SR$Code<- paste(SR$SPECIES, SR$VAR, SR$row.names, sep = "_")
TRI_NDVI_gg$Code<- paste(TRI_NDVI_gg$SPECIES, TRI_NDVI_gg$VAR, TRI_NDVI_gg$YEAR, sep = "_")

TRI_NDVI_gg<- merge(TRI_NDVI_gg, SR[,c(6,5)], by.x = "Code", by.y = "Code", all.x = T)

ggplot(TRI_NDVI_gg, aes(x=as.numeric(YEAR), y=VAL))+
  geom_line(aes(group=LOK, colour=TOPO))+
  geom_line(aes(x=as.numeric(YEAR), y=SR_1, color="SR"), size=1)+
  facet_grid2(VAR~SPECIES, axes = "all", scales = "free_y")+
  labs(x= "", y= "")+
  scale_color_manual(values = c("#f5a925", "#cccccc", "#004ca8", "#828282", "#005ce6"), breaks = c("SR", "NS", "PL", "SS", "VA"))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  geom_vline(xintercept = c(1994, 2003, 2006, 2015), color="#f5a925")+
  geom_text(x=1993, y=0.6, label="1994", size=3)+
  geom_text(x=2002, y=0.6, label="2003", size=3)+
  geom_text(x=2007, y=0.6, label="2006", size=3)+
  geom_text(x=2014, y=0.6, label="2015", size=3)+
  theme(strip.text  = element_text(size = 14))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legendy
  guides(col = guide_legend(nrow = 1)) ## počet řádků legendy

ggsave("Grafy/Fig. 2 (Variables pointers).tiff", height = 130, width = 220, units = "mm", dpi = 300)

##############################################################
#             7. Superposed epoch analysis (SEA)
##############################################################

PISY_CPY<- CPY[CPY$SPECIES=="PISY",]
PCAB_CPY<- CPY[CPY$SPECIES=="PCAB",]

Levels<- as.data.frame(unique(CPY$CLIM)); colnames(Levels)<- "Code"

# SEA for PISY TRI

PISY_TRI_SEA<- data.frame(lag=NA,	se=NA, se.unscaled=NA,	p=NA,	ci.95.lower=NA,	ci.95.upper=NA,	ci.99.lower=NA,	ci.99.upper=NA,	LOK=NA)

for (e in c(1:3)) {
  
  Set<- subset(PISY_CPY, subset = PISY_CPY$CLIM==Levels[e,"Code"])
  
  P_YEARS<- as.vector(Set$YEAR)
  
  for (f in c(1:ncol(PISY_TRI))) {
    
    Chrono<- as.data.frame(PISY_TRI[,f]); rownames(Chrono)<- rownames(PISY_TRI); colnames(Chrono)<- colnames(PISY_TRI)[f]
    
    SEA<- sea(Chrono, P_YEARS, 4)# possible to change number of considered years
    
    counter<- nrow(PISY_TRI_SEA)+1
    
    for (g in c(1:nrow(SEA))) {
      
      PISY_TRI_SEA[counter,"lag"]<- SEA[g,"lag"]
      PISY_TRI_SEA[counter,"se"]<- SEA[g,"se"]
      PISY_TRI_SEA[counter,"se.unscaled"]<- SEA[g,"se.unscaled"]
      PISY_TRI_SEA[counter,"p"]<- SEA[g,"p"]
      PISY_TRI_SEA[counter,"ci.95.lower"]<- SEA[g,"ci.95.lower"]
      PISY_TRI_SEA[counter,"ci.95.upper"]<- SEA[g,"ci.95.upper"]
      PISY_TRI_SEA[counter,"ci.99.lower"]<- SEA[g,"ci.99.lower"]
      PISY_TRI_SEA[counter,"ci.99.upper"]<- SEA[g,"ci.99.upper"]
      PISY_TRI_SEA[counter,"LOK"]<- colnames(Chrono)
      PISY_TRI_SEA[counter,"CLIMA"]<- Levels[e, "Code"]
      
      counter<- counter+1
      
    }
  }
}

PISY_TRI_SEA<- PISY_TRI_SEA[c(-1),]
PISY_TRI_SEA$SPECIES<- "PISY"
PISY_TRI_SEA$VAR<- "TRI"

# SEA for PCAB TRI

PCAB_TRI_SEA<- data.frame(lag=NA,	se=NA, se.unscaled=NA,	p=NA,	ci.95.lower=NA,	ci.95.upper=NA,	ci.99.lower=NA,	ci.99.upper=NA,	LOK=NA)

for (e in c(1:3)) {
  
  Set<- subset(PCAB_CPY, subset = PCAB_CPY$CLIM==Levels[e,"Code"])
  
  P_YEARS<- as.vector(Set$YEAR)
  
  for (f in c(1:ncol(PCAB_TRI))) {
    
    Chrono<- as.data.frame(PCAB_TRI[,f]); rownames(Chrono)<- rownames(PCAB_TRI); colnames(Chrono)<- colnames(PCAB_TRI)[f]
    
    SEA<- sea(Chrono, P_YEARS, 4)# possible to change number of considered years
    
    counter<- nrow(PCAB_TRI_SEA)+1
    
    for (g in c(1:nrow(SEA))) {
      
      PCAB_TRI_SEA[counter,"lag"]<- SEA[g,"lag"]
      PCAB_TRI_SEA[counter,"se"]<- SEA[g,"se"]
      PCAB_TRI_SEA[counter,"se.unscaled"]<- SEA[g,"se.unscaled"]
      PCAB_TRI_SEA[counter,"p"]<- SEA[g,"p"]
      PCAB_TRI_SEA[counter,"ci.95.lower"]<- SEA[g,"ci.95.lower"]
      PCAB_TRI_SEA[counter,"ci.95.upper"]<- SEA[g,"ci.95.upper"]
      PCAB_TRI_SEA[counter,"ci.99.lower"]<- SEA[g,"ci.99.lower"]
      PCAB_TRI_SEA[counter,"ci.99.upper"]<- SEA[g,"ci.99.upper"]
      PCAB_TRI_SEA[counter,"LOK"]<- colnames(Chrono)
      PCAB_TRI_SEA[counter,"CLIMA"]<- Levels[e, "Code"]
      
      counter<- counter+1
      
    }
  }
}

PCAB_TRI_SEA<- PCAB_TRI_SEA[c(-1),]
PCAB_TRI_SEA$SPECIES<- "PCAB"
PCAB_TRI_SEA$VAR<- "TRI"

# SEA for PISY NDVI

PISY_NDVI_SEA<- data.frame(lag=NA,	se=NA, se.unscaled=NA,	p=NA,	ci.95.lower=NA,	ci.95.upper=NA,	ci.99.lower=NA,	ci.99.upper=NA,	LOK=NA)

for (e in c(1:3)) {
  
  Set<- subset(PISY_CPY, subset = PISY_CPY$CLIM==Levels[e,"Code"])
  
  P_YEARS<- as.vector(Set$YEAR)
  
  for (f in c(1:ncol(PISY_NDVI))) {
    
    Chrono<- as.data.frame(PISY_NDVI[,f]); rownames(Chrono)<- rownames(PISY_NDVI); colnames(Chrono)<- colnames(PISY_NDVI)[f]
    
    SEA<- sea(Chrono, P_YEARS, 4)# possible to change number of considered years
    
    counter<- nrow(PISY_NDVI_SEA)+1
    
    for (g in c(1:nrow(SEA))) {
      
      PISY_NDVI_SEA[counter,"lag"]<- SEA[g,"lag"]
      PISY_NDVI_SEA[counter,"se"]<- SEA[g,"se"]
      PISY_NDVI_SEA[counter,"se.unscaled"]<- SEA[g,"se.unscaled"]
      PISY_NDVI_SEA[counter,"p"]<- SEA[g,"p"]
      PISY_NDVI_SEA[counter,"ci.95.lower"]<- SEA[g,"ci.95.lower"]
      PISY_NDVI_SEA[counter,"ci.95.upper"]<- SEA[g,"ci.95.upper"]
      PISY_NDVI_SEA[counter,"ci.99.lower"]<- SEA[g,"ci.99.lower"]
      PISY_NDVI_SEA[counter,"ci.99.upper"]<- SEA[g,"ci.99.upper"]
      PISY_NDVI_SEA[counter,"LOK"]<- colnames(Chrono)
      PISY_NDVI_SEA[counter,"CLIMA"]<- Levels[e, "Code"]
      
      counter<- counter+1
      
    }
  }
}

PISY_NDVI_SEA<- PISY_NDVI_SEA[c(-1),]
PISY_NDVI_SEA$SPECIES<- "PISY"
PISY_NDVI_SEA$VAR<- "NDVI"

# SEA for PCAB NDVI

PCAB_NDVI_SEA<- data.frame(lag=NA,	se=NA, se.unscaled=NA,	p=NA,	ci.95.lower=NA,	ci.95.upper=NA,	ci.99.lower=NA,	ci.99.upper=NA,	LOK=NA)

for (e in c(1:3)) {
  
  Set<- subset(PCAB_CPY, subset = PCAB_CPY$CLIM==Levels[e,"Code"])
  
  P_YEARS<- as.vector(Set$YEAR)
  
  for (f in c(1:ncol(PCAB_NDVI))) {
    
    Chrono<- as.data.frame(PCAB_NDVI[,f]); rownames(Chrono)<- rownames(PCAB_NDVI); colnames(Chrono)<- colnames(PCAB_NDVI)[f]
    
    SEA<- sea(Chrono, P_YEARS, 4) # possible to change number of considered years
    
    counter<- nrow(PCAB_NDVI_SEA)+1
    
    for (g in c(1:nrow(SEA))) {
      
      PCAB_NDVI_SEA[counter,"lag"]<- SEA[g,"lag"]
      PCAB_NDVI_SEA[counter,"se"]<- SEA[g,"se"]
      PCAB_NDVI_SEA[counter,"se.unscaled"]<- SEA[g,"se.unscaled"]
      PCAB_NDVI_SEA[counter,"p"]<- SEA[g,"p"]
      PCAB_NDVI_SEA[counter,"ci.95.lower"]<- SEA[g,"ci.95.lower"]
      PCAB_NDVI_SEA[counter,"ci.95.upper"]<- SEA[g,"ci.95.upper"]
      PCAB_NDVI_SEA[counter,"ci.99.lower"]<- SEA[g,"ci.99.lower"]
      PCAB_NDVI_SEA[counter,"ci.99.upper"]<- SEA[g,"ci.99.upper"]
      PCAB_NDVI_SEA[counter,"LOK"]<- colnames(Chrono)
      PCAB_NDVI_SEA[counter,"CLIMA"]<- Levels[e, "Code"]
      
      counter<- counter+1
      
    }
  }
}

PCAB_NDVI_SEA<- PCAB_NDVI_SEA[c(-1),]
PCAB_NDVI_SEA$SPECIES<- "PCAB"
PCAB_NDVI_SEA$VAR<- "NDVI"

##########################################################
#             Fig. 3 (SEA)
##########################################################

# gathering SEA together
SEA<- rbind(PISY_TRI_SEA, PCAB_TRI_SEA, PISY_NDVI_SEA, PCAB_NDVI_SEA)
SEA$TOPO<- substring(SEA$LOK,3,4)

# preparing for plot
SEA_SR<- SEA[SEA$CLIMA=="SR",]

ggplot(SEA_SR, aes(x=lag, y=se, colour=TOPO))+geom_point(aes(alpha = p<0.05), size = 3)+
  facet_grid2(VAR~SPECIES, axes = "all")+
  scale_color_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
  scale_x_continuous(breaks = seq(-4, 4, 1))+
  coord_cartesian(ylim= c(-2, 2))+
  scale_alpha_discrete(range = c(0.05, 1))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t= -10))+ ## umístění legendy dolu
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+ ## odstranění šedého podbaevení legend
  guides(alpha = "none")

ggsave("Grafy/Fig. 3 (SEA).tiff", height = 110, width = 170, units = "mm", dpi = 300)


##############################################################
#         8. Selection of TRI and NDVI for LMEM
##############################################################

counter<-1

List_var<- as.data.frame(c("TRI", "NDVI")); colnames(List_var)<- "Code"

YEARS<- as.data.frame(c(1994, 1995, 1996, 1997, 1998, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010)); colnames(YEARS)<- "YEAR"

Model_data<- data.frame(SPECIES=NA, CODE=NA, VAL=NA, LOC=NA, YEAR=NA, EXT=NA, EXT_4=NA)

for (radek in c(1:nrow(List_var))) {
  
  VAR<- read.table(paste("C:/Users/jirka/Desktop/Pointers/Vystupy/", List_var[radek,"Code"], ".txt", sep = ""), check.names=FALSE, dec=".", header = T)
  
  PISY_var<- VAR[,1:20]
  PCAB_var<- VAR[,21:40]
  
  for (y in c(1:nrow(YEARS))) {
    
    year<- YEARS[y, "YEAR"]
    
    for (a in c(1:ncol(PISY_var))) {
      
      lag1<-mean(PISY_JJA$SR)-PISY_JJA[rownames(PISY_JJA) %in% as.numeric(year+1),"SR"]
      lag2<-mean(PISY_JJA$SR)-PISY_JJA[rownames(PISY_JJA) %in% as.numeric(year+2),"SR"]
      lag3<-mean(PISY_JJA$SR)-PISY_JJA[rownames(PISY_JJA) %in% as.numeric(year+3),"SR"]
      lag4<-mean(PISY_JJA$SR)-PISY_JJA[rownames(PISY_JJA) %in% as.numeric(year+4),"SR"]
      
      Model_data[counter, "SPECIES"]<- "PISY"
      Model_data[counter, "CODE"]<-  List_var[radek,"Code"]
      Model_data[counter, "LOC"]<- colnames(PISY_var)[a]
      Model_data[counter, "YEAR"]<- as.numeric(year)
      Model_data[counter, "VAL"]<- PISY_var[rownames(PISY_var) %in% year,a]
      Model_data[counter, "EXT"]<- mean(PISY_JJA$SR)-PISY_JJA[rownames(PISY_JJA)%in% year,"SR"]
      Model_data[counter, "EXT_4"]<- mean(lag1, lag2, lag3, lag4)
      
      counter<- counter+1
      
    }
    
    for (a in c(1:ncol(PCAB_var))) {
      
      lag1<-mean(PCAB_JJA$SR)-PCAB_JJA[rownames(PCAB_JJA) %in% as.numeric(year+1),"SR"]
      lag2<-mean(PCAB_JJA$SR)-PCAB_JJA[rownames(PCAB_JJA) %in% as.numeric(year+2),"SR"]
      lag3<-mean(PCAB_JJA$SR)-PCAB_JJA[rownames(PCAB_JJA) %in% as.numeric(year+3),"SR"]
      lag4<-mean(PCAB_JJA$SR)-PCAB_JJA[rownames(PCAB_JJA) %in% as.numeric(year+4),"SR"]
      
      Model_data[counter, "SPECIES"]<- "PCAB"
      Model_data[counter, "CODE"]<-  List_var[radek,"Code"]
      Model_data[counter, "LOC"]<- colnames(PCAB_var)[a]
      Model_data[counter, "YEAR"]<- as.numeric(year)
      Model_data[counter, "VAL"]<- PCAB_var[rownames(PCAB_var) %in% year,a]
      Model_data[counter, "EXT"]<- mean(PCAB_JJA$SR)-PCAB_JJA[rownames(PCAB_JJA)%in% year,"SR"]
      Model_data[counter, "EXT_4"]<- mean(lag1, lag2, lag3, lag4)
      
      counter<- counter+1
      
    }
  }
}


Model_data$TOPO<- substring(Model_data$LOC,3,4)
Model_data<- merge(Model_data, List_loc[,c(2,7)], by.x = "LOC", by.y = "Code", all.x = T)
Model_data$CODE<- paste(Model_data$CODE, Model_data$SPECIES, sep = "_")

##############################################################
#         9. Calculation of linear mixed effect models
##############################################################

Levels3<- as.data.frame(unique(Model_data$CODE)); colnames(Levels3)<- "Code"

LMEM<- data.frame(DATA=NA, Model=NA, Rm=NA, Rc=NA, pR=NA, t_Intercept=NA, t_EXT=NA, t_EXT_4=NA, t_DBH=NA, p_EXT=NA, p_EXT_4=NA, p_DBH=NA, RAND=NA)

options(na.action = "na.fail")

counter<- 1

for (m in c(1:nrow(Levels3))) {
  
  Set4<- Model_data[Model_data$CODE==Levels3[m, "Code"],]
  Set4<- na.omit(Set4)
  
  ## full model
  Mixed_model<-lmer(VAL~EXT+EXT_4+DBH+(1+EXT+EXT_4+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model))
  dr<- dredge(Mixed_model)
  Sum<- summary(Mixed_model)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "Full"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model), Set4$VAL)^2
  LMEM[counter, c(6:9)]<- Sum
  LMEM[counter, c(11, 10, 12)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  ## just variable intercept
  Mixed_model1<-lmer(VAL~EXT+EXT_4+DBH+(1|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model1))
  dr<- dredge(Mixed_model1)
  Sum<- summary(Mixed_model1)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model1, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "Var intercept"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model1), Set4$VAL)^2
  LMEM[counter, c(6:9)]<- Sum
  LMEM[counter, c(11, 10, 12)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  ## just variable slope
  Mixed_model2<-lmer(VAL~EXT+(1+EXT+EXT_4+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model2))
  dr<- dredge(Mixed_model2)
  Sum<- summary(Mixed_model2)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model2, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "Var slope"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model2), Set4$VAL)^2
  LMEM[counter, c(6:9)]<- Sum
  LMEM[counter, c(11, 10, 12)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  ## without EXT_4
  Mixed_model3<-lmer(VAL~EXT+DBH+(1+EXT+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model3))
  dr<- dredge(Mixed_model3)
  Sum<- summary(Mixed_model3)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model3, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "without EXT_4"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model3), Set4$VAL)^2
  LMEM[counter, c(6,7,9)]<- Sum
  LMEM[counter, c(10, 12)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  ## without EXT
  Mixed_model4<-lmer(VAL~EXT_4+DBH+(1+EXT_4+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model4))
  dr<- dredge(Mixed_model4)
  Sum<- summary(Mixed_model4)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model4, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "without EXT"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model4), Set4$VAL)^2
  LMEM[counter, c(6,8,9)]<- Sum
  LMEM[counter, c(11, 12)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  
  ## without DBH
  Mixed_model5<-lmer(VAL~EXT+EXT_4+(1+EXT+EXT_4|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model5))
  dr<- dredge(Mixed_model5)
  Sum<- summary(Mixed_model5)
  Sum<- Sum[["coefficients"]]
  Sum<-Sum[,4]
  im<- sw(dr)
  rand<- ranova(Mixed_model5, reduce.terms = F)
  
  LMEM[counter, "DATA"]<- Levels3[m, "Code"]
  LMEM[counter, "Model"]<- "without DBH"
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model5), Set4$VAL)^2
  LMEM[counter, c(6:8)]<- Sum
  LMEM[counter, c(11, 10)]<- im
  LMEM[counter, "RAND"]<- rand[2,6]
  
  counter<- counter+1
  
  
  GG_EXTR<- ggplot(Set4, aes(EXT, VAL, color=TOPO))+geom_point()+
    geom_smooth(method=lm, aes(fill=TOPO), alpha=0.2)+
    scale_color_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
    scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
    theme(legend.position = "none")
  
  GG_DBH<- ggplot(Set4, aes(DBH, VAL, color=TOPO))+geom_point()+
    geom_smooth(method=lm, aes(fill=TOPO), alpha=0.2)+
    scale_color_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
    scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
    theme(legend.position = "none")
  
  plot_grid(GG_EXTR, GG_DBH, nrow=2, ncol=1, labels = c("", ""), label_size = 12)
  #ggsave(paste("C:/Users/jirka/Desktop/Pointers/Grafy/Model_", Levels3[m, "Code"], ".tiff", sep = ""), height = 150, width = 150, units = "mm", dpi = 300)
  
}

LMEM$t_Intercept<- NULL

write.table(LMEM, "Vystupy/3. LMEM.txt", row.names = F, col.names = T, sep = "\t")


##############################################################
#             7. Resilience indices calculation
##############################################################

# for TRI
TRI_res<- res.comp(TRI)
TRI_RS<- as.data.frame(TRI_res$resil)
TRI_RT<- as.data.frame(TRI_res$resist)
TRI_RC<- as.data.frame(TRI_res$recov)

write.table(TRI_RS, "Vystupy/TRI_RS.txt", sep="\t")
write.table(TRI_RT, "Vystupy/TRI_RT.txt", sep="\t")
write.table(TRI_RC, "Vystupy/TRI_RC.txt", sep="\t")

# for NDVI
NDVI_res<- res.comp(NDVI)
NDVI_RS<- as.data.frame(NDVI_res$resil)
NDVI_RT<- as.data.frame(NDVI_res$resist)
NDVI_RC<- as.data.frame(NDVI_res$recov)

write.table(NDVI_RS, "Vystupy/NDVI_RS.txt", sep="\t")
write.table(NDVI_RT, "Vystupy/NDVI_RT.txt", sep="\t")
write.table(NDVI_RC, "Vystupy/NDVI_RC.txt", sep="\t")

## reorganization of data
List_indices<- as.data.frame(c("TRI_RS", "TRI_RT", "TRI_RC", "NDVI_RS", "NDVI_RT", "NDVI_RC")); colnames(List_indices)<- "Code"

for (ind in c(1:nrow(List_indices))) {
  
  Index <- read.table(paste("C:/Users/jirka/Desktop/Pointers/Vystupy/", List_indices[ind, "Code"], ".txt", sep = ""), dec=".", check.names = F)
  
  Index$YEAR<- rownames(Index)
  Index$Index <- List_indices[ind, "Code"]
  
  Index<- gather(Index, "LOC", "VAL", 1:40)
  Index$TOPO<- substring(Index$LOC, 3,4)
  
  Index$VAR<- lapply(strsplit(as.character(Index$Index), "\\_"), "[",1)
  Index$Index<- lapply(strsplit(as.character(Index$Index), "\\_"), "[",2)
  Index <- as.data.frame(lapply(Index, unlist))
  
  Index[1:660, "SPECIES"]<- "PISY"
  Index[661:1320, "SPECIES"]<- "PCAB"
  Index<- na.omit(Index)
  
  write.table(Index, paste("C:/Users/jirka/Desktop/Pointers/Vystupy/", List_indices[ind, "Code"], ".txt", sep = ""), dec=".", sep = "\t")
  
}

TRI_RS<- read.table("Vystupy/TRI_RS.txt", sep="\t")
TRI_RT<- read.table("Vystupy/TRI_RT.txt", sep="\t")
TRI_RC<- read.table("Vystupy/TRI_RC.txt", sep="\t")

NDVI_RS<- read.table("Vystupy/NDVI_RS.txt", sep="\t")
NDVI_RT<- read.table("Vystupy/NDVI_RT.txt", sep="\t")
NDVI_RC<- read.table("Vystupy/NDVI_RC.txt", sep="\t")

RES_data<- rbind(TRI_RS, TRI_RT, TRI_RC, NDVI_RS, NDVI_RT, NDVI_RC)

### selection of drought spells
RES_data<- RES_data[RES_data$YEAR %in% c(1994, 2003, 2006),]

##########################################################
#             Fig. S3 (Resil boxplots)
##########################################################

RS_gg<- ggplot(RES_data[RES_data$Index=="RS",], aes(x=as.character(YEAR), y=VAL, fill=TOPO))+geom_boxplot()+
  facet_grid2(VAR~SPECIES, axes = "all")+
  scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
  labs(x= "", y= "Resilience")+
  coord_cartesian(ylim= c(0.5, 1.8))+
  theme(legend.position = "none", legend.margin=margin(t= -25))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black")) ## barva os

RT_gg<- ggplot(RES_data[RES_data$Index=="RT",], aes(x=as.character(YEAR), y=VAL, fill=TOPO))+geom_boxplot()+
  facet_grid2(VAR~SPECIES, axes = "all")+
  scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
  labs(x= "", y= "Resistance")+
  coord_cartesian(ylim= c(0.5, 1.8))+
  theme(legend.position = "none", legend.margin=margin(t= -25))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black")) ## barva os

RC_gg<- ggplot(RES_data[RES_data$Index=="RC",], aes(x=as.character(YEAR), y=VAL, fill=TOPO))+geom_boxplot()+
  facet_grid2(VAR~SPECIES, axes = "all")+
  scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
  labs(x= "", y= "Recovery")+
  coord_cartesian(ylim= c(0.5, 1.8))+
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black")) ## barva os

plot_grid(RS_gg, RT_gg, RC_gg, nrow=3, ncol=1, labels = c("A", "B", "C"), label_size = 16, label_x = 0.08)


ggsave("Grafy/Fig. S3 (Resil boxplots) .tiff", height = 350, width = 200, units = "mm", dpi = 300)

##############################################################
#       11. Kruskal-Wallis test of resilience for TOPO
##############################################################

RES_data$CODE<- paste(RES_data$SPECIES, RES_data$VAR, RES_data$Index, RES_data$YEAR, sep = "_")

Levels2<- as.data.frame(unique(RES_data$CODE)); colnames(Levels2)<- "CODE"

KW_res_topo<- as.data.frame(matrix(nrow = nrow(Levels2), ncol = 9)); rownames(KW_res_topo)<- Levels2$CODE; colnames(KW_res_topo)<- c("NS-PL", "NS-SS", "NS-VA", "PL-PL", "PL-SS", "PL-VA", "SS-PL", "SS-SS", "SS-VA")

for (h in c(1:nrow(Levels2))) {
  
  Set3<- subset(RES_data, subset = RES_data$CODE==Levels2[h, "CODE"])
  
  PostHoc<- pairwise.wilcox.test(Set3$VAL,Set3$TOPO)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  KW_res_topo[h,] <-PH_P_val
  
}

KW_res_topo<- KW_res_topo[,-c(4,7,8)]

write.table(KW_res_topo, "Vystupy/5. KW_topo_res.txt", sep="\t")


##############################################################
#       11. Kruskal-Wallis test of resilience for years
##############################################################

RES_data$Code<- paste(RES_data$SPECIES, RES_data$VAR, RES_data$Index, RES_data$TOPO, sep = "_")

Levels_res<- as.data.frame(unique(RES_data$Code)); colnames(Levels_res)<- "Code"

KW_res<- as.data.frame(matrix(ncol = 4, nrow = nrow(Levels_res))); rownames(KW_res)<- Levels_res$Code; colnames(KW_res)<- c("1994-2003", "1994-2006", "2003-2003", "2003-2006")

for (d in c(1:nrow(Levels_res))) {
  
  Set2<- subset(RES_data, subset = RES_data$Code==Levels_res[d, "Code"])
  
  PostHoc<- pairwise.wilcox.test(Set2$VAL,Set2$YEAR)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  KW_res[d,] <-PH_P_val
  
}

KW_res<- KW_res[,-c(3)]

write.table(KW_res, "Vystupy/4. KW_year_res.txt", sep="\t")

##############################################################
#                 Fig. S4 (Biomass)
##############################################################

TRI<- TRI[rownames(TRI) %in% c(1994, 1995, 1996, 1997, 1998, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010),]
NDVI<- NDVI[rownames(NDVI) %in% c(1994, 1995, 1996, 1997, 1998, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010),]

TRI$YEAR<- rownames(TRI)
NDVI$YEAR<- rownames(NDVI)

TRI<- gather(TRI, "LOC", "VAL", 1:40)
NDVI<- gather(NDVI, "LOC", "VAL", 1:40)

TRI$VAR<- "TRI"
NDVI$VAR<- "NDVI"

TRI[1:260, "SPECIES"]<- "PISY"
TRI[261:520, "SPECIES"]<- "PCAB"

NDVI[1:260, "SPECIES"]<- "PISY"
NDVI[261:520, "SPECIES"]<- "PCAB"

Data<- rbind(TRI, NDVI)
Data$TOPO<- substring(Data$LOC, 3,4)

ggplot(Data, aes(x=TOPO, y=VAL, fill=TOPO))+geom_boxplot()+
  facet_grid2(VAR~SPECIES, scales = "free", axes = "all")+
  scale_fill_manual(values = c("#cccccc", "#004ca8", "#828282", "#005ce6"))+
  labs(x= "", y= "")+
  theme(legend.position = "bottom", legend.margin=margin(t= -15))+
  theme(strip.text  = element_text(size = 12))+ ## velikost nadpisů
  theme(axis.text.x = element_text(size = 12))+ ## velikost názvu osy x
  theme(axis.text.y = element_text(size = 12))+ ## velikost názvu osy y
  theme(axis.title.x = element_text(size = 15))+ ## velikost popisků osy x
  theme(axis.title.y = element_text(size = 15))+ ## velikost popisků osy y
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(legend.title = element_blank())+ ## odstranění názvu legendy
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))+
  theme(axis.line = element_line(colour = "black")) ## barva os

ggsave("Grafy/Fig. S4 (Biomass).tiff", height = 110, width = 170, units = "mm", dpi = 300)

Data$Code<- paste(Data$SPECIES, Data$VAR, sep = "_")

Levels4<- as.data.frame(unique(Data$Code)); colnames(Levels4)<- "Code"

for (i in c(1:4)) {
  
  Set5<- subset(Data, subset = Data$Code==Levels4[i, "Code"])
  
  AOV<- aov(Set5$VAL~Set5$TOPO)
  
  PH<- TukeyHSD(AOV)
  PH<- as.data.frame(PH[["Set5$TOPO"]])
  
  Levels4[i, c(2:7)]<- PH$`p adj`
  
}

colnames(Levels4)<- c("Code", "PL-NS", "SS-NS", "VA-NS", "SS-PL", "VA-PL", "VA-SS")

write.table(Levels4, "Vystupy/6. AOV_biomass.txt", row.names = F, col.names = T, sep = "\t")


##############################################################
#                 Fig. S1 (NDVI sensitivity)
##############################################################

NDVI_MS<- read.table("Data/NDVI_MS.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8", check.names = F, row.names = 1)
NDVI_274<- read.table("Data/NDVI_274.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8", check.names = F, row.names = 1)
NDVI_JA<- read.table("Data/NDVI_JJA.txt", header=TRUE, sep="\t", na.strings="NA", dec=",", strip.white=TRUE, encoding = "UTF-8", check.names = F, row.names = 1)

COR<- as.data.frame(matrix(nrow = 40, ncol = 3)); rownames(COR)<- colnames(NDVI_274); colnames(COR)<- c("MS-274", "MS-JA", "274-JA")

for (i in c(1:40)) {
  
  COR[i, 1]<- cor(NDVI_MS[,i], NDVI_274[,i])
  COR[i, 2]<- cor(NDVI_MS[,i], NDVI_JA[,i])
  COR[i, 3]<- cor(NDVI_274[,i], NDVI_JA[,i])
  
}



COR[1:20, "SPECIES"]<- "PISY"
COR[21:40, "SPECIES"]<- "PCAB"

COR$LOC<- rownames(COR)

COR<- gather(COR, "Pair", "COR", 1:3)

PISY<- ggplot(COR[COR$SPECIES=="PISY",], aes(x=LOC, y=COR, fill=Pair))+
  geom_bar(stat='identity', position=position_dodge())+
  labs(x= "", y= "Correlation")+
  coord_cartesian(ylim= c(0, 1))+
  scale_fill_manual(values= c("#80d45f", "#759ffa", "#f27777"))+
  theme(strip.text  = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 14))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")+
  theme(legend.title = element_blank())+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))


PCAB<- ggplot(COR[COR$SPECIES=="PCAB",], aes(x=LOC, y=COR, fill=Pair))+
  geom_bar(stat='identity', position=position_dodge())+
  labs(x= "", y= "Correlation")+
  coord_cartesian(ylim= c(0, 1))+
  scale_fill_manual(values= c("#80d45f", "#759ffa", "#f27777"))+
  theme(strip.text  = element_text(size = 16))+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(legend.text = element_text(size = 14))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+
  theme(legend.title = element_blank())+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

plot_grid(PISY, PCAB, nrow=2, ncol=1, labels = c("PISY", "PCAB"), label_size = 12, label_x = 0.08)

ggsave("Grafy/Fig. S1 (NDVI sensitivity).tiff", height = 150, width = 170, units = "mm", dpi = 300)


##############################################################
#                       Climadiagrams
##############################################################

Climadiagramy<-read.table("Data/Climadiagrams.txt", check.names=FALSE, dec=",", header = T)

### přepočty pro druhou osu
ylim.prim<- c(-2, 55)
ylim.sec<- c(0, 110)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

Climadiagramy$Species_f<- factor(Climadiagramy$Species, levels = c("PISY", "PCAB"))

order<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")

ggplot(Climadiagramy, aes(x=factor(Month, level=order), group=1))+
  geom_bar(aes(y=a+precipitation*b), fill="#527af2", stat='identity')+
  geom_line(aes(y=temperatures), color="#fa051e",  stat='identity', size=2)+
  scale_y_continuous("Temperatures (°C)", sec.axis = sec_axis(~ (. - a)/b, name = "Precipitation (mm)"))+
  facet_wrap(~Species_f)+
  labs(x="")+
  coord_cartesian(ylim= c(-2,60))+
  theme(strip.text  = element_text(size = 17))+ ## velikost nadpisů
  theme(axis.text.y = element_text(size = 14))+ ## velikost názvu osy y
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1))+ ## velikost názvu osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black")) ## barva os

ggsave("Grafy/Climadiagrams.png", height = 100, width = 200, units = "mm", dpi = 300)
