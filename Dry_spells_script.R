#############################################################################################
#           Responses of stem and leaf biomass of temperate conifers to dry spells 
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
library(multcompView)
library(car)

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

########################################################################
#    Table S2 (Wilcox of correlations between site categories)
########################################################################

COR$Code<- paste(COR$SPECIES, COR$VAR, COR$CLIMA, sep = "_")
COR <- as.data.frame(lapply(COR, unlist))

Levels1<- as.data.frame(unique(COR$Code)); colnames(Levels1)<- "Code"

Wilcox1<- data.frame(matrix(nrow = nrow(Levels1), ncol = 9)); rownames(Wilcox1)<- Levels1$Code; colnames(Wilcox1)<- c("PL_NS", "SS_NS", "VA_NS", "PL_PL", "SS_PL", "VA_PL", "PL_SS", "SS_SS", "VA_SS")

for (q in c(1:nrow(Levels1))) {
  
  Set1<- COR[COR$Code == Levels1[q, "Code"],]
  
  PostHoc<- pairwise.wilcox.test(Set1$COR,Set1$TOPO)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  Wilcox1[q, c(2:10)]<- PH_P_val
  
}

Wilcox1<- Wilcox1[, -c(4,7,8)]

write.table(Wilcox1, "Vystupy/Table S2.txt", row.names = T, col.names = T, sep = "\t")

##############################################################
#    5. Selection of nonconsecutive dry events
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
  scale_color_manual(values = c("#f5a925", "#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"), breaks = c("SR", "NS", "PL", "SS", "VA"))+
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
#             6. Superposed epoch analysis (SEA)
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
  labs(x= "Lag (years)", y= " Scaled anomaly")+
  scale_color_manual(values = c("#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"))+
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
#         7. Selection of TRI and NDVI for LMEM
##############################################################

Event1<- c(1994, 1995, 1996, 1997, 1998)
Event2<- c(2003, 2004, 2005, 2006, 2007)
Event3<- c(2006, 2007, 2008, 2009, 2010)
Event4<- c(2015, 2016, 2017, NA, NA)

## for Supplementary variant Lag 1, 2, 3
#Event1<- c(1994, 1995, 1996)
#Event2<- c(2003, 2004, 2005)
#Event3<- c(2006, 2007, 2008)
#Event4<- c(2015, 2016, 2017)

Events<-as.data.frame(cbind(Event1, Event2, Event3, Event4))

Model_data<- data.frame(SPECIES=NA, VAR=NA, LOC=NA ,YEAR=NA, VAL=NA, SEV0=NA, SEV=NA, Event=NA, Lag=NA)

counter<-1

for (i in c(1:4)) {
  
  years<- Events[,i]
  years<- na.omit(years)
  
  for (j in c(1:length(years))) {
    
    year<- years[j]
    
    for (a in c(1:20)) {
      
      Model_data[counter, "SPECIES"]<- "PISY"
      Model_data[counter, "VAR"]<- "TRI"
      Model_data[counter, "LOC"]<- colnames(PISY_TRI)[a]
      Model_data[counter, "YEAR"]<- year
      Model_data[counter, "VAL"]<- PISY_TRI[rownames(PISY_TRI) %in% year,a]
      Model_data[counter, "SEV0"]<- PISY_JJA[rownames(PISY_JJA)%in% years[1],"SR"]-mean(PISY_JJA$SR)
      Model_data[counter, "SEV"]<- PISY_JJA[rownames(PISY_JJA)%in% year[1],"SR"]-mean(PISY_JJA$SR)
      Model_data[counter, "Event"]<- colnames(Events)[i]
      Model_data[counter, "Lag"]<- paste("Lag",j-1, sep = "")
      
      counter<- counter+1
      
      Model_data[counter, "SPECIES"]<- "PISY"
      Model_data[counter, "VAR"]<- "NDVI"
      Model_data[counter, "LOC"]<- colnames(PISY_NDVI)[a]
      Model_data[counter, "YEAR"]<- year
      Model_data[counter, "VAL"]<- PISY_NDVI[rownames(PISY_NDVI) %in% year,a]
      Model_data[counter, "SEV0"]<- PISY_JJA[rownames(PISY_JJA)%in% years[1],"SR"]-mean(PISY_JJA$SR)
      Model_data[counter, "SEV"]<- PISY_JJA[rownames(PISY_JJA)%in% year[1],"SR"]-mean(PISY_JJA$SR)
      Model_data[counter, "Event"]<- colnames(Events)[i]
      Model_data[counter, "Lag"]<- paste("Lag",j-1, sep = "")
      
      counter<- counter+1
      
      Model_data[counter, "SPECIES"]<- "PCAB"
      Model_data[counter, "VAR"]<- "TRI"
      Model_data[counter, "LOC"]<- colnames(PCAB_TRI)[a]
      Model_data[counter, "YEAR"]<- year
      Model_data[counter, "VAL"]<- PCAB_TRI[rownames(PCAB_TRI) %in% year,a]
      Model_data[counter, "SEV0"]<- PCAB_JJA[rownames(PCAB_JJA)%in% years[1],"SR"]-mean(PCAB_JJA$SR)
      Model_data[counter, "SEV"]<- PCAB_JJA[rownames(PCAB_JJA)%in% year[1],"SR"]-mean(PCAB_JJA$SR)
      Model_data[counter, "Event"]<- colnames(Events)[i]
      Model_data[counter, "Lag"]<- paste("Lag",j-1, sep = "")
      
      counter<- counter+1
      
      Model_data[counter, "SPECIES"]<- "PCAB"
      Model_data[counter, "VAR"]<- "NDVI"
      Model_data[counter, "LOC"]<- colnames(PCAB_NDVI)[a]
      Model_data[counter, "YEAR"]<- year
      Model_data[counter, "VAL"]<- PCAB_NDVI[rownames(PCAB_NDVI) %in% year,a]
      Model_data[counter, "SEV0"]<- PCAB_JJA[rownames(PCAB_JJA)%in% years[1],"SR"]-mean(PCAB_JJA$SR)
      Model_data[counter, "SEV"]<- PCAB_JJA[rownames(PCAB_JJA)%in% year[1],"SR"]-mean(PCAB_JJA$SR)
      Model_data[counter, "Event"]<- colnames(Events)[i]
      Model_data[counter, "Lag"]<- paste("Lag",j-1, sep = "")
      
      counter<- counter+1
      
    }
  }
}

Model_data<- merge(Model_data, List_loc[,c(2,7)], by.x = "LOC", by.y = "Code", all.x = T)
Model_data$TOPO<- substring(Model_data$LOC, 3,4)
Model_data$CODE<- paste(Model_data$SPECIES,Model_data$VAR, sep = " ")

##############################################################
#         8. Calculation of linear mixed effect models
##############################################################

Levels3<- as.data.frame(unique(Model_data$CODE)); colnames(Levels3)<- "CODE"

LMEM<- data.frame(DATA=NA, Model=NA, AIC=NA, Rm=NA, Rc=NA, pR=NA, t_SEV0=NA, t_SEV=NA, t_DBH=NA, p_SEV0=NA, p_SEV=NA, p_DBH=NA, TOPO=NA, v_SEV0=NA, v_SEV=NA, v_DBH=NA)

options(na.action = "na.fail")

counter<-1

for (i in c(1:nrow(Levels3))) {
  
  Set4<- subset(Model_data, subset = Model_data$CODE==Levels3[i, "CODE"])
  Set4$SEV0 <- ordered(Set4$SEV0)
  
  #######################################################################################################
  ##                            Full models
  #######################################################################################################
  
  Mixed_model<-lmer(VAL~SEV0+SEV+DBH+(1+SEV0+SEV+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model))
  #dr<- dredge(Mixed_model)
  #im<- as.data.frame(sw(dr))
  Sum<- summary(Mixed_model)
  Sum<- as.data.frame(Sum[["coefficients"]])
  rand<- ranova(Mixed_model, reduce.terms = F)
  VIF<- as.data.frame(vif(Mixed_model))
  qqnorm(resid(Mixed_model))
  
  LMEM[counter, "DATA"]<- Levels3[i, "CODE"]
  LMEM[counter, "Model"]<- "Full model"
  LMEM[counter, "AIC"]<- AIC(Mixed_model)
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model), Set4$VAL)^2
  LMEM[counter, "t_SEV0"]<- Sum["SEV0.L",4]
  LMEM[counter, "t_SEV"]<- Sum["SEV",4]
  LMEM[counter, "t_DBH"]<- Sum["DBH",4]
  LMEM[counter, "p_SEV0"]<- Sum["SEV0.L", 5]
  LMEM[counter, "p_SEV"]<- Sum["SEV", 5]
  LMEM[counter, "p_DBH"]<- Sum["DBH", 5]
  LMEM[counter, "TOPO"]<- rand[2,6]
  LMEM[counter, "v_SEV0"]<- VIF["SEV0", 1]
  LMEM[counter, "v_SEV"]<- VIF["SEV", 1]
  LMEM[counter, "v_DBH"]<- VIF["DBH", 1]
  
  counter<- counter+1
  
  #######################################################################################################
  ##                            Without SEV0
  #######################################################################################################
  
  Mixed_model1<-lmer(VAL~SEV+DBH+(1+SEV+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model1))
  #dr<- dredge(Mixed_model1)
  #im<- as.data.frame(sw(dr))
  Sum<- summary(Mixed_model1)
  Sum<- as.data.frame(Sum[["coefficients"]])
  rand<- ranova(Mixed_model1, reduce.terms = F)
  VIF<- as.data.frame(vif(Mixed_model1))
  qqnorm(resid(Mixed_model1))
  
  LMEM[counter, "DATA"]<- Levels3[i, "CODE"]
  LMEM[counter, "Model"]<- "without SEV0"
  LMEM[counter, "AIC"]<- AIC(Mixed_model1)
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model1), Set4$VAL)^2
  #LMEM[counter, "t_SEV0"]<- Sum["SEV0.L",4]
  LMEM[counter, "t_SEV"]<- Sum["SEV",4]
  LMEM[counter, "t_DBH"]<- Sum["DBH",4]
  #LMEM[counter, "p_SEV0"]<- Sum["SEV0.L", 5]
  LMEM[counter, "p_SEV"]<- Sum["SEV", 5]
  LMEM[counter, "p_DBH"]<- Sum["DBH", 5]
  LMEM[counter, "TOPO"]<- rand[2,6]
  #LMEM[counter, "v_SEV0"]<- VIF["SEV0", 1]
  LMEM[counter, "v_SEV"]<- VIF["SEV", 1]
  LMEM[counter, "v_DBH"]<- VIF["DBH", 1]
  
  counter<- counter+1
  
  #######################################################################################################
  ##                            Without SEV
  #######################################################################################################
  
  Mixed_model2<-lmer(VAL~SEV0+DBH+(1+SEV0+DBH|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model2))
  #dr<- dredge(Mixed_model2)
  #im<- as.data.frame(sw(dr))
  Sum<- summary(Mixed_model2)
  Sum<- as.data.frame(Sum[["coefficients"]])
  rand<- ranova(Mixed_model2, reduce.terms = F)
  VIF<- as.data.frame(vif(Mixed_model2))
  qqnorm(resid(Mixed_model2))
  
  LMEM[counter, "DATA"]<- Levels3[i, "CODE"]
  LMEM[counter, "Model"]<- "without SEV"
  LMEM[counter, "AIC"]<- AIC(Mixed_model2)
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model2), Set4$VAL)^2
  LMEM[counter, "t_SEV0"]<- Sum["SEV0.L",4]
  #LMEM[counter, "t_SEV"]<- Sum["EXT",4]
  LMEM[counter, "t_DBH"]<- Sum["DBH",]
  LMEM[counter, "p_SEV0"]<- Sum["SEV0.L", 5]
  #LMEM[counter, "p_SEV"]<- Sum["EXT", 5]
  LMEM[counter, "p_DBH"]<- Sum["DBH", 5]
  LMEM[counter, "TOPO"]<- rand[2,6]
  LMEM[counter, "v_SEV0"]<- VIF["SEV0", 1]
  #LMEM[counter, "v_SEV"]<- VIF["SEV", 1]
  LMEM[counter, "v_DBH"]<- VIF["DBH", 1]
  
  counter<- counter+1
  
  #######################################################################################################
  ##                            Without DBH
  #######################################################################################################
  
  Mixed_model3<-lmer(VAL~SEV0+SEV+(1+SEV0+SEV|TOPO), data=Set4)
  
  Rmc<- as.data.frame(r.squaredGLMM(Mixed_model3))
  #dr<- dredge(Mixed_model3)
  #im<- as.data.frame(sw(dr))
  Sum<- summary(Mixed_model3)
  Sum<- as.data.frame(Sum[["coefficients"]])
  rand<- ranova(Mixed_model3, reduce.terms = F)
  VIF<- as.data.frame(vif(Mixed_model3))
  qqnorm(resid(Mixed_model3))
  
  LMEM[counter, "DATA"]<- Levels3[i, "CODE"]
  LMEM[counter, "Model"]<- "without DBH"
  LMEM[counter, "AIC"]<- AIC(Mixed_model3)
  LMEM[counter, "Rm"]<- Rmc[1,1]
  LMEM[counter, "Rc"]<- Rmc[1,2]
  LMEM[counter, "pR"]<- cor(predict(Mixed_model3), Set4$VAL)^2
  LMEM[counter, "t_SEV0"]<- Sum["SEV0.L",4]
  LMEM[counter, "t_SEV"]<- Sum["SEV",4]
  #LMEM[counter, "t_DBH"]<- Sum["DBH",4]
  LMEM[counter, "p_SEV0"]<- Sum["SEV0.L", 5]
  LMEM[counter, "p_SEV"]<- Sum["SEV", 5]
  #LMEM[counter, "p_DBH"]<- Sum["DBH", 5]
  LMEM[counter, "TOPO"]<- rand[2,6]
  LMEM[counter, "v_SEV0"]<- VIF["SEV0", 1]
  LMEM[counter, "v_SEV"]<- VIF["SEV", 1]
  #LMEM[counter, "v_DBH"]<- VIF["DBH", 1]
  
  counter<- counter+1
  
}

write.table(LMEM, "Vystupy/Table 1.txt", row.names = F, col.names = T, sep = "\t")
#write.table(LMEM, "Vystupy/Table S3txt", row.names = F, col.names = T, sep = "\t")

##############################################################
#             9. Resilience indices calculation
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
  scale_fill_manual(values = c("#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"))+
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
  scale_fill_manual(values = c("#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"))+
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
  scale_fill_manual(values = c("#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"))+
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
#          Wilcox test of resilience for TOPO
##############################################################

RES_data$CODE<- paste(RES_data$SPECIES, RES_data$VAR, RES_data$Index, RES_data$YEAR, sep = "_")

Levels2<- as.data.frame(unique(RES_data$CODE)); colnames(Levels2)<- "CODE"

W_res_topo<- as.data.frame(matrix(nrow = nrow(Levels2), ncol = 9)); rownames(W_res_topo)<- Levels2$CODE; colnames(W_res_topo)<- c("NS-PL", "NS-SS", "NS-VA", "PL-PL", "PL-SS", "PL-VA", "SS-PL", "SS-SS", "SS-VA")

for (h in c(1:nrow(Levels2))) {
  
  Set3<- subset(RES_data, subset = RES_data$CODE==Levels2[h, "CODE"])
  
  PostHoc<- pairwise.wilcox.test(Set3$VAL,Set3$TOPO)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  W_res_topo[h,] <-PH_P_val
  
}

W_res_topo<- W_res_topo[,-c(4,7,8)]

#write.table(W_res_topo, "Vystupy/5. W_topo_res.txt", sep="\t")

##############################################################
#       Table S4 (Wilcox test of resilience for years)
##############################################################

RES_data$Code<- paste(RES_data$SPECIES, RES_data$VAR, RES_data$Index, RES_data$TOPO, sep = "_")

Levels_res<- as.data.frame(unique(RES_data$Code)); colnames(Levels_res)<- "Code"

W_res<- as.data.frame(matrix(ncol = 4, nrow = nrow(Levels_res))); rownames(W_res)<- Levels_res$Code; colnames(W_res)<- c("1994-2003", "1994-2006", "2003-2003", "2003-2006")

for (d in c(1:nrow(Levels_res))) {
  
  Set2<- subset(RES_data, subset = RES_data$Code==Levels_res[d, "Code"])
  
  PostHoc<- pairwise.wilcox.test(Set2$VAL,Set2$YEAR)
  
  PH_P_val<- as.vector(PostHoc[["p.value"]])
  
  W_res[d,] <-PH_P_val
  
}

W_res<- W_res[,-c(3)]

write.table(W_res, "Vystupy/Table S4.txt", sep="\t")

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
  theme(strip.text  = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
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
  theme(strip.text  = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 12, angle=45, hjust = 1))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.title.x = element_text(size = 12))+
  theme(axis.title.y = element_text(size = 12))+
  theme(legend.text = element_text(size = 12))+
  theme(panel.background = element_blank())+
  theme(panel.grid = element_blank())+
  theme(axis.line = element_line(colour = "black"))+
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+
  theme(legend.title = element_blank())+
  theme(legend.key = element_rect(colour = "transparent", fill = "white"))

plot_grid(PISY, PCAB, nrow=2, ncol=1, labels = c("PISY", "PCAB"), label_size = 12, label_x = 0.08)

ggsave("Grafy/Fig. S1 (NDVI sensitivity).tiff", height = 150, width = 170, units = "mm", dpi = 300)

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

Data$Code<- paste(Data$SPECIES, Data$VAR, sep = "_")

Levels4<- as.data.frame(unique(Data$Code)); colnames(Levels4)<- "Code"

Dif_let<- as.data.frame(matrix(ncol = 5, nrow = 4)); colnames(Dif_let)<- c("CODE", "NS", "PL", "SS", "VA")

for (i in c(1:4)) {
  
  Set5<- subset(Data, subset = Data$Code==Levels4[i, "Code"])
  
  AOV<- aov(Set5$VAL~Set5$TOPO)
  Tukey<- TukeyHSD(AOV)
  
  cld<- multcompLetters4(AOV, Tukey)
  cld<- as.data.frame.list(cld$`Set5$TOPO`)
  
  Dif_let[i, "CODE"]<- Levels4[i, "Code"]
  Dif_let[i, "NS"]<- cld[rownames(cld)=="NS", "Letters"]
  Dif_let[i, "PL"]<- cld[rownames(cld)=="PL", "Letters"]
  Dif_let[i, "SS"]<- cld[rownames(cld)=="SS", "Letters"]
  Dif_let[i, "VA"]<- cld[rownames(cld)=="VA", "Letters"]
  
}

Dif_let<- gather(Dif_let, "TOPO", "Letters", 2:5)
Dif_let$VAR<- lapply(strsplit(as.character(Dif_let$CODE), "\\_"), "[",2)
Dif_let$SPECIES<- lapply(strsplit(as.character(Dif_let$CODE), "\\_"), "[",1)
Dif_let <- as.data.frame(lapply(Dif_let, unlist))

ggplot(Data, aes(x=TOPO, y=VAL, fill=TOPO))+geom_boxplot()+
  facet_grid2(VAR~SPECIES, scales = "free", axes = "all")+
  scale_fill_manual(values = c("#bcc0c4", "#0d63ba", "#656d75", "#4d9ef0"))+
  geom_text(data = Dif_let, aes(label=Letters, x=TOPO, y=0.8), hjust= -0.3, size=5)+
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

##############################################################
#                       Climadiagrams
##############################################################

PCAB_P<-read.table("C:/Users/jirka/Desktop/Pointers/Data/PCAB_P.txt", dec=",")
PCAB_T<-read.table("C:/Users/jirka/Desktop/Pointers/Data/PCAB_T.txt", dec=",")
PISY_P<-read.table("C:/Users/jirka/Desktop/Pointers/Data/PISY_P.txt", dec=",")
PISY_T<-read.table("C:/Users/jirka/Desktop/Pointers/Data/PISY_T.txt", dec=",")

PCAB_P<- as.data.frame(t(PCAB_P))
PCAB_T<- as.data.frame(t(PCAB_T))
PISY_P<- as.data.frame(t(PISY_P))
PISY_T<- as.data.frame(t(PISY_T))

PCAB_P$Normal<- rowMeans(PCAB_P)
PCAB_T$Normal<- rowMeans(PCAB_T)
PISY_P$Normal<- rowMeans(PISY_P)
PISY_T$Normal<- rowMeans(PISY_T)

PCAB_P$Dry<- rowMeans(PCAB_P[,c(10,19,22,31)])
PCAB_T$Dry<- rowMeans(PCAB_T[,c(10,19,22,31)])
PISY_P$Dry<- rowMeans(PISY_P[,c(10,19,22,31)])
PISY_T$Dry<- rowMeans(PISY_T[,c(10,19,22,31)])

PCAB_P<- as.data.frame(PCAB_P[,c(34,35)])
PCAB_T<- as.data.frame(PCAB_T[,c(34,35)])
PISY_P<- as.data.frame(PISY_P[,c(34,35)])
PISY_T<- as.data.frame(PISY_T[,c(34,35)])

PCAB_P$MONTH<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")
PCAB_T$MONTH<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")
PISY_P$MONTH<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")
PISY_T$MONTH<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")

PCAB_P$SPECIES<- "PCAB"
PCAB_T$SPECIES<- "PCAB"
PISY_P$SPECIES<- "PISY"
PISY_T$SPECIES<- "PISY"

PCAB_P<- gather(PCAB_P, "Type", "Prec", 1:2)
PCAB_T<- gather(PCAB_T, "Type", "Temp", 1:2)
PISY_P<- gather(PISY_P, "Type", "Prec", 1:2)
PISY_T<- gather(PISY_T, "Type", "Temp", 1:2)

Climadiagrams<- rbind(PCAB_P, PISY_P)
Temp<- rbind(PCAB_T, PISY_T)

Climadiagrams$Temp<- as.numeric(Temp$Temp)

Climadiagrams$CODE<- paste(Climadiagrams$SPECIES, Climadiagrams$Type)

aggregate(Climadiagrams$Prec, list(Climadiagrams$CODE), "sum")
aggregate(Climadiagrams$Temp, list(Climadiagrams$CODE), "mean")

### přepočty pro druhou osu
ylim.prim<- c(-5, 65)
ylim.sec<- c(0, 130)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

Climadiagrams$Species_f<- factor(Climadiagrams$SPECIES, levels = c("PISY", "PCAB"))

order<- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",	"Jul",	"Aug",	"Sep",	"Oct",	"Nov",	"Dec")

ggplot(Climadiagrams, aes(x=factor(MONTH, level=order), group=Type))+
  geom_bar(aes(y=a+Prec*b, fill=as.character(Type)), stat='identity', position=position_dodge(.9))+
  geom_line(aes(y=Temp, color=as.character(Type)),  stat='identity', linewidth=2)+
  scale_y_continuous("Temperatures (°C)", sec.axis = sec_axis(~ (. - a)/b, name = "Precipitation (mm)"))+
  facet_wrap(~Species_f)+
  scale_color_manual(values= c("#ad0a0a", "#f28383"))+
  scale_fill_manual(values= c("#a5cef0", "#085696"))+
  
  labs(x="")+
  coord_cartesian(ylim= c(-2,65))+
  theme(strip.text  = element_text(size = 17))+ ## velikost nadpisů
  theme(axis.text.y = element_text(size = 14))+ ## velikost názvu osy y
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1))+ ## velikost názvu osy x
  theme(axis.title.y = element_text(size = 14))+ ## velikost popisků osy y
  theme(axis.title.x = element_text(size = 14))+ ## velikost popisků osy x
  theme(legend.text = element_text(size = 14))+ ## velikost popisků legendy
  theme(panel.background = element_blank())+ ## nabarvení pozadí bíle
  theme(panel.grid = element_blank())+ ## odstranění sítě
  theme(axis.line = element_line(colour = "black"))+ ## barva os
  theme(legend.position = "bottom", legend.margin=margin(t= -25))+
  theme(legend.title = element_blank())

ggsave("Grafy/Climadiagrams.png", height = 100, width = 200, units = "mm", dpi = 300)

