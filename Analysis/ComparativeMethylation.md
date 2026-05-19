# Comparative Methylation
Along with differential methylation I compared overall methylation between populations.
Much of this code is borrowed from 


##### Import Libraries
`library(ggplot2)
library("FactoMineR")
library("factoextra")
library(dplyr)
library(egg)
library(ggsignif)`

##### Import methylation files produced by nanopolish/f5c
###### Make sure the headeds for the file are REMVOED!
lst <- list.files(path="./Projects/West_vs_East/BothCoastData/FreqTbl_2", 
                  pattern = "*.tsv", full.names=TRUE)`
                  
filelist <- lapply(lst, read.table)

names(filelist) <- list.files(path="./Projects/West_vs_East/BothCoastData/FreqTbl_2",
                              pattern = "*.tsv", full.names=FALSE)

mergedData <- Reduce(function(...) merge(..., by = c("V1","V2","V3")), filelist)

###### Just get columns of interest and provide headers
mergedData2 <- mergedData[,c(1:7,9:12,14:17,19:22,24:27,29:32,34:37,39:42,44:47,49:52,54:57,59:62,64:67,69:72,74:77,79:82,84:87,89:92,94:97,99:102,104:107,109:112,114:117,119:122,124:127,129:132,134:137,139:142,144:147,149:152,154:157,159:162,164:167,169:172,174:177,179:182,184:187,189:192,194:197,199:202,204:207,209:212,214:217,219:222,224:227,229:232,234:237,239:242)]

colnames(mergedData2) <- c("chrom", "start", "end",	
                           "num_cpgs_CoosA", "called_CoosA",	"called_meth_CoosA", "meth_freq_CoosA",
                           "num_cpgs_CoosB", "called_CoosB",	"called_meth_CoosB", "meth_freq_CoosB",
                           "num_cpgs_CoosE", "called_CoosE",	"called_meth_CoosE", "meth_freq_CoosE",
                           "num_cpgs_CoosF", "called_CoosF",	"called_meth_CoosF", "meth_freq_CoosF",
                           "num_cpgs_CoosI", "called_CoosI",	"called_meth_CoosI", "meth_freq_CoosI",
                           "num_cpgs_CoosJ", "called_CoosJ",	"called_meth_CoosJ", "meth_freq_CoosJ",
                           "num_cpgs_CRB1", "called_CRB1",	"called_meth_CRB1",	"meth_freq_CRB1",
                           "num_cpgs_CRB2", "called_CRB2",	"called_meth_CRB2",	"meth_freq_CRB2",
                           "num_cpgs_CRB3", "called_CRB3",	"called_meth_CRB3",	"meth_freq_CRB3",
                           "num_cpgs_CRB5", "called_CRB5",	"called_meth_CRB5",	"meth_freq_CRB5",
                           "num_cpgs_CRB6", "called_CRB6",	"called_meth_CRB6",	"meth_freq_CRB6",
                           "num_cpgs_CRB8", "called_CRB8",	"called_meth_CRB8",	"meth_freq_CRB8",
                           "num_cpgs_MA1", "called_MA1",	"called_meth_MA1", "meth_freq_MA1",
                           "num_cpgs_MA2", "called_MA2",	"called_meth_MA2", "meth_freq_MA2",
                           "num_cpgs_MA3", "called_MA3",	"called_meth_MA3", "meth_freq_MA3",
                           "num_cpgs_MA4", "called_MA4",	"called_meth_MA4", "meth_freq_MA4",
                           "num_cpgs_MA5", "called_MA5",	"called_meth_MA5",	"meth_freq_MA5",
                           "num_cpgs_MA6", "called_MA6",	"called_meth_MA6",	"meth_freq_MA6",
                           "num_cpgs_MA7", "called_MA7",	"called_meth_MA7",	"meth_freq_MA7",
                           "num_cpgs_MA8", "called_MA8",	"called_meth_MA8",	"meth_freq_MA8",
                           "num_cpgs_ME1", "called_ME1",	"called_meth_ME1",	"meth_freq_ME1",
                           "num_cpgs_ME2", "called_ME2",	"called_meth_ME2",	"meth_freq_ME2",
                           "num_cpgs_ME3", "called_ME3",	"called_meth_ME3",	"meth_freq_ME3",
                           "num_cpgs_ME4", "called_ME4",	"called_meth_ME4",	"meth_freq_ME4",
                           "num_cpgs_ME5", "called_ME5",	"called_meth_ME5",	"meth_freq_ME5",
                           "num_cpgs_ME6", "called_ME6",	"called_meth_ME6",	"meth_freq_ME6",
                           "num_cpgs_ME7", "called_ME7",	"called_meth_ME7",	"meth_freq_ME7",
                           "num_cpgs_ME8", "called_ME8",	"called_meth_ME8",	"meth_freq_ME8",
                           "num_cpgs_NH1", "called_NH1",	"called_meth_NH1",	"meth_freq_NH1",
                           "num_cpgs_NH2", "called_NH2",	"called_meth_NH2",	"meth_freq_NH2",
                           "num_cpgs_NH3", "called_NH3",	"called_meth_NH3",	"meth_freq_NH3",
                           "num_cpgs_NH4", "called_NH4",	"called_meth_NH4",	"meth_freq_NH4",
                           "num_cpgs_NH5", "called_NH5",	"called_meth_NH5",	"meth_freq_NH5",
                           "num_cpgs_NH6", "called_NH6",	"called_meth_NH6",	"meth_freq_NH6",
                           "num_cpgs_NH7", "called_NH7",	"called_meth_NH7",	"meth_freq_NH7",
                           "num_cpgs_NH8", "called_NH8",	"called_meth_NH8",	"meth_freq_NH8",
                           "num_cpgs_NW2", "called_NW2",	"called_meth_NW2",	"meth_freq_NW2",
                           "num_cpgs_NW3", "called_NW3",	"called_meth_NW3",	"meth_freq_NW3",
                           "num_cpgs_NW4", "called_NW4",	"called_meth_NW4",	"meth_freq_NW4",
                           "num_cpgs_NW6", "called_NW6",	"called_meth_NW6",	"meth_freq_NW6",
                           "num_cpgs_NW7", "called_NW7",	"called_meth_NW7",	"meth_freq_NW7",
                           "num_cpgs_NWA", "called_NWA",	"called_meth_NWA",	"meth_freq_NWA",
                           "num_cpgs_SSA", "called_SSA",	"called_meth_SSA",	"meth_freq_SSA",
                           "num_cpgs_SSB", "called_SSB",	"called_meth_SSB",	"meth_freq_SSB",
                           "num_cpgs_SSC", "called_SSC",	"called_meth_SSC",	"meth_freq_SSC",
                           "num_cpgs_SSD", "called_SSD",	"called_meth_SSD",	"meth_freq_SSD",
                           "num_cpgs_SSE", "called_SSE",	"called_meth_SSE",	"meth_freq_SSE",
                           "num_cpgs_SSF", "called_SSF",	"called_meth_SSF",	"meth_freq_SSF")

##### Calculate and plot means
meth_means <- colMeans(subset(mergedData2, select = c(meth_freq_CoosA,meth_freq_CoosB,meth_freq_CoosE,meth_freq_CoosF,meth_freq_CoosI,meth_freq_CoosJ,
                                                      meth_freq_CRB1,meth_freq_CRB2,meth_freq_CRB3,meth_freq_CRB5,meth_freq_CRB6,meth_freq_CRB8,
                                                      meth_freq_MA1,meth_freq_MA2,meth_freq_MA3,meth_freq_MA4,meth_freq_MA5,meth_freq_MA6,meth_freq_MA7,meth_freq_MA8,
                                                      meth_freq_ME1,meth_freq_ME2,meth_freq_ME3,meth_freq_ME4,meth_freq_ME5,meth_freq_ME6,meth_freq_ME7,meth_freq_ME8,
                                                      meth_freq_NH1,meth_freq_NH2,meth_freq_NH3,meth_freq_NH4,meth_freq_NH5,meth_freq_NH6,meth_freq_NH7,meth_freq_NH8,
                                                      meth_freq_NW2,meth_freq_NW3,meth_freq_NW4,meth_freq_NW6,meth_freq_NW7,meth_freq_NWA,
                                                      meth_freq_SSA,meth_freq_SSB,meth_freq_SSC,meth_freq_SSD,meth_freq_SSE,meth_freq_SSF)))


plot(meth_means)

meth_means2 <- as.data.frame(cbind(meth_means, c("Coos","Coos","Coos","Coos","Coos","Coos",
                                                 "CRB","CRB","CRB", "CRB", "CRB","CRB",
                                                 "Mass","Mass","Mass","Mass","Mass","Mass","Mass","Mass",
                                                 "Maine","Maine","Maine","Maine","Maine","Maine","Maine","Maine",
                                                 "N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp",
                                                 "NW","NW", "NW","NW","NW","NW",
                                                 "SS","SS","SS","SS","SS","SS")))

meth_means2$meth_means <- as.numeric(as.character((meth_means2$meth_means)))

##### Same but with just West and East
wvw_methmeans <- as.data.frame(cbind(meth_means, c("West","West","West","West","West","West",
                                                   "West","West","West","West","West","West",
                                                   "East","East","East","East","East","East","East","East",
                                                   "East","East","East","East","East","East","East","East",
                                                   "East","East","East","East","East","East","East","East",
                                                   "West","West","West","West","West","West",
                                                   "West","West","West","West","West","West")))
wvw_methmeans$meth_means <- as.numeric(as.character((wvw_methmeans$meth_means)))

##### Make Plots
meth_means2 %>%
  arrange(meth_means) %>%
  mutate(V2=factor(V2,levels=c("Coos","CRB","NW","SS",
                               "Mass","Maine","N.Hamp")))%>% 
  ggplot(aes(x=V2, y=meth_means, color=V2))+
  geom_boxplot(linewidth=0.8)+
  geom_point(size=1, shape=1)+
  scale_color_manual(values=c("#A40122","#Dc3440" , "red", "#ff7079",  "#01366c", "#0071b3", "#7aaae1")) +
  theme_bw()+
  labs(x="Location",y="Methylation frequency") +
  theme(legend.position = "none") 

wvw_methmeans%>%
  arrange(meth_means)%>%
  mutate(V2=factor(V2,levels=c("West","East")))%>%
  ggplot(aes(x=V2, y=meth_means, color=V2))+
  geom_boxplot(linewidth=0.8)+
  geom_point(size=2)+
  scale_color_manual(values=c("#Dc3440","#0071b3")) +
  theme_minimal()+
  labs(x="Coast",y="Methylation Frequency") +
  theme(legend.position = "none")
  
##### Stats tests 

###### Overall variance homogeneity test and ttest

var.test(meth_means[c(1:13,37:48)], meth_means[13:36])
t.test(meth_means[c(1:13,37:48)], meth_means[13:36], paired = FALSE, var.equal = TRUE)

west_means<- c(0.3489694,
               0.3437167,
               0.3439507,
               0.3356036,
               0.3417704,
               0.3481926,
               0.3413838,
               0.3397191,
               0.3476416,
               0.3399276,
               0.3307446,
               0.3431359,
               0.343127,
               0.3374195,
               0.3395435,
               0.3356902,
               0.3497697,
               0.3252493,
               0.3418478,
               0.342266,
               0.3354772,
               0.3504246,
               0.3387371,
               0.3439276)

east_means<- c(0.3433764,
               0.3426709,
               0.3355774,
               0.3384365,
               0.3473978,
               0.3408298,
               0.3363567,
               0.3418084,
               0.344749,
               0.3330713,
               0.3295194,
               0.3278431,
               0.3381264,
               0.3422575,
               0.3407639,
               0.3410703,
               0.3414416,
               0.3393886,
               0.3477096,
               0.3475109,
               0.3420333,
               0.3379718,
               0.3406484,
               0.3380496)

#Overall variance homogeneity test and ttest
var.test(west_means,east_means)

t.test(west_means,east_means, paired = FALSE, var.equal = TRUE)

res_aov<-aov(meth_means~V2,meth_means2)
summary(res_aov)
plot(res_aov)

kruskal.test(meth_means~V2,meth_means2)

##### Make a PCA
meth_cols <- mergedData2[,c(7,11,15,19,23,27,31,35,39,43,47,51,55,59,63,67,71,75,79,83,87,91,95,99,103,107,111,115,119,123,127,131,135,139,143,147,151,155,159,163,167,171,175,179,183,187,191,195)]
row.names(meth_cols) <- mergedData2$chrom

pca <- prcomp(meth_cols, scale=T)
scores <- data.frame(pca$rotation)
scores2 <- cbind(c("Coos","Coos","Coos","Coos","Coos","Coos",
                   "CRB","CRB","CRB", "CRB", "CRB","CRB",
                   "Mass","Mass","Mass","Mass","Mass","Mass","Mass","Mass",
                   "Maine","Maine","Maine","Maine","Maine","Maine","Maine","Maine",
                   "N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp","N.Hamp",
                   "NW","NW", "NW","NW","NW","NW",
                   "SS","SS","SS","SS","SS","SS"), scores)
colnames(scores2)[1] <- "Location"


eig <- fviz_eig(pca, geom = "bar", barfill = "white", barcolor = "black",
                ggtheme = theme_classic(), main= "", ylab= "% Variance", addlabels = TRUE)
eig

scores2 %>%
  arrange(scores) %>%
  mutate(Location=factor(Location,levels=c("Coos","CRB","NW","SS",
                               "Mass","Maine","N.Hamp")))%>%
  ggplot(aes(PC1, PC2, colour = Location)) + 
  theme_bw() +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#A40122", "#cd1c2d", "red", "#ff7079",  "#01366c", "#0071b3", "#7aaae1")) +
  xlab("PC1 54.5%") +
  ylab("PC2 1.9%")+
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) +
  theme(text = element_text(size = 12))


##### PCA for coast verses coast

scores_w_e<-cbind(c("West","West","West","West","West","West",
                    "West","West","West", "West","West","West",
                    "East","East","East","East","East","East","East","East",
                    "East","East","East","East","East","East","East","East",
                    "East","East","East","East","East","East","East","East",
                    "West","West","West","West","West","West",
                    "West","West","West","West","West","West"),c("Coos","Coos","Coos","Coos","Coos","Coos",
                                                                 "CA","CA","CA","CA","CA","CA",
                                                                 "MA","MA","MA","MA","MA","MA","MA","MA",
                                                                 "ME","ME","ME","ME","ME","ME","ME","ME",
                                                                 "NH","NH","NH","NH","NH","NH","NH","NH",
                                                                 "NW","NW","NW","NW","NW","NW",
                                                                 "SS","SS","SS","SS","SS","SS"), scores)
colnames(scores_w_e)[1] <- "Location"
colnames(scores_w_e)[2] <- "Pop"

eig_we <- fviz_eig(pca, geom = "bar", barfill = "white", barcolor = "black",
                ggtheme = theme_classic(), main= "", ylab= "% Variance", addlabels = TRUE)
eig_we

ggplot(scores_w_e, aes(PC1, PC2)) + 
  theme_minimal() +
  geom_point(size=2.5,aes(shape= Pop, colour = Location)) +
  scale_shape_manual(values=c(0,1,15,16,17,2,5))+
  #scale_shape_manual(values=c(16,16,15,15,15,16,16))+
  scale_colour_manual(values = c("#0071b3","#cd1c2d")) +
  xlab("PC1 54.5%") +
  ylab("PC2 1.9%")+
  theme(axis.text = element_text(size = 12, color = "black")) +
  theme(axis.title = element_text(size = 12, color = "black")) +
  theme(axis.ticks = element_line(color = "black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(size = 1)) +
  theme(text = element_text(size = 14))


#THis is more like it!
#Take this out and combine first three columns so I know what loci are which
write.csv(mergedData2, "mergedData2_dimondstyle.csv", row.names = FALSE)

merge_freqtable<-read.table("./Projects/West_vs_East/BothCoastData/merge_methfreq.tsv", header = TRUE)

merge_freqtable_2<-merge_freqtable[,-1]
rownames(merge_freqtable_2)<-merge_freqtable[,1]

head(merge_freqtable_2)


#calculate and plot means
coos_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_CoosA,meth_freq_CoosB,meth_freq_CoosE,meth_freq_CoosF,meth_freq_CoosI,meth_freq_CoosJ)))
ca_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_CRB1,meth_freq_CRB2,meth_freq_CRB3,meth_freq_CRB5,meth_freq_CRB6,meth_freq_CRB8)))
ss_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_SSA,meth_freq_SSB,meth_freq_SSC,meth_freq_SSD,meth_freq_SSE,meth_freq_SSF)))
nw_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_NW2,meth_freq_NW3,meth_freq_NW4,meth_freq_NW6,meth_freq_NW7,meth_freq_NWA)))

ma_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_MA1,meth_freq_MA2,meth_freq_MA3,meth_freq_MA4,meth_freq_MA5,meth_freq_MA6,meth_freq_MA7,meth_freq_MA8)))
me_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_ME1,meth_freq_ME2,meth_freq_ME3,meth_freq_ME4,meth_freq_ME5,meth_freq_ME6,meth_freq_ME7,meth_freq_ME8)))
nh_means <- rowMeans(subset(merge_freqtable_2, select = c(meth_freq_NH1,meth_freq_NH2,meth_freq_NH3,meth_freq_NH4,meth_freq_NH5,meth_freq_NH6,meth_freq_NH7,meth_freq_NH8)))

merge_freq_means<- as.data.frame(cbind(coos_means,ca_means,ss_means,nw_means,me_means,ma_means,nh_means))

#Have to rename rows AFTER making numeric
merge_freq_means2<-sapply(merge_freq_means,as.numeric)
colnames(merge_freq_means2)<-c("Coos","CA","S.Slough","Newport","Maine","Mass","New Hamp")
rownames(merge_freq_means2)<-merge_freqtable[,1]
heatmap.2(merge_freq_means2, scale = "none", trace = "none", Colv = F, col=rgb.palette)
heatmap.2(merge_freq_means2, trace= "none", Colv = F, col=rgb.palette)

write.csv(merge_freq_means2, "merged_freqmeans2_dimondstyle.csv", row.names = TRUE)


