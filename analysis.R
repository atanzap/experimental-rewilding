##################################################################################
# Trade-offs between passive and trophic rewilding for biodiversity and ecosystem functioning
# File contains the following sections:
# 0) Prepare workspace
# 1) Tree diversity
# 2) Arthropod diversity and biomass
# 3) Arthropod composition and function
# 4) Vegetation structure
# 5) Does soil C differ between plots?
# 6) Does above-ground C differ between plots?
# 7) Does total C differ between plots?
# 8) Historical baseline comparison
# 9) Prepare display items (Figs 1, S1-S3)

######################################################################################
## 0) Prepare workspace
######################################################################################
# setwd and loack packages
setwd('./data/')
library(lme4)
library(car)
library(coin)
library(vegan)
library(ggplot2)
library(pals)
library(reshape2)
library(dabestr)
library(vioplot)
library(LDM)
library(MASS)
library(emmeans)

# function to calculate stem C for oaks (inc blackthorn, willow) from FC Carbon Assessment Protocol (https://www.woodlandcarboncode.org.uk/images/PDFs/WCC_CarbonAssessmentProtocol_V2.0_March2018.pdf)
# could consider just calculating ABG biomass
broadleaf_C <- function(DBH_cm,Ht_m){
                                      Tariff = 5.88300 + (2.01230 * Ht_m) + (-0.0054780 * DBH_cm) + (-0.0057397 * DBH_cm * Ht_m)
                                      quadratic_mean_dbh = sqrt(DBH_cm^2)
                                      BA = pi * quadratic_mean_dbh^2 / 40000
                                      BA = pi * (DBH_cm/2/100)^2
                                      Vol = (0.0360541 * Tariff) - ((0.315049301 * (Tariff - 0.138763302)) * 0.118288) + ((0.315049301 * (Tariff - 0.138763302)) * BA)
                                      stem_biomss = Vol * merch_vol_to_stem_vol_FC[match(round(DBH_cm),merch_vol_to_stem_vol_FC$dbh_cm),'factor'] * 0.56 #oak-specific gravity
                                      crown_biomass = 0.0000168513 * DBH_cm^2.4767
                                      root_biomass =  0.000022700 * DBH_cm^2.5
                                      treeC = (stem_biomss+crown_biomass+root_biomass)*.5*1000
                                      return(as.numeric(treeC))
                                    }
# scaling function
scale.max.1 <- function(x) x / max(x)

# carbon conversions for saplings from Table 6.1.3 in FC Carbon Assessment Protocol (https://www.woodlandcarboncode.org.uk/images/PDFs/WCC_CarbonAssessmentProtocol_V2.0_March2018.pdf)
sapling_carbon_FC <- read.csv('sapling_carbon_FC.csv',header=T,fileEncoding="UTF-8-BOM")
merch_vol_to_stem_vol_FC <- read.csv('merch_vol_to_stem_vol_FC.csv',header=T,fileEncoding="UTF-8-BOM")

# function to calculate overdispersion (HT: https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)
overdisp_fun <- function(model) {
    rdf <- df.residual(model)
    rp <- residuals(model,type="pearson")
    Pearson.chisq <- sum(rp^2)
    prat <- Pearson.chisq/rdf
    pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
    c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# function to plot dimensionality vs stress
NMDS.scree<-function(x) { 
  mm1 <- metaMDS(x,distance="bray",k=1)$stress
  plot(1,mm1, xlim=c(1,10),ylim=c(0,0.6),xlab="# of Dimensions",ylab="Stress",main="NMDS stress plot")
  for (i in 2:10) {
      points(i,metaMDS(x,distance="bray",k=i)$stress)
  }
}

# read in datafiles
soilC <- read.csv('soils.csv',header=T,fileEncoding="UTF-8-BOM")
treeC <- read.csv('trees.csv',header=T,fileEncoding="UTF-8-BOM")
grndC <- read.csv('groundveg.csv',header=T,fileEncoding="UTF-8-BOM")
gd_CC <- read.csv('groundveg_conversions.csv',header=T,fileEncoding="UTF-8-BOM")
veg_struc <- read.csv('RawVeg.csv',header=T)
bugs <- read.csv('arthropods.csv',header=T,fileEncoding="UTF-8-BOM")
hist_biom <- read.csv('historical_biomass.csv',header=T)
hist_swht <- read.csv('historical_swardht.csv',header=T)
hist_gcov <- read.csv('historical_cover.csv',header=T)
hist_spp_list <- read.csv('tree_shrub_list.csv',header=T)
hist_soils <- read.csv('historical_soils.csv',header=T)
art_guilds <- read.csv('arthropod_guilds.csv',header=T)
animals <- read.csv('animals.csv',header=T,fileEncoding="UTF-8-BOM")
animals_by_class_v1 <- read.csv('animals_by_class.csv',header=T,fileEncoding="UTF-8-BOM")
animals_by_class_v2 <- read.csv('animals_by_class_v2.csv',header=T,fileEncoding="UTF-8-BOM")
height_by_field <- read.csv('knepp_height2.csv')


######################################################################################
## 1) Woody species diversity classed based on http://ecoflora.org.uk/
######################################################################################
# prepare current data
treeD <- with(treeC, table(Species,Site,Treatment))
treeD[rownames(treeD) == 'Blackthorn',,] <- treeD[rownames(treeD) == 'Blackthorn',,] + treeD[rownames(treeD) == 'Blackthorn_patch',,]
treeD[rownames(treeD) == 'Rose',,] <- treeD[rownames(treeD) == 'Rose',,] + treeD[rownames(treeD) == 'Rose_patch',,]
#treeD <- treeD[-which(rownames(treeD) %in% c('Blackthorn_patch','Rose_patch','Rose')),,]
treeD <- treeD[-which(rownames(treeD) %in% c('Blackthorn_patch','Rose_patch')),,]
treeD <- cbind(rowSums(t(treeD[,,1]) > 0),rowSums(t(treeD[,,2]) > 0))

# create community matrix
tree_spp <- dcast(treeC[,c('Site','Treatment','Species','DBH_cm')], Site+Treatment ~ Species)
tree_spp$Blackthorn <- tree_spp$Blackthorn + tree_spp$Blackthorn_patch
tree_spp$Rose <- tree_spp$Rose + tree_spp$Rose_patch
tree_spp <- tree_spp[,which(colnames(tree_spp) %in% c('Blackthorn_patch','Rose_patch') == F)]

# add missing rows
tree_spp_emp1 <- cbind(names(which(table(tree_spp$Site) == 1)),rep('C',sum(table(tree_spp$Site) == 1)),as.data.frame(matrix(0,sum(table(tree_spp$Site) == 1),4)))
tree_spp_emp2 <- cbind(rep(unique(soilC$Site)[which(unique(soilC$Site) %in% unique(tree_spp$Site) == F)],2),rep(c('C','E'),each=2),as.data.frame(matrix(0,length(which(unique(soilC$Site) %in% unique(tree_spp$Site) == F))*2,4)))
colnames(tree_spp_emp1) <- colnames(tree_spp)
colnames(tree_spp_emp2) <- colnames(tree_spp)
tree_spp_all <- rbind(tree_spp,tree_spp_emp1,tree_spp_emp2)

# subset historical vegetation composition only to woody taxa
hist_spp_list_sub <- as.character(hist_spp_list$species)#[which(hist_spp_list$species %in% c('Lonicera_periclymenum','Ribes_idaeus','Rosa_canina','Rubus_frut._Agg','Rubus_frut._Agg_shade_canopy','Ulex_europaeus') == F),])
wood_cov <- hist_gcov[,which(colnames(hist_gcov) %in% c('FieldName','Plot','Quadrat',hist_spp_list_sub))]

# replace p with minimum value of 1
wood_cov[,which(sapply(wood_cov, 'class') == 'factor')[-c(1:3)]] <- apply(wood_cov[,which(sapply(wood_cov, 'class') == 'factor')[-c(1:3)]],2,function(x){as.numeric(gsub('p','1',x))})

# sum same species between understorey and canopy
mnames1 <- sapply(colnames(wood_cov[-c(1:3)]),function(x){ grep(x,colnames(wood_cov)) })
mnames_rep <- apply(sapply(which(sapply(mnames1,function(x){length(x>1)})==2), function(x){mnames1[[x]]}),2,function(y){ rowSums(wood_cov[,y]) })
wood_cov_sub <- wood_cov[,-as.vector(sapply(which(sapply(mnames1,function(x){length(x>1)})==2), function(x){mnames1[[x]]}))]
wood_cov_sub <- cbind(wood_cov_sub,mnames_rep)

# Salix_caprea and Salix_cinerea hybridise - little Salix_caprea (1% in only 2 plots), so merge
wood_cov_sub$Sallow <- with(wood_cov_sub,Salix_caprea+Salix_cinerea)
wood_cov_sub <- subset(wood_cov_sub,select=-c(Salix_caprea,Salix_cinerea))

# get plot-level diversity
wood_cov_sub$Treat <- substr(wood_cov_sub$Plot,1,1)
wood_cov_sub$Plot <- with(wood_cov_sub,paste(FieldName,Plot,sep='_'))
wood_rich <- sapply(unique(wood_cov_sub$Plot), function(x) {
                                                         foo1 <- wood_cov_sub[which(wood_cov_sub$Plot==x),];
                                                         foo2 <- colSums(foo1[,which(sapply(foo1, 'class') %in% c('factor','character')==F)]);
                                                         return(sum(foo2>0)) })
wood_rich <- as.data.frame(cbind(names(wood_rich),wood_rich))
rownames(wood_rich) <- NULL
colnames(wood_rich)[1] <- 'Plot'
wood_rich$Site <- sapply(as.character(wood_rich$Plot),function(x){strsplit(x,'_')[[1]][1]})
wood_rich$Treat <- as.factor(substr(sapply(as.character(wood_rich$Plot),function(x){strsplit(x,'_')[[1]][2]}),1,1))
levels(wood_rich$Treat)[2] <- 'C'
wood_rich$BA <- c('before')

# merge historical and current data
colnames(treeD) <- c('C','E')
treeD <- melt(treeD)
colnames(treeD) <- c('Site','Treat','wood_rich')
treeD$BA <- c('later')
# add in plots with 0 in present-day as somehow these get dropped from aggregation
treeD_mis <- as.data.frame(cbind(rep(unique(wood_rich$Site[which(wood_rich$Site %in% treeD$Site == F)]),2),
                                 rep(c('C','E'),each=length(unique(wood_rich$Site[which(wood_rich$Site %in% treeD$Site == F)]))),
                                 rep(0,2*length(unique(wood_rich$Site[which(wood_rich$Site %in% treeD$Site == F)]))),
                                 rep('later',2*length(unique(wood_rich$Site[which(wood_rich$Site %in% treeD$Site == F)])))  
                              ))
colnames(treeD_mis) <- colnames(treeD)
treeD_mis$wood_rich <- as.numeric(as.character(treeD_mis$wood_rich))
treeD_all <- rbind(treeD, treeD_mis)
# and merge everything together
wood_rich_merge <- rbind(treeD_all,wood_rich[,match(colnames(treeD_all),colnames(wood_rich))])
wood_rich_merge$wood_rich <- as.numeric(wood_rich_merge$wood_rich)
wood_rich_merge$TreatC <- as.factor(wood_rich_merge$Treat=='C')

# fit model with historical comparison
trich_m1 <- glmer(wood_rich ~ TreatC*BA + (1|Site),family='poisson',data=wood_rich_merge,control=glmerControl(optimizer="bobyqa"))
overdisp_fun(trich_m1)

# get values for main text
trich_m1_emmeans <- emmeans(trich_m1,pairwise ~ TreatC|BA, type = "response")$emmeans		


######################################################################################
## 2) Arthropod diversity and biomass
######################################################################################
# drop 27 acres day 1 as sample was misplaced
bugs <- bugs[which(with(bugs,Site=='27acres' & Day=='1')==F),]

#change Canacidae to Carnidae as these were obviously mixed up (former lives only on beaches)
bugs$Carnidae <- bugs$Canacidae + bugs$Carnidae
bugs <- subset(bugs, select = -Canacidae)

# world spider catalogue now lists the Zoridae under the Miturgidae
bugs$Miturgidae <- bugs$Miturgidae + bugs$Zoridae
bugs <- subset(bugs, select = -Zoridae)

# create community matrix and calculate diversity
bugcom <- aggregate(subset(bugs,select=-c(Site,Exclosure,Day,Trap,Trap.biomass)),by=list(bugs$Site,bugs$Exclosure), FUN=sum)
bugdiv <- as.data.frame(cbind(diversity(subset(bugcom,select=-c(Group.1,Group.2))),specnumber(subset(bugcom,select=-c(Group.1,Group.2)))))
colnames(bugdiv) <- c('shannon','sppnum')
bugdiv$Site <- bugcom$Group.1
bugdiv$Exclosure <- bugcom$Group.2
bugdiv$Abund <- rowSums(subset(bugcom,select=-c(Group.1,Group.2)))

###########################################################
## herbivory does not directly enhance arthropod diversity
# this analysis is not in the text
#shapiro.test(bugdiv$sppnum)
#shapiro.test(bugdiv$shannon)
#t.test(bugdiv[which(bugdiv$Exclosure=='Control'),'shannon'],bugdiv[which(bugdiv$Exclosure=='Exclosure'),'shannon'],paired=T)
#t.test(bugdiv[which(bugdiv$Exclosure=='Control'),'sppnum'],bugdiv[which(bugdiv$Exclosure=='Exclosure'),'sppnum'],paired=T)

## drop site with just 1 sampling night from diversity analyses
#t.test(bugdiv[which(bugdiv$Exclosure=='Control' & bugdiv$Site!='27acres'),'shannon'],bugdiv[which(bugdiv$Exclosure=='Exclosure' & bugdiv$Site!='27acres'),'shannon'],paired=T)
#t.test(bugdiv[which(bugdiv$Exclosure=='Control' & bugdiv$Site!='27acres'),'sppnum'],bugdiv[which(bugdiv$Exclosure=='Exclosure' & bugdiv$Site!='27acres'),'sppnum'],paired=T)
###########################################################

# but herbivory does enhance arthropod diversity by enhancing structural diversity when considering direct/indirect effects
# create new matrix of diversity data
bugdiv_pa <- bugdiv
levels(bugdiv_pa$Site)[which(levels(bugdiv_pa$Site) %in% c('Broomerscorner','Oaklandslegg'))] <- c('Broomer','Oaklands')

# calculate vegetation structural diversity
veg_struc$Treatment <- sapply(as.character(veg_struc$Site),function(x){strsplit(x,'_')[[1]][2]})
veg_struc$Site <- sapply(as.character(veg_struc$Site),function(x){strsplit(x,'_')[[1]][1]})

# calculate diversity (http://www.countrysideinfo.co.uk/simpsons.htm)
# following approaches here: https://doi.org/10.1007/s10021-013-9694-8; https://link.springer.com/content/pdf/10.1007/BF00384260.pdf
# methods papers: https://doi.org/10.1016/j.foreco.2016.09.003
# use simpson rather than shannon because assumes finite community - which is the case here with specific number of height tiers
Si_D <- with(veg_struc, sapply(sort(unique(paste(Site,Treatment))), function(x){ foo1 <- veg_struc[which(paste(veg_struc$Site,veg_struc$Treatment)==x),];
                                                                                   foo2 <- rowSums(sapply(1:5, function(y){ foo1[foo1$Subsample.number==y,'Hits'] }));
                                                                                   return(sum(foo2*(foo2-1))/(sum(foo2)*(sum(foo2)-1))) }))
Si_D <- cbind(Si_D[seq(2,length(Si_D),by=2)],Si_D[seq(1,length(Si_D),by=2)])
colnames(Si_D) <- c('Control','Exclosure')
rownames(Si_D) <- as.vector(sapply(rownames(Si_D),function(x){strsplit(x,' ')[[1]][1]}))
bugdiv_pa$strucD <- 1/c(Si_D[match(bugdiv_pa[which(bugdiv$Exclosure=='Control'),'Site'],rownames(Si_D)),'Control'],Si_D[match(bugdiv_pa[which(bugdiv$Exclosure=='Exclosure'),'Site'],rownames(Si_D)),'Exclosure'])

# build dataframe for analysis that accounts for day effects
spp_num_day <- aggregate(subset(bugs,select=-c(Site,Exclosure,Day,Trap,Trap.biomass)),by=list(bugs$Site,bugs$Exclosure,bugs$Day), FUN=sum)
spp_num_day <- as.data.frame(cbind(specnumber(subset(spp_num_day,select=-c(Group.1,Group.2,Group.3))),spp_num_day[,1:3]))
colnames(spp_num_day) <- c('sppdiv','Site','Plot','Day')
levels(spp_num_day$Site)[which(levels(spp_num_day$Site) %in% c('Broomerscorner','Oaklandslegg'))] <- c('Broomer','Oaklands')
spp_num_day$strucD <- bugdiv_pa[match(with(spp_num_day,paste(Site,Plot)),with(bugdiv_pa,paste(Site,Exclosure))),]$strucD
# calculate differences in metrics between exclosures
sppnum_diff <- log(spp_num_day[spp_num_day$Plot=='Exclosure','sppdiv']) - log(spp_num_day[spp_num_day$Plot=='Control','sppdiv'])
strdiv_diff <- as.numeric(scale(log(spp_num_day[spp_num_day$Plot=='Exclosure','strucD']) - log(spp_num_day[spp_num_day$Plot=='Control','strucD'])))
# fit model intercept = direct effect of fencing, structural diversity = indirect effect
a_div_m1 <- lm(sppnum_diff ~ as.factor(spp_num_day[which(spp_num_day$Plot=='Exclosure'),]$Day) + strdiv_diff)
# model siplification
a_div_m2 <- lm(sppnum_diff ~ strdiv_diff)
# values to report in the text
tapply(spp_num_day$sppdiv,spp_num_day$Plot,mean)
tapply(spp_num_day$sppdiv,spp_num_day$Plot,quantile,c(0.025,0.975))


# repeat analysis for biomass
# collate data nicely
bugmass_day <- as.data.frame(with(bugs, tapply(Trap.biomass,list(Site,Exclosure,Day),sum)))
rownames(bugmass_day)[which(rownames(bugmass_day) %in% c('Broomerscorner','Oaklandslegg'))] <- c('Broomer','Oaklands')
bugmass_day_wide <- melt(bugmass_day)
colnames(bugmass_day_wide) <- c('trap','biomass')
bugmass_day_wide$Site <- rep(rownames(bugmass_day),ncol(bugmass_day))
bugmass_day_wide <- bugmass_day_wide[!is.na(bugmass_day_wide$biomass),]
bugmass_day_wide$Plot <- sapply(strsplit(as.character(bugmass_day_wide$trap),"[.]"),function(x){x[1]})
bugmass_day_wide$Day <- sapply(strsplit(as.character(bugmass_day_wide$trap),"[.]"),function(x){x[2]})
spp_num_day_biomass <-  merge(subset(bugmass_day_wide, select = -trap),spp_num_day,by=c('Site','Plot','Day'))

# repeat collation for abundance
bugs_abund <- cbind(aggregate(subset(bugs,select=-c(Site,Exclosure,Day,Trap,Trap.biomass)),by=list(bugs$Site,bugs$Exclosure,bugs$Day),FUN=sum)[,1:3],
                     rowSums(aggregate(subset(bugs,select=-c(Site,Exclosure,Day,Trap,Trap.biomass)),by=list(bugs$Site,bugs$Exclosure,bugs$Day),FUN=sum)[,-c(1:3)]))
colnames(bugs_abund) <- c('Site','Plot','Day','abund')
levels(bugs_abund$Site)[which(levels(bugs_abund$Site) %in% c('Broomerscorner','Oaklandslegg'))] <- c('Broomer','Oaklands')
spp_num_day_biomass <-  merge(bugs_abund,spp_num_day_biomass,by=c('Site','Plot','Day'))

# fit models
biomass_diff <- log(spp_num_day_biomass[spp_num_day_biomass$Plot=='Exclosure','biomass']) - log(spp_num_day_biomass[spp_num_day_biomass$Plot=='Control','biomass'])
abundan_diff <- log(spp_num_day_biomass[spp_num_day_biomass$Plot=='Exclosure','abund']) - log(spp_num_day_biomass[spp_num_day_biomass$Plot=='Control','abund'])
a_mas_m1 <- lm(biomass_diff ~ as.factor(spp_num_day_biomass[which(spp_num_day_biomass$Plot=='Exclosure'),]$Day) + strdiv_diff)
# model siplification
a_mas_m2 <- lm(biomass_diff ~ strdiv_diff)
# values to report in the text
tapply(spp_num_day_biomass$biomass,spp_num_day_biomass$Plot,mean)
tapply(spp_num_day_biomass$biomass,spp_num_day_biomass$Plot,quantile,c(0.025,0.975))



######################################################################################
## 3) Arthropod composition and function
######################################################################################
# are exclosures more different among each other than controls?
bugcom <- aggregate(subset(bugs,select=-c(Site,Exclosure,Day,Trap,Trap.biomass)),by=list(bugs$Site,bugs$Exclosure,bugs$Day), FUN=sum)
colnames(bugcom)[1:3] <- c('Site','Exclosure','Day')
bugcom_mat_only <- bugcom[,-c(1:3)]

# drop empty species
bugcom_mat_only <- bugcom_mat_only[,which(colSums(bugcom_mat_only) != 0)]

# compare community composition
# first test null hypothesis either centroids and/or dispersion of the groups are equivalent for all groups
# rejection of null = either the centroid and/or the spread of samples differs between groups
adonis2(raupcrick(bugcom_mat_only)~bugcom$Exclosure+bugcom$Site+bugcom$Day,by='margin')
permutest(betadisper(raupcrick(bugcom_mat_only), bugcom$Exclosure))

# look at functional data
art_guilds$Decomposer <- with(art_guilds, Detritivore+Saprophyte+Scavenger)
art_guilds$Herbivore <- with(art_guilds, Fungivore+Herbivore)
art_guilds$Predator <- with(art_guilds, Predator+Parasitoid)
art_guilds <- subset(art_guilds,select=-c(Detritivore,Saprophyte,Scavenger,Parasitoid,Fungivore))

# calculate guild diversity
bugcom$herbS <- diversity(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Herbivore>0),1])])
bugcom$decoS <- diversity(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Decomposer>0),1])])
bugcom$predS <- diversity(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Predator>0),1])])
bugcom$herbA <- rowSums(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Herbivore>0),1])])
bugcom$predA <- rowSums(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Predator>0),1])])
bugcom$decoA <- rowSums(bugcom_mat_only[,which(colnames(bugcom_mat_only) %in% art_guilds[which(art_guilds$Decomposer>0),1])])
bugcom$he_pr <-  bugcom$herbA/bugcom$predA
bugcom$he_de <-  bugcom$herbA/bugcom$decoA
bugcom$de_pr <-  bugcom$decoA/bugcom$predA
bug_trop <- bugcom[,c(1:3,ncol(bugcom):(ncol(bugcom)-8))]
bug_trop$id <- with(bug_trop, paste(Site,Day))          
bug_trop <- reshape(bug_trop, idvar = "id", timevar = "Exclosure", direction = "wide")         
          
# just look at abundances of herbivores, predators, and decomposers
shapiro.test(with(bug_trop,he_pr.Exclosure-he_pr.Control))
shapiro.test(with(bug_trop,he_de.Exclosure-he_de.Control))
shapiro.test(with(bug_trop,de_pr.Exclosure-de_pr.Control))
 
summary(lm(with(bug_trop,he_pr.Exclosure-he_pr.Control)~bug_trop$Day.Control))
summary(lm(with(bug_trop,he_de.Exclosure-he_de.Control)~bug_trop$Day.Control))
summary(lm(with(bug_trop,de_pr.Exclosure-de_pr.Control)~bug_trop$Day.Control))

summary(lm(with(bug_trop,he_pr.Exclosure-he_pr.Control)~1))
summary(lm(with(bug_trop,he_de.Exclosure-he_de.Control)~1))
summary(lm(with(bug_trop,de_pr.Exclosure-de_pr.Control)~1))



######################################################################################
## 4) Vegetation structure (max height, number of unique and non-unique 20 cm intervals hit, structural diversity)
######################################################################################
# calculate max heights per transect
Max_ht_m <- with(veg_struc[which(veg_struc$Hits>0),],tapply(Height,list(Site,Treatment,Subsample.number),max))#/100+.1

# calculate proportion of unique heights that are hit (= index of diversity)
N_unique_hits <- with(veg_struc[which(veg_struc$Hits>0),],tapply(Height,list(Site,Treatment),function(x){length(unique(x))}))/with(veg_struc,tapply(Height,list(Site,Treatment),function(x){length(unique(x))}))

# calculate total proportion of all heights that are hit (= index of biomass)
N_hits <- with(veg_struc[which(veg_struc$Hits>0),],tapply(Height,list(Site,Treatment),length))/with(veg_struc,tapply(Height,list(Site,Treatment),length))

# merge everything together
veg_struc_plot <- as.data.frame(cbind(N_unique_hits,N_hits,Si_D))
colnames(veg_struc_plot) <- as.vector(sapply(c('N_u_hits','N_hits','Si_D'), function(x){c(paste(x,'_Con',sep=''),paste(x,'_Exc',sep=''))}))

# check normality of variables
with(veg_struc_plot, shapiro.test(c(N_u_hits_Con,N_u_hits_Exc)))
with(veg_struc_plot, shapiro.test(log(c(N_u_hits_Con,N_u_hits_Exc))))
with(veg_struc_plot, shapiro.test(c(N_hits_Con,N_hits_Exc)))

shapiro.test(1/Si_D)
shapiro.test(log(1/Si_D))

# perform t-tests
#with(veg_struc_plot, t.test(log(c(N_u_hits_Con,N_u_hits_Exc))~rep(c('C','E'),each=nrow(veg_struc_plot)),paired=T))
#with(veg_struc_plot, t.test(c(N_hits_Con,N_hits_Exc)~rep(c('C','E'),each=nrow(veg_struc_plot)),paired=T))
with(veg_struc_plot, t.test(log(c(1/Si_D[,1],1/Si_D[,2]))~rep(c('C','E'),each=nrow(Si_D)),paired=T))
apply(1/Si_D, 2, quantile, c(0.025,0.975))
colMeans(1/Si_D)


######################################################################################
## 5) Does soil C differ between plots?
######################################################################################
# find kg soil C in top 10 cm of each 7.2 x 7.2m plot by multiplying bulk density by volume
# and converting LOI to %C using standard ratio https://www.tandfonline.com/doi/full/10.1080/00103620500306080
soilC$kgC_plot <- with(soilC,drymass_g/vol_cm * 10 * 7.2*7.2 * 10000 * LOI/100 * 0.58 / 1000)

# reconvert using UK equation from https://doi.org/10.1007/BF00634106 
soilC$TC <- (soilC$LOI-3.627)/1.670 
# and convert to kg / m2 in top 10-cm 
soilC$kgC_plot2 <- with(soilC,drymass_g/vol_cm * 10 * 7.2*7.2 * 10000 * TC/100 / 1000)

# merge with historical soils
allsoilC <- rbind(soilC[,c('Site','Treatment','TC')], hist_soils[,c('Site','Treatment','TC')])
allsoilC$BA <- c(rep('later',nrow(soilC)),rep('before',nrow(hist_soils)))
allsoilC <- allsoilC[which(allsoilC$Site != 'Oakland'),]

# fit lmer on TC as historically didn't measure BD
shapiro.test(allsoilC$TC)
soilC_m1 <- lmer(TC ~ Treatment*BA + (1|Site), data=allsoilC)

# this carries on with modern (i.e. after-treatment) C - pooling samples
soilC_perplot <- with(soilC,tapply(kgC_plot,list(Site,Treatment),mean))
soilC_perplot2 <- with(soilC,tapply(kgC_plot2,list(Site,Treatment),mean))
shapiro.test(soilC_perplot)
shapiro.test(log(soilC_perplot))
t.test(log(soilC_perplot[,1]),log(soilC_perplot[,2]),paired=T)

# get kg/m2 estimate
apply(soilC_perplot/(7.2*7.2),2,quantile,c(0.025,0.975))
colMeans(soilC_perplot/(7.2*7.2))

# compare bulk density
bulk_d <- with(soilC,tapply(drymass_g,list(Site,Treatment),sum))/with(soilC,tapply(vol_cm,list(Site,Treatment),sum))
shapiro.test(bulk_d)
t.test(bulk_d[,1],bulk_d[,2],paired=T)


######################################################################################
## 6) Does above-ground C differ between plots?
######################################################################################
### calculate C content of trees using FC practices
# for saplings, defined as <7 cm DBH
# consider oak, blackthorn and sallow together as in FC Carbon Assessment Protocol (e.g. table 5.2.6)
treeC$C_kg <- rep(NA,nrow(treeC))
treeC[which(treeC$Species %in% c('Oak','Blackthorn','Sallow') & treeC$DBH_cm < 7),]$C_kg <- sapling_carbon_FC[match(treeC[which(treeC$Species %in% c('Oak','Blackthorn','Sallow') & treeC$DBH_cm < 7),'Ht_m'],sapling_carbon_FC$height_m),'C_stem_t']*1000

# for trees follow FC Carbon Assessment Protocol
treeC[which(treeC$Species %in% c('Oak','Blackthorn','Sallow') & treeC$DBH_cm >= 7),]$C_kg <- apply(treeC[which(treeC$Species %in% c('Oak','Blackthorn','Sallow') & treeC$DBH_cm >= 7),c('DBH_cm','Ht_m')],1,function(y){ broadleaf_C(DBH_cm=y[1], Ht_m=y[2]) })

# for blackthorn patches follow table 2 https://www.jstor.org/stable/42568965
treeC[which(treeC$Species == 'Blackthorn_patch'),'C_kg'] <- 696.7 * with(treeC[which(treeC$Species == 'Blackthorn_patch'),],Ht_m*Width_m*Length_m) * 0.5 /1000

# derive biomass density for rose patches
# here Georgie subsampled three bushes in much smaller volumes than observed in the field
# and assume that these relationships hold more generally
# convert biomass to C at 50%
rose_allometry <- as.data.frame(cbind(c(0.19125,0.07968,0.091875),c(116.64,76.26,85.07)))
colnames(rose_allometry) <- c('volume_m3','biomass_g')
treeC[treeC$Species == 'Rose_patch','C_kg'] <- with(treeC[treeC$Species == 'Rose_patch',], Ht_m*Width_m*Length_m * mean(rose_allometry[,2]/rose_allometry[,1]*.5/1000) )

# for individual stems take value of 0.20 mg/mm3 from doi:10.1007/s004420050563
treeC[treeC$Species == 'Rose','C_kg'] <- with(treeC[treeC$Species == 'Rose',], pi*(DBH_cm/100/2)^2*Ht_m * 0.20*1000*.5 )
                                                                                                 
# calculate total tree C per plot
# first add sites where no trees recorded
levels(treeC$Site) <- c(levels(treeC$Site),as.character(unique(soilC$Site[which(soilC$Site %in% treeC$Site == F)])))

treeC_perplot <- with(treeC, tapply(C_kg,list(Site,Treatment),sum))
treeC_perplot[is.na(treeC_perplot)] <- c(0)

                                                     
### fit different models to C content of herbaceous vegetation
herb_conv <- vector(mode = "list", length = 11)
herb_conv[[1]] <- lm(log(drymass_g) ~ percent_cover + ht_cm, data = gd_CC)
herb_conv[[2]] <- lm(log(drymass_g) ~ log(percent_cover) + ht_cm, data = gd_CC)
herb_conv[[3]] <- lm(log(drymass_g) ~ log(percent_cover) + log(ht_cm), data = gd_CC)
herb_conv[[4]] <- lm(log(drymass_g) ~ percent_cover + log(ht_cm), data = gd_CC)
herb_conv[[5]] <- lm(log(drymass_g) ~ log(percent_cover) * log(ht_cm), data = gd_CC)
herb_conv[[6]] <- lm(log(drymass_g) ~ percent_cover * log(ht_cm), data = gd_CC)
herb_conv[[7]] <- lm(log(drymass_g) ~ log(percent_cover) * ht_cm, data = gd_CC)
herb_conv[[8]] <- lm(log(drymass_g) ~ percent_cover * ht_cm, data = gd_CC)
herb_conv[[9]] <- lm(log(drymass_g) ~ I(percent_cover * ht_cm), data = gd_CC)
herb_conv[[10]] <- lm(log(drymass_g) ~ I(log(percent_cover * ht_cm)), data = gd_CC)
herb_conv[[11]] <- lm(log(drymass_g) ~ I(log(percent_cover * ht_cm)) + I(log(percent_cover * ht_cm)^2), data = gd_CC)
                                                                                     
# find best fitting model
which.min(sapply(herb_conv,function(x){AIC(x) + (2*(length(coefficients(x))^2)+2*length(coefficients(x)))/(15-length(coefficients(x))-1)}))

# extract table s2 information
AICc_vals <- round(sapply(herb_conv,function(x){AIC(x) + (2*(length(coefficients(x))^2)+2*length(coefficients(x)))/(15-length(coefficients(x))-1)}),2)
rel_like <- exp(-0.5 * AICc_vals)
AICc_wts <- round(rel_like/sum(rel_like),2)
R2vals <- round(sapply(herb_conv,function(x){summary(x)$adj.r.squared}),2)
write.csv(cbind(AICc_vals,AICc_wts,R2vals),file='tables2.csv')

# predict values in plots
grndC$C_kg <- exp(predict(herb_conv[[10]], grndC))*.5/1000

# calculate total herbaceous C per plot
grndC_perplot <- with(grndC, tapply(C_kg,list(Site,Treatment),sum))*7.2*7.2

# assess representativeness of sampling, i.e. does variance plateau?
herbc_plotsd <- sapply(unique(with(grndC,paste(Site,Treatment))),function(x){
					sapply(2:5,function(y){ mean(apply(combn(1:5,y),2,
							   function(z){ sd(grndC[which(with(grndC,paste(Site,Treatment)) == x),'C_kg'][z])/sqrt(y) }
							   )) 
							 })
							})   
plot(2:5,seq(min(herbc_plotsd),max(herbc_plotsd),length.out=4),type='n',xlab='sampled area (m2)', ylab = 'standard error',las=1,xaxt='n')
axis(1,at=c(2:5))
apply(herbc_plotsd,2,function(x){lines(2:5,x,col='gray80')})				
lines(2:5,rowMeans(herbc_plotsd),lwd=2)	


### compare total above-ground
tot_ag_C_perplot <- grndC_perplot[order(rownames(grndC_perplot)),] + treeC_perplot[order(rownames(treeC_perplot)),]
shapiro.test(log(tot_ag_C_perplot))
t.test(tot_ag_C_perplot[,1],tot_ag_C_perplot[,2],paired=T)
colMeans(tot_ag_C_perplot/(7.2*7.2))
apply(tot_ag_C_perplot/(7.2*7.2),2,quantile,c(0.025,0.975))

 
 
######################################################################################
## 7) Does total C differ between plots?
###################################################################################### 
tot_agbg_C_perplot <- soilC_perplot[order(rownames(soilC_perplot)),] + tot_ag_C_perplot[order(rownames(tot_ag_C_perplot)),]
shapiro.test(tot_agbg_C_perplot)
t.test(tot_agbg_C_perplot[,1],tot_agbg_C_perplot[,2],paired=T)
colMeans(tot_agbg_C_perplot/(7.2*7.2))
apply(tot_agbg_C_perplot/(7.2*7.2),2,quantile,c(0.025,0.975))



######################################################################################
## 8) Historical baseline comparison
######################################################################################

#################### woody species diversity
summary(trich_m1)


#################### vegetation biomass
# set <0.01 to minimum observed value - which is 0.01
hist_biom[,which(sapply(hist_biom, 'class') == 'factor')[-c(1:2)]] <- apply(hist_biom[,which(sapply(hist_biom, 'class') == 'factor')[-c(1:2)]],2,function(x){as.numeric(gsub('<0.01','0.01',x))})

# sum values per plot
hist_biom$summed_veg_dry_wt_g <- rowSums(hist_biom[,-c(1:5)])
hist_biom$allveg_biomass <- with(hist_biom, mixed_grasses_dry_wt_g + Unsorted_residue_dry_wt_g + summed_veg_dry_wt_g)
hist_biom$total_abg_biomass <- with(hist_biom, allveg_biomass + Dead_Organic_Matter_dry_wt_g)

# fit LMM beacuse we have 1 exclosure plot vs 2 controls, so somewhat pseudo-replicated 
# can't use straight t-test here because would be comparing each of the 2 open plots in a site to the same exclosure 
hist_biom$type <- substr(hist_biom$Plot,1,1)
abg_base_m1 <- with(hist_biom, lmer(log(total_abg_biomass) ~ type + (1|Site)))
summary(abg_base_m1)


#################### sward height
# NewBarn	Open_east	East has been changed from >1 to max value of 1.7
# convert from wide to long
hist_swht_long <- melt(hist_swht, id=c("FieldName","Plot","Quadrat"))
hist_swht_long <- hist_swht_long[,-which(colnames(hist_swht_long)%in%'variable')]
hist_swht_long$Plot <- with(hist_swht_long,paste(FieldName,Plot,sep='_'))
hist_swht_long$Treat <- sapply(strsplit(hist_swht_long$Plot,'_'), function(x){paste(x[2],x[3])})
hist_swht_long$Treat <- as.factor(hist_swht_long$Treat)
levels(hist_swht_long$Treat)[which(levels(hist_swht_long$Treat) %in% c('Open north','Open south'))] <- c('Open east','Open west')

# aggregate at plot-level         
mean_max_ht_m_old <- tapply(hist_swht_long$value, list(hist_swht_long$FieldName, hist_swht_long$Treat), mean, na.rm=T)
std_max_ht_m_old <- tapply(hist_swht_long$value, list(hist_swht_long$FieldName, hist_swht_long$Treat), sd, na.rm=T)

Max_ht_m[,1,][which(is.na(Max_ht_m[,1,]))] <- c(0)
mean_Max_ht_m <- cbind(rowMeans(Max_ht_m[,1,]),rowMeans(Max_ht_m[,2,]))
std_Max_ht_m <- cbind(apply(Max_ht_m[,1,],1,sd),apply(Max_ht_m[,2,],1,sd))
colnames(mean_Max_ht_m) <- colnames(std_Max_ht_m) <- colnames(Max_ht_m[,,1])

rownames(mean_max_ht_m_old) <- rownames(mean_Max_ht_m)
colnames(mean_max_ht_m_old) <- colnames(std_max_ht_m_old) <- c('Exclosure','Control','Control')

# combine before and after surveys
Max_ht_m_all <- melt(cbind(mean_max_ht_m_old*100,mean_Max_ht_m), id.vars=c('Exclosure','Control'))
Max_ht_m_all$stdev <- melt(cbind(std_max_ht_m_old*100,std_Max_ht_m), id.vars=c('Exclosure','Control'))$value
Max_ht_m_all$stdev[which(Max_ht_m_all$stdev==0)] <- min(Max_ht_m_all$stdev[which(Max_ht_m_all$stdev!=0)]) # one value = 0 and can't have that as an inverse weight so set to smallest observed value
Max_ht_m_all$invar <- 1/sqrt(Max_ht_m_all$stdev)
Max_ht_m_all$BA <- c(rep('before',nrow(mean_max_ht_m_old)*3),rep('later',nrow(mean_Max_ht_m)*2))
colnames(Max_ht_m_all)[1:2] <- c('Plot','Treat')
Max_ht_m_all$Treat2 <- as.factor(Max_ht_m_all$Treat == "Control")

# fit model and calculate values for main text
maxht_m1 <- lmer(log(value) ~ Treat2 * BA + (1|Plot), weights=invar, data=Max_ht_m_all)   # this is singular because RE ~ 0 (not something to worry about)
#maxht_m1 <- blmer(log(value) ~ Treat2 * BA + (1|Plot), weights=invar,data=Max_ht_m_all)  # so use blmer and find that parameter estimates are unchanged
maxht_m1 <- lm(log(value) ~ Treat2 * BA , weights=invar, data=Max_ht_m_all)   # drop Plot because of the singularity and to avoid problem with emmeans
maxht_m1_emmeans <- emmeans(maxht_m1,pairwise ~ Treat2|BA, type = "response")$emmeans			# save for plotting
	
# compare with changes estimated across fields from Knepp WildVeg Geodatabase (https://doi.org/10.5281/zenodo.5556944)
# drop site that is both exclosure and perm pasture so is confounded
height_by_field <- height_by_field[which(height_by_field$Field_name != 'Wild Flower Meadow'),]
height_by_field$troprew <- as.factor(height_by_field$permpast==0)
# save emmeans for plotting purposes
lidar_lm_emmeans <- emmeans(lm(log(X_mean*100) ~ troprew,weights=1/sqrt(X_stdev),data=height_by_field),'troprew',type='response')
     
#################### vegetation community composition
# subset just the community data for the 5 quadrats per plot 
hist_veg_comm <- hist_gcov[which(hist_gcov$Quadrat != 'rest_20by20'),-c(1:4)]

# replace p with minimum value of 0.1 as per what USDA does for 'trace' cover (https://www.fs.fed.us/pnw/pubs/pnw_gtr781.pdf; and in AUS https://onlinelibrary.wiley.com/doi/full/10.1111/avsc.12437) 
# and others that modify simple presence on Braun-Blanquet scale (https://doi.org/10.1111/j.1654-1103.2004.tb02321.x)  
hist_veg_comm[,which(sapply(hist_veg_comm, 'class') == 'factor')] <- apply(hist_veg_comm[,which(sapply(hist_veg_comm, 'class') == 'factor')],2,function(x){as.numeric(gsub('p','0.1',x))})

# sum same species between understorey and canopy
mnames2 <- sapply(colnames(hist_veg_comm),function(x){ grep(x,colnames(hist_veg_comm)) })
mnames_rep2 <- apply(sapply(which(sapply(mnames2,function(x){length(x>1)})==2), function(x){mnames2[[x]]}),2,function(y){ rowSums(hist_veg_comm[,y]) })
hist_veg_comm <- hist_veg_comm[,-as.vector(sapply(which(sapply(mnames2,function(x){length(x>1)})==2), function(x){mnames2[[x]]}))]
hist_veg_comm <- cbind(hist_veg_comm,mnames_rep2)

# aggregate quadrats in plots to reduce pseudoreplication driving dispersion results
hist_veg_comm_plot <- apply(hist_veg_comm, 2, function(x){ tapply(x, paste0(hist_gcov[which(hist_gcov$Quadrat != 'rest_20by20'),2], hist_gcov[which(hist_gcov$Quadrat != 'rest_20by20'),1]), sum)/4 })
hist_veg_comm_plot <- hist_veg_comm_plot[,which(colSums(hist_veg_comm_plot) != 0)]

# test for differences in dispersion (i.e. variance) of communities 
#permutest(betadisper(vegdist(hist_veg_comm), substr(hist_gcov[which(hist_gcov$Quadrat != 'rest_20by20'),]$Plot,1,1)))
permutest(betadisper(raupcrick(hist_veg_comm_plot), substr(rownames(hist_veg_comm_plot),1,1)))

# check for difference in location (i.e. mean) of communities
#adonis2(vegdist(hist_veg_comm) ~ substr(hist_gcov[which(hist_gcov$Quadrat != 'rest_20by20'),]$Plot,1,1), perm = 999)
adonis2(raupcrick(hist_veg_comm_plot) ~ substr(rownames(hist_veg_comm_plot),1,1), perm = 999)
			
# run screeplot                  
NMDS.scree(hist_veg_comm_plot)
dev.off()

# fit k-dimensional NMDS 
veg_NMDS <- metaMDS(raupcrick(hist_veg_comm_plot),k=3)        
          

#################### historical soils
# Three samples taken per location (either side of exclosure and another)
# Soil cores taken using dutch auger into dry soils to depth 0.1-0.15m  I.E. SAME AS YOUR STUDENTS
soil_pH <- lmer(pH ~ Treatment+(1|Site), data=hist_soils)

# only use TOTAL COMBUSTION USING DUMAS MACHINE as problem with OC based on emails with M Heard
soilC_old <- lmer(TC ~ Treatment+(1|Site), data=hist_soils)

# other macronutrients
soilTN_old <- lmer(TN ~ Treatment+(1|Site), data=hist_soils)
soilK_old <- lmer(log(K) ~ Treatment+(1|Site), data=hist_soils)
soilP_old <- lmer(log(P) ~ Treatment+(1|Site), data=hist_soils)

                  
######################################################################################
## 9) Prepare display items
######################################################################################
# build tables for control
control_metrics <- as.data.frame(cbind(soilC_perplot[,1],tot_ag_C_perplot[,1],tot_agbg_C_perplot[,1]))/(7.2*7.2)                   
colnames(control_metrics) <- c('BG_solC', 'AG_vegC', 'totC')
control_metrics$Site <- rownames(control_metrics)

control_metrics2 <- spp_num_day_biomass[spp_num_day_biomass$Plot=='Control',c('Site','sppdiv','biomass')]
control_metrics2[,3] <- control_metrics2[,3]*100 # do this so scaling works properly later, i.e. don't have negative values
control_metrics2[,2:3] <- log(control_metrics2[,2:3])
colnames(control_metrics2)[2:3] <- c('ArthroRich','ArthroMass')
levels(control_metrics2$Site) <- c("27Acres","Broomers","Oakland","Rainbow" ,"Wildflower")

control_metrics3 <- bug_trop[,c('de_pr.Control','he_de.Control','he_pr.Control')]
control_metrics3$Site <- as.factor(bug_trop$Site.Control)
levels(control_metrics3$Site) <- c("27Acres","Broomers","Oakland","Rainbow" ,"Wildflower")

control_metrics4 <- as.data.frame(veg_struc_plot[,'Si_D_Con'])
control_metrics4[,1] <- log(1/control_metrics4[,1])
control_metrics4$Site <- as.factor(rownames(control_metrics4))
levels(control_metrics4$Site) <- c("27Acres",'Brook5','Brook6','Brookbuild',"Broomers",'NewBarn',"Oakland","Rainbow" ,"Wildflower",'WaterWorks')

control_metrics5 <- wood_rich_merge[wood_rich_merge$BA=='later'&wood_rich_merge$Treat=='C',c('Site','wood_rich')]

# build tables for exclosures
exclosure_metrics <- as.data.frame(cbind(soilC_perplot[,2],tot_ag_C_perplot[,2],tot_agbg_C_perplot[,2]))/(7.2*7.2)                   
colnames(exclosure_metrics) <- c('BG_solC', 'AG_vegC', 'totC')
exclosure_metrics$Site <- rownames(exclosure_metrics)

exclosure_metrics2 <- spp_num_day_biomass[spp_num_day_biomass$Plot=='Exclosure',c('Site','sppdiv','biomass')]
exclosure_metrics2[,3] <- exclosure_metrics2[,3]*100 
exclosure_metrics2[,2:3] <- log(exclosure_metrics2[,2:3])
colnames(exclosure_metrics2)[2:3] <- c('ArthroRich','ArthroMass')
levels(exclosure_metrics2$Site) <- c("27Acres","Broomers","Oakland","Rainbow" ,"Wildflower")

exclosure_metrics3 <- bug_trop[,c('de_pr.Exclosure','he_de.Exclosure','he_pr.Exclosure')]
exclosure_metrics3$Site <- as.factor(bug_trop$Site.Exclosure)
levels(exclosure_metrics3$Site) <- c("27Acres","Broomers","Oakland","Rainbow" ,"Wildflower")

exclosure_metrics4 <- as.data.frame(veg_struc_plot[,'Si_D_Exc'])
exclosure_metrics4[,1] <- log(1/exclosure_metrics4[,1])
exclosure_metrics4$Site <- as.factor(rownames(exclosure_metrics4))
levels(exclosure_metrics4$Site) <- c("27Acres",'Brook5','Brook6','Brookbuild',"Broomers",'NewBarn',"Oakland","Rainbow" ,"Wildflower",'WaterWorks')

exclosure_metrics5 <- wood_rich_merge[wood_rich_merge$BA=='later'&wood_rich_merge$Treat=='E',c('Site','wood_rich')]

# create master table of means for each treatment
control_means <- unlist(sapply(list(control_metrics,control_metrics2,control_metrics3,control_metrics4,control_metrics5),function(x){colMeans(subset(x,select=-c(Site)))}))
exclosure_means <- unlist(sapply(list(exclosure_metrics,exclosure_metrics2,exclosure_metrics3,exclosure_metrics4,exclosure_metrics5),function(x){colMeans(subset(x,select=-c(Site)),na.rm=T)}))
all_mets_treat_means <- as.data.frame(cbind(control_means,exclosure_means))
all_mets_treat_means[(1+nrow(all_mets_treat_means)),] <- c(mean(Max_ht_m_all[Max_ht_m_all$Treat=='Control' & Max_ht_m_all$BA=='later','value']),mean(Max_ht_m_all[Max_ht_m_all$Treat=='Exclosure' & Max_ht_m_all$BA=='later','value']))
rownames(all_mets_treat_means)[c(9,11)] <- c('Si_D_Con','Max_ht_Con')

# collapse into means for each treatment
all_mets_treat_means_scaled <- apply(t(all_mets_treat_means), 2, function(x) {x/sum(x)})
all_mets_treat_means_scaled <- all_mets_treat_means_scaled[,c('wood_rich','ArthroRich','ArthroMass','he_pr.Control','he_de.Control','de_pr.Control','Si_D_Con','Max_ht_Con','BG_solC','AG_vegC')]
all_mets_treat_means_scaleddf <- data.frame(melt(all_mets_treat_means_scaled, varnames=c("treatment", "metric"), value.name = "rate"))


######################### Figure 3                                                  
colors <- pals::parula(10)[c(1,4,7,9)]
spot.theme <- list(
  theme_classic(),
  theme(axis.ticks.x=element_blank(), axis.text.x=element_text(size = 19, angle = 90, hjust = 0)),
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size = 19)),
  theme(axis.line=element_blank()),
  theme(text = element_text(size = 22)),
  theme(legend.position = "none"),
  theme(plot.margin = unit(c(10,10,10,10), "mm")),
  scale_size_continuous(range = c(-0.3, 15)),
  scale_x_discrete(position = "top"))

#pdf("fig3_j23.pdf", height = 5.767*1.2, width = 8.549*1.2, useDingbats = F)
ggplot(all_mets_treat_means_scaleddf, aes(metric, treatment)) + spot.theme + geom_point(colour = colors[[1]], aes(size = rate))
#dev.off()

######################### Figure 1         
with(animals, plot(as.numeric(Date),LonghornCattle/4.5,type='b',pch=19,las=1,xlab='Year',ylab='Density (# km2)',bty='n',ylim=c(0,60),xaxt='n',cex=1.2))
with(animals, points(as.numeric(Date),ExmoorPonies/4.5,type='b',pch=19,cex=1.2))
with(animals, points(as.numeric(Date),TamworthPigs/4.5,type='b',pch=19,cex=1.2,col='gray40'))
with(animals, points(as.numeric(Date),RedDeer/4.5,type='b',pch=17,col='gray',cex=1.2))
with(animals, points(as.numeric(Date),FallowDeer/4.5,type='b',pch=17,col='gray',cex=1.2))
axis(1,at=1:9,lab=2010:2018)

# add overall LUs from Knepp Feasability Assessment
# use first set of LUs from pg 136 of the Assessment
LU_v1 <- as.vector(as.matrix(animals_by_class_v1[-10,-1]) %*% c(0.8,0.3,0.6,0.8,0.4,0.15,0.3,0.4,0.2,0.05,0.15,0.15,0.2,0.05,0.15,0.2,0.2,0.05,0.15,0.2))
points(as.numeric(animals$Date),LU_v1/4.5,type='l',lwd=2)

# use second set of LUs from pg 151
LU_v2 <- as.vector(as.matrix(animals_by_class_v2[-10,-1]) %*% c(0.9,0.8,0.65,0.2,0.65,0.8,0.8,0.7,0.8,0.1,0.3,0.35,0.5,0.2,0.2,0.15,0.15))
points(as.numeric(animals$Date),LU_v2/4.5,type='l',lwd=2,col='red')

# calculate CH4 emissions (t) using production per individual per year (kg)
as.matrix(animals[,-1]) %*% c(15,15,1.5,55,7.2)/1000
# convert to kg CO2e and then convert to C and divide by study area
(as.matrix(animals[,-1]) %*% c(15,15,1.5,55,7.2)) * 25 / 450
# convert to kg CO2e and then divide by study area and add to difference between treatments in total C in ha-1 yr-1
(sum((as.matrix(animals[,-1]) %*% c(15,15,1.5,55,7.2)) * 25 / 450) + 1.8*10000*44/12)/9 	# with CH4
(1.8*10000*44/12)/9																			# without CH4


# calculate animal removals
range(c(LU_v2,LU_v1))/450
as.vector(as.matrix(animals_by_class_v1[10,-1]) %*% c(0.8,0.3,0.6,0.8,0.4,0.15,0.3,0.4,0.2,0.05,0.15,0.15,0.2,0.05,0.15,0.2,0.2,0.05,0.15,0.2))/450
as.vector(as.matrix(animals_by_class_v2[10,-1]) %*% c(0.9,0.8,0.65,0.2,0.65,0.8,0.8,0.7,0.8,0.1,0.3,0.35,0.5,0.2,0.2,0.15,0.15))/450

######################### Figure 4 
bug_ts <- bugcom[,c('Exclosure','Day','decoA','herbA','predA')]
bug_ts$Treat_Day <- with(bug_ts, paste(Exclosure,Day))
bug_ts <- cbind(with(bug_ts, tapply(decoA,Treat_Day,mean)),with(bug_ts, tapply(herbA,Treat_Day,mean)),with(bug_ts, tapply(predA,Treat_Day,mean)))
barplot(t(bug_ts)[,c(1,3,2,4)],
        col=colors()[c(23,89,12)] ,
        border="white",
        space=0.04,
        las = 1,
        ylab='Individuals per trap',
        xlab="group")


                  
######################### Figure S1
par(mfrow = c(1,3), mar =c(4.5,4.5,.5,.5))
sapply(1:3, function(x){
vegMDS.plot <- plot(veg_NMDS, type = "none", 
                    choices=as.numeric(expand.grid(c(1,2,3), c(1,2,3))[c(2,3,6),2:1][x,]), 
                    xlim = c(-0.6,0.6), ylim = c(-0.6,0.6), cex.axis = 1.2, cex.lab = 1.2 )
points(vegMDS.plot, "sites", cex = 3,pch=21,
       bg = c("gray20", "gray90")[as.numeric(as.factor(substr(rownames(hist_veg_comm_plot),1,1)))] )
ordiellipse(vegMDS.plot, substr(rownames(hist_veg_comm_plot),1,1), display="sites", kind = "ehull", draw = "polygon", border = "black", alpha = 0, lwd = 1 )
})
    
    

######################### Figure 5
rescale_fun <- function(x,a,b){(b-a)*((x-min(x))/(max(x)-min(x)))+a}

#pdf("fig5.pdf", height = 0.951*4, width = 0.829*4, useDingbats = F)
plot(runif(nrow(height_by_field[which(height_by_field$permpast == 1),]),.8,1.2),height_by_field[which(height_by_field$permpast == 1),]$X_mean*100,
		xlim=c(0,5),ylim=c(.1,1000),ylab='height (cm)',xaxt='n', xlab='', log='y', col='gray90',pch=19,bty='n',las=1,
		cex=rescale_fun(1/sqrt(height_by_field$X_stdev),0.5,2)[which(height_by_field$permpast == 1)])
axis(1,at=1:4, lab=c('Control','Trophic','Trophic','Passive'),las=2)
points(runif(nrow(height_by_field[which(height_by_field$permpast == 0),]),1.8,2.2),height_by_field[which(height_by_field$permpast == 0),]$X_mean*100,pch=19,col='gray75',
		cex=rescale_fun(1/sqrt(height_by_field$X_stdev),0.5,2)[which(height_by_field$permpast == 0)])
points(runif(length(Max_ht_m_all[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Control'),'value']),2.8,3.2),Max_ht_m_all[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Control'),'value'],pch=19,col='gray60',
		cex=rescale_fun(Max_ht_m_all$invar,0.5,2)[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Control')])
points(runif(length(Max_ht_m_all[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Exclosure'),'value']),3.8,4.2),Max_ht_m_all[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Exclosure'),'value'],pch=19,col='gray75',
		cex=rescale_fun(Max_ht_m_all$invar,0.5,2)[which(Max_ht_m_all$BA=='later' & Max_ht_m_all$Treat=='Exclosure')])
sapply(1:2, function(x){ lines(rep(x,2),c(summary(lidar_lm_emmeans)$lower.CL[x],summary(lidar_lm_emmeans)$upper.CL[x]),lwd=2)})
points(1:2,summary(lidar_lm_emmeans)$response,pch=15,cex=1.4,col='black')
sapply(3:4, function(x){ lines(rep(x,2),c(summary(maxht_m1_emmeans)$lower.CL[4:3][x-2],summary(maxht_m1_emmeans)$upper.CL[4:3][x-2]),lwd=2)})
points(3:4,summary(maxht_m1_emmeans)$response[4:3],pch=15,cex=1.4,col='black')
#dev.off()
