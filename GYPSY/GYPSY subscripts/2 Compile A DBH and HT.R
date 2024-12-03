## Dead trees are deleted
## If trees have DBH and height < 1.3m, height is estimated
## If trees have a missing DBH and height < 1.3m, sent to regen tally table
#/**********************************************************/
# Append measurement info to trees data and make adjustments
#/**********************************************************/
trees_merge <- fread("GYPSY data/intermediate/i_trees_merge.csv")
plot_measurement <- fread("GYPSY data/intermediate/i_plot_measurement.csv")
trees <- trees_merge[plot_measurement,on=c("company","company_plot_number","measurement_number")]
plot <- fread("GYPSY data/intermediate/i_plot.csv")[,.(company,company_plot_number,natural_subregion)]
plot <- plot[fread("GYPSY data/lookup/natsub.csv"), on="natural_subregion"]
trees <- trees[plot,on=c("company","company_plot_number")]
rm(trees_merge,plot_measurement);gc()

trees[,species_og:=toupper(species)]
trees[,dbh_og:=dbh]
trees[,height_og:=height]
trees[,unique:=paste(company,company_plot_number,tree_number,measurement_number,sep="_")]

trees[height<1.3 & !is.na(dbh),height:=NA] # Predict new height where DBH is measured but height is less than 1.3

# Fix suspected 0 filling
trees[dbh==0,dbh:=NA]
trees[height==0,height:=NA]

trees[,species:=NULL]
trees <- trees[fread("GYPSY data/lookup/species.csv"),on="species_og"] # Joins species and species group from a lookup table

trees[is.na(species),species:=species_og] # Use original species if no replacement
trees[,dbh_stat:="XX"]
trees[,ht_stat:="XX"]
trees[!is.na(dbh),dbh_stat:="MM"]
trees[height<1.3,dbh_stat:="NA"]
trees[!is.na(height),ht_stat:="MM"]
trees <- trees[!(species%in%c("NO","MS"))]
#trees <- trees[!(tree_type=="B" & is.na(height))] # Excludes age trees without heights
trees <- trees[!(!is.na(height) & height<0.3)]

#trees <- trees[!(species%in%c("AW","BW","PB") & (!is.na(height) & height<1.3)),] # Excludes deciduous under 1.3m

trees <- trees[!(condition_code1%in%c(1,2,13,14,15)),]
trees <- trees[!(species%in%c("DC","DD","DU"))] # could be combined with MS and NO
trees <- trees[tree_type!=""]
trees <- trees[,!c("measurement_month","measurement_day","stand_type",
                   "tree_plot_shape","sapling_plot_shape","regen_plot_shape",
                   "contractor","cruiser_1_name","cruiser_2_name","shrub_cover", 
                   "herb_forb_cover","grass_cover","moss_lichen_cover",
                   "avi_field_call","plot_measurement_comment")];gc()

#### Fill missing diameters

## Predict missing DBH from height and add dbh_stat=
#'MM' (measured), 
#'PA' (predicted, adjusted), 
#'PU' (predicted unadjusted, 
#'IP' (inferred from previous mmt), 
#'IS' (inferred from subsequent mmt), 
#'IX' (inferred using growth rate derived from previous & subsequent), 
#'IA' (inferred average by species / tree_type),
#'IF' (inferred manual fill ht=1.3 set dbh=0.1), 
#'NA' (not applicable under 1.3 m), 
#'IM' (inferred using DBH from minimum tagging limit)
# IA should only be assigned in cases where stem was only measured once, and either both height and DBH was missing, or stem was under 1.3 m tall and height was missing
trees[fread("GYPSY data/lookup/ht_to_dbh_2016r1.csv"),on=c("species","natsub"),dbh_p:=b1*(height-1.3)^b2*exp(-b3*(height-1.3))]

trees[!is.na(dbh) & !is.na(dbh_p) & condition_code1!=3 & !(condition_code1%in%c(5,10) & severity1==3), 
      rate:=mean(dbh/dbh_p),
      by=.(company,company_plot_number,measurement_number,species)]

trees[!is.na(dbh_p) & dbh_stat=="XX" & height>1.3 & !is.na(rate),dbh_PA:=dbh_p*rate]
trees[!is.na(dbh_p) & dbh_stat=="XX" & height>1.3 & is.na(rate),dbh_PU:=dbh_p]
trees[,c("dbh_p","rate"):=NULL]

## Try to fill in missing DBHs based on one of: previous DBH, subsequent DBH, or previous + growth rate (if available)
trees[,tree_id:=paste(company,company_plot_number,tree_number,sep="_")]

loop <- trees %>%
        group_by(tree_id) %>%
        filter(length(unique)>1 &
               any(dbh_stat == "XX", na.rm=T) &
               (any(dbh_stat == "MM", na.rm=T) |
                any(dbh_stat == "NA", na.rm=T))) %>%
  data.table

loop <- loop[,unique(tree_id)]

temp <- trees[tree_id%in%loop,c("tree_id","unique","measurement_year","dbh_stat","dbh")] # Smaller subset to improve processing time
temp[,c("lm_yr","nm_yr","l_dbh","n_dbh"):=as.numeric(NA)]

dbh_check <- function(df,i){
  dat <- df[tree_id==i][order(measurement_year)]
  for(j in 1:dat[,.N]){
    if(dat[j,dbh_stat]=="XX"){
      m_yr <- dat[j,measurement_year]
      p_meas <- dat[1:j]
      f_meas <- dat[j:.N]
      if("NA"%in%f_meas[,dbh_stat]){
        dat[j,dbh_stat:="NA"]
      } 
      if("MM"%in%p_meas[,dbh_stat]){
        lm.yr <- tail(p_meas[dbh_stat=="MM"],1)
        dat[j,lm_yr:=lm.yr[,measurement_year]]
        dat[j,l_dbh:=lm.yr[,dbh]]
      } 
      if("MM"%in%f_meas[,dbh_stat]){
        nm.yr <- head(f_meas[dbh_stat=="MM"],1)
        dat[j,nm_yr:=nm.yr[,measurement_year]]
        dat[j,n_dbh:=nm.yr[,dbh]]
      }
    }
  }
  dat
}

cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(data.table))
clusterExport(cl, c("temp","loop","dbh_check"))
temp <- rbindlist(parLapply(cl,loop,function(i) dbh_check(temp,i)))
stopCluster(cl)

temp[!is.na(l_dbh),
     dbh_IP:=l_dbh]
temp[!is.na(n_dbh),
     dbh_IS:=n_dbh]
temp[,dbh_IX:=if_else(l_dbh > n_dbh, NA, round((n_dbh-l_dbh)/(nm_yr-lm_yr)*(measurement_year-lm_yr)+l_dbh,1))]

trees[temp, on="unique", c("dbh_stat","dbh_IP","dbh_IS","dbh_IX"):=list(i.dbh_stat,dbh_IP,dbh_IS,dbh_IX)]

## A number of observations with missing height and diameter, no previous or subsequent measurement from which to populate, use average observed for that plot/mmt/species/tree_type
trees[,dbh_IA:=mean(dbh,na.rm=T),by=.(company,company_plot_number,measurement_number,tree_type,spp_grp)]

rm("temp");gc()

## Pick which method hierarchy to use to assign missing DBH values
trees[dbh_stat=="XX" & height<1.3, c("dbh","dbh_stat"):=list(0,NA)]
trees[dbh_stat=="XX" & height==1.3 & !is.na(dbh_IP),c("dbh","dbh_stat"):=list(dbh_IP,"IP")]
trees[dbh_stat=="XX" & height==1.3 & species%in%c("AW","BW","PB","LT","PJ","PL"),c("dbh","dbh_stat"):=list(0.2,"IF")]
trees[dbh_stat=="XX" & height==1.3 & !species%in%c("AW","BW","PB","LT","PJ","PL"),c("dbh","dbh_stat"):=list(0.3,"IF")]
trees[dbh_stat=="XX" & height>1.3 & !is.na(dbh_IX),c("dbh","dbh_stat"):=list(dbh_IX,"IX")]
trees[dbh_stat=="XX" & height>1.3 & !is.na(dbh_PA) & ((dbh_PA>=dbh_IP | is.na(dbh_IP)) & (dbh_PA<=dbh_IS | is.na(dbh_IS))),c("dbh","dbh_stat"):=list(dbh_PA,"PA")]
trees[dbh_stat=="XX" & height>1.3 & !is.na(dbh_IP),c("dbh","dbh_stat"):=list(dbh_IP,"IP")]
trees[dbh_stat=="XX" & height>1.3 & !is.na(dbh_IS),c("dbh","dbh_stat"):=list(dbh_IS,"IS")]
trees[dbh_stat=="XX" & height>1.3 & !is.na(dbh_PU),c("dbh","dbh_stat"):=list(dbh_PU,"PU")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","ES1") & !is.na(trees$dbh_IP),c("dbh","dbh_stat"):=list(dbh_IP,"IP")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","ES1") & is.na(trees$dbh_IP),c("dbh","dbh_stat"):=list(0,"NA")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("T","S1","ET") & !is.na(dbh_IX),c("dbh","dbh_stat"):=list(dbh_IX,"IX")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("T","S1","ET") & !is.na(dbh_IP),c("dbh","dbh_stat"):=list(dbh_IP,"IP")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("T","S1","ET") & !is.na(dbh_IS),c("dbh","dbh_stat"):=list(dbh_IS,"IS")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("T","S1","ET") & !is.na(dbh_IA),c("dbh","dbh_stat"):=list(dbh_IA,"IA")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type%in%c("S1","ET") & !is.na(sapling_tagging_limit_dbh),c("dbh","dbh_stat"):=list(sapling_tagging_limit_dbh,"IM")]
trees[dbh_stat=="XX" & ht_stat=="XX" & tree_type=="T" & !is.na(tree_tagging_limit),c("dbh","dbh_stat"):=list(tree_tagging_limit,"IM")]

trees[,dbh:=round(dbh,1)]
trees[dbh_stat=="PU" & dbh==0 & height>1.3,c("dbh","dbh_stat"):=list(0.1,"IF")]

trees <- trees[,!c("dbh_PA","dbh_PU","dbh_IP","dbh_IS","dbh_IA","dbh_IX")]
table(trees[,c("tree_type","dbh_stat")])

#### Fill missing heights

## Predict missing height from DBH and add ht_stat=
#'MM' (measured) 
#'PA' (predicted, adjusted)
#'PU' (predicted unadjusted
#'IP' (inferred from previous mmt)
#'IS' (inferred from subsequent mmt)
#'IX' (inferred using growth rate derived from previous & subsequent)
#'IA' (inferred average by species / tree_type) 
#'IM' (inferred from minimum tagging limit)
# IA and IM should only be assigned in cases where stem was under 1.3 m tall and height was missing

# Predict heights for all trees including ones that already have a height
trees[fread("GYPSY data/lookup/dbh_to_ht_2016r4.csv"),on=c("natsub","species"),
      ht_PU:=ifelse(model==1,1.3+b1*(1-exp(-1*b2*dbh))^b3,1.3+(b1/(1+exp(b2+b3*log(dbh+0.01)))))]

## Calculate ratio of means to create adjustment factor for height
# Don't want trees with missing actual heights to be used in averaging
# Also don't want trees with unusual (i.e. data issues) height-DBH relationships to affect ratios
trees[!is.na(height) & !is.na(dbh),hdr:=height/dbh]

# Attempt 1 - ratio by species, plot and measurement
trees[!is.na(ht_PU) &
      !(height>5 & hdr>3) & 
      !(height>5 & hdr<0.3) & 
      !((height<=5 & height>=2) & hdr>5) &
      !((height<=5 & height>=2) & hdr<0.3) &
      !(height<2 & hdr>15) &
      !(height<2 & hdr<0.3) &
      !(condition_code1%in%c(1:3,13:15)) & # Dead, broken or dead top, etc.
      !(condition_code1==10 & severity1==3), # Severe lean
      c("hr","freq","ratio"):=list(T,
                                   .N,
                                   mean(height)/mean(ht_PU)),
      by=.(company,company_plot_number,measurement_number,species)]
trees[,c("freq","ratio"):=list(mean(freq,na.rm=T),mean(ratio,na.rm=T)),by=.(company,company_plot_number,measurement_number,species)]
trees[freq>3 & ht_stat=="XX", ht_PA1:= ht_PU*ratio]
trees <- trees[,!c("hdr","ratio","freq")]

# Attempt 2 - ratio by species and plot only (to fill in missing after attempt #1)
trees[hr==T,
      c("freq","ratio"):=list(.N,
                              mean(height)/mean(ht_PU)),
      by=.(company,company_plot_number,species)]
trees[,c("freq","ratio"):=list(mean(freq,na.rm=T),mean(ratio,na.rm=T)),by=.(company,company_plot_number,species)]
trees[freq>3 & ht_stat=="XX",ht_PA2:=ht_PU*ratio]

trees <- trees[,!c("hr","ratio","freq")]

# Height from previous measurement
loop <- trees %>%
        group_by(tree_id) %>%
        filter(any(ht_stat == "XX", na.rm=T) &
               any(ht_stat == "MM", na.rm=T))
loop <- unique(loop$tree_id)
temp <- trees[tree_id%in%loop,c("tree_id","unique","measurement_year","ht_stat","height")] # Smaller subset to improve processing time
temp[,c("lm_yr","nm_yr","l_ht","n_ht"):=as.numeric(NA)]

ht_check <- function(df,i){
  dat <- df[tree_id==i][order(measurement_year)]
  for(j in 1:dat[,.N]){
    if(dat[j,ht_stat]=="XX"){
      m_yr <- dat[j,measurement_year]
      p_meas <- dat[1:j]
      f_meas <- dat[j:.N]
      if("MM"%in%p_meas[,ht_stat]){
        lm.yr <- tail(p_meas[ht_stat=="MM"],1)
        dat[j,lm_yr:=lm.yr[,measurement_year]]
        dat[j,l_ht:=lm.yr[,height]]
      } 
      if("MM"%in%f_meas[,ht_stat]){
        nm.yr <- head(f_meas[ht_stat=="MM"],1)
        dat[j,nm_yr:=nm.yr[,measurement_year]]
        dat[j,n_ht:=nm.yr[,height]]
      }
    }
  }
  dat
}

cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(data.table))
clusterExport(cl, c("temp","loop","ht_check"))
temp <- rbindlist(parLapply(cl,loop,function(i) ht_check(temp,i)))
stopCluster(cl)

# calculating height based on previous and next year
temp[!is.na(l_ht) & is.na(n_ht),
     ht_IP:=l_ht]
temp[!is.na(n_ht) & is.na(l_ht),
     ht_IS:=n_ht]
temp[,ht_IX:=if_else(l_ht > n_ht, NA, round((n_ht-l_ht)/(nm_yr-lm_yr)*(measurement_year-lm_yr)+l_ht,1))]

trees <- temp[,.(unique,ht_IP,ht_IS,ht_IX)][trees, on="unique"]

# Average height by species / measurement / tree type
trees[,ht_IA:=mean(height,na.rm=T), by=.(company,company_plot_number,measurement_number,tree_type,spp_grp)]
trees[ht_stat=="MM",ht_IA:=NA]

# Assignment rules for which height to select
trees[ht_stat=="XX" & dbh>0 & ht_IX>=1.3,c("height","ht_stat"):=list(ht_IX,"IX")]

trees[ht_stat=="XX" & dbh>0 & ht_PA1 >=1.3,c("height","ht_stat"):=list(ht_PA1,"PA")]

trees[ht_stat=="XX" & dbh>0 & ht_PA2>=1.3,c("height","ht_stat"):=list(ht_PA2,"Pa")]

trees[ht_stat=="XX" & dbh>0 & ht_IP>=1.3,c("height","ht_stat"):=list(ht_IP,"IP")]

trees[ht_stat=="XX" & dbh>0 & ht_PU>=1.3,c("height","ht_stat"):=list(ht_PU,"PU")]

trees[ht_stat=="XX" & dbh>0 & ht_IS>=1.3,c("height","ht_stat"):=list(ht_IS,"IS")]

trees[ht_stat=="XX" & dbh==0 & ht_IX<1.3,c("height","ht_stat"):=list(ht_IX,"IX")]

trees[ht_stat=="XX" & dbh==0 & ht_IX<1.3,c("height","ht_stat"):=list(ht_IP,"IP")]

trees[ht_stat=="XX" & dbh==0 & ht_IX<1.3,c("height","ht_stat"):=list(ht_IA,"IA")]

trees[ht_stat=="XX" & dbh==0 & tree_type%in%c(paste0("R",1:10),"ES1") & !species%in%c("AW","BW","PB") & !is.na(regen_tagging_limit_conifer),c("height","ht_stat"):=list(regen_tagging_limit_conifer,"IM")]

trees[ht_stat=="XX" & dbh==0 & tree_type%in%c(paste0("R",1:10),"ES1") & species%in%c("AW","BW","PB") & !is.na(regen_tagging_limit_decid),c("height","ht_stat"):=list(regen_tagging_limit_decid,"IM")]
  
trees[ht_stat=="XX" & dbh==0 & tree_type%in%c(paste0("R",1:10),"ES1") & !species%in%c("AW","BW","PB") & is.na(regen_tagging_limit_conifer),c("height","ht_stat"):=list(0.3,"IM")]

trees[ht_stat=="XX" & dbh==0 & tree_type%in%c(paste0("R",1:10),"ES1") & species%in%c("AW","BW","PB") & is.na(regen_tagging_limit_decid),c("height","ht_stat"):=list(0.3,"IM")]

trees[,height:=round(height,2)]

trees[dbh==0 & height>=1.3, dbh:=0.1]

trees <- trees[,!c("ht_PU","ht_PA1","ht_PA2","ht_IP","ht_IS","ht_IA","ht_IX","dbh_og","height_og")]
  
fwrite(trees,"GYPSY data/intermediate/i_trees_measurement_filled.csv")