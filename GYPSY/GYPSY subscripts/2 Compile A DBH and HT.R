## Dead trees are deleted
## If trees have DBH and height < 1.3m, height is estimated
## If trees have a missing DBH and height < 1.3m, sent to regen tally table
#/**********************************************************/
# Append measurement info to trees data and make adjustments
#/**********************************************************/
start_time <- Sys.time()

trees_merge <- fread("GYPSY data/intermediate/i_trees_merge.csv")
plot_measurement <- fread("GYPSY data/intermediate/i_plot_measurement.csv")
trees <- left_join(trees_merge,plot_measurement,by=c("company","company_plot_number","measurement_number"))
plot <- fread("GYPSY data/intermediate/i_plot.csv")[,.(company,company_plot_number,natural_subregion)]
plot <- left_join(plot,
                  fread("GYPSY data/lookup/natsub.csv"),
                  by="natural_subregion")
trees <- left_join(trees,plot,by=c("company","company_plot_number"))
rm(trees_merge,plot_measurement)

trees$species_og <- toupper(trees$species)
trees$dbh_og <- trees$dbh
trees$height_og <- trees$height
trees$unique <- paste(trees$company,trees$company_plot_number,trees$tree_number,trees$measurement_number,sep="_")

trees[height<1.3 & !is.na(dbh),"height"] <- NA # Predict new height where DBH is measured but height is less than 1.3

# Fix suspected 0 filling
trees[dbh==0,"dbh"] <- NA
trees[height==0,"height"] <- NA

trees <- left_join(trees[,!"species"],
                   fread("GYPSY data/lookup/species.csv"),
                   by="species_og") # Joins species and species group from a lookup table

trees[is.na(species),"species"] <- trees[is.na(trees$species),"species_og"] # Use original species if no replacement
trees$dbh_stat <- "XX"
trees$ht_stat <- "XX"
trees[!is.na(dbh),"dbh_stat"] <- "MM"
trees[height<1.3 & !is.na(height),"dbh_stat"] <- "NA"
trees[!is.na(height), "ht_stat"] <- "MM"
trees <- trees[!(species%in%c("NO","MS")),]
#trees <- trees[!(tree_type=="B" & is.na(height)),] # Excludes age trees without heights
trees <- trees[!(!is.na(height) & height<0.3),]

#trees <- trees[!(species%in%c("AW","BW","PB") & (!is.na(height) & height<1.3)),] # Excludes deciduous under 1.3m

trees <- trees[!(condition_code1%in%c(1,2,13,14,15)),]
trees <- trees[!(species%in%c("DC","DD","DU")),] # could be combined with MS and NO
trees <- trees[tree_type!="",]
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
temp <- left_join(trees,
                  fread("GYPSY data/lookup/ht_to_dbh_2016r1.csv"),
                  by=c("species","natsub"))
temp[,dbh_p:=b1*(height-1.3)^b2*exp(-b3*(height-1.3))]
trees[height>1.3, "dbh_p"] <- temp[height>1.3, "dbh_p"]
rm(temp);gc()

trees[!is.na(dbh) & !is.na(dbh_p) & condition_code1!=3 & !(condition_code1%in%c(5,10) & severity1==3), 
      rate0:=dbh/dbh_p]

rate <- trees[!is.na(rate0)] %>%
        group_by(company,company_plot_number,measurement_number,species) %>%
        summarize(rate=mean(rate0))

## Append revised dbhs to dataset and add category for tracking type of DBH
temp <- left_join(trees,rate,by=c("company","company_plot_number","measurement_number","species"))[!is.na(dbh_p) & dbh_stat=="XX" & height>1.3]

temp[!is.na(rate), dbh_PA:=dbh_p * rate]
temp[is.na(rate), dbh_PU:=dbh_p]

trees <- left_join(trees[,!c("dbh_p","rate0")], temp[,c("unique","dbh_PA","dbh_PU")],by="unique")
rm(temp,rate);gc()

## Try to fill in missing DBHs based on one of: previous DBH, subsequent DBH, or previous + growth rate (if available)
trees$tree_id <- paste(trees$company,trees$company_plot_number,trees$tree_number,sep="_")

trees[,c("dbh_IP","dbh_IS","dbh_IX"):=as.numeric(NA)]

loop <- trees %>%
        group_by(tree_id) %>%
        filter(length(unique)>1 &
               any(dbh_stat == "XX", na.rm=T) &
               (any(dbh_stat == "MM", na.rm=T) |
                any(dbh_stat == "NA", na.rm=T)))
loop <- unique(loop$tree_id)
temp <- trees[trees$tree_id%in%loop,c("tree_id","unique","measurement_year","dbh_stat","dbh")] # Smaller subset to improve processing time
temp[,c("lm_yr","nm_yr","l_dbh","n_dbh"):=as.numeric(NA)]

dbh_check <- function(df,i){
  dat <- df[tree_id==i,]
  for(j in dat[dbh_stat=="XX",unique]){
    k <- dat$unique==j
    m_yr <- dat[k,measurement_year]
    p_meas <- dat[measurement_year<m_yr,]
    f_meas <- dat[measurement_year>m_yr,]
    if("NA"%in%f_meas[,dbh_stat]){
      dat$dbh_stat[k] <- "NA"
    }
    if("MM"%in%p_meas[,dbh_stat]){
      dat$lm_yr[k] <- max(p_meas[dbh_stat=="MM",measurement_year])
      dat$l_dbh[k] <- p_meas[measurement_year==max(measurement_year[dbh_stat=="MM"]),dbh]
    } 
    if("MM"%in%f_meas[,dbh_stat]){
      nm_yr <- min(f_meas[dbh_stat=="MM",measurement_year])
      dat$nm_yr[k] <- min(f_meas[dbh_stat=="MM",measurement_year])
      dat$n_dbh[k] <- f_meas[measurement_year==min(measurement_year[dbh_stat=="MM"]),dbh]
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
temp <- temp %>%
  mutate(
    dbh_IX = if_else(l_dbh > n_dbh, NA, round((n_dbh-l_dbh)/(nm_yr-lm_yr)*(measurement_year-lm_yr)+l_dbh,1))
  )


trees <- rows_update(trees,temp[,.(unique,dbh_IP,dbh_IS,dbh_IX,dbh_stat)],by="unique")

## A number of observations with missing height and diameter, no previous or subsequent measurement from which to populate, use average observed for that plot/mmt/species/tree_type
temp <- trees[trees$dbh_stat=="MM",] 
temp <- temp %>%
  group_by(company,company_plot_number,measurement_number,tree_type,spp_grp) %>%
  summarize(dbh_IA = mean(dbh))
trees <- left_join(trees,temp,by=c("company","company_plot_number","measurement_number","tree_type","spp_grp"))
rm("temp");gc()

## Pick which method hierarchy to use to assign missing DBH values
i <- trees$dbh_stat=="XX" & trees$height<1.3 & !is.na(trees$height)
trees$dbh[i] <- 0
trees$dbh_stat[i] <- "NA"

i <- trees$dbh_stat=="XX" & trees$height==1.3 & !is.na(trees$height) & !is.na(trees$dbh_IP)
trees$dbh[i] <- trees$dbh_IP[i]
trees$dbh_stat[i] <- "IP"

i <- trees$dbh_stat=="XX" & trees$height==1.3 & !is.na(trees$height) & trees$species%in%c("AW","BW","PB","LT","PJ","PL")
trees$dbh[i] <- 0.2
trees$dbh_stat[i] <- "IF"

i <- trees$dbh_stat=="XX" & trees$height==1.3 & !is.na(trees$height) & !trees$species%in%c("AW","BW","PB","LT","PJ","PL")
trees$dbh[i] <- 0.3
trees$dbh_stat[i] <- "IF"

i <- trees$dbh_stat=="XX" & trees$height>1.3 & !is.na(trees$height) & !is.na(trees$dbh_IX)
trees$dbh[i] <- trees$dbh_IX[i]
trees$dbh_stat[i] <- "IX"

i <- trees$dbh_stat=="XX" & trees$height>1.3 & !is.na(trees$height) & !is.na(trees$dbh_PA) & ((!is.na(trees$dbh_IP) & trees$dbh_PA>=trees$dbh_IP) | is.na(trees$dbh_PA)) & ((is.na(trees$dbh_IS & trees$dbh_PA<=trees$dbh_IS) | is.na(trees$dbh_IS)))
trees$dbh[i] <- trees$dbh_PA[i]
trees$dbh_stat[i] <- "PA"

i <- trees$dbh_stat=="XX" & trees$height>1.3 & !is.na(trees$height) & !is.na(trees$dbh_IP)
trees$dbh[i] <- trees$dbh_IP[i]
trees$dbh_stat[i] <- "IP"

i <- trees$dbh_stat=="XX" & trees$height>1.3 & !is.na(trees$height) & !is.na(trees$dbh_IS)
trees$dbh[i] <- trees$dbh_IS[i]
trees$dbh_stat[i] <- "IS"

i <- trees$dbh_stat=="XX" & trees$height>1.3 & !is.na(trees$height) & !is.na(trees$dbh_PU)
trees$dbh[i] <- trees$dbh_PU[i]
trees$dbh_stat[i] <- "PU"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","ES1") & !is.na(trees$dbh_IP)
trees$dbh[i] <- trees$dbh_IP[i]
trees$dbh_stat[i] <- "IP"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","ES1") & is.na(trees$dbh_IP)
trees$dbh[i] <- 0
trees$dbh_stat[i] <- "NA"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("T","S1","ET") & !is.na(trees$dbh_IX)
trees$dbh[i] <- trees$dbh_IX[i]
trees$dbh_stat[i] <- "IX"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("T","S1","ET") & !is.na(trees$dbh_IP)
trees$dbh[i] <- trees$dbh_IP[i]
trees$dbh_stat[i] <- "IP"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("T","S1","ET") & !is.na(trees$dbh_IS)
trees$dbh[i] <- trees$dbh_IS[i]
trees$dbh_stat[i] <- "IS"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("T","S1","ET") & !is.na(trees$dbh_IA)
trees$dbh[i] <- trees$dbh_IA[i]
trees$dbh_stat[i] <- "IA"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type%in%c("S1","ET")
trees$dbh[i] <- trees$sapling_tagging_limit_dbh[i]
trees$dbh_stat[i] <- "IM"

i <- trees$dbh_stat=="XX" & trees$ht_stat=="XX" & trees$tree_type=="T"
trees$dbh[i] <- trees$tree_tagging_limit[i]
trees$dbh_stat[i] <- "IM"

trees$dbh <- round(trees$dbh,1)

i <- trees$dbh_stat=="PU" & trees$dbh==0 & trees$height>1.3
trees$dbh[i] <- 0.1
trees$dbh_stat[i] <- "IF"

trees_test <- trees %>%
  select(
    company, company_plot_number, tree_number, measurement_year, unique, species, dbh, height, dbh_stat, dbh_PU, dbh_IP, dbh_IS, dbh_IX, dbh_IA
  )

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
temp <- left_join(trees,fread("GYPSY data/lookup/dbh_to_ht_2016r4.csv"),by=c("natsub","species"))

# Predict heights for all trees including ones that already have a height
temp[,ht_PU:=ifelse(model==1,1.3+b1*(1-exp(-1*b2*dbh))^b3,1.3+(b1/(1+exp(b2+b3*log(dbh+0.01)))))]
trees$ht_PU <- temp$ht_PU

## Calculate ratio of means to create adjustment factor for height
# Don't want trees with missing actual heights to be used in averaging
# Also don't want trees with unusual (i.e. data issues) height-DBH relationships to affect ratios
temp <- trees[!is.na(height) & !is.na(dbh)]
temp[,hdr:= height/dbh]
temp <- temp[!(height>5 & hdr>3) & 
             !(height>5 & hdr<0.3) & 
             !((height<=5 & height>=2) & hdr>5) &
             !((height<=5 & height>=2) & hdr<0.3) &
             !(height<2 & hdr>15) &
             !(height<2 & hdr<0.3) &
             !(condition_code1%in%c(1:3,13:15)) & # Dead, broken or dead top, etc.
             !(condition_code1==10 & severity1==3),] # Severe lean


# Attempt 1 - ratio by species, plot and measurement
height_ratio <- temp %>% 
  group_by(company,company_plot_number,measurement_number,species) %>%
  summarize(freq=sum(!is.na(height)&!is.na(ht_PU)),
            height_r=mean(height,na.rm=T),
            ht_PU=mean(ht_PU,na.rm=T)) %>%
  mutate(ratio=height_r/ht_PU,.keep="unused")

trees <- left_join(trees,height_ratio,by=c("company","company_plot_number","measurement_number","species"))
trees[freq>3 & ht_stat=="XX", ht_PA1:= ht_PU*ratio]
trees <- trees[,!c("ratio","freq")]

# Attempt 2 - ratio by species and plot only (to fill in missing after attempt #1)
height_ratio_2 <- temp %>% 
  group_by(company,company_plot_number,species) %>%
  summarize(freq=sum(!is.na(height)&!is.na(ht_PU)),
            height_r=mean(height,na.rm=T),
            ht_PU=mean(ht_PU,na.rm=T)) %>%
  mutate(
    ratio=height_r/ht_PU
  )
trees <- left_join(trees,height_ratio_2,by=c("company","company_plot_number","species"))
#rm(height_ratio)
i <- trees$freq>3 & !is.na(trees$freq) & trees$ht_stat=="XX"
trees$ht_PA2[i] <- trees$ht_PU.x[i]*trees$ratio[i]
trees <- trees[,!c("ratio","freq")]

# Height from previous measurement
trees[,c("ht_IP","ht_IS","ht_IX"):=as.numeric(NA)]

loop <- trees %>%
        group_by(tree_id) %>%
        filter(any(ht_stat == "XX", na.rm=T) &
               any(ht_stat == "MM", na.rm=T))
loop <- unique(loop$tree_id)
temp <- trees[trees$tree_id%in%loop,c("tree_id","unique","measurement_year","ht_stat","height","ht_IP","ht_IS","ht_IX")] # Smaller subset to improve processing time
temp[,c("lm_yr","nm_yr","l_ht","n_ht"):=as.numeric(NA)]

ht_check <- function(df,i){
  dat <- df[tree_id==i,]
  for(j in dat[ht_stat=="XX",unique]){
    k <- dat$unique==j
    m_yr <- dat[k,measurement_year]
    p_meas <- dat[measurement_year<m_yr,]
    f_meas <- dat[measurement_year>m_yr,]
    if("MM"%in%p_meas[,ht_stat]){
      dat$lm_yr[k] <- max(p_meas[ht_stat=="MM",measurement_year])
      dat$l_ht[k] <- p_meas[measurement_year==max(measurement_year[ht_stat=="MM"]),height] # This formatting is a workaround for what seems to be a bug
    } 
    if("MM"%in%f_meas[,ht_stat]){
      dat$nm_yr[k] <- min(f_meas[ht_stat=="MM",measurement_year])
      dat$n_ht[k] <- f_meas[measurement_year==min(measurement_year[ht_stat=="MM"]),height]  
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
## need to add condition if height n_ht > l_ht
temp[!is.na(l_ht)&!is.na(n_ht) & nm_yr > lm_yr,
     ht_IX:=round((n_ht-l_ht)/(nm_yr-lm_yr)*(measurement_year-lm_yr)+l_ht,1)]

########################
# temp[!is.na(l_ht)&is.na(n_ht),
#      ht_IP:=l_ht]
# temp[is.na(l_ht)&!is.na(n_ht),
#      ht_IS:=n_ht]

temp[!is.na(l_ht),
     ht_IP:=l_ht]
temp[!is.na(n_ht),
     ht_IS:=n_ht]
temp[l_ht > n_ht,
     ht_IX:= NA]

temp <- temp %>%
  mutate(
    ht_IX = if_else(l_ht > n_ht, NA, round((n_ht-l_ht)/(nm_yr-lm_yr)*(measurement_year-lm_yr)+l_ht,1))
  )

trees_1 <- trees %>%
  select(-c(ht_IP, ht_IS, ht_IX,))

trees_2 <- trees_1 %>%
  left_join(temp %>% select(ht_IP, ht_IS, ht_IX, ht_stat, unique), by = "unique") %>%
  mutate(ht_PU = coalesce(ht_PU.x, ht_PU, ht_PU.y)) %>%
  mutate(ht_stat = coalesce(ht_stat.x, ht_stat.y))

trees_origin <- trees

trees <- trees_2
################
# Average height by species / measurement / tree type
i <- trees$ht_stat=="MM"
htavg <- trees[i,] %>% 
  group_by(company,company_plot_number,measurement_number,tree_type,spp_grp) %>%
  summarize(ht=mean(height,na.rm=T))

trees <- left_join(trees,htavg,by=c("company","company_plot_number","measurement_number","tree_type","spp_grp"))
#rm(list=c("htavg","temp"));gc()


i<- trees$ht_stat=="XX" & !is.na(trees$ht)
trees$ht_IA[i] <- trees$ht[i]
trees <- trees[,!"ht"]

# Assignment rules for which height to select
i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_IX) & trees$ht_IX>=1.3
trees$height[i] <- trees$ht_IX[i]
trees$ht_stat[i] <- "IX"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_PA1) & trees$ht_PA1>=1.3
trees$height[i] <- trees$ht_PA1[i]
trees$ht_stat[i] <- "PA"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_PA2) & trees$ht_PA2>=1.3
trees$height[i] <- trees$ht_PA2[i]
trees$ht_stat[i] <- "Pa"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_IP) & trees$ht_IP>=1.3
trees$height[i] <- trees$ht_IP[i]
trees$ht_stat[i] <- "IP"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_PU) & trees$ht_PU>=1.3
trees$height[i] <- trees$ht_PU[i]
trees$ht_stat[i] <- "PU"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh>0 & !is.na(trees$ht_IS) & trees$ht_IS>=1.3
trees$height[i] <- trees$ht_IS[i]
trees$ht_stat[i] <- "IS"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & !is.na(trees$ht_IX) & trees$ht_IX<1.3
trees$height[i] <- trees$ht_IX[i]
trees$ht_stat[i] <- "IX"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & !is.na(trees$ht_IP) & trees$ht_IP<1.3
trees$height[i] <- trees$ht_IP[i]
trees$ht_stat[i] <- "IP"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & !is.na(trees$ht_IA) & trees$ht_IA<1.3
trees$height[i] <- trees$ht_IA[i]
trees$ht_stat[i] <- "IA"
  
i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & trees$tree_type%in%c(paste0("R",1:10),"ES1") & !trees$species%in%c("AW","BW","PB") & !is.na(trees$regen_tagging_limit_conifer)  
trees$height[i] <- trees$regen_tagging_limit_conifer[i]
trees$ht_stat[i] <- "IM"
  
i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & trees$tree_type%in%c(paste0("R",1:10),"ES1") & trees$species%in%c("AW","BW","PB") & !is.na(trees$regen_tagging_limit_decid)    
trees$height[i] <- trees$regen_tagging_limit_decid[i]
trees$ht_stat[i] <- "IM"
  
i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & trees$tree_type%in%c(paste0("R",1:10),"ES1") & !trees$species%in%c("AW","BW","PB") & is.na(trees$regen_tagging_limit_conifer)  
trees$height[i] <- 0.3
trees$ht_stat[i] <- "IM"

i <- trees$ht_stat=="XX" & !is.na(trees$dbh) & trees$dbh==0 & trees$tree_type%in%c(paste0("R",1:10),"ES1") & trees$species%in%c("AW","BW","PB") & is.na(trees$regen_tagging_limit_decid)     
trees$height[i] <- 0.3
trees$ht_stat[i] <- "IM"

trees$height <- round(trees$height,2)

i <- !is.na(trees$dbh) & trees$dbh==0 & !is.na(trees$height) & trees$height>=1.3
trees$dbh[i] <- 0.1

trees_3 <- trees[!is.na(trees$height) & trees$height>0.3,]

trees_test <- trees_3 %>%
  select(
    company, company_plot_number, tree_number, measurement_year, unique, species, dbh, height, dbh_stat, ht_PU, ht_PA1, ht_PA2, ht_IP, ht_IS, ht_IX, ht_IA
  )

trees <- trees[,!c("ht_PU","ht_PA1","ht_PA2","ht_IP","ht_IS","ht_IA","ht_IX","dbh_og","height_og")]
  
fwrite(trees,"GYPSY data/intermediate/i_trees_measurement_filled.csv")
end_time <- Sys.time()-start_time
minutes <- round(end_time)
seconds <- round((end_time-minutes)*60)
message(paste(minutes,"minutes and",seconds,"seconds for DBH and HT fill"))
