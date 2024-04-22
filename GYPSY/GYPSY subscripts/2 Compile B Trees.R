#/**********************************************************/
# Compile individual tree basal area, volumes, sph
#/**********************************************************/
start_time <- Sys.time()

trees <- fread("GYPSY data/intermediate/i_trees_measurement_filled.csv")

trees <- trees[!(trees$tree_type%in%c("T","ET") & (trees$tree_plot_area==0 | is.na(trees$tree_plot_area))) &
               !(trees$tree_type%in%c("S","ES") & (trees$sapling_plot_area==0 | is.na(trees$sapling_plot_area))) &
               !(trees$tree_type%in%c(paste0("R",1:10)) & (trees$regen_plot_area==0 | is.na(trees$regen_plot_area))),]

trees$index_1 <- paste(trees$company, trees$company_plot_number, trees$measurement_number, trees$tree_number, sep = "_")
 
trees$ba <- ((trees$dbh/200)^2)*3.14159

i <- trees$tree_type%in%c("T","ET")
trees$sph[i] <- 10000/trees$tree_plot_area[i]

i <- trees$tree_type%in%c("S1","ES1")
trees$sph[i] <- 10000/trees$sapling_plot_area[i]

i <- trees$tree_type%in%c(paste0("R",1:10))
trees$sph[i] <- 10000/(trees$regen_plot_area[i]*trees$number_regen_plots[i])

i <- trees$tree_type=="B"
trees$sph[i] <- 0
trees$ba[i] <- 0

# assign new species names 
trees <- trees %>%
  mutate(species = case_when(
    species_og %in% c("AW", "AX") ~ "AW",
    species_og == "BW" ~ "BW",
    species_og == "PB" ~ "PB",
    species_og %in% c("LA", "LT", "LW") ~ "LT",
    species_og %in% c("P", "PX", "PF", "PL", "PW") ~ "PL",
    species_og == "PJ" ~ "PJ",
    species_og %in% c("FA", "FB") ~ "FB",
    species_og == "FD" ~ "FD",
    species_og %in% c("SE", "SW", "SX") ~ "SW",
    species_og == "SB" ~ "SB"
  ))

trees.tap <- left_join(trees,fread("GYPSY data/lookup/natsub.csv"),by="natural_subregion")

# changed species to species_og for testing purpose
trees.tap$natsub[trees$species%in%c("BW","FD","LT","PJ","SE")] <- 0

trees.tap <- left_join(trees.tap,fread("GYPSY data/lookup/taper.csv"),by=c("species","natsub"))
trees.tap <- trees.tap[,.(company,company_plot_number,tree_number,measurement_number,unique,height,dbh,b1,b2,b3,b4,b5,a0,a1,a2,k7,k8)];gc()

## Size parameters
p <- 0.2250
stumpH <- 0.3
topdib <- list("T"=0.0125,
               "M1"=7,
               "M2"=10)
stumpD <- list(13,15)
min_mlen <- 3.66

## Volume formulae
# Height ratio formula
acc <- 0.00000001
r1 <- 1
r0 <- 0.8

hr <- function(df,i,d,diff=1){
  hr1 <- function(df,i){
    df$b1[i]*r0^2+df$b2[i]*log(r0+0.001)+df$b3[i]*sqrt(r0)+df$b4[i]*exp(r0)+df$b5[i]*(df$dbh[i]/df$height[i])
  }
  hr2 <- function(df,i,d){
    (1-((d/(df$a0[i]*df$dbh[i]^df$a1[i]*df$a2[i]^df$dbh[i]))^(1/hr1(df,i)))*(1-sqrt(p)))^2
  }
  tryCatch({
    withTimeout({
      while(diff>acc){
        r1 <- hr2(df,i,d)
        r0 <- (r1+r0)/2
        diff <- abs(r1-r0)
      }
    },timeout=1)
  }, TimeoutException=function(ex){
    warning(paste0("Timeout. Skipping row'",i,"'"))
  })
  r0
}

# Diameter inside bark formula
dib <- function(df,i,x){
  (df$a0[i]*df$dbh[i]^df$a1[i])*(df$a2[i]^df$dbh[i])*((1-sqrt(x/df$height[i]))/(1-sqrt(p)))^
    (df$b1[i]*(x/df$height[i])^2+df$b2[i]*log(x/df$height[i]+0.001)+df$b3[i]*sqrt(x/df$height[i])+df$b4[i]*exp(x/df$height[i])+df$b5[i]*df$dbh[i]/df$height[i])  
}

# Volume formula
fvol <- function(df,i){
  svol <- function(df,i){(2*df$sl[i]/6)*((1/200)^2*pi)*(d0^2+4*d1^2+d2^2)}
  v <- c()
  h0 <- stumpH
  d0 <- df$dibs[i]
  h1 <- h0 + df$sl[i]
  d1 <- dib(df,i,h1)
  h2 <- h1 + df$sl[i]
  d2 <- dib(df,i,h2)
  v[1] <- svol(df,i)
  for(n in 2:10){
    h0 <- h2
    d0 <- d2
    h1 <- h0 + df$sl[i]
    d1 <- dib(df,i,h1)
    h2 <- h1 + df$sl[i]
    d2 <- dib(df,i,h2)
    v[n] <- svol(df,i)
  }
  sum(v)
}

#### Calculate merchantable volume 13/7
temp <- trees.tap[dbh>topdib[["M1"]],]

## Merchantable height ratio
len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl,c("hr","acc","topdib","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl,1:len,function(i) hr(temp,i,topdib[["M1"]])))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(temp,1:nrow(temp),stumpH)
temp$dobs <- temp$k7+temp$k8*temp$dibs

## Merchantable volume
temp <- temp[temp$dobs>=stumpD[[1]],]
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH
temp$sl <- temp$ml/20
temp <- temp[ml>min_mlen,]

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl,c("temp","p","stumpH","dib","fvol"))
temp$vol_1307 <- unlist(parLapply(cl,1:len,function(x) fvol(temp,x)))
stopCluster(cl)

trees.vol <- left_join(trees,temp[,.(company,company_plot_number,tree_number,measurement_number,vol_1307)],by=c("company","company_plot_number","tree_number","measurement_number"))
trees.vol[is.na(vol_1307),"vol_1307"] <- 0
rm(temp);gc()

#### Calculate merchantable volume 15/10
temp <- trees.tap[dbh>topdib[["M2"]]]

## Merchantable height ratio
len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl,c("hr","acc","topdib","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl,1:len,function(i) hr(temp,i,topdib[["M2"]])))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(temp,1:nrow(temp),stumpH)
temp$dobs <- temp$k7+temp$k8*temp$dibs

## Merchantable volume
temp <- temp[temp$dobs>=stumpD[[2]],]
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH
temp$sl <- temp$ml/20
temp <- temp[ml>min_mlen,]

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl,c("temp","p","stumpH","dib","fvol"))
temp$vol_1510 <- unlist(parLapply(cl,1:len,function(x) fvol(temp,x)))
stopCluster(cl)

trees.vol <- left_join(trees.vol,temp[,.(company,company_plot_number,tree_number,measurement_number,vol_1510)],by=c("company","company_plot_number","tree_number","measurement_number"))
trees.vol[is.na(vol_1510),"vol_1510"] <- 0
rm(temp);gc()

#### Calculate total volume
temp <- trees.tap[dbh>topdib[["T"]],]

## Height ratio
len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl,c("hr","acc","topdib","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl,1:len,function(i) hr(temp,i,topdib[["T"]])))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(temp,1:nrow(temp),stumpH)
temp$dobs <- temp$k7+temp$k8*temp$dibs

## Total volume
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH
temp$sl <- temp$ml/20

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl,c("temp","p","stumpH","dib","fvol"))
tvol <- unlist(parLapply(cl,1:len,function(x) fvol(temp,x)))
stopCluster(cl)
temp$vol_0000 <- tvol+pi*(topdib[["T"]]/200)^2*(temp$height-temp$mh)/3+pi*(temp$dibs/200)^2*stumpH

trees.vol <- left_join(trees.vol,temp[,.(company,company_plot_number,tree_number,measurement_number,vol_0000)],by=c("company","company_plot_number","tree_number","measurement_number"))
trees.vol[is.na(vol_0000),"vol_0000"] <- 0

##### Left off at the at the end of the volume calculation section (before biomass)

# biomass and carbon calculation

trees.bio <- trees.vol %>%
  mutate(
    biomass_tree = if_else(dbh < 15 | vol_1307 == 0,
                      case_when(
                        species == 'AW' ~ 0.26738 + 0.01917 * dbh * dbh * height,
                        species == 'BW' ~ 2.47035 + 0.02454 * dbh * dbh * height,
                        species == 'FA' ~ 7.03447 + 0.01477 * dbh * dbh * height,
                        species %in% c('FB', 'FD') ~ 7.94339 + 0.01465 * dbh * dbh * height,
                        species %in% c('LT', 'LA', 'LW') ~ 4.48372 + 0.01768 * dbh * dbh * height,
                        species == 'PB' ~ 10.74706 + 0.01350 * dbh * dbh * height,
                        species == 'PJ' ~ 2.91931 + 0.01678 * dbh * dbh * height,
                        species %in% c('PL', 'PW', 'PF', 'PX') ~ 8.18267 + 0.01597 * dbh * dbh * height,
                        species == 'SB' ~ 2.79552 + 0.01698 * dbh * dbh * height,
                        species %in% c('SW', 'SE', 'SX') ~ 6.03377 + 0.01500 * dbh * dbh * height,
                        TRUE ~ NA_real_
                      ),
                      case_when(
                        species == 'AW' ~ 499.508 * vol_1307 ^ 0.980765,
                        species == 'BW' ~ 703.360 * vol_1307 ^ 0.946751,
                        species == 'FA' ~ 434.694 * vol_1307 ^ 0.903315,
                        species %in% c('FB', 'FD') ~ 444.532 * vol_1307 ^ 0.873007,
                        species == 'PJ' ~ 477.288 * vol_1307 ^ 0.983019,
                        species %in% c('PL', 'PW', 'PF', 'PX') ~ 436.564 * vol_1307 ^ 0.962308,
                        species %in% c('LT', 'LA', 'LW') ~ 530.347 * vol_1307 ^ 0.922289,
                        species == 'SB' ~ 516.226 * vol_1307 ^ 1.001660,
                        species %in% c('SW', 'SE', 'SX') ~ 451.544 * vol_1307 ^ 0.958852,
                        TRUE ~ NA_real_
                      )
    ),
    biomass = biomass_tree * sph,
    carbon = biomass * 0.5,
    vol_1307 = vol_1307 * sph,
    vol_1510 = vol_1510 * sph,
    vol_0000 = vol_0000 * sph
  )

treelist <- trees.bio %>%
  mutate(
    vol_1510 = ifelse(tree_type == 'B' | height < 1.3, 0, vol_1510),
    vol_1307 = ifelse(tree_type == 'B' | height < 1.3, 0, vol_1307),
    vol_0000 = ifelse(tree_type == 'B' | height < 1.3, 0, vol_0000),
    biomass = ifelse(tree_type == 'B' | height < 1.3, 0, biomass),
    carbon = ifelse(tree_type == 'B' | height < 1.3, 0, carbon)
  ) %>%
  filter(
    tree_type != "B"
  ) %>%
  select(company, company_plot_number, measurement_number, measurement_year, tree_number, 
         tree_location, tree_origin, tree_type, species, distance, azimuth, dbh, height, 
         crown_class, condition_code1, severity1, condition_code2, severity2, condition_code3, 
         severity3, ht_stat, dbh_stat, ba, sph, vol_0000, vol_1307, vol_1510, biomass, carbon, spp_grp, condec) 

fwrite(treelist,"GYPSY data/intermediate/i_tree_list1.csv") 

tree_ages <- trees.bio %>%
  filter(!(total_age %in% c(0, NA) & stump_age %in% c(0, NA) & dbh_age %in% c(0, NA))) %>%
  filter(!(tree_origin %in% c(9, 10))) %>%
  filter(ht_stat %in% c('MM', 'IX', 'PU')) 

fwrite(tree_ages,"GYPSY data/intermediate/i_tree_age.csv")

minutes <- round(Sys.time()-start_time)
seconds <- round((Sys.time()-start_time-minutes)*60)
message(paste(minutes,"minutes and",seconds,"seconds for compile trees"))