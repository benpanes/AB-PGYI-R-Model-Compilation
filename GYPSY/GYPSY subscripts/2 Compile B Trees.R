#/**********************************************************/
# Compile individual tree basal area, volumes, sph
#/**********************************************************/
start_time <- Sys.time()

trees <- fread("GYPSY data/intermediate/i_trees_measurement_filled.csv")

trees <- trees[!(trees$tree_type%in%c("T","ET") & (trees$tree_plot_area==0 | is.na(trees$tree_plot_area))) &
               !(trees$tree_type%in%c("S","ES") & (trees$sapling_plot_area==0 | is.na(trees$sapling_plot_area))) &
               !(trees$tree_type%in%c(paste0("R",1:10)) & (trees$regen_plot_area==0 | is.na(trees$regen_plot_area))),]

trees$ba <- ((trees$dbh/200)^2)*pi

i <- trees$tree_type%in%c("T","ET")
trees$sph[i] <- 10000/trees$tree_plot_area[i]

i <- trees$tree_type%in%c("S","ES")
trees$sph[i] <- 10000/trees$sapling_plot_area[i]

i <- trees$tree_type%in%c(paste0("R",1:10))
trees$sph[i] <- 10000/(trees$regen_plot_area[i]*trees$number_regen_plots[i])

i <- trees$tree_type=="B"
trees$sph[i] <- 0
trees$ba[i] <- 0

trees.tap <- left_join(trees,fread("GYPSY data/lookup/natsub.csv"),by="natural_subregion")

trees.tap$natsub[trees$species%in%c("BW","FD","LT","PJ","SE")] <- 0

trees.tap <- left_join(trees.tap,fread(paste0(work.dir,"lookup/taper.csv")),by=c("species","natsub"));gc()

## Size parameters
# All
p <- 0.2250
stumpH <- 0.3
topdib <- 0.0125
  
#Merchantable
stumpD <- 13
min_mlen <- 3.66
topdibM <- 7

#### Merchantable volume
## Merchantable height / height ratio (h/H)
temp <- trees[!is.na(trees$height) & trees$height>=1.3 & !is.na(trees$dbh) & topdibM<trees$dbh,
              c("unique","height","dbh","b1","b2","b3","b4","b5","a0","a1","a2","k7","k8")]

# Diameter inside bark formula
dib <- function(x,z){
  (temp$a0[z]*temp$dbh[z]^temp$a1[z])*(temp$a2[z]^temp$dbh[z])*((1-sqrt(x/temp$height[z]))/(1-sqrt(p)))^
  (temp$b1[z]*(x/temp$height[z])^2+temp$b2[z]*log(x/temp$height[z]+0.001)+temp$b3[z]*sqrt(x/temp$height[z])+temp$b4[z]*exp(x/temp$height[z])+temp$b5[z]*temp$dbh[z]/temp$height[z])
}

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
    message(paste0("Timeout. Skipping row'",i,"'"))
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
  svol <- function(df,i){(2*df$sl[i]/6)*0.00007854*(d0^2+4*d1^2+d2^2)}
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

#### Calculate merchantable volume
temp <- trees.tap[dbh>topdibM,]

## Merchantable height ratio
len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl,c("hr","acc","topdibM","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl,1:len,function(i) hr(temp,i,topdibM)))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(temp,1:nrow(temp),stumpH)
temp$dobs <- temp$k7+temp$k8*temp$dibs

## Merchantable volume
temp <- temp[temp$dobs>=stumpD,]
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH
temp$sl <- temp$ml/20
temp <- temp[ml>min_mlen,]

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl,c("temp","p","stumpH","dib","fvol"))
temp$mvol <- unlist(parLapply(cl,1:len,function(x) fvol(temp,x)))
stopCluster(cl)

trees.vol <- left_join(trees,temp[,.(company,company_plot_number,tree_number,measurement_number,mvol)],by=c("company","company_plot_number","tree_number","measurement_number"))
trees.vol[is.na(mvol),"mvol"] <- 0
rm(temp);gc()

#### Calculate total volume
temp <- trees.tap[dbh>topdib,]

## Height ratio
len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl,c("hr","acc","topdib","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl,1:len,function(i) hr(temp,i,topdib)))
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
temp$tvol <- tvol+pi*(topdib/200)^2*(temp$height-temp$mh)/3+pi*(temp$dibs/200)^2*stumpH

trees.vol <- left_join(trees.vol,temp[,.(company,company_plot_number,tree_number,measurement_number,tvol)],by=c("company","company_plot_number","tree_number","measurement_number"))
trees.vol[is.na(tvol),"tvol"] <- 0
rm(temp);gc()

fwrite(trees.vol,"GYPSY data/intermediate/i_tree_ages.csv") # wrong output

minutes <- round(Sys.time()-start_time)
seconds <- round((Sys.time()-start_time-minutes)*60)
message(paste(minutes,"minutes and",seconds,"seconds for compile trees"))