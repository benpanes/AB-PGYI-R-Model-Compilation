#/**********************************************************/
#Append measurement info to trees data and make adjustments
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

natsub <- fread("GYPSY data/lookup/natsub.csv")
trees <- left_join(trees,natsub,by="natural_subregion")
rm(natsub)

trees$natsub[trees$species%in%c("BW","FD","LT","PJ","SE")] <- 0

i <- trees$species=="AW" & trees$dbh==7.1 & trees$height==12.8 & trees$natsub==10 # This combination hangs up the height equation for some reason
trees$dbh[i] <- 7.2

coef <- fread("GYPSY data/lookup/taper.csv")
trees <- left_join(trees,coef,by=c("species","natsub"))
rm(coef)

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

# Volume formula
fvol <- function(x){
  svol <- function(x){(2*temp$sl[x]/6)*0.00007854*(d0^2+4*d1^2+d2^2)}
  v <- c()
  h0 <- stumpH
  d0 <- temp$dibs[x]
  h1 <- h0 + temp$sl[x]
  d1 <- dib(h1,x)
  h2 <- h1 + temp$sl[x]
  d2 <- dib(h2,x)
  v[1] <- svol(x)
  for(n in 2:10){
    h0 <- h2
    h1 <- h0 + temp$sl[x]
    d1 <- dib(h1,x)
    h2 <- h1 + temp$sl[x]
    d2 <- dib(h2,x)
    v[n] <- svol(x)
  }
  sum(v)
}

# Height ratio formula
acc <- 0.00000001 # Desired difference between penultimate and last h/H
r1 <- 1
r0 <- 0.8

hr <- function(i,d,diff=1){
  hr1 <- function(i){
    temp$b1[i]*r0^2+temp$b2[i]*log(r0+0.001)+temp$b3[i]*sqrt(r0)+temp$b4[i]*exp(r0)+temp$b5[i]*(temp$dbh[i]/temp$height[i])
  }
  hr2 <- function(i,td){
    (1-(d/(temp$a0[i]*temp$dbh[i]^temp$a1[i]*temp$a2[i]^temp$dbh[i]))^(1/hr1(i))*(1-sqrt(p)))^2
  }
  tryCatch({ # There was a specific combination where the function was getting stuck (noted in Ken's code), so I added this to account for any calculations taking too long
    withTimeout({
      while(diff>acc){
        r1 <- hr2(i)
        r0 <- (r1+r0)/2
        diff <- abs(r1-r0)
      }
    }, timeout=1)
  }, TimeoutException=function(ex){
    message(paste0("Timeout. Skipping row'",i,"'"))
    r0 <- NULL
  })
r0
}

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl, c("hr","acc","topdibM","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl, 1:len, function(i) hr(i,topdibM)))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(stumpH,1:nrow(temp))
temp$dobs <- temp$k7+temp$k8*temp$dibs

## Volume
temp <- temp[temp$dobs>=stumpD,]
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH
temp$sl <- temp$ml/20
temp <- temp[temp$ml>min_mlen,]

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl, c("temp","p","stumpH","dib","fvol"))
temp$mvol <- unlist(parLapply(cl, 1:len, function(x) fvol(x)))
stopCluster(cl)

volume <- temp[,c("unique","mvol")]
trees <- left_join(trees,volume,by="unique")

#### Gross total volume
## Merchantable height / height ratio (h/H)
temp <- trees[!is.na(trees$height) & trees$height>=1.3 & !is.na(trees$dbh) & topdib<trees$dbh,
              c("unique","height","dbh","b1","b2","b3","b4","b5","a0","a1","a2","k7","k8")]

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils))
clusterExport(cl, c("hr","acc","topdib","temp","r0","r1","p"))
temp$r0 <- unlist(parLapply(cl, 1:len, function(i) hr(i,topdib)))
stopCluster(cl)

## Stump diameters
temp$dibs <- dib(stumpH,1:nrow(temp))

temp$dobs <- temp$k7+temp$k8*temp$dibs

## Volume
temp$mh <- temp$height*temp$r0
temp$ml <- temp$mh-stumpH[3]
temp$sl <- temp$ml/20

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterExport(cl, c("temp","p","stumpH","dib","fvol"))
temp$vol <- unlist(parLapply(cl, 1:len, function(x) fvol(x)))
stopCluster(cl)

temp$vol <- temp$vol+pi*(topdib/200)^2*(temp$height-temp$mh)/3+pi*(temp$dibs/200)^2*stumpH
vol_0000 <- temp[,c("unique","vol")]
trees <- left_join(trees,vol_0000,by="unique")

fwrite(trees,"GYPSY data/intermediate/i_tree_ages.csv")

minutes <- round(Sys.time()-start_time)
seconds <- round((Sys.time()-start_time-minutes)*60)
message(paste(minutes,"minutes and",seconds,"seconds for compile trees"))