treeage <- read.csv(""GYPSY data/intermediate/i_tree_age.csv")
plot_level2 <- read.csv("GYPSY data/intermediate/i_plot_level2.csv")
plot_mmt <- fread(""GYPSY data/intermediate/i_plot_measurement.csv")
tree_list1 <- fread(""GYPSY data/intermediate/tree_list1.csv")
  
################################################################################
#  STEP 1. fill in missing ages & remove duplicates - goal to have an age for both total and breast height, where available
################################################################################

treeage$stump_age <- as.numeric(treeage$stump_age)
treeage$total_age <- as.numeric(treeage$total_age)
treeage$dbh_age <- as.numeric(treeage$dbh_age)

site1 <- treeage %>% 
  filter(!(tree_origin %in% c(9, 10))) %>%
  mutate(
    totage = ifelse(!is.na(total_age), total_age, NA),
    bhage =  ifelse(!is.na(dbh_age), dbh_age, NA)) %>%              
    mutate(
    totage = ifelse(!is.na(bhage) & totage %in% c(0, NA),
                                case_when(
                                  spp_grp == 'AW' ~ dbh_age + 4,
                                  spp_grp == 'PL' ~ dbh_age + 8,
                                  spp_grp == 'SB' ~ dbh_age + 15,
                                  spp_grp == 'SW' ~ dbh_age + 12
                                ), totage),
    bhage = ifelse(!is.na(totage) & bhage %in% c(0, NA),
                               case_when(
                                 spp_grp == 'AW' & totage > 4 ~ totage - 4,
                                 spp_grp == 'PL' & totage > 8 ~ totage - 8,
                                 spp_grp == 'SB' & totage > 15 ~ totage - 15,
                                 spp_grp == 'SW' & totage > 12 ~ totage - 12
                               ), bhage)) 

site1 <- site1%>%
  mutate(
    totage = ifelse(!(stump_age %in% c(0, NA)) & totage %in% c(0, NA) & bhage %in% c(0, NA),
                         case_when(
                           spp_grp == 'AW' ~ stump_age + 0.5,
                           spp_grp == 'PL' ~ stump_age + 3,
                           spp_grp == 'SB' ~ stump_age + 5,
                           spp_grp == 'SW' ~ stump_age + 4
                         ), totage),
    bhage = ifelse(!(stump_age %in% c(0, NA)) & totage %in% c(0, NA) & bhage %in% c(0, NA),
                        case_when(
                          spp_grp == 'AW' & totage > 4 ~ totage - 4,
                          spp_grp == 'PL' & totage > 8 ~ totage - 8,
                          spp_grp == 'SB' & totage > 15 ~ totage - 15,
                          spp_grp == 'SW' & totage > 12 ~ totage - 12
                        ), bhage)
  )%>% 
  arrange(desc(measurement_number))

managed2 <- site1 %>% 
  group_by(unique) %>%
  filter(row_number() == 1)%>%
  filter(stand_origin != 'N') %>%
  arrange(company, company_plot_number)

natural2 <- site1 %>% 
  group_by(unique) %>%
  filter(row_number() == 1) %>%
  filter(stand_origin == 'N')%>%
  arrange(company, company_plot_number)

managed_natural <- bind_rows(managed2, natural2)

################################################################################
#  STEP 2. calculate site index from measured age trees and calculate average by species
################################################################################

# solve for si by totage and toph

managed_natural$height <- as.numeric(managed_natural$height)
site2 <- managed_natural %>%
  rename(
    topht = height
  ) %>%
  filter(
    !(spp_grp %in% c('AW', 'PL') &  totage < 7.5),
    !(spp_grp %in% c('SB', 'SW') &  totage < 9.5)
  ) %>%
  filter(!((condition_code1 %in% c(3, 8, 10) & severity1 %in% c(2, 3, 9)) |
             (condition_code2 %in% c(3, 8, 10) & severity2 %in% c(2, 3)) |
             (condition_code3 %in% c(3, 8, 10) & severity3 %in% c(2, 3)))) 

################################################################################
# calculate si_t, si_bh, y2bh
covar <- data.frame(spp_grp = c("AW","PL","SB","SW"), 
                    b1 = c(9.908888,12.84571,14.56236,12.14943),
                    b2 = c(-3.92451,-5.73936,-6.04705,-3.77051),
                    b3 = c(-0.32778,-0.91312,-1.53715,-0.28534),
                    b4 = c(0.134376,0.150668,0.240174,0.165483))

temp <- left_join(site2[!is.na(site2$spp_grp)&
                         !is.na(site2$topht)& 
                         !is.na(site2$totage),],
                  covar,by="spp_grp")

si_solve <- function(i,df){
  b1 = df$b1[i]
  b2 = df$b2[i]
  b3 = df$b3[i]
  b4 = df$b4[i]
  age_tot = df$totage[i]
  topht = df$topht[i]
  si0 = 10
  si1 = 100
  tryCatch({ # These functions add a max time and an error message if a row takes too long
    withTimeout({
      while(abs(si0-si1)>0.00000001){
        if(df$spp_grp[i]=="AW"){
          k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^2+b4*sqrt(50))
          k2 = si0^(b3*log(si0))
          k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
          Y2BH = exp((log(k3)/b2)^2)-1
          x10 = (1+exp(b1+b2*sqrt(log(0+age_tot+1.0))+b3*log(si0)^2+b4*sqrt(50)))
          x20 = (1+exp(b1+b2*sqrt(log(50+1.0))+b3*log(si0)^2+b4*sqrt(50)))
          si1 = topht*x10/x20
          si0 = (si0+si1)/2
          si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^2+b4*sqrt(50)))/
            (1+exp(b1+b2*sqrt(log(50+Y2BH+1))+b3*(log(si1))^2+b4*sqrt(50)))
          si_t = si0
        }
        if(df$spp_grp[i]=="SB"){
          k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^1+b4*sqrt(50))
          k2 = si0^b3
          k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
          Y2BH = exp((log(k3)/b2)^2)-1
          x10 = (1+exp(b1+b2*sqrt(log(0+age_tot+1))+b3*log(si0)^1+b4*sqrt(50)))
          x20 = (1+exp(b1+b2*sqrt(log(50+1))+b3*log(si0)^1+b4*sqrt(50)))
          si1 = topht*x10/x20
          si0 = (si0+si1)/2
          si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^1+b4*sqrt(50)))/
            (1+exp(b1+b2*sqrt(log(50+Y2BH+1))+b3*(log(si1))^1+b4*sqrt(50)))
          si_t = si0
        }
        if(df$spp_grp[i]=="SW"){
          k1 = exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50))
          k2 = si0^(b3*log(si0))
          k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
          Y2BH = sqrt(exp((log(k3)/b2)^2)-1)
          x10 = (1+exp(b1+b2*sqrt(log(0+age_tot^2+1))+b3*log(si0)^2+b4*sqrt(50)))
          x20 = (1+exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50)))
          si1 = topht*x10/x20
          si0 = (si0+si1)/2
          si_bh = si1*(1+exp(b1+b2*sqrt(log(50^2+1))+b3*(log(si1))^2+b4*sqrt(50)))/
            (1+exp(b1+b2*sqrt(log((50+Y2BH)^2+1))+b3*(log(si1))^2+b4*sqrt(50)))
          si_t = si0
        }
        if(df$spp_grp[i]=="PL"){
          k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^1+b4*sqrt(50))
          k2 = si0^b3
          k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
          Y2BH = exp((log(k3)/b2)^2)-1
          x10 = (1+exp(b1+b2*sqrt(log(0+age_tot+1))+b3*log(si0)^1+b4*sqrt(50)))
          x20 = (1+exp(b1+b2*sqrt(log(50+1))+b3*log(si0)^1+b4*sqrt(50)))
          si1 = topht*x10/x20
          si0 = (si0+si1)/2
          si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^1+b4*sqrt(50)))/
            (1+exp(b1+b2*sqrt(log((50+Y2BH)+1))+b3*(log(si1))^1+b4*sqrt(50)))
          si_t = si0
        }
      }
    }, timeout=1)
  }, TimeoutException=function(ex){
    message(paste0("Timeout. Skipping row'",i,"'")) # Error message and returned value on timeout
    Y2BH <- NULL
    si_t <- NULL
    si_bh <- NULL
  })
  return(list(Y2BH = Y2BH, si_t = si_t, si_bh = si_bh))
}

len <- nrow(temp)
cl <- makeCluster(cluster.num)
clusterEvalQ(cl,library(R.utils)) # packages used in function
clusterExport(cl, c("temp","len","si_solve")) # function and variables used in function
result <- parLapply(cl, 1:len, function(i) si_solve(i, temp))

Y2BH <- unlist(lapply(result, function(x) x$Y2BH))
si_bh <- unlist(lapply(result, function(x) x$si_bh))
si_t <- unlist(lapply(result, function(x) x$si_t))

temp$Y2BH <- Y2BH
temp$si_bh <- si_bh
temp$si_t <- si_t

stopCluster(cl)
################################################################################
# solve for total age 
# 
# 
# agesolve <- function(i,df){
#   b1 = df$b1[i]
#   b2 = df$b2[i]
#   b3 = df$b3[i]
#   b4 = df$b4[i]
#   si_bh = df$si_bh[i]
#   topht = df$topht[i]
#   si0 = 10
#   si1 = 100
#   tryCatch({ # These functions add a max time and an error message if a row takes too long
#     withTimeout({
#       while(abs(si0-si1)>0.00000001){
#         if(df$species[i]=="AW"){
#           k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^2+b4*sqrt(50))
#           k2 = si0^(b3*log(si0))
#           k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
#           Y2BH = exp((log(k3)/b2)^2)-1
#           x10 = (1+exp(b1+b2*sqrt(log(50+Y2BH+1.0))+b3*log(si0)^2+b4*sqrt(50)))
#           x20 = (1+exp(b1+b2*sqrt(log(50+1.0))+b3*log(si0)^2+b4*sqrt(50)))
#         }
#         if(df$species[i]=="SW"){
#           k1 = exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50))
#           k2 = si0^(b3*log(si0))
#           k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
#           Y2BH = sqrt(exp((log(k3)/b2)^2)-1)
#           x10 = (1+exp(b1+b2*sqrt(log((50+Y2BH^2+1)))+b3*log(si0)^2+b4*sqrt(50)))
#           x20 = (1+exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50)))
#         }
#         if(df$species[i]%in%c("PL","SB")){
#           k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^1+b4*sqrt(50))
#           k2 = si0^b3
#           k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
#           Y2BH = exp((log(k3)/b2)^2)-1
#           x10 = (1+exp(b1+b2*sqrt(log(50+Y2BH+1))+b3*log(si0)^1+b4*sqrt(50)))
#           x20 = (1+exp(b1+b2*sqrt(log(50+1))+b3*log(si0)^1+b4*sqrt(50)))
#         }
#         si1 = si_bh*x10/x20
#         si0 = (si0+si1)/2
#       }
#       k3x <- ((si0)*(1+k1)/topht-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
#       age_tot <- exp((log(k3x)/b2)^2)-1
#     }, timeout=1)
#   }, TimeoutException=function(ex){
#     message(paste0("Timeout. Skipping row'",i,"'")) # Error message and returned value on timeout
#     age_tot <- NULL
#   })
#   age_tot
# }
# 
# len <- nrow(temp)
# cl <- makeCluster(cluster.num)
# clusterEvalQ(cl,library(R.utils)) # packages used in function
# clusterExport(cl, c("temp","len","agesolve")) # function and variables used in function
# temp$age_tot <- unlist(parLapply(cl, 1:len, function(i) agesolve(i,temp)))
# stopCluster(cl)


# Compute bhage
site2 <- temp %>%
  mutate(
    bhage = ifelse(Y2BH > 0, totage - Y2BH, NA_real_))%>%
  select(-c(b1, b2, b3, b4))%>% 
  arrange(company, company_plot_number, species)

# average si by species
# si_n has a problem
site3 <- site2 %>%
  group_by(company, company_plot_number, species) %>%
  summarise(si_bh = mean(si_bh, na.rm = TRUE)) %>%
  mutate(scale = 'species', si_n = n()) 

plot_level2 <- plot_level2 %>% 
  arrange(company, company_plot_number, scale, species)

site4 <- site3 %>% 
  arrange(company, company_plot_number, scale, species)
 
plot_level2 <- plot_level_2 %>%
  arrange(company, company_plot_number, scale, species)

site5 <- left_join(plot_level2, site3, by = c("company", "company_plot_number", "scale", "species")) 

################################################################################
#  STEP 3. calculate top height for each plot/measurement/species group
################################################################################

trees <- left_join(tree_list1, plot_mmt, by = c("company", "company_plot_number", "measurement_number"))

trees <- trees %>%
  filter(!(tree_origin %in% c(9, 10)), height >= 1.3) 

###################
#Round 1 selection
###################
trees2 <- trees %>%
  mutate(target = ifelse(round(tree_plot_area / 100, 0) == 1, 2, round(tree_plot_area / 100, 0)),
         unique2 = paste(company, company_plot_number, measurement_number, species, sep = "_")) %>%
  select(unique2, company, company_plot_number, measurement_number, tree_number, species, dbh, height, condition_code1, severity1, condition_code2, severity2, target) %>%
  arrange(company, company_plot_number, measurement_number, species)

th <- trees2 %>%
  group_by(unique2) %>%
  arrange(desc(dbh)) %>%
  mutate(num = row_number()) %>%
  filter(num <= target) %>%
  ungroup()%>%
  arrange(company, company_plot_number, measurement_number, species)

th3 <- th %>%
  group_by(unique2) %>%
  summarise( 
    company = first(company),
    company_plot_number = first(company_plot_number),
    measurement_number = first(measurement_number),
    species = first(species),
    FREQ = first(num),
          minht = min(height),
          maxht = max(height),
          th_aa = mean(height),
          th_n_a = n()) %>%
  ungroup()%>%
  arrange(company, company_plot_number, measurement_number, species)
###################
#Round 2 selection
#select top height trees from measured OR predicted/adjusted heights, excluding leaning or swept trees and broken tops;
###################
trees2a <- trees %>%
  filter(!(condition_code1 == 3 & severity1 != 1),
          !(condition_code2 %in% c(8, 10) & severity2 %in% c(2, 3)),
           ht_stat %in% c('MM', 'PA', 'Pa', 'IX')) %>%
  mutate(target = ifelse(round(tree_plot_area / 100, 0) == 1, 2, round(tree_plot_area / 100, 0)),
         unique2 = paste(company, company_plot_number, measurement_number, species, sep = "_")) %>%
  select(unique2, company, company_plot_number, measurement_number, tree_number, species, dbh, height, condition_code1, severity1, condition_code2, severity2, target) 

tha <- trees2a %>%
  group_by(unique2) %>%
  arrange(desc(dbh)) %>%
  mutate(num = row_number()) %>%
  filter(num <= target) %>%
  ungroup()

th3a <- tha %>%
  group_by(unique2) %>%
  summarise( 
    company = first(company),
    company_plot_number = first(company_plot_number),
    measurement_number = first(measurement_number),
    species = first(species),
    FREQ = first(num),
    minht = min(height),
    maxht = max(height),
    th_mp = mean(height),
    th_n_mp = n()) %>%
  ungroup()%>%
  arrange(company, company_plot_number, measurement_number, species)

###################
#Round 3 selection
#select top height trees from measured heights, excluding leaning or swept trees and broken tops;
###################
trees2b <- trees %>%
  filter(!(condition_code1 == 3 & severity1 != 1),
         !(condition_code2 %in% c(8, 10) & severity2 %in% c(2, 3)),
         !(condition_code3 %in% c(8, 10) & severity3 %in% c(2, 3)),
         ht_stat %in% 'MM') %>%
  mutate(target = ifelse(round(tree_plot_area / 100, 0) == 1, 2, round(tree_plot_area / 100, 0)),
         unique2 = paste(company, company_plot_number, measurement_number, species, sep = "_")) %>%
  select(unique2, company, company_plot_number, measurement_number, tree_number, species, spp_grp, dbh, height, condition_code1, severity1, condition_code2, severity2, target)

thb <- trees2b %>%
  group_by(unique2) %>%
  arrange(desc(dbh)) %>%
  mutate(num = row_number()) %>%
  filter(num <= target) %>%
  ungroup()

th3b <- thb %>%
  group_by(unique2) %>%
  summarise( 
    company = first(company),
    company_plot_number = first(company_plot_number),
    measurement_number = first(measurement_number),
    species = first(species),
    FREQ = first(num),
    minht = min(height),
    maxht = max(height),
    th_mm = mean(height),
    th_n_mm = n()) %>%
  ungroup()%>%
  arrange(company, company_plot_number, measurement_number, species)

###################
#Round 4 selection(strictest) 
#average the measured heights of the thickest 100 sph of each species regardless of whether they have a height mmt or not
###################
trees2c <- trees %>%
  filter(!(condition_code1 == 3 & severity1 != 1) &
           !(condition_code1 %in% c(8, 10) & severity1 %in% c(2, 3)) &
           !(condition_code2 %in% c(8, 10) & severity2 %in% c(2, 3)) &
           !(condition_code3 %in% c(8, 10) & severity3 %in% c(2, 3))) %>%
  mutate(
    height = if_else(!(ht_stat %in% 'MM'), NA, height )
  ) %>%
  mutate(target = ifelse(round(tree_plot_area / 100, 0) == 1, 2, round(tree_plot_area / 100, 0)),
         unique2 = paste(company, company_plot_number, measurement_number, species, sep = "_")) %>%
  select(unique2, company, company_plot_number, measurement_number, tree_number, species, spp_grp, dbh, height, ht_stat, condition_code1, severity1, condition_code2, severity2, condition_code3, severity3, target)

thc <- trees2c %>%
  group_by(unique2) %>%
  arrange(desc(dbh)) %>%
  mutate(num = row_number()) %>%
  filter(num <= target) %>%
  ungroup()

th3c <- thc %>%
  group_by(unique2) %>%
  summarise( 
    company = first(company),
    company_plot_number = first(company_plot_number),
    measurement_number = first(measurement_number),
    species = first(species),
    FREQ = first(num),
    minht = min(height),
    maxht = max(height),
    th_ms = mean(height, na.rm = TRUE),
    th_n_ms = sum(!is.na(height))) %>%
  ungroup()%>%
  arrange(company, company_plot_number, measurement_number, species)

###################
#Round 5 selection
#Exclude lower cohorts by the empirical range of 6m. 
#if the current stem is less than the max C and D crown class ht (veterans removed) minus this range variable, do not include this in the top ht population;
###################
trees2d <- trees %>%
  mutate(
    target = round(tree_plot_area / 100, 1),
    target = ifelse(target == 1, 2, target),
    unique2 = paste(company, "_", company_plot_number, "_", measurement_number, "_", species, sep = ""),
    keep_flag = ifelse(target <= 1, TRUE, FALSE)
  ) %>%
  select(unique2, company, company_plot_number, measurement_number, tree_number, species, spp_grp, dbh, height, ht_stat, condition_code1, severity1, condition_code2, severity2, condition_code3, severity3, crown_class, target) %>%
  arrange(desc(unique2), desc(dbh))

thd <- trees2d %>%
  group_by(unique2) %>%
  mutate(
    maxht = first(height)
  ) %>%
  filter(row_number() <= target & height >= maxht - 6) %>%
  ungroup() %>%
  arrange(company, company_plot_number, measurement_number, species)

th3d <- thd %>%
  group_by(company, company_plot_number, measurement_number, species) %>%
  summarise(
    th_mpe = mean(height, na.rm = TRUE),
    th_n_mpe = n()
  )

# Merge all datasets
th4 <- left_join(th3, th3a, by = c("company", "company_plot_number", "measurement_number", "species"))
th4aaa <- left_join(th4, th3b, by = c("company", "company_plot_number", "measurement_number", "species"))
th4aa <- left_join(th4aaa, th3c, by = c("company", "company_plot_number", "measurement_number", "species"))
th4a <- left_join(th4aa, th3d, by = c("company", "company_plot_number", "measurement_number", "species"))

# th_aa, th_ms do not match with sas output
th4a <- th4a %>%
  select(company, company_plot_number, measurement_number, species, FREQ.x, th_aa, th_n_a, th_mp, th_n_mp, th_mm, th_n_mm, th_ms, th_n_ms, th_mpe, th_n_mpe)

#################################
trees_target2 <- trees2 %>% 
  group_by(unique2) %>% 
  filter(row_number() == 1) %>%
  ungroup() %>%
  select(
    company, company_plot_number, measurement_number, species, target
  )

trees_target3 <- left_join(trees_target2, th4a, by = c("company", "company_plot_number", "measurement_number", "species"))

# assign based on meeting target number of top height trees;
trees_target4 <- trees_target3 %>%
  mutate(
    th_m = 0,
    th_n_m = 0,
    th_n_aa = 0
  ) %>%
  mutate(                                         #what are th_n_m and th_m and th_n_aa? why are the names different
    topht = case_when(
      th_n_m >= (target - 1) & th_n_m >= 2 ~ th_m,
      th_n_mp == th_n_m & th_n_m >= 2 ~ th_m,
      th_n_mp >= 2 ~ th_mp,
      th_n_m == 1 & th_n_aa == 1 ~ th_m,
      th_n_mp == 1 & th_n_aa == 1 ~ th_mp,
      TRUE ~ th_aa
    ),
    topht_n = case_when(
      th_n_m >= (target - 1) & th_n_m >= 2 ~ th_n_m,
      th_n_mp == th_n_m & th_n_m >= 2 ~ th_n_m,
      th_n_mp >= 2 ~ th_n_mp,
      th_n_m == 1 & th_n_aa == 1 ~ th_n_m,
      th_n_mp == 1 & th_n_aa == 1 ~ th_n_mp,
      TRUE ~ th_n_aa
    ),
    topht_stat = case_when(
      th_n_m >= (target - 1) & th_n_m >= 2 ~ 'MM',
      th_n_mp == th_n_m & th_n_m >= 2 ~ 'MM',
      th_n_mp >= 2 ~ 'MP',
      th_n_m == 1 & th_n_aa == 1 ~ 'MM',
      th_n_mp == 1 & th_n_aa == 1 ~ 'MP',
      TRUE ~ 'AA'
    ),
    scale = 'species'
  ) %>%
  mutate(
    th_m = ifelse(th_m == "0", "", th_m),
    th_n_m = ifelse(th_n_m == "0", "", th_n_m),
    th_n_aa = ifelse(th_n_aa == "0", "", th_n_aa),
    topht_n = ifelse(topht_n == "0", "", topht_n)
  )

site6 <- site5 %>%
  left_join(
    trees_target4, by = c("company", "company_plot_number", "measurement_number", "scale", "species")
  )%>%
  arrange(company, company_plot_number, measurement_number, scale, species)

# check to make sure that everything with a measurement has a topheight;
check <- site6 %>%
  filter(
    scale == 'species' & sph > 0 & (is.na(topht) | topht == 0),
    scale == 'species' & sph == 0 & topht > 0
  )


################################################################################
#  STEP 4. calculate average age by species as well as maximum age;
################################################################################

age1 <- managed_natural %>%
  mutate(
    birth_yr_tot = if_else(totage > 0, measurement_year - totage, NA_real_),
    birth_yr_bh = if_else(bhage > 0, measurement_year - bhage, NA_real_)
  )

age3 <- age1 %>%
  group_by(company, company_plot_number, species) %>%
  summarise(
    birth_yr_tot = mean(birth_yr_tot, na.rm = TRUE),
    birth_yr_bh = mean(birth_yr_bh, na.rm = TRUE),
    FREQ.x = n() ) %>%
  mutate(scale = 'species') %>%
  filter(
    !(species %in% c(50))
  )

site7 <- site6 %>%
  left_join(age3, by = c("company", "company_plot_number", "scale", "species"))%>%
  arrange(company, company_plot_number, measurement_number, scale, species)

site8 <- site7 %>%
  mutate(age_tot = measurement_year - birth_yr_tot,
         age_bh = measurement_year - birth_yr_bh) %>%
  mutate(age_n = ifelse(is.na(age_tot), NA, age_tot)) %>%
  select(-birth_yr_tot, -birth_yr_bh)%>%
  arrange(company, company_plot_number, measurement_number, scale, species)

site9 <- site8 %>%
  filter(scale == 'species') %>%
  mutate(unique_age = paste(company, company_plot_number, measurement_number, sep = '_'))%>%
  arrange(desc(unique_age), desc(age_tot))

site10 <- site9 %>%
  group_by(unique_age) %>%
  slice(1) %>%
  rename(age_max = age_tot) %>%
  select(company, company_plot_number, measurement_number, age_max)

site11 <- left_join(site8, site10, by = c("company", "company_plot_number", "measurement_number"))

site11$standage <- ifelse(!is.na(site11$age_max), site11$age_max,
                          ifelse(!is.na(site11$age_harv), site11$age_harv,
                                 ifelse(!is.na(site11$age_fire), site11$age_fire, NA)))

pgyi_compiled <- site11 %>%
  mutate(
    standage = as.integer(standage)
  ) %>%
  select(
    company, company_plot_number, establishment_year, measurement_number, measurement_year, measurement_month,
    standage, stand_type, scale, species, sphRegen, sphBH, sphD15, sphD91, ba, vol_0000, vol_1307, vol_1510,
    biomass, carbon, topht_n, topht, topht_stat, age_n, age_tot, age_bh, age_max, age_harv, age_fire, age_avi,
    si_n, si_bh, fmu, fma, opening_number, sampling_unit_number, topographic_position, elevation, slope, aspect,
    x_coord, y_coord, utm_zone, datum, latitude, longitude, natural_subregion, ecosite_guide, ecosite,
    ecosite_phase, shared, stand_origin, plot_type, tree_plot_area, tree_tagging_limit, sapling_plot_area,
    sapling_tagging_limit_dbh, sapling_tagging_limit_height, regen_plot_area, regen_tagging_limit_conifer,
    regen_tagging_limit_decid, number_regen_plots
  ) %>%
  rename(
    plot = company_plot_number,
    est_yr = establishment_year,
    mmt_num = measurement_number,
    mmt_yr = measurement_year,
    mmt_mo = measurement_month,
    age = standage,
    sphRegenOnly = sphRegen,
    opening = opening_number,
    su = sampling_unit_number,
    topo = topographic_position,
    elev = elevation,
    utm = utm_zone,
    natsub = natural_subregion,
    ecoguide = ecosite_guide,
    eco = ecosite,
    ecophs = ecosite_phase,
    tree_area = tree_plot_area,
    tree_tag = tree_tagging_limit,
    sap_area = sapling_plot_area,
    sap_dbh = sapling_tagging_limit_dbh,
    sap_ht = sapling_tagging_limit_height,
    regen_area = regen_plot_area,
    regen_conht = regen_tagging_limit_conifer,
    regen_decht = regen_tagging_limit_decid,
    regen_n = number_regen_plots
  )  


fwrite(pgyi_compiled, ""GYPSY data/intermediate/i_pgyi_compiled.csv")
