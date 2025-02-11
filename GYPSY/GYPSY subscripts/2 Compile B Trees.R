#/**********************************************************/
# Compile individual tree basal area, volumes, sph
#/**********************************************************/
start_time <- Sys.time()

trees <- fread("GYPSY data/intermediate/i_trees_measurement_filled.csv")

trees <- trees[!(trees$tree_type%in%c("T","ET") & (trees$tree_plot_area==0 | is.na(trees$tree_plot_area))) &
               !(trees$tree_type%in%c("S","ES") & (trees$sapling_plot_area==0 | is.na(trees$sapling_plot_area))) &
               !(trees$tree_type%in%c(paste0("R",1:10)) & (trees$regen_plot_area==0 | is.na(trees$regen_plot_area))),]

trees[,index_1:=paste(company, company_plot_number, measurement_number, tree_number, sep = "_")]
 
trees[,ba:=((dbh/200)^2)*3.14159]


trees[tree_type%in%c("T","ET"),sph:= 10000/tree_plot_area]

trees[tree_type%in%c("S1","ES1"),sph:= ifelse(sapling_plot_area==0,10000/tree_plot_area,10000/sapling_plot_area)] # Assume all saplings were measured if saplings present and sapling plot area = 0

trees[tree_type%in%c(paste0("R",1:10)),sph:= 10000/(regen_plot_area*number_regen_plots)]

trees[tree_type=="B",c("sph","ba"):=0] # test first
#trees[tree_type=="B",c("sph","ba"):=list(0,0)]

# assign new species names 
trees[,species:=fcase(
  species_og %in% c("AW","AX"), "AW",
  species_og %in% c("BW","PB","PJ","FD","SB"), species,
  species_og %in% c("LA", "LT", "LW"), "LT",
  species_og %in% c("P", "PX", "PF", "PL", "PW"), "PL",
  species_og %in% c("FA", "FB"), "FB",
  species_og %in% c("SE", "SW", "SX"), "SW")]

# changed species to species_og for testing purpose
trees[species%in%c("BW","FD","LT","PJ","SE"), natsub:=0]

trees[,table(natural_subregion,natsub)]

# Volume function
calc_volume <- function(table,spp="species",ht="height",dbh="dbh",natsub="natural_subregion",merch="blank",max.i=1000){
  hr <- function(df,topdib,max.i=1000){
    df[,diff:=ifelse(dbh>topdib,1,NA)]
    df[,r0:=ifelse(dbh>topdib,1,NA)]
    df[,r1:=ifelse(dbh>topdib,1,NA)]
    hr1 <- function(b1,b2,b3,b4,b5,dbh,ht,r0){
      b1*r0^2+b2*log(r0+0.001)+b3*sqrt(r0)+b4*exp(r0)+b5*(dbh/ht)
    }
    hr2 <- function(a0,a1,a2,b1,b2,b3,b4,b5,dbh,ht,r0,topdib){
      (1-((topdib/(a0*dbh^a1*a2^dbh))^(1/hr1(b1,b2,b3,b4,b5,dbh,ht,r0)))*(1-sqrt(0.2250)))^2
    }
    n=0
    while(df[,max(diff,na.rm=T)]>0.00000001){
      df[dbh>topdib & diff>0.00000001,r1:=hr2(a0,a1,a2,b1,b2,b3,b4,b5,dbh,ht,r0,topdib)]
      df[dbh>topdib & diff>0.00000001,r0:=(r1+r0)/2]
      df[dbh>topdib & diff>0.00000001,diff:=abs(r1-r0)]
      n=n+1
      if(n>max.i){
        message("Max height ratio iterations reached.")
        message(paste("company","company_plot_number","tree_number","measurement_number",natsub,spp,ht,dbh,"diff"))
        for(i in 1:df[diff>0.00000001,.N]){
          message(df[diff>0.00000001,paste(company,company_plot_number,tree_number,measurement_number,natsub,spp,ht,dbh,diff)][i])
        }
        break
      }
    }
  }
  dib <- function(a0,a1,a2,b1,b2,b3,b4,b5,dbh,dheight,ht,r=(dheight/ht)){
    (a0*dbh^a1)*(a2^dbh)*((1-sqrt(0.2250)))^(b1*(r)^2+b2*log(r+0.001)+b3*sqrt(r)+b4*exp(r)+b5*dbh/ht)
  }
  fvol <- function(sl,stumpH,dbh,ht,a0,a1,a2,b1,b2,b3,b4,b5){
    svol <- function(sl,d0,d1,d2){(2*sl/6)*((1/200)^2*pi)*(d0^2+4*d1^2+d2^2)}
    v <- data.table()
    h0 <- stumpH
    d0 <- dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,stumpH,ht)
    h1 <- h0+sl
    d1 <- dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,h1,ht)
    h2 <- h1+sl
    d2 <- dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,h2,ht)
    v[,v1:=svol(sl,d0,d1,d2)]
    for(n in 2:10){
      h0 <- h2
      d0 <- d2
      h1 <- h0+sl
      d1 <- dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,h1,ht)
      h2 <- h1+sl
      d2 <- dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,h2,ht)
      v[,paste0("v",n):=svol(sl,d0,d1,d2)]
    }
    rowSums(v)
  }
  if(merch=="blank"){
    message("Merchantability not specified. Calculating total volume by default")
    topdib = 0.0125
    merch="total"
  } else if(merch=="total"){
    topdib = 0.0125
  } else if(merch=="1307"){
    topdib = 7
    stumpD = 13
    minlen = 3.66
  } else if(merch=="1510"){
    topdib = 10
    stumpD = 15
    minlen = 3.66
  } else {
    stop("Invalid merch parameter. Use 'total', '1307', or '1510'")
  }
  df <- copy(data.table(table))
  setnames(df,c(dbh,ht,spp,natsub),c("dbh","ht","spp","natsub"),skip_absent=T)
  if(any(is.na(dbh) | is.na(ht))){
    stop("Table contains missing height or DBH.")
  }
  stumpH = 0.3
  df[,spp:=toupper(spp)]
  if(any(!df[,spp]%in%c("AW","AX","BW","PB",
                        "LA","LT","LW",
                        "P","PX","PF","PL","PW","PJ",
                        "FA","FB","FD",
                        "SE","SW","SX",
                        "SB"))){
    stop("Invalid species. Accepted codes are: AW, AX, BW, PB, LA, LT, LW, P, PX, PF, PL, PW, PJ, FA, FB, FD, SE, SW, SX, or SB")
  }
  if(any(!df[,natsub]%in%c("CM","DM","NM","BSA","PAD","UBH","LBH","ALP","SA","MT","UF","LF","AP","KU","FP","PRP","CP","DMG","FF","NF","MG"))){
    stop("Invalid natural subregion. Accepted codes are: CM, DM, NM, BSA, PAD, UBH, LBH, ALP, SA, MT, UF, LF, AP, KU, FP, PRP, CP, DMG, FF, NF, MG")
  }
  df[,spp:=fcase(spp%in%c("AW","AX"),"AW",
                 spp%in%c("LA","LT","LW"),"LT",
                 spp%in%c("P","PX","PF","PL","PW"),"PL",
                 spp%in%c("FA","FB"),"FB",
                 spp%in%c("SE","SW","SX"),"SW",
                 spp%in%c("BW","PB","PJ","FD","SB"),spp)]
  df[,c(                                                                                                                   "a0",     "a1",     "a2",     "b1",      "b2",      "b3",      "b4",      "b5",      "k7",      "k8"):=
       rbindlist(fcase(
         spp=="AW" &           natsub%in%c(NA,"DMG","FF","NF","MG"),.(.(                                                   0.790406, 1.026943, 0.997524, 0.600584,  -0.065681, -0.173812, 0.121363,  0.063253,  0.091433,  1.077082)),
         spp=="AW" &           natsub%in%c("CM","NM","BSA","PAD","UBH","LBH","AP","KU"),.(.(                               0.841897, 0.997064, 0.998713, 0.536865,  -0.06402,  -0.234471, 0.179963,  0.03155,   0.061755,  1.079512)),
         spp=="AW" &           natsub%in%c("DM","FP","PRP","CP"),.(.(                                                      0.944522, 0.93803,  1.001644, 0.695363,  -0.067849, 0.050603,  -0.01633,  0.116432,  0.024127,  1.081189)),
         spp=="AW" &           natsub%in%c("ALP","SA","UF"),.(.(                                                           0.588838, 1.161895, 0.992096, 0.7093,    -0.075446, -0.116041, 0.040949,  0.113638,  0.21134,   1.073573)),
         spp=="AW" &           natsub%in%c("MT","LF"),.(.(                                                                 0.905615, 0.964894, 1.000054, 0.553236,  -0.049737, -0.280768, 0.170687,  0.075789,  0.135261,  1.072734)),
         spp=="BW",.(.(                                                                                                    0.894358, 1.007721, 0.991384, -0.483072, 0.155593,  -2.273122, 1.326501,  0.168897,  0.077197,  1.062798)),
         spp%in%c("FA","FB") & natsub%in%c(NA,"DMG","FF","NF","MG"),.(.(                                                   1.002016, 0.944076, 0.999921, 1.33633,   -0.320352, 2.839497,  -1.324815, 0.077452,  0.28994,   1.051003)),
         spp%in%c("FA","FB") & natsub%in%c("CM","DM","NM","BSA","PAD","UBH","LBH","MT","LF","AP","KU","FP","PRP","CP"),.(.(0.918647, 0.990225, 0.997292, 1.568514,  -0.384262, 3.503466,  -1.677185, 0.128169,  0.249035,  1.050236)),
         spp%in%c("FA","FB") & natsub%in%c("ALP","SA","UF"),.(.(                                                           1.108006, 0.89838,  1.001816, 1.338336,  -0.30463,  2.694363,  -1.277617, 0.087438,  0.323616,  1.050716)),
         spp=="FD",.(.(                                                                                                    0.913153, 0.964386, 0.998391, 1.386315,  -0.286495, 1.783899,  -0.916932, 0.05883,   -0.095253, 1.123153)),
         spp=="LT",.(.(                                                                                                    0.933517, 0.965471, 0.998393, 2.079455,  -0.462028, 3.732057,  -1.950194, 0.190425,  0.37887,   1.034008)),
         spp=="PB" &           natsub%in%c(NA,"DMG","FF","NF","MG"),.(.(                                                   0.861179, 0.951483, 1.000957, 0.752581,  -0.167305, 0.693611,  -0.224137, 0.008214,  0.149322,  1.117988)),
         spp=="PB" &           natsub%in%c("CM","DM","NM","BSA","PAD","UBH","LBH","AP","KU","PRP","CP"),.(.(               0.80437,  0.982874, 0.999527, 0.996958,  -0.223248, 1.106731,  -0.459817, -0.003392, 0.109731,  1.12512)),
         spp=="PB" &           natsub%in%c("ALP","SA","MT","UF","LF","FP"),.(.(                                            0.913329, 0.92259,  1.002574, 0.308448,  -0.06567,  -0.10213,  0.226336,  0.023148,  0.257085,  1.10315)),
         spp=="PJ",.(.(                                                                                                    0.940832, 0.955575, 0.999333, 0.116311,  -0.028172, -0.384427, 0.304055,  0.072192,  0.161727,  1.045672)),
         spp=="PL" &           natsub%in%c(NA,"DMG","FF","NF","MG"),.(.(                                                   0.897617, 0.988518, 0.998735, 0.675759,  -0.130313, 0.570634,  -0.275457, 0.105403,  0.283173,  1.025305)),
         spp=="PL" &           natsub%in%c("CM","DM","NM","PAD","AP","KU","PRP","CP"),.(.(                                 1.033572, 0.913621, 1.000765, 0.256633,  -0.049091, -0.252118, 0.174267,  0.123722,  0.189744,  1.04681)),
         spp=="PL" &           natsub%in%c("BSA","UF"),.(.(                                                                0.828665, 1.024196, 0.997492, 0.596193,  -0.118777, 0.465591,  -0.196176, 0.083094,  0.308258,  1.024549)),
         spp=="PL" &           natsub%in%c("UBH","LBH","MT","LF","FP"),.(.(                                                0.957164, 0.959992, 0.999774, 0.766747,  -0.140758, 0.666037,  -0.35505,  0.13214,   0.294015,  1.024582)),
         spp=="PL" &           natsub%in%c("ALP","SA"),.(.(                                                                0.800648, 1.053544, 0.995568, 0.568347,  -0.125114, 0.610085,  -0.238442, 0.045398,  0.240347,  1.020105)),
         spp=="SB" &           natsub%in%c(NA,"DMG","FF","MF","MG"),.(.(                                                   0.940695, 0.957211, 0.99964,  1.395784,  -0.344672, 2.835917,  -1.39646,  0.152487,  0.382746,  1.033405)),
         spp=="SB" &           natsub%in%c("CM","DM","NM","BSA","PAD","UBH","LBH","AP","KU","FP","PRP","CP"),.(.(          0.929037, 0.967718, 0.998511, 1.236597,  -0.308204, 2.535507,  -1.22206,  0.146243,  0.349765,  1.036689)),
         spp=="SB" &           natsub%in%c("ALP","SA","MT","UF","LF"),.(.(                                                 0.957624, 0.94674,  1.000452, 1.430462,  -0.356702, 2.950725,  -1.455471, 0.154263,  0.414614,  1.030781)),
         spp=="SE",.(.(                                                                                                    1.072576, 0.897766, 1.001919, 1.301834,  -0.305439, 2.265717,  -1.119671, 0.123519,  0.46127,   1.023752)),
         spp=="SW" &           natsub%in%c(NA,"DMG","FF","NF","MG"),.(.(                                                   0.860438, 0.995406, 0.998493, 1.040218,  -0.252387, 1.842818,  -0.852227, 0.110359,  0.484768,  1.024893)),
         spp=="SW" &           natsub%in%c("CM","DM","NM","BSA","PAD","UBH","LBH","AP","KU","PRP","CP"),.(.(               0.903528, 0.975136, 0.999018, 0.846981,  -0.244969, 1.783097,  -0.730236, 0.040997,  0.413577,  1.028342)),
         spp=="SW" &           natsub%in%c("ALP","SA","UF"),.(.(                                                           0.713393, 1.071533, 0.996067, 1.153679,  -0.283807, 2.022713,  -0.953783, 0.101608,  0.521645,  1.024172)),
         spp=="SW" &           natsub%in%c("MT","LF","FP"),.(.(                                                            0.862685, 0.993148, 0.998773, 1.135018,  -0.252377, 1.885321,  -0.921437, 0.150228,  0.536767,  1.0227))))]
  hr(df,topdib)
  df[,dibs:=dib(a0,a1,a2,b1,b2,b3,b4,b5,dbh,stumpH,ht)]
  df[,dobs:=k7+k8*dibs]
  if(merch=="total"){
    df[,mh:=ht*r0]
    df[,ml:=mh-stumpH]
    df[,sl:=ml/20]
    df[,vol:=fvol(sl,stumpH,dbh,ht,a0,a1,a2,b1,b2,b3,b4,b5)]
    df[,tvol:=vol+pi*(topdib/200)^2*(ht-mh)/3+pi*(dibs/200)^2*stumpH]
    df[is.na(tvol),tvol:=0]
    return(df[,tvol ])
  } else {
    df[dobs>stumpD,mh:=ht*r0]
    df[dobs>stumpD,ml:=mh-stumpH]
    df[ml>minlen,sl:=ml/20]
    df[ml>minlen,vol:=fvol(sl,stumpH,dbh,ht,a0,a1,a2,b1,b2,b3,b4,b5)]
    df[is.na(vol),vol:=0]
    return(df[,vol])
  }
}

trees[,vol_0000:=calc_volume(trees,merch="total")]
trees[,vol_1307:=calc_volume(trees,merch="1307")]
trees[,vol_1510:=calc_volume(trees,merch="1510")]

# biomass and carbon calculation

trees.bio <- trees %>%
  mutate(
    biomass = if_else(dbh < 15 | vol_1307 == 0,
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
    carbon = biomass * 0.5,
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

minutes <- floor(Sys.time()-start_time)
seconds <- round((Sys.time()-start_time-minutes)*60)
message(paste(minutes,"minutes and",seconds,"seconds for compile trees"))