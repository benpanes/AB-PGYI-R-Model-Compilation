## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply "data smoothing" rules
# Read .csv, convert "." to NA
trees <- fread(paste0(input.dir,"/trees.csv"), na.string=".")
trees_measurement <- fread(paste0(input.dir,"/trees_measurement.csv"), na.string=".")
trees$id <- paste(trees$company,trees$company_plot_number,trees$tree_number,sep="_")
trees_measurement$id <- paste(trees_measurement$company,trees_measurement$company_plot_number,trees_measurement$tree_number,sep="_")

################################################################################
trees_mmt_chk <- read.csv("H:/Shared drives/Growth & Yield Lab/Data Sets/PGYI/2023-09-08 PGYI Compiled/2_interim/tree_list1.csv") 

# Report trees with no measurements or vice versa
id_list <- unique(append(trees$id,trees_measurement$id))
trees_report <- data.frame(id = id_list, error = as.character(NA))
trees_report$error[!(id_list%in%unique(trees$id))] <- "No record"
trees_report$error[!(id_list%in%unique(trees_measurement$id))] <- "No measurement"
trees_report <- na.omit(trees_report)
if(nrow(trees_report>0)){View(trees_report)} # Show report if errors are found

################################################################################
# Merge and filter duplicate columns
trees_merge <- inner_join(trees,trees_measurement[,!(c("company","company_plot_number","tree_number"))],by="id")
rm(list=c("trees","trees_measurement"));gc()
# Fix suspected zero filling
trees_merge[azimuth==0 & is.na(distance),"azimuth"] <- NA
trees_merge[is.na(azimuth) & distance==0,"distance"] <- NA
# Reset 0 azimuth to 360
trees_merge[azimuth==0 & distance>0,"azimuth"] <- 360
# Assumptions in absence of other info
trees_merge[rcd>0 & is.na(rcd_height),"rcd_height"] <- 0.1
trees_merge[dbh>0 & is.na(dbh_height),"dbh_height"] <- 1.3
trees_merge[stump_age>0 & is.na(stump_height), "stump_height"] <- 0.3
i <- trees_merge$condition_code1==0
trees_merge[condition_code1==0,"cause1"] <- 99
trees_merge[condition_code1==0,"severity1"] <- 9
trees_merge[condition_code1%in%c(1,2,13,14,15,16), c("severity1","severity2","severity3")] <- 9
# Fix non-buffer trees (tree_type="B") outside of plot (tree_location="BP")
# Remove assumed trees outside plot
trees_merge <- trees_merge[!(tree_type!="B" & 
                               is.na(dbh_age) &
                               is.na(stump_age) &
                               is.na(total_age) &
                               tree_location=="BP"),]
# Reclassify assumed age trees outside of plot misclassified as T or S
trees_merge[tree_type!="B" & 
              (!is.na(dbh_age)|
                 !is.na(stump_age)|
                 !is.na(total_age))&
                  tree_location=="BP","tree_type"] <- "B"

# Fix improper use of species="No"
trees_merge[species=="No" &
              tree_origin!=0 &
              tree_number!=0 &
              tree_location!="U","species"] <- "Ms"
# Fix broadleafs classified as artificial seed or planted origin
trees_merge[species%in%c("Aw","Bw","Pb") & tree_origin%in%c(5,6,7,8),"tree_origin"] <- 0
# Fix needleleafs classified as coppice or sucker origin
trees_merge[species%in%c("Fa","Fb","Pl","Sb","Sw") & tree_origin==2,"tree_origin"] <- 0
# Remove nil plots
trees_merge$m_id <- paste(trees_merge$company,trees_merge$company_plot_number,trees_merge$measurement_number,sep="_")
m_count <- as.data.frame(table(trees_merge$m_id))
m_count <- data.frame(m_id=as.character(m_count[,1]),m_count=m_count[,2])
trees_merge <- merge(trees_merge,m_count,by="m_id")
rm("m_count")
trees_merge <- trees_merge[!(species=="No" & m_count==1),
                           !c("m_id","m_count","id")]
# Check species="No" individually to see if fix is possible
#View(trees_merge[trees_merge$species=="No",])

fwrite(trees_merge,"GYPSY data/intermediate/i_trees_merge.csv")
