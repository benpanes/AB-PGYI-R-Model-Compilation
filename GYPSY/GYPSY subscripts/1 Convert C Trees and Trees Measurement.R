## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply "data smoothing" rules
# Read .csv, convert "." to NA
trees <- fread(paste0(input.dir,"/trees.csv"), na.string=".")
trees_measurement <- fread(paste0(input.dir,"/trees_measurement.csv"), na.string=".")
trees$id <- paste(trees$company,trees$company_plot_number,trees$tree_number,sep="_")
trees_measurement$id <- paste(trees_measurement$company,trees_measurement$company_plot_number,trees_measurement$tree_number,sep="_")

# Report trees with no measurements or vice versa
id_list <- unique(append(trees$id,trees_measurement$id))
trees_report <- data.frame(id = id_list, error = as.character(NA))
trees_report$error[!(id_list%in%unique(trees$id))] <- "No record"
trees_report$error[!(id_list%in%unique(trees_measurement$id))] <- "No measurement"
trees_report <- na.omit(trees_report)
fwrite(trees_report,"GYPSY data/reports/unmatched tree records.csv")

# Merge and filter duplicate columns
trees_merge <- inner_join(trees,trees_measurement[,!(c("company","company_plot_number","tree_number"))],by="id")
# Fix suspected zero filling
trees_merge$azimuth[trees_merge$azimuth==0 & 
                    is.na(trees_merge$distance)] <- NA
trees_merge$distance[is.na(trees_merge$azimuth) & 
                     trees_merge$distance==0] <- NA
# Reset 0 azimuth to 360
trees_merge$azimuth[trees_merge$azimuth==0 & 
                    trees_merge$distance>0] <- 360
# Assumptions in absence of other info
trees_merge$rcd_height[trees_merge$rcd>0 & 
                       is.na(trees_merge$rcd_height)] <- 0.1
trees_merge$dbh_height[trees_merge$dbh>0 & 
                       is.na(trees_merge$dbh_height)] <- 1.3
trees_merge$stump_height[trees_merge$stump_age>0 & 
                         is.na(trees_merge$stump_height)] <- 0.3
i <- trees_merge$condition_code1==0
trees_merge$cause1[i] <- 99
trees_merge$severity1[i] <- 9
trees_merge[trees_merge$condition_code1%in%c(1,2,13,14,15,16),
            c("severity1","severity2","severity3")] <- 9
# Fix non-buffer trees (tree_type="B") outside of plot (tree_location="BP")
# Remove assumed trees outside plot
trees_merge <- trees_merge[!(trees_merge$tree_type!="B" & 
                           is.na(trees_merge$dbh_age) &
                           is.na(trees_merge$stump_age) &
                           is.na(trees_merge$total_age) &
                           trees_merge$tree_location=="BP"),]
# Reclassify assumed age trees outside of plot misclassified as T or S
trees_merge$tree_type[trees_merge$tree_type!="B" & 
                      (!is.na(trees_merge$dbh_age)|
                       !is.na(trees_merge$stump_age)|
                       !is.na(trees_merge$total_age))] <- "B"

# Fix improper use of species="No"
trees_merge$species[trees_merge$species=="No" &
                    trees_merge$tree_origin!=0 &
                    trees_merge$tree_number!=0 &
                    trees_merge$tree_location!="U"] <- "Ms"
# Fix broadleafs classified as artificial seed or planted origin
trees_merge$tree_origin[trees_merge$species%in%c("Aw","Bw","Pb") &
                        trees_merge$tree_origin%in%c(5,6,7,8)] <- 0
# Fix needleleafs classified as coppice or sucker origin
trees_merge$tree_origin[trees_merge$species%in%c("Fa","Fb","Pl","Sb","Sw") &
                        trees_merge$tree_origin==2] <- 0
# Remove nil plots
trees_merge$m_id <- paste(trees_merge$company,trees_merge$company_plot_number,trees_merge$measurement_number,sep="_")
m_count <- as.data.frame(table(trees_merge$m_id))
m_count <- data.frame(m_id=as.character(m_count[,1]),m_count=m_count[,2])
trees_merge <- merge(trees_merge,m_count,by="m_id")
trees_merge <- trees_merge[!(trees_merge$species=="No" &
                           trees_merge$m_count==1),
                           !c("m_id","m_count","id")]
# Check species="No" individually to see if fix is possible
fwrite(trees_merge[trees_merge$species=="No",],"GYPSY data/reports/No species.csv")

fwrite(trees_merge,"GYPSY data/intermediate/i_trees_merge.csv")