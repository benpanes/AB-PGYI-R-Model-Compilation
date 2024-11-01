## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply 'data smoothing' rules
# Read .csv, convert "." to NA
plot <- fread(paste0(input.dir,"/plot.csv"), na.string=".")

# Remove suspected zero filling
plot$elevation[plot$elevation==0] <- NA
plot$sampling_unit_number[plot$sampling_unit_number==0] <- ""
plot$aspect[plot$aspect=="" & 
            plot$slope==0] <- NA

# Fill in missing values where possible
plot$ecosite[plot$ecosite_phase!="" & 
             plot$ecosite==""] <- substr(plot$ecosite_phase,1,1)

# Check NSR/ecosite combinations and toss if not eligible
plot$ecosite[plot$natural_subregion%in%c("ALP","DMG","FF","MG","NF","CP","FP","PRP") & 
             plot$ecosite!=""] <- ""
plot$ecosite[plot$natural_subregion%in%c("CM","DM","NM","PAD") &
             plot$ecosite%in%c("m","n")] <- ""
plot$ecosite[plot$natural_subregion%in%c("LBH","UBH") & 
             plot$ecosite%in%c("k","l","m","n")] <- ""
plot$ecosite[plot$natural_subregion%in%c("BSA","AP","KU") & 
             plot$ecosite%in%c("i","j","k","l","m","n")] <- ""
plot$ecosite[plot$natural_subregion%in%c("LF","UF") &
             plot$ecosite%in%c("l","m","n") &
             plot$ecosite_guide=="SW"] <- ""
plot$ecosite[plot$natural_subregion=="UF" &
             plot$ecosite=="n" &
             plot$ecosite_guide=="WC"] <- ""
plot$ecosite[plot$natural_subregion=="MT" &
             plot$ecosite%in%c("h","i","j","k","l","m","n") &
             plot$ecosite_guide=="SW"] <- ""
plot$ecosite[plot$natural_subregion=="MT" &
             plot$ecosite%in%c("i","j","k","l","m","n") &
             plot$ecosite_guide=="WC"] <- ""
plot$ecosite[plot$natural_subregion=="SA" &
             plot$ecosite%in%c("i","j","k","l","m","n") &
             plot$ecosite_guide=="SW"] <- ""
plot$ecosite[plot$natural_subregion=="SA" &
             plot$ecosite%in%c("j","k","l","m","n") &
             plot$ecosite_guide=="WC"] <- ""

# Check NSR/ecosite guide combinations and adjust if not eligible - fore example if WC or SW and BM, reset to N
plot$ecosite_guide[plot$natural_subregion%in%c("LF","UF","MT","SA") & 
                   plot$ecosite_guide=="N"] <- "" 
plot$ecosite_guide[!(plot$natural_subregion%in%c("LF","UF","MT","SA")) & 
                   plot$ecosite_guide%in%c("SW","WC")] <- "N" 

fwrite(plot,"GYPSY data/intermediate/i_plot.csv")