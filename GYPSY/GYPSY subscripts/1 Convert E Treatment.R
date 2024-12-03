## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply "data smoothing" rules
# Read .csv, convert "." to NA
treatment <- fread(paste0(input.dir,"/treatment.csv"), na.string=".") %>%
  distinct

fwrite(treatment,"GYPSY data/intermediate/i_treatment.csv",row.names=F)