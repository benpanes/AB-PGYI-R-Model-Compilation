## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply "data smoothing" rules

# Read .csv, convert "." to NA
regeneration <- fread(paste0(input.dir,"/regeneration.csv"), na.string=".")

# Delete entries without regeneration counts
regeneration <- regeneration[!is.na(regeneration$regeneration_count),]

fwrite(regeneration,"GYPSY data/intermediate/i_regeneration.csv",row.names=F)