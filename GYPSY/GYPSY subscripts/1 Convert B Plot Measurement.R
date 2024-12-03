## Import .csv, convert PGYI specs, drop unwanted fields added by Tesera, and apply "data smoothing" rules
# Read .csv, convert "." to NA
plot_measurement <- fread(paste0(input.dir,"/plot_measurement.csv"), na.string=".") %>%
  distinct

fwrite(plot_measurement,"GYPSY data/intermediate/i_plot_measurement.csv",row.names=F)