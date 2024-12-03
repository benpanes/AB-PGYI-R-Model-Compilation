## Import .csv, convert to correct PGYI specs, drop unwanted fields added by Tesera, and apply 'data smoothing' rules
# Read .csv, convert "." to NA
plot <- fread(paste0(input.dir,"/plot.csv"), na.string=".") %>%
  distinct

# Remove suspected zero filling
plot[elevation==0,elevation:=NA]
plot[sampling_unit_number==0,sampling_unit_number:=""]
plot[aspect=="" & slope==0,aspect:=NA]

# Fill in missing values where possible
plot[ecosite_phase!="" & ecosite=="",ecosite:=substr(ecosite_phase,1,1)]

# Check NSR/ecosite combinations and toss if not eligible
plot[natural_subregion%in%c("ALP","DMG","FF","MG","NF","CP","FP","PRP") &
     ecosite!="",
     ecosite:=""]
plot[natural_subregion%in%c("CM","DM","NM","PAD") &
     ecosite%in%c("m","n"),
     ecosite:=""]
plot[natural_subregion%in%c("LBH","UBH") & 
     ecosite%in%c("k","l","m","n"),
     ecosite:=""]
plot[natural_subregion%in%c("BSA","AP","KU") & 
     ecosite%in%c("i","j","k","l","m","n"),
     ecosite:=""]
plot[natural_subregion%in%c("LF","UF") &
     ecosite%in%c("l","m","n") &
     ecosite_guide=="SW",
     ecosite:=""]
plot[natural_subregion=="UF" &
     ecosite=="n" &
     ecosite_guide=="WC",
     eccosite:=""]
plot[natural_subregion=="MT" &
     ecosite%in%c("h","i","j","k","l","m","n") &
     ecosite_guide=="SW",
     ecosite:=""]
plot[natural_subregion=="MT" &
     ecosite%in%c("i","j","k","l","m","n") &
     ecosite_guide=="WC",
     ecosite:=""]
plot[natural_subregion=="SA" &
     ecosite%in%c("i","j","k","l","m","n") &
     ecosite_guide=="SW",
     ecosite:=""]
plot[natural_subregion=="SA" &
     ecosite%in%c("j","k","l","m","n") &
     ecosite_guide=="WC",
     ecosite:=""]

# Check NSR/ecosite guide combinations and adjust if not eligible - fore example if WC or SW and BM, reset to N
plot[natural_subregion%in%c("LF","UF","MT","SA") & 
     ecosite_guide=="N",
     ecosite_guide:=""]
plot[!(natural_subregion%in%c("LF","UF","MT","SA")) & 
     ecosite_guide%in%c("SW","WC"),
     ecosite_guide:="N"]

fwrite(plot,"GYPSY data/intermediate/i_plot.csv")
