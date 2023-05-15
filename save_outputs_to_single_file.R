# COLLATE ALL SCENARIO OUTPUTS INTO SINGLE CSV 

library(moments)

setwd("/home/alison/Dropbox/03_IPAC/P1_paper_scenarios")

path = "/home/alison/Dropbox/03_IPAC/P1_paper_scenarios/csv_output/summary"

df <- list.files(path, full.names = TRUE) %>%
  map_dfr(read_csv)

write.csv(df, "all_output2.csv")



