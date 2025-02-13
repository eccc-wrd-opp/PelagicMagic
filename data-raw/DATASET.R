## code to prepare `DATASET` dataset goes here

cms_table <- read.csv('data-raw/cms.csv', na.strings = '')
cms_table$min_date <- as.Date(cms_table$min_date, format = c('%m/%d/%Y'))
cms_table$max_date <- as.Date(cms_table$max_date, format = c('%m/%d/%Y'))

usethis::use_data(cms_table, overwrite = TRUE)
usethis::use_data(cms_table, overwrite = TRUE, internal = TRUE)
