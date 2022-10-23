# CRISPR Conjoint Analysis 
setwd("~/Desktop/CRISPR FOOD/CRISPR_R")

# install packages for conjoint 
# install.packages(c("cluster", "survival"))
# install.packages("conjoint")
# install.packages("tidyverse")
# install.packages("radiant")
# install.packages('gtools')
library(tidyverse)
library(tidyselect)
library(conjoint)
library(readxl)
library(radiant)
library(data.table)
library(tidyr)
library('fastDummies')
library(gtools)


################################################################################
# Cleaning the data 

df <- read_csv("AT_raw.csv") # get data from csv 
exp_design <- read.csv(file = "AppleTomato_tomato1_Design.csv", header = TRUE, stringsAsFactors = FALSE)

# explain variable in design_exp: 
#Concept is the choice 4 profile, none option is not included 

df <-rename(df, Version = sys_CBCVersion_tomato1)
df <-rename(df, id = sys_SequentialRespNum)

# drop unrelated variables, now only keep tomato first

#df <- df[-c(2:21,28:33)]
#df <- df[-c(1:5,7:31),]

# reshape data to wide format one id with 6 task

df <- pivot_longer(df,
                   cols = starts_with("tomato1_Random"),
                   names_to = "task",
                   names_prefix = "tomato1_Random",
                   values_to = "choice",
                   values_drop_na = TRUE
)
# cast task from character to numeric
df$task <- as.numeric(df$task)

# head(select(df, id, task),4)

############ Change to standard wide ############

N = 1
#Every N rows after which empty rows should be inserted
after_rows = 4
exp_design <- do.call(rbind, lapply(split(exp_design, ceiling(1:NROW(exp_design)/after_rows)),
                                    function(a) rbind(a, replace(a[1:N,], TRUE, ""))))

insert <- matrix(5,1800)
exp_design$Concept[c(FALSE, FALSE, FALSE, FALSE, TRUE)] <- insert[, 1]

for (row in 2:length(exp_design$Task)){ # 2 so you don't affect column names
  if(exp_design$Task[row] == "") {    # if its empty...
    exp_design$Task[row] = exp_design$Task[row-1] # ...replace with previous row's value
  }
}

for (row in 2:length(exp_design$Version)){ # 2 so you don't affect column names
  if(exp_design$Version[row] == "") {    # if its empty...
    exp_design$Version[row] = exp_design$Version[row-1] # ...replace with previous row's value
  }
}


# merge 

####### ?? why ?? why id 6 is behind ??????????????????

df$id[df$id==6] <- 06

df_merge = merge(df, exp_design, by.x=c('id','task', 'choice'), by.y=c('Version','Task','Concept'),sort = T)
df_merge2<- df_merge[ ,c(1:3,62:65)]

###### create Dummy for each level of attributes for each profile  
# First change the exp_design_ori from long to wide  
exp_design_ori <- read.csv(file = "AppleTomato_tomato1_Design.csv", header = TRUE, stringsAsFactors = FALSE)

exp_wide <- exp_design_ori %>% 
  pivot_wider(names_from = Concept, 
              values_from = c("Att.1...Production.Method", "Att.3...Purchase.Method", "Att.2...Price"))

# merge again  
df_merge3 = merge(df_merge, exp_wide, by.x=c('id','task'), by.y=c('Version','Task'))
df_merge4<- df_merge3[ ,c(1:3,62,66:77)]

# create dummy 
#install.packages('fastDummies')

df_merge4 <- df_merge4 %>% 
  rename(treatment = sys_block_set_1)

df_merge5 <- dummy_cols(df_merge4,
                        select_columns = c( 'Att.1...Production.Method_1', 
                                            'Att.1...Production.Method_2',
                                            'Att.1...Production.Method_3',
                                            'Att.1...Production.Method_4',
                                            'Att.1...Production.Method_5',
                                            'treatment',
                                            
                                            'Att.3...Purchase.Method_1',
                                            'Att.3...Purchase.Method_2',
                                            'Att.3...Purchase.Method_3',
                                            'Att.3...Purchase.Method_4',
                                            'Att.3...Purchase.Method_5'
                                            ))

#'Att.2...Price_1','Att.2...Price_2',
#'Att.2...Price_3', 'Att.2...Price_4',




# rename for profile 1 
df_merge6 <- df_merge5 %>% 
  rename(
    organic_1 = Att.1...Production.Method_1_1, nongmo_1 = Att.1...Production.Method_1_2,
    conventional_1 = Att.1...Production.Method_1_3, CRISPR_1 = Att.1...Production.Method_1_4,
    GMO_1 = Att.1...Production.Method_1_5, price_1 = Att.2...Price_1, 
    instore_1 = Att.3...Purchase.Method_1_1, online_1 = Att.3...Purchase.Method_1_2,
  )

# rename for profile 2

df_merge6 <- df_merge6 %>% 
  rename(
    organic_2 = Att.1...Production.Method_2_1, nongmo_2 = Att.1...Production.Method_2_2,
    conventional_2 = Att.1...Production.Method_2_3, CRISPR_2 = Att.1...Production.Method_2_4,
    GMO_2 = Att.1...Production.Method_2_5, price_2 = Att.2...Price_2, 
    instore_2 = Att.3...Purchase.Method_2_1, online_2 = Att.3...Purchase.Method_2_2,
  )

# rename profile 3&4
df_merge7 <- df_merge6 %>% 
  rename(
    organic_3 = Att.1...Production.Method_3_1, nongmo_3 = Att.1...Production.Method_3_2,
    conventional_3 = Att.1...Production.Method_3_3, CRISPR_3 = Att.1...Production.Method_3_4,
    GMO_3 = Att.1...Production.Method_3_5, price_3 = Att.2...Price_3, 
    instore_3 = Att.3...Purchase.Method_3_1, online_3 = Att.3...Purchase.Method_3_2,
    
    organic_4 = Att.1...Production.Method_4_1, nongmo_4 = Att.1...Production.Method_4_2,
    conventional_4 = Att.1...Production.Method_4_3, CRISPR_4 = Att.1...Production.Method_4_4,
    GMO_4 = Att.1...Production.Method_4_5, price_4 = Att.2...Price_4, 
    instore_4 = Att.3...Purchase.Method_4_1, online_4 = Att.3...Purchase.Method_4_2,
  )

drop <- c( 'Att.1...Production.Method_1', 'Att.1...Production.Method_2',
           'Att.1...Production.Method_3', 'Att.1...Production.Method_4',
           'Att.3...Purchase.Method_1',  'Att.3...Purchase.Method_2',
           'Att.3...Purchase.Method_3', 'Att.3...Purchase.Method_4',
           'Att.2...Price_1','Att.2...Price_2',
           'Att.2...Price_3', 'Att.2...Price_4',
          'task'
           ) 

# 'GMO_1','GMO_2','GMO_3','GMO_4'

df_clean3 = df_merge7[,!(names(df_merge7) %in% drop)]



######select participants
df_clean3 <- df_clean3[df_clean3$id >= 1, ]


# export clean data in data_clean
OutPath<- "~/Desktop/CRISPR FOOD/CRISPR_R/data_clean/df_clean3.csv"

OutTbl <- df_clean3
write.csv(OutTbl, file = OutPath)

#end 


