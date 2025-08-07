getwd()
setwd("C:\Users\aanis\OneDrive\Documents\general\afu_1\term_2_afu\data_1\AI_Omics_Internship_2025\Module_I.R")

# Create common subfolders for project organization
dir.create("raw_data")     # for storing raw files
dir.create("clean_data")     # for storing cleaned data files
dir.create("scripts")   # for saving R scripts
dir.create("results")  # for saving analysis outputs
dir.create("plots")  # for saving plots
dir.create("tasks")  # for saving tasks

data <- read.csv(file.choose())
View(data)

str(data)

# convert character to factor for categorical variables

data$gender <- as.factor(data$gender)
data$diagnosis <- as.factor(data$diagnosis)
data$smoker <- as.factor(data$smoker)

str(data)

# create binary variable for smokers
data$smoking_status <- ifelse(data$smoker == "Yes", 1, 0)
class(data$smoking_status)

data=data[,-7]

# save cleaned data as csv file in clean data folder
write.csv(data, "clean_data/patient_info_clean.csv")

# save the entire R workspace
save.image(file = "full_workspace_1.RData")

save.image(file = "AmnaMAnis_Class_Ib_Assignment")
