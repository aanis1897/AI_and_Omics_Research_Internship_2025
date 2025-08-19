gene <- "TP53"

# Explanation:
# gene is the variable name
# <- is an assignment operator that stores the gene name "TP53" into variable
# "TP53" is now saved in variable gene

# Retrieve the value in the console
 # or use print () function


# Store numeric values in one variable (vector)
expression_levels <- c(2.3, 4.6, 3.6, 7.2, 4.7)

# Import a CSV file as a variable
raw_data <- read.csv(file.choose())
raw_data$patient_id <- NULL

data <- raw_data
clean_data <- data[, -1]
# Sort 'age' from largest to smallest
sorted_age <- sort(raw_data$age, decreasing = TRUE)
sorted_age

# Sort 'age' from smallest to largest
sorted_age2 <- sort(raw_data$age, decreasing = FALSE)
sorted_age2

gene_expression <- 30

# if TRUE, prints the message "Gene is highly expressed"
if (gene_expression > 30) {
  print("Gene is highly expressed")
}

if (gene_expression > 50) {
  print("Gene is highly expressed")
} else {
  print("Gene expression is low")
}

# check structure
str(raw_data)

clean_data1 <- raw_data
str(clean_data1)

# convert all character columns to factors
# loop through each column of clean_data using its index (i)

for (i in 1:ncol(clean_data1)) {
  if (is.character(clean_data1[[i]])) {
    clean_data1[[i]] <- as.factor(clean_data1[[i]])
  }
}

str(clean_data1)

# practice exercises 
# 1. check cholesterol level (using if) 
cholesterol <- 230
if (cholesterol > 240) {
  print("high cholesterol")
} 


# 2. blood pressure status (using if else)

systolic_bp <- 131

if (systolic_bp < 120) {
  print("blood pressure is normal")
} else {
print("blood pressure is high")
}

# 3. Automating Data Type Conversion with for loop

# use patient_info.csv data and metadata.csv

# patient_info.csv
raw_data2 <- read.csv(file.choose())

clean_data3 <- raw_data2
str(clean_data3)

factor_cols <- c("gender", "diagnosis", "smoker" )


for (col in factor_cols) {
  clean_data3[[col]] <- as.factor(clean_data3[[col]]) 
}
str(clean_data3)


# metadata.csv
raw_data_meta <- read.csv(file.choose())
clean_data_meta2 <- raw_data_meta
str(clean_data_meta2)

factor_cols_m <- c("gender", "height" )
for (col in factor_cols_m) {
  clean_data_meta2[[col]] <- as.factor(clean_data_meta2[[col]]) 
}
str(clean_data_meta2)


# 4. Converting Factors to Numeric Codes

# create binary variable for gender
binary_cols <- c("gender")   # store column names in a vector


for (col in binary_cols) {
  clean_data_meta[[col]] <-ifelse (clean_data_meta[[col]] == "Female", 1, 0)
}
str(clean_data_meta)


clean_data_meta3 <- raw_data_meta

factor_cols_m3 <- c("gender", "height" )
for (col in factor_cols_m3) {
  clean_data_meta3[[col]] <- as.factor(clean_data_meta3[[col]])
}
# create binary variable for height (0,1,2)
library(dplyr)
clean_data_meta3$height_num <- recode(clean_data_meta3$height,
                                      "Tall" = 0,
                                      "Medium" = 1,
                                      "Short" = 2)

class(clean_data_meta3$height_num)

str(clean_data_meta2)

clean_data_meta3$height_num <- as.integer(clean_data_meta3$height_num)
str(clean_data_meta3)
str(clean_data_meta)