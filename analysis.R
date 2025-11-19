library(dplyr)
library(openxlsx)

# Load the dataset
data <- read.xlsx("data/Assignment Data Nov 2025.xlsx")

# remove the rows with missing values
df <- na.omit(data)

