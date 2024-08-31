library(dplyr)
library(tidyr)

# Function to check if one EC number is less complete than another
is_less_complete <- function(ec1, ec2) {
  parts1 <- strsplit(ec1, "\\.")[[1]]
  parts2 <- strsplit(ec2, "\\.")[[1]]
  for (i in seq_along(parts1)) {
    if (parts1[i] == "-" && parts2[i] != "-") {
      return(TRUE)
    } else if (parts1[i] != parts2[i]) {
      return(FALSE)
    }
  }
  return(FALSE)
}

# Function to filter less complete EC numbers within each group
filter_less_complete <- function(group) {
  ec_numbers <- group$EC_numbers
  to_remove <- logical(length(ec_numbers))
  
  for (i in seq_along(ec_numbers)) {
    for (j in seq_along(ec_numbers)) {
      if (i != j && is_less_complete(ec_numbers[i], ec_numbers[j])) {
        to_remove[i] <- TRUE
      }
    }
  }
  
  group[!to_remove, ]
}

# Integrate EC numbers
EC_number <- annotation_results_aggregated |>
  select(locus_tag_SC, starts_with("EC"))

aaa <- EC_number |>
  pivot_longer(cols = starts_with("EC"), names_to = "EC_column", values_to = "EC_numbers") |>
  separate_rows(EC_numbers, sep = "; ") |>
  filter(!is.na(EC_numbers)) # Optional: Remove rows with NA EC numbers

# Group by locus_tag and EC_column, then filter less complete EC numbers
filtered_aaa <- aaa |>
  group_by(locus_tag_SC, EC_column) |>
  group_modify(~ filter_less_complete(.x)) |>
  distinct() |>
  ungroup()

# View the result
print(filtered_aaa)