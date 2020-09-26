# Lint the RNAeditr Package
# Gabriel Odom
# 2020-09-21

# We have completed the BiocCheck of the package, and we found this warning:
#   " Use TRUE/FALSE instead of T/F"
# For the life of me, I cannot find where this WARNING is coming from. If it's
#   real, the lintr:: package should find it.

# install.packages("lintr")
library(tidyverse)

rnaEditr_lint <-
  lintr::lint_package() %>%
  as.data.frame %>%
  group_by(linter)

linterGroups_char <- 
  rnaEditr_lint %>% 
  pull(linter) %>% 
  unique()

map(
  .x = linterGroups_char,
  .f = ~{
    rnaEditr_lint %>% 
      filter(linter == .x) %>% 
      pull(message) %>% 
      unique()
  }
)

