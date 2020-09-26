# Time comparison of psych::alpha() with Lanyu's CreateRdrop()
# Gabriel Odom
# 2019-05-07

system.time(
  a <- CreateRdrop(data = data_df, method = "pearson")
)
a[, 2]
# Lanyu for() loop: 5.647 seconds for 942 x 723

system.time(
  a <- CreateRdrop(data = data_df, method = "pearson")
)
# Gabriel lapply + tibble: 5.831 for 942 x 723
a[, 2]

system.time(
  a <- CreateRdrop(data = data_df, method = "pearson")
)
# Gabriel lapply + data.frame: 5.536 for 942 x 723
a[, 2]

system.time(
  b <- psych::alpha(data_df, warnings = FALSE)
)
# 1184.579 seconds for 942 x 723
b$item.stats$r.drop
