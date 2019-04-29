###############################################################################
# Use test_that to test mina functions / methods.
# Methods "norm_tab" in norm_tab.R are tested.
# Methods "fit_tabs" in fit_tabs.R are tested.
###############################################################################

library("mina"); library("testthat")

# read in data
maize <- new("mina", tab = maize_asv, des_tab = maize_des)

test_that("Test that when method is missed, error message is printed", {
              expect_error(norm_tab(maize),
                  "Must specify a `method`. See `? norm_tab_method_list`")
})
