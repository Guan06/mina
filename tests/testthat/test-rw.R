###############################################################################
# Use test_that to test mina functions / methods.
# Functions in add.R are tested, including add_tab and add_des.
# Functions in check.R are tested, including check_mina, check_mina_qu and
# check_mina_de.
###############################################################################

library("mina"); library("testthat")

# add a new object of the class `mina` and add raw data into the object

maize1 <- new ("mina")
test_that("Test that this object does not have @tab defined / added", {
             expect_error(check_mina_qu(maize1),
                          "The @tab of this object does not exist!")
})

asv_file <- system.file("extdata", "maize_asv.rds", package = "mina")
des_file <- system.file("extdata", "maize_des.txt", package = "mina")

maize1 <- add_tab(maize1, asv_file)
test_that("Test that this is an object of class mina with @tab added", {
             expect_true(check_mina_qu(maize1))
})

maize1 <- add_des(maize1, des_file)
test_that("Test that this is an object of class mina with @des_tab added", {
            expect_true(class(maize1@des_tab) == "data.frame")
})

test_that("Test that this object does not have same samples in @tab and
         @des_tab, should output error message", {
         expect_error(check_mina_de(maize1),
                     "The samples in @tab and @des_tab are different!")
})

# add a new object of the class `mina` with @tab and @des_tab

maize2 <- new("mina", tab = maize_asv, des_tab = maize_des)

test_that("Test that this is an object of class mina with @tab and @des_tab and
         they have the same samples", {
         expect_true(class(maize2@tab) == "matrix")
         expect_true(class(maize2@des_tab) == "data.frame")
         expect_true(check_mina(maize2))
})
