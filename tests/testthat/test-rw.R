###############################################################################
# Use test_that to test mina functions / methods.
# Functions in check.R are tested, including check_mina, check_mina_qu and
# check_mina_de.
###############################################################################

library("mina"); library("testthat")

# add a new object of the class `mina` and add raw data into the object

maize1 <- new("mina", tab = maize_asv, des = maize_des)

#test_that("Test that this object does not have @tab defined / added", {
#             expect_error(check_mina_qu(maize1),
#                          "The @tab of this object does not exist!")
#})


test_that("Test that this is an object of class mina with @tab added", {
             expect_true(check_mina_qu(maize1))
})

test_that("Test that this is an object of class mina with @des added", {
            expect_true(class(maize1@des) == "data.frame")
})

test_that("Test that this object has the same samples in @tab and @des", {
            expect_true(check_mina_de(maize1))
})

maize1@des <- maize1@des[1:100, ]
test_that("Test that this object does not have same samples in @tab and
         @des, should output error message", {
         expect_error(check_mina_de(maize1),
                     "The samples in @tab and @des are different!")
})

# add a new object of the class `mina` with @tab and @des

maize2 <- new("mina", tab = maize_asv, des = maize_des)

test_that("Test that this is an object of class mina with @tab and @des and
         they have the same samples", {
         expect_true(class(maize2@tab) == "matrix")
         expect_true(class(maize2@des) == "data.frame")
         expect_true(check_mina(maize2))
})
