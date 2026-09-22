test_that("src/RcppExports.cpp is synchronized with Stan models and Rcpp attributes", {
  skip_if_not_installed("Rcpp")

  pkg_dir <- tryCatch(
    rprojroot::find_package_root_file(),
    error = function(e) "."
  )
  cpp_file <- file.path(pkg_dir, "src", "RcppExports.cpp")
  skip_if(
    !file.exists(cpp_file),
    "src/RcppExports.cpp not found (installed package)"
  )

  before <- readLines(cpp_file)
  Rcpp::compileAttributes(pkg_dir)
  after <- readLines(cpp_file)

  expect_identical(
    before,
    after,
    info = paste(
      "src/RcppExports.cpp is out of date.",
      "Run 'just docs' or 'Rcpp::compileAttributes()' and commit the changes."
    )
  )
})
