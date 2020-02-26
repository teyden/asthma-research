context("test-translate-sql-helpers.r")

old <- NULL
setup(old <<- set_current_con(simulate_dbi()))
teardown(set_current_con(old))

test_that("aggregation functions warn if na.rm = FALSE", {
  sql_mean <- sql_aggregate("MEAN")

  expect_warning(sql_mean("x"), "Missing values")
  expect_warning(sql_mean("x", na.rm = TRUE), NA)
})

test_that("missing window functions create a warning", {
  sim_scalar <- sql_translator()
  sim_agg <- sql_translator(`+` = sql_infix("+"))
  sim_win <- sql_translator()

  expect_warning(
    sql_variant(sim_scalar, sim_agg, sim_win),
    "Translator is missing"
  )
})

test_that("missing aggregate functions filled in", {
  sim_scalar <- sql_translator()
  sim_agg <- sql_translator()
  sim_win <- sql_translator(mean = function() {})

  trans <- sql_variant(sim_scalar, sim_agg, sim_win)
  expect_error(trans$aggregate$mean(), "only available in a window")
})

test_that("output of print method for sql_variant is correct", {
  sim_trans <- sql_translator(`+` = sql_infix("+"))
  expect_known_output(
    sql_variant(sim_trans, sim_trans, sim_trans),
    test_path("test-sql-variant.txt"),
    print = TRUE
  )
})

test_that("win_rank() is accepted by the sql_translator", {
  expect_known_output(
    print(sql_variant(
      sql_translator(
        test = win_rank("test")
      )
    )),
    test_path("test-sql-translator.txt"))
})
