test_that("table_to_rdm works with basic similarity table", {
  # Create a simple similarity table
  sim_table <- data.frame(
    label1 = c("A", "A", "B"),
    label2 = c("B", "C", "C"),
    similarity = c(0.8, 0.5, 0.6)
  )
  
  labels <- c("A", "B", "C")
  
  # Test as dist object
  rdm_dist <- table_to_rdm(sim_table, labels)
  expect_s3_class(rdm_dist, "dist")
  expect_equal(attr(rdm_dist, "Size"), 3)
  
  # Test as matrix
  rdm_mat <- table_to_rdm(sim_table, labels, as_dist = FALSE)
  expect_true(is.matrix(rdm_mat))
  expect_equal(dim(rdm_mat), c(3, 3))
  
  # Check values (converting similarity to dissimilarity)
  expect_equal(rdm_mat[1, 2], 1 - 0.8)  # A-B
  expect_equal(rdm_mat[2, 1], 1 - 0.8)  # B-A (symmetric)
  expect_equal(rdm_mat[1, 3], 1 - 0.5)  # A-C
  expect_equal(rdm_mat[2, 3], 1 - 0.6)  # B-C
  
  # Check diagonal
  expect_equal(as.numeric(diag(rdm_mat)), rep(0, 3))  # Self-dissimilarity is 0
})

test_that("table_to_rdm handles missing labels correctly", {
  sim_table <- data.frame(
    label1 = c("A", "B"),
    label2 = c("B", "C"),
    similarity = c(0.7, 0.9)
  )
  
  # Include labels not in the table
  labels <- c("A", "B", "C", "D")
  
  rdm_mat <- table_to_rdm(sim_table, labels, 
                          default_similarity = 0,
                          as_dist = FALSE)
  
  # Check that missing pairs use default similarity
  expect_equal(rdm_mat[1, 4], 1 - 0)  # A-D (not in table, default = 0)
  expect_equal(rdm_mat[3, 4], 1 - 0)  # C-D (not in table)
  expect_equal(rdm_mat[1, 3], 1 - 0)  # A-C (not in table)
  
  # Check specified pairs
  expect_equal(rdm_mat[1, 2], 1 - 0.7)  # A-B
  expect_equal(rdm_mat[2, 3], 1 - 0.9)  # B-C
})

test_that("table_to_rdm handles asymmetric tables", {
  # Create asymmetric similarity table
  sim_table <- data.frame(
    label1 = c("A", "B", "A"),
    label2 = c("B", "A", "C"),
    similarity = c(0.7, 0.3, 0.5)  # A→B ≠ B→A
  )
  
  labels <- c("A", "B", "C")
  
  rdm_mat <- table_to_rdm(sim_table, labels,
                          symmetric = FALSE,
                          as_dist = FALSE)
  
  # Check asymmetric values
  expect_equal(rdm_mat[1, 2], 1 - 0.7)  # A→B
  expect_equal(rdm_mat[2, 1], 1 - 0.3)  # B→A (different!)
  expect_equal(rdm_mat[1, 3], 1 - 0.5)  # A→C
  expect_equal(rdm_mat[3, 1], 1 - 0)    # C→A (default)
})

test_that("table_to_rdm handles duplicate entries", {
  sim_table <- data.frame(
    label1 = c("A", "A", "B"),
    label2 = c("B", "B", "C"),
    similarity = c(0.7, 0.9, 0.5)  # Duplicate A-B pair
  )
  
  labels <- c("A", "B", "C")
  
  expect_warning(
    rdm_mat <- table_to_rdm(sim_table, labels, as_dist = FALSE),
    "Duplicate label pairs"
  )
  
  # Should use first occurrence
  expect_equal(rdm_mat[1, 2], 1 - 0.7)
})

test_that("table_to_rdm handles self-similarity options", {
  sim_table <- data.frame(
    label1 = c("A", "A", "B"),
    label2 = c("A", "B", "C"),
    similarity = c(0.95, 0.7, 0.5)  # Includes self-similarity for A
  )
  
  labels <- c("A", "B", "C")
  
  # Test with default self-similarity
  rdm_mat1 <- table_to_rdm(sim_table, labels, 
                           self_similarity = 1,
                           as_dist = FALSE)
  expect_equal(as.numeric(diag(rdm_mat1)), rep(0, 3))  # All diagonal = 0
  
  # Test with NA self-similarity (look up in table)
  rdm_mat2 <- table_to_rdm(sim_table, labels,
                           self_similarity = NA,
                           as_dist = FALSE)
  expect_equal(rdm_mat2[1, 1], 1 - 0.95)  # A-A from table
  expect_equal(rdm_mat2[2, 2], 0)         # B-B not in table, default
  expect_equal(rdm_mat2[3, 3], 0)         # C-C not in table, default
})

test_that("table_to_rdm validates input", {
  # Not a data frame
  expect_error(
    table_to_rdm(matrix(1:9, 3, 3), c("A", "B", "C")),
    "must be a data frame"
  )
  
  # Missing columns
  bad_table <- data.frame(x = 1:3, y = 4:6)
  expect_error(
    table_to_rdm(bad_table, c("A", "B")),
    "must contain columns"
  )
  
  # Invalid similarity values
  sim_table <- data.frame(
    label1 = c("A", "B"),
    label2 = c("B", "C"),
    similarity = c(1.5, -0.2)  # Out of range
  )
  
  expect_warning(
    table_to_rdm(sim_table, c("A", "B", "C")),
    "should be between 0 and 1"
  )
})

test_that("category_rdm creates correct category-based RDM", {
  # Create category structure
  categories <- c(cat = "animal", dog = "animal", bird = "animal",
                  car = "vehicle", plane = "vehicle")
  
  rdm_mat <- category_rdm(categories, 
                          within_category_sim = 0.8,
                          between_category_sim = 0.2,
                          as_dist = FALSE)
  
  # Check dimensions
  expect_equal(dim(rdm_mat), c(5, 5))
  expect_equal(rownames(rdm_mat), names(categories))
  
  # Check within-category similarity
  expect_equal(rdm_mat["cat", "dog"], 1 - 0.8)    # Same category
  expect_equal(rdm_mat["car", "plane"], 1 - 0.8)   # Same category
  
  # Check between-category similarity
  expect_equal(rdm_mat["cat", "car"], 1 - 0.2)     # Different categories
  expect_equal(rdm_mat["bird", "plane"], 1 - 0.2)  # Different categories
  
  # Check diagonal
  expect_equal(as.numeric(diag(rdm_mat)), rep(0, 5))
})

test_that("category_rdm handles different input formats", {
  # Test with named vector
  categories <- c(cat = "animal", dog = "animal", car = "vehicle", plane = "vehicle")
  
  rdm_mat <- category_rdm(categories, as_dist = FALSE)
  
  # Check it uses names as labels
  expect_equal(rownames(rdm_mat), c("cat", "dog", "car", "plane"))
  
  # Test error with unnamed vector
  expect_error(
    category_rdm(c("animal", "animal", "vehicle")),
    "must be a named vector"
  )
  
  # Test error with factor
  expect_error(
    category_rdm(factor(c("animal", "vehicle"))),
    "For factors"
  )
})

test_that("table_to_rdm works with custom column names", {
  sim_table <- data.frame(
    condition1 = c("A", "B"),
    condition2 = c("B", "C"),
    correlation = c(0.7, 0.5)
  )
  
  labels <- c("A", "B", "C")
  
  rdm_mat <- table_to_rdm(sim_table, labels,
                          label1_col = "condition1",
                          label2_col = "condition2",
                          similarity_col = "correlation",
                          as_dist = FALSE)
  
  expect_equal(rdm_mat[1, 2], 1 - 0.7)
  expect_equal(rdm_mat[2, 3], 1 - 0.5)
})

test_that("table_to_rdm integration with RSA", {
  # This tests that the output format is compatible with RSA functions
  sim_table <- data.frame(
    label1 = c("A", "B", "A"),
    label2 = c("B", "C", "C"),
    similarity = c(0.7, 0.8, 0.5)
  )
  
  labels <- c("A", "B", "C", "D")
  
  # Create RDM as dist object (typical RSA input)
  rdm <- table_to_rdm(sim_table, labels, default_similarity = 0.1)
  
  # Check it's a proper dist object
  expect_s3_class(rdm, "dist")
  expect_equal(attr(rdm, "Size"), 4)
  
  # Check it can be converted to matrix and back
  mat <- as.matrix(rdm)
  expect_equal(dim(mat), c(4, 4))
  expect_true(isSymmetric(mat))
  
  # Check values are in expected range [0, 1]
  expect_true(all(rdm >= 0 & rdm <= 1))
})