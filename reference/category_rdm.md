# Create Hypothesis RDM from Category Structure

Creates an RDM based on category membership, where items in the same
category have high similarity and items in different categories have low
similarity.

## Usage

``` r
category_rdm(
  categories,
  within_category_sim = 0.8,
  between_category_sim = 0.2,
  as_dist = TRUE
)
```

## Arguments

- categories:

  A named vector or factor where names/levels are labels and values are
  category assignments

- within_category_sim:

  Similarity value for items within the same category (default: 0.8)

- between_category_sim:

  Similarity value for items in different categories (default: 0.2)

- as_dist:

  Logical; if TRUE return a dist object, otherwise return matrix
  (default: TRUE)

## Value

A dist object or matrix representing the category-based RDM

## Examples

``` r
# Create category structure
categories <- c(cat = "animal", dog = "animal", bird = "animal",
               car = "vehicle", plane = "vehicle", boat = "vehicle")

# Create category-based RDM
rdm <- category_rdm(categories)

# Custom similarity values
rdm <- category_rdm(categories, 
                   within_category_sim = 0.9,
                   between_category_sim = 0.1)
```
