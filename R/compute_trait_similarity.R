#' Compute Trait Similarity for Numeric and Categorical Variables
#'
#' This function calculates within-trait similarity as a percentage for each column
#' in a data frame. Numeric traits use scaled mean pairwise similarity; categorical
#' traits use the proportion of identical pairs.
#'
#' @param df A data.frame or tibble where each column is a trait vector (numeric or factor/character).
#' @return A tibble with two columns:
#'   \describe{
#'     \item{Trait}{The original column name}
#'     \item{Similarity}{Percentage similarity (0-100) for that trait}
#'   }
#' @examples
#' df <- data.frame(
#'   height = c(10, 15, 15, 20),
#'   color  = c("red", "blue", "red", "red")
#' )
#' compute_trait_similarity(df)
#' @export
compute_trait_similarity <- function(df) {
  # ensure required packages are available
  if (!requireNamespace("purrr", quietly = TRUE)) stop("purrr required")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("tibble required")

  #---- Helper: numeric similarity ----
  sim_numeric <- function(x) {
    # remove NA values
    x2 <- x[!is.na(x)]
    # if fewer than 2 values, return perfect similarity
    if (length(x2) <= 1) {
      return(1)
    }
    # compute range
    r <- diff(range(x2))
    # if no variation, perfect similarity
    if (r == 0) {
      return(1)
    }
    # scale to [0,1]
    x_scaled <- (x2 - min(x2)) / r
    # mean pairwise distance → convert to similarity
    1 - mean(dist(x_scaled))
  }

  #---- Helper: categorical similarity ----
  sim_categ <- function(x) {
    x <- as.character(x)
    freqs <- table(x)
    # count same-level pairs over total pairs
    sum(choose(freqs, 2)) / choose(sum(freqs), 2)
  }

  #---- Apply similarity by column ----
  trait_sim <- purrr::map_dbl(
    df,
    ~ {
      col <- .x
      if (is.numeric(col)) {
        sim_numeric(col)
      } else {
        sim_categ(col)
      }
    }
  )

  #---- Convert to percentage and return tibble ----
  tibble::enframe(
    trait_sim * 100,
    name  = "Trait",
    value = "Similarity"
  )
}
