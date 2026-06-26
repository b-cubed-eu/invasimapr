# ======================================================================
# 7) ONE-LINER CONVENIENCE (optional)
# ======================================================================

#' @title Run the full invasimapr pipeline
#' @description Convenience wrapper that chains the stages with sensible defaults.
#' @param traits_inv Invader traits (rows = invaders).
#' @param option invasion-fitness option; see \code{predict_establishment()}.
#' @param prob_method NULL or "probit"/"logit"/"hard".
#' @param prob_args list forwarded to \code{compute_establishment_probability()}.
#' @param ... arguments passed to \code{assemble_matrices()} in the very first stage.
#' @return An object of class \code{invasimapr_fit} containing the results of
#'   the full pipeline (inputs, trait space, resident model, sensitivities,
#'   invader predictions, invasion fitness and summaries).
#' @examples
#' \dontrun{
#' fit <- invasimapr(traits_inv = my_invader_traits, long_df = my_long_table)
#' }
#' @export
#' @keywords internal
invasimapr = function(traits_inv,
                      option = "D",
                      prob_method = "probit",
                      prob_args = list(sigma = 1),
                      ...){
  fit = prepare_inputs(...)
  fit = prepare_trait_space(fit, traits_inv = traits_inv)
  fit = model_residents(fit, saturation_mode = "evenness_deficit")
  fit = learn_sensitivities(fit, use_site_random_slopes = TRUE, lrt = TRUE)
  fit = predict_invaders(fit, traits_inv = traits_inv)
  fit = predict_establishment(fit, option = option,
                               calibrate_kappa = TRUE,
                               prob_method = prob_method,
                               prob_args = prob_args)
  fit = summarise_results(fit)
  fit
}
