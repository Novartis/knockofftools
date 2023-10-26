#' Heatmap of multiple variable selections ordered by importance
#'
#' @param S data.frame of variable selections from multiple knockoffs (each entry is either 1 if variable is selected and 0 otherwise). Columns correspond to different knockoffs and rows correspond to the underlying variables. row.names(S) records the variable names.
#' @param nbcocluster bivariate vector c(number of variable clusters, number of selection clusters).
#' The former number must be specified less than nrow(S) and the latter must be less than ncol(S).
#'
#' @details To help visualize most important variables we perform clustering both selections and variables.
#'
#' @return plot of heatmap
#' @export
#'
#' @examples
#' library(knockofftools)
#'
#' set.seed(1)
#'
#' # Simulate 8 Gaussian covariate predictors and 2 binary factors:
#' X <- generate_X(n=100, p=10, p_b=2, cov_type="cov_equi", rho=0.2)
#'
#' # create linear predictor with first 5 beta-coefficients = 1 (all other zero)
#' lp <- generate_lp(X, p_nn = 5, a=1)
#'
#' # Gaussian
#'
#' # Simulate response from a linear model y = lp + epsilon, where epsilon ~ N(0,1):
#' y <- lp + rnorm(100)
#'
#' # Calculate M independent knockoff feature statistics:
#' W <- knockoff.statistics(y=y, X=X, type="regression", M=5)
#'
#' S = variable.selections(W, error.type = "pfer", level = 1)
#'
#' # plot heatmap of knockoff selections:
#' plot(S)
plot.variable.selections <- function(S, nbcocluster=c(7,7)) {

  if (class(S)[1]!="variable.selections") {
    stop("Input S must be of class \'variable.selections\'. Please see ?variable.selections.")
  }

  stable.vars <- S$stable.variables

  S <- S$selected

  # Make sure user doesn't specify more clusters than data to support:
  nbcocluster <- pmin(nbcocluster, dim(S))

  selections <- data.frame(draw = factor(rep(rep(1:ncol(S)),each=nrow(S))),
                           variable = factor(rownames(S)),
                           selected = as.numeric(as.matrix(S)))

  `%>%` <- dplyr::`%>%`

  sel.mat <- matrix(selections$selected,nrow=nrow(S))
  hclust.row <- hclust(dist(sel.mat, method="binary"), method="ward.D")
  hclust.col <- hclust(dist(t(sel.mat), method="binary"), method="ward.D")

  selections$varclass = factor(cutree(hclust.row, k=nbcocluster[1]), ordered=TRUE)
  selections$drawclass = factor(rep(cutree(hclust.col, k=nbcocluster[2]), each=nrow(S)), ordered=TRUE)

  # Calculate means per draw cluster (to order heatmap)
  meta.draw <- selections %>%
    dplyr::group_by(drawclass) %>%
    dplyr::summarise(selected = mean(selected)) %>%
    dplyr::arrange(desc(selected))

  # Calculate means per variable block (to order heatmap)
  meta.var <- selections %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(selected = mean(selected))%>%
    dplyr::arrange(selected)

  # Order selections according to block cluster means:
  selections <- selections %>%
    dplyr::mutate(drawclass = factor(drawclass, levels = as.character(meta.draw$drawclass)),
                  variable = factor(variable, levels = as.character(meta.var$variable))) %>%
    dplyr::arrange(drawclass) %>%
    dplyr::mutate(draw = factor(draw, levels = unique(draw))) %>%
    dplyr::arrange(variable)

  selections$selected <- factor(selections$selected)

  # vectorized input to ggplot2::element_text() will likely be deprecated:
  #axis.text.y.color <- ifelse(unique(selections$variable) %in% stable.vars, "red", "black")

  ggplot2::ggplot(data=selections, mapping=ggplot2::aes(x = draw, y = variable)) +
    ggplot2::geom_tile(ggplot2::aes(fill = selected)) +
    ggplot2::xlab("knockoff repetition") +
    ggplot2::ylab("variable (stable selections colored in red)") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   #axis.text.y =  ggplot2::element_text(colour = axis.text.y.color),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values=c("0"="#132B43","1"="#56B1F7"))

}
