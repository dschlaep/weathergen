#' Adjust Markov Chain Transition Matrix
#'
#' @param transition markov transition matrix
#' @param dry_spell change factor for dry spell (default=1, no change)
#' @param wet_spell change factor for wet spell (default=1, no change)
#' @export
#' @examples
#' transitions <- mc_fit(x=sample(c('d', 'w', 'e'), size=720, replace=TRUE, prob=c(0.5, 0.3, 0.2)), months=rep(rep(seq(1, 12), each=30), times=2))
#' mc_adjust_transition(transition[[1]], dry_spell=1, wet_spell=1)
#'
mc_adjust_transition <- function(transition, dry_spell=1, wet_spell=1) {
  states <- rownames(transition)

  p <- transition
  p_new <- p

  # adjust dry spell
  # states changed order somewhere during execution, considering a change to
  p_dw <- (p['d', 'w'] + p['d', 'e'])/dry_spell -  p['d', 'e']
  p_dd <-  p['d', 'd'] + (p['d', 'w'] - p_dw)
  p_new['d','w'] <- p_dw
  p_new['d','d'] <- p_dd
  # p_dw <- (p[states[1], states[3]] + p[states[1], states[2]])/dry_spell - p[states[1], states[2]]
  # p_dd <- p[states[1], states[1]] + (p[states[1], states[3]] - p_dw)
  # p_new[states[1],states[3]] <- p_dw
  # p_new[states[1],states[1]] <- p_dd

  # adjust wet spell
  # as above considering change to
  p_wd <- (p['w', 'd'] + p['w', 'e'])/wet_spell - p['w', 'e']
  p_ww <- p['w', 'w'] + (p['w', 'd'] - p_wd)
  p_new['w','d'] <- p_wd
  p_new['w','w'] <- p_ww
  # p_wd <- (p[states[3], states[1]] + p[states[3], states[2]])/wet_spell - p[states[3], states[2]]
  # p_ww <- p[states[3], states[3]] + (p[states[3], states[1]] - p_wd)
  # p_new[states[3],states[1]] <- p_wd
  # p_new[states[3],states[3]] <- p_ww
  
  stopifnot(all(p_new >= 0))

  p_new
}
