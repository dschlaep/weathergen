#' Fit Markov Chain Transition Matrices to State Sequence
#'
#' @param states character array of states as ordered factor
#' @param months numeric array of months for each daily time step
#' @export
#' @return monthly list of transition matrices
#' @examples
#' transitions <- mc_fit(x=sample(c('d', 'w', 'e'), size=720, replace=TRUE, prob=c(0.5, 0.3, 0.2)), months=rep(rep(seq(1, 12), each=30), times=2))
#'
mc_fit <- function(states, months) {
  stopifnot(length(states) == length(months))

  states_next <- dplyr::lead(states)

  transitions <- lapply(seq(1, 12), function(m) {
    idx <- which(months==m)
    trans <- prop.table(table(states[idx], states_next[idx]), 1)
    # Added in case prop.table returns non square matrix
    if (dim(trans)[1] < 3 | dim(trans)[2] < 3) {
      temptrans <- matrix(0,3,3, dimnames = list(c('d', 'e', 'w'),c('d', 'e', 'w')))
      temptrans[,1] <- 1
      if(dim(trans)[1] < dim(trans)[2]) {
        temptrans[rownames(trans),] <- trans
      } else if (dim(trans)[1] > dim(trans)[2]) {
        temptrans[,colnames(trans)] <- trans
      } else {
        temptrans[rownames(trans),colnames(trans)] <- trans
      }
      trans <- temptrans 
    }
    trans
  })

  transitions
}
