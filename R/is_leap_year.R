#' Check if a year is a leap year
#' Since this package uses lubridate, consider using lubridate::leap_year instead
#' @param year to check
#' @return TRUE in case year is a leap year
is_leap_year=function(year){
  #http://en.wikipedia.org/wiki/Leap_year
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}
