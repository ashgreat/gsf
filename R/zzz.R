#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "gsf: Generalized Stochastic Frontier Model (version ",
    utils::packageVersion("gsf"), ")\n",
    "Based on Kumbhakar and Tsionas (2021)\n",
    "Type ?gsf for help, and citation('gsf') for citation info."
  )
}

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "TI", "firm", "input", "lower", "type", "upper", "value", "variable"
  ))
}
