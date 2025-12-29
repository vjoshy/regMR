.onAttach <- function(libname, pkgname){
  packageStartupMessage(
    "\n",
    "── regMR ──────────────────────────────\n",
    "Version: ", utils::packageVersion(pkgname), "\n",
    "Type ?regMR for help\n",
    "────────────────────────────────────────"
  )
}
