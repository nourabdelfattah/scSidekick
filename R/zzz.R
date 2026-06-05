# =============================================================================
# scSidekick — package startup message  (zzz.R)
# Named zzz.R by convention so it is sourced last during package build.
# =============================================================================

.onAttach <- function(libname, pkgname) {
  ver <- tryCatch(
    as.character(utils::packageVersion("scSidekick")),
    error = function(e) "?"
  )

  packageStartupMessage(
    "\n",
    "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
    "                                                         \n",
    "             ___   _      _         _  __  _        _   \n",
    "   ___  __  / __| (_)  __| |  ___  | |/ / (_)  __  | |__\n",
    "  (_-< / _| \\__ \\ | | / _` | / -_) | ' <  | | / _| | / /\n",
    "  /__/ \\__| |___/ |_| \\__,_| \\___| |_|\\_\\ |_| \\__| |_\\_\\\n",
    "                                                         \n",
    "   v", ver, "  ·  Your New Best Friend in Visualization\n",
    "                                                         \n",
    "   ✨  It's a good day to make pretty figures!  ✨\n",
    "  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
  )
}
