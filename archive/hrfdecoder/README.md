This directory preserves the retired `hrfdecoder` integration outside the live
package surface.

Contents:
- `R/`: archived adapter code
- `man/`: archived Rd files
- `tests/testthat/`: archived tests
- `vignettes/`: archived tutorial source
- `docs/`: archived notes and comparisons

`archive/` is excluded from package builds via `.Rbuildignore`, so these files
stay in the repository without making `hrfdecoder` part of the shipped package
or its active documentation set.
