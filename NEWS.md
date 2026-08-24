# scatools (development version)

* Fixed an unbalanced `{}` expression in the `identify_normal()` warning about
  `n_normal_clusts`. cli could not parse it, so that branch raised
  `Could not parse cli {} expression` instead of warning and clamping the value.
* Debug-level diagnostics in `scale_mat()` and `read_vartrix()` are hidden again
  unless `options(scatools.debug = TRUE)` is set. They sat below `logger`'s
  default threshold before the cli migration and had become unconditional.

# scatools 0.1.2

* Replace `add_gc_cor`, `segment_cnv`, and `merge_segments` parallel backend to use `BiocParallel`
* Added `segment=FALSE` as default
* `get_label_centers` works with `Milo` objects
* Add TCGA sample vignette
* Minor bugfixes


# scatools 0.1.1

* Fixed H5AD conversion
* Cleanup build warnings
* Reworked `bin_atac_frags` to efficiently utilize `data.table` for loading

# SCAtools 0.1

* Initial release
