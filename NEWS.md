# scatools (development version)

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
