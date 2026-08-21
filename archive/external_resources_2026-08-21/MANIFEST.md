# External resource archive -- 2026-08-21

Frozen copies of every version-sensitive external resource used by the
TF-activity / SNF breast cancer subtyping pipeline. Deposited so that the
analysis can be reproduced exactly, including by readers running it after
these resources have been updated at source.

## Why this archive exists

CollecTRI, MSigDB and the Bioconductor annotation packages are curated
resources that change between releases. Code that re-fetches them at run
time will not fail when they change -- it will silently produce different
numbers. These frozen copies, together with the checksums below, remove
that ambiguity.

## Environment

- R version: 4.4.3
- Bioconductor release: 3.20
- Archive created: 2026-08-21
- Full package environment: `sessionInfo.txt`
- Package lockfile: `renv.lock` (repository root)

## Contents

| Resource | Rows | Package | Version | MD5 (tsv) |
|---|---|---|---|---|
| collectri_network | 62411 | OmnipathR | NA | `846a128211ef6f4b` |
| dorothea_hs_regulons | 454504 | dorothea | 1.18.0 | `7ee081a85bddda8d` |
| msigdb_hallmark | 7331 | msigdbr | 26.1.0 | `02c543a205addc63` |
| msigdb_reactome | 113653 | msigdbr | 26.1.0 | `3d6414295bec86ba` |
| ensembl_symbol_map_used | 20303 | org.Hs.eg.db | 3.20.0 | `6246eeb6ce0df7ca` |

Each resource is provided twice: as `.rds` (R-native, fast to reload) and
as `.tsv.gz` (plain text, readable without R). Full MD5 checksums for both
formats are in `manifest.csv`.

## Verifying a file

```r
tools::md5sum("collectri_network.tsv.gz")
# compare against manifest.csv
```

## Using these instead of live downloads

Place `collectri_network.rds` at `data/processed/network_collectri.rds`
before running `05_tf_activity.R`. The script uses the cache when present
and only downloads when it is absent, so the archived version will be used
and no network access will occur.

## Primary data (not redistributed here)

- **TCGA-BRCA**: downloaded via TCGAbiolinks from the GDC. See
  `logs/01_download_tcga.log` for the access date and
  `logs/sessioninfo/01_download_tcga_package_versions.csv` for the
  TCGAbiolinks version used.
- **METABRIC**: obtained from cBioPortal (`brca_metabric`). Redistribution
  is restricted by the original data use terms; the access date is recorded
  in `logs/13_external_validation.log`.

## Citation

If you use this archive, please cite the project report and the original
resources (CollecTRI, DoRothEA, MSigDB) in their own right.
