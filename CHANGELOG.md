# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-07-25
### Change
- Split the single monolithic `runMavis` task into one task per MAVIS stage (`setup`,
  `convert`, `cluster`, `validate`, `annotate`, `pairing`, `mavisSummary`) so that Cromwell
  tracks, sizes and retries each stage individually. `validate` and `annotate` are fused
  into one task because they are one-to-one per batch.
- MAVIS is no longer run through snakemake; Cromwell owns the dependency graph. `setup`
  computes the per-library `total_batches` that `mavis cluster` requires, and `pairing` and
  `mavisSummary` create their own output directories. Snakemake used to handle both.
- Memory is now requested per stage (16 GB, 32 GB for `validate` and `annotate`) instead of one
  120 GB reservation for the whole pipeline. The per-stage limits are enforced by the
  scheduler, which the Snakefile's own declarations never were.
- `SINGULARITY_TMPDIR` is set to node-local disk in every stage. Singularity cannot
  loop-mount the image on this cluster, and extracting it to the NFS-backed default TMPDIR
  cost about 140 seconds per invocation instead of 3.
- Set the default value of cluster.min_clusters_per_file to 30 instead of 100

Workflow inputs, outputs and vidarr labels are unchanged.

## [1.5.2] - 2026-06-09
### Fix
- Regression: exclude non-deterministic supplementary `_NA` drawings from
  svg/json md5 validation in `calculate.sh`.

### Change
- Regression directory updated with new tag

## [1.5.1] - 2026-06-03
### Change
- Regression directory updated with new tag
- Regression: exclude non-deterministic columns from md5 checksums in calculate.sh to match mavis 3.1.2 columns output.

## [1.5.0] - 2026-05-25
### Change
- Updated MAVIS from 3.1.0 to 3.1.2

### Added
- Added option for reference genome hg38_noAlt

## [1.4.1] - 2026-02-11
### Change
- Updated hg38 annotation file path to ensembl v110 (HG38V110_MAVIS_ROOT)
- Adjusted regression test outputs and workflow commands
- Updated README

## [1.4.0] - 2025-05-09
### Change
- [GRD-931] (https://jira.oicr.on.ca/browse/GRD-931) - Update mavis ensembl json file to nsembl v110

## [Unreleased] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only)

## [1.3.0] - 2024-05-03
### Added
- Add conditional logic to re-attempt mavis with adjusted sample_bin_size if the config is not successfully generated

## [1.2.0] - 2024-04-09
### Added
- Add GRIDSS as an allowable SV input

## [1.1.0] - 2024-02-21
### Added
- Add delly filtering task to handle large delly files

## [1.0.3] - 2023-09-18
### Changed
- Update 'for loop' to append different library types to the "libraries" section of config rather than overwriting them

## [1.0.2] - 2023-09-15
### Changed
- Update 'for loop' to append different SV files to "convert" section of config rather than overwriting them

## [1.0.1] - 2023-08-22
### Added
- Initial Release
  
## [1.0.0] - 2023-07-26
### Added
- Premature tag
