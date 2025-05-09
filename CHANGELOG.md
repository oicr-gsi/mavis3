# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.4.0] - 2024-05-03
### Change
- Update mavis ensembl json file to Ensembl 110

## [Unreleased] - 2024-06-25
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
