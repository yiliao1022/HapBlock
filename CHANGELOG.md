# Changelog

## 17.0.0 — 2026-10-20

- Reframed block detection around the cultivar population; wild accessions are optional and used only for post-call projection.
- Replaced connected-component grouping with complete-linkage clustering and direct maximum within-cluster distance validation, preventing A–B–C chain merging.
- Added common-SNP coarse parent blocks, fine child blocks, 10 kb boundary scanning, and explicit short-gap annotations.
- Added recent subhaplotypes within each backbone haplotype using group-specific polymorphic variants.
- Removed the fixed maximum haplotype count and added association-ready filtering based on minimum 10% frequency per major haplotype and 50% combined assignment.
- Added one-pass compact genotype caching, cache fingerprint validation, sample call-rate/heterozygosity QC, PCA QC, wild projection, reports, and output indexes.
- Added cultivar-only, cultivar+wild, chain-merging, short-gap, and sustained-boundary regression tests.
