Changelog
==========

<!--

Newest changes should be on top.

This document is user facing. Please word the changes in such a way
that users understand how the changes affect the new version.
-->

version 0.1.0-dev
---------------------------
+ Create a long-read variant calling pipeline with:
  + Mapping with either minimap2 or pbmm2
  + Variant calling with clair3 and/or DeepVariant (clair3 is enabled by 
    default, both are optional.)
  + Methylation analysis with modkit (optional).
  + Phasing with whatshap (optional)
  + Statistics:
    + Sequali
    + Mosdepth
    + bcftools stats (when clair3 or deepvariant are run)

