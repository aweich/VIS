### Limitations

#### Sequencing depth

The pipeline detects insertions based on the available data: While the number of detected insertions may occasionally appear low, this often reflects the nature of the underlying sequencing data, which typically contains a limited number of insertions. Long-read sequencing is well-suited for detecting and precisely locating insertions at the base level. However, if the base-level probability of the insertion sequences is low, the number of detectable insertions may also be restricted.

#### Input file inconsistencies

The workflow strives to follow best practices for all data types whenever possible. However, this may not always work optimally if the input reference FASTAs, BED files, or BAM files do not strictly adhere to these practices, or if they contain unexpected symbols, which can lead to downstream issues following some kind of malfunctioning string split. Not all such cases can be anticipated in advance. If you encounter problems and suspect that input file names might be the cause, please let us know, ideally with a minimal reproducible example. 

#### Base-level accuracy for Split/Join splitmode

Currently, only the `Buffer` splitmode is fully implemented to allow insertions to be detected with base-level accuracy. In contrast, the `Split` and `Join` modes only report the coordinates of the aligned read(s) that contain the insertion(s). This functionality is expected to be improved in future versions.
