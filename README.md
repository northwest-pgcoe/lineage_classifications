# lineage_classifications

The lineage_classifications.csv file is maintained by the WA DOH Molecular Epidemiology Program in order to group the lineages for the Sequencing & Variants Report.

Variables include:

| Variable  | Definition | Source |
| ------------- | ------------- | ------------- |
| lineage_extracted  | SC2 lineage | extracted from lineage_notes.txt from cov-lineages  |
| description  | description of lineage_extracted | extracted from lineage_notes.txt from cov-lineages  |
| status  | "Active" or "Withdrawn" | created from "description" column from lineage_notes.txt  |
| ~~cdc_class~~  | ~~Indicates CDC lineage class~~ | ~~created in R code~~  |(Commented out of script 9/17/2024 due to CDC no longer maintaining this information)
| who_class  | Indicates WHO lineage class | created in R code  |
| doh_variant_name  | Indicates grouping of lineages into parent group used in the DOH Sequencing & Variants Report | created in R code  |
| hex_code  | Hex code for color of doh_variant_name; follows CDC Nowast colors | created in R code  |
| lineage_reporting_group  | 1: currently monitored variants; 2: formerly monitored variants; 3: formerly circulating variants, not monitored | created in R code  |
| doh_variant_name_tables  | Indicates variant name in numerical/pango form for report tables | created in R code  |
