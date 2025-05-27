# 📂 R Utility Scripts

This folder contains reusable R functions developed for the manuscript:

**Raghavan et al., _"Multiple mechanisms suppress the initiation of meiotic recombination near chromosome ends"_**  
bioRxiv 2025.02.27.640173 — https://doi.org/10.1101/2025.02.27.640173

Each script defines a core function used in signal normalization, telomeric profiling, bootstrapped genome-wide signal estimation, or genomic region annotation.

---

## 📋 Contents

| File                        | Function                         | Description                                                   |
|-----------------------------|-----------------------------------|---------------------------------------------------------------|
| `normalize_signal.R`        | `gendiv()`                        | Normalizes signal values to the genome-wide average           |
| `telo_signal.R`             | `teloSeqSignal()`                 | Extracts and smooths signal near telomeres                    |
| `bootstrap_distribution.R`  | `AxisDistributionBootstrapped()`  | Performs genome-wide signal bootstrapping                    |
| `annotate_ranges.R`         | `annotate_ranges()`               | Labels bins as telomeric, centromeric, or other              |

---

## 🔧 Usage

These functions are commonly sourced in figure-generating scripts. To load them:

```r
# Load shared utility functions
source("R/normalize_signal.R")
source("R/telo_signal.R")
source("R/bootstrap_distribution.R")
source("R/annotate_ranges.R")
