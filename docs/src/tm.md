# Melting Temperature Calculation

`SeqFold.jl` provides functions for calculating DNA melting temperatures using nearest-neighbor thermodynamics.

## Basic Usage

```@docs
tm
SeqFold.tm_cache
MeltingConditions
SeqFold.gc_cache
```

## Implementation Details

The calculation uses nearest-neighbor thermodynamic parameters from related literature, accounting for initialization terms, nearest-neighbor pairs, and terminal mismatches when present.

!!! note
    For implementation details, see the source code documentation of the internal constants 
    `DNA_NN`, `DNA_INTERNAL_MM`, and `DNA_TERMINAL_MM`.