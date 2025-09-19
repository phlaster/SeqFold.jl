# Sequence Folding

[`SeqFold.jl`](https://github.com/phlaster/SeqFold.jl) implements the Zuker and Stiegler (1981) dynamic programming algorithm for predicting nucleic acid secondary structures.

## Showcase

```julia
julia> seq = "GGGAGGTCAGCAAACCTGAACCTGTTGAGATGTTGACGTCAGGAAACCCT";

julia> # Get the list of possible secondary structures for a sequence

julia> folded = fold(seq) # calculated at 37°C by default, Julia uses 1-based indexing
12-element Vector{SeqFold.Structure}:
    3   40   -1.3  STACK:GA/CT    
    4   39   -1.4  STACK:AGG/TGC  
    6   37   -1.5  STACK:GT/CA    
    7   36   -1.3  STACK:TC/AG    
    8   35   -1.5  STACK:CA/GT    
    9   34   -2.0  STACK:AGC/TTG  
   11   32   -1.5  STACK:CA/GT    
   12   31    2.8  INTERIOR_LOOP:4/2
   16   29   -1.3  STACK:CT/GA    
   17   28   -0.1  STACK:TGA/AGT  
   19   26   -1.0  STACK:AA/TT    
   20   25    3.5  HAIRPIN:AC/TG 

julia> dot_bracket(seq, folded) # get the dot-bracket notation
"..((.((((.((...((.((....)).)).)).)))).)).........."

julia> folded_hot = fold(seq, temp=80)
4-element Vector{SeqFold.Structure}:
    7   19   -0.4  STACK:TC/AG    
    8   18   -0.5  STACK:CA/GT    
    9   17   -0.4  STACK:AG/TC    
   10   16    3.6  HAIRPIN:GC/CC 

julia> dot_bracket(seq, folded_hot) # at a higher temperature number of secondary structures is decreased
"......((((.....))))..............................."

julia> dg(folded) # get the total free energy of the predicted secondary structures
-6.6

julia> dg(folded_hot)
2.3

julia> dg(seq, temp=4) # this method calls the computationally intensive `fold` internally
-16.2
```

## Exported names
```@docs
fold
dot_bracket
dg
```

## Public names
```@docs
SeqFold.Structure
SeqFold.dg_cache
```