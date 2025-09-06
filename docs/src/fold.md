# Sequence Folding

`SeqFold.jl` implements the Zuker and Stiegler (1981) dynamic programming algorithm for predicting nucleic acid secondary structures.

## Basic Usage

```@docs
fold
```

## Structure Representation

```@docs
SeqFold.Structure
```

## Dot-Bracket Notation

```@docs
dot_bracket
```

## Free Energy Calculation

```@docs
dg
SeqFold.dg_cache
```

## GC Content Cache

```@docs
SeqFold.gc_cache
```