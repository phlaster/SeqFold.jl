# Calling from Python scripts
SeqFold can be called from Python scripts via `juliacall package`:

## Showcase
1. Install [`juliacall`](https://pypi.org/project/juliacall) package through you favourite package manager (pip, Conda, uv):
```bash
$ pip install juliacall
```
2. In Python REPL:
```py
Python 3.13.7 | packaged by conda-forge | (main, Sep  3 2025, 14:30:35) [GCC 14.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from juliacall import Main as jl # this will download Julia binaries
>>> jl.seval("""using Pkg; Pkg.add("SeqFold")""")
>>> jl.seval("using SeqFold") # this exports tm, fold, dg, ...
>>> jl.fold("ACACGCAGCGGGTGTGGCTGTGGTGC")
Julia:
4-element Vector{SeqFold.Structure}:
    1   13   -1.9  STACK_DE:AC/TG 
    2   12   -2.3  STACK:CAC/GGG  
    4   10   -2.2  STACK:CG/GC    
    5    9    3.5  HAIRPIN:GC/CG  
>>> jl.dg("ACACGCAGCGGGTGTGGCTGTGGTGC")
-2.9
```

!!! note
    Remember, that calling a function in Julia for the first time takes more time due to JIT compiliation.