# Correlate Bed 

### Assumptions:
This program can only correlate values at identical chromosome, start and end positions with the 5th column being the region the correlation is calculated at.  

### Installation:

#### Install julia >= 1.6 on your machine:
[https://julialang.org/downloads/](https://julialang.org/downloads/)

#### Install argparse library using the julia REPL:
```
$ julia
$ julia> ]
$ (@v1.8) pkg> add ArgParse 
 
```
You are now ready to run!

### Usage:
```
julia /path/to/CorrelateBed.jl --bed1 /path/to/bed1.bed --bed2 /path/to/bed2.bed
```

### Output:



