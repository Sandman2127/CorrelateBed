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
Provides standard output and err on the run operations, comparisons made and the correlation between all matching regions in bed1 & bed2.  

```
(base) sandman2127@Unix4Life CorrelateBed % julia CorrelateBed.jl --bed1 ./testData/bed4_low_magnitude_correlation.bed --bed2 ./testData/bed3_neg_correlation.bed
[STDOUT]: out of the total of 9 lines in bed1, we saw perfect matches at:9 positions or a total of: 100.0%
[STDOUT]: printing first 5 values of bed1 correlations
[0.437, 0.501, 0.109, 0.314, 0.203]
[STDOUT]: printing first 5 values of bed2 correlations
[-43.7, -50.1, -10.9, -31.4, -20.3]
[STDOUT]: the correlation between bed1:./testData/bed4_low_magnitude_correlation.bed and bed2: ./testData/bed3_neg_correlation.bed is: -0.9946758470074051
  0.549360 seconds (1.95 M allocations: 102.559 MiB, 2.93% gc time, 96.99% compilation time)
```
