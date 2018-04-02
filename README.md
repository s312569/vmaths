# vmaths


Some common maths functions for vectors. Works on normal vectors and various core.matrix vectors.

## Usage

Add `[vmaths "0.1.4"]` to your project.clj file.

Functions include:

```clj
;; Summary stats
mean, variance, pop-variance, ssd, psd, mad, meanad, quantile, median, mode, vrange,
midrange

;; Distance measures
euclidean, euclidean-squared, manhattan, mae, canberra

;; Similarity measures
pearson, cosine

;; normalise vectors using `normalise` and a normalisation function:
vmaths.core> (normalise [4 15 16 18 32 0 39 32 15 19] standard-score)
[[-1.2175510066851167 -0.3246802684493645 -0.24351020133702336 -0.08117006711234112
 1.0552108724604345 -1.5422312751344813 1.6234013422468225 1.0552108724604345
 -0.3246802684493645 0.0] [19.0 12.31981240838422]]
vmaths.core> (:doc (meta #'normalise))
"Normalises an object using function (nf) which can be one of:

   positional-unitization - ((x-median)/range)
   unit - (x/norm)
   positional-norm-in-range - ((x-median)/max(abs(x-median)))
   quotient-trans-sum - (x/sum)
   quotient-trans - (x/sd)
   quotient-trans-max - (x/max)
   positional-norm - ((x-median)/sqrt(sum((x-median)^2)))
   unitization - ((x-mean)/range)
   positional-standardization - ((x-median)/mad)
   norm-central-zero - ((x-midrange)/(range/2))
   quotient-trans-mean - (x/mean)
   quotient-trans-range - (x/range)
   positional-quotient-trans-median - (x/median)
   standard-score - ((x-mean)/sd)
   normalization - ((x-midrange)/(range/2))
   norm-in-range - ((x-mean)/max(abs(x-mean)))
   min-max - (x-minimum)/(maximum-minimum)
   unitization-zero-min - ((x-min)/range)
   positional-quotient-trans - (x/mad)

Returns a vector containing normalised object and a vector of values
used to normalize the object."
vmaths.core>
```

;; Online functions

The following functions have online versions that work on lazy
sequnces. They return a lazy sequence of their respective values after
the corresponding value in the input list:

mean, variance, pop-variance, ssd, psd, mad, meanad, quantile, median,
vrange, midrange

## License

Copyright © 2017 Jason Mulvenna

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
