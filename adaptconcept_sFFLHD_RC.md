adaptconcept_sFFLHD_RC
========================================================
author: Collin Erickson
date: 7/13/2016
autosize: true
width: 1920
height: 1080

What I've worked on
=====================

- Implemented sFFLHD
- Adaptive sampling concept
 - Sample, then focus on subregion or return up level
- Combined these two to get adaptive sFFLHD sampling

sFFLHD
================




```r
s <- sFFLHD.seq(D = 3, L = 5)
s$get.batch()
```

```
          [,1]      [,2]       [,3]
[1,] 0.8382688 0.5280977 0.74472828
[2,] 0.5814938 0.6391916 0.02985124
[3,] 0.7542784 0.1638133 0.80409794
[4,] 0.1813992 0.9460096 0.40016803
[5,] 0.3191431 0.2637954 0.23064120
```

sFFLHD plot test
=============== 

```r
s <- sFFLHD.seq(D = 2, L = 5)
plot(NULL, xlim=0:1, ylim=0:1)
abline(h=(0:5)/5, v=(0:5)/5)
for(i in 1:5) points(s$get.batch(), col=i, pch=19)
```

![plot of chunk unnamed-chunk-3](adaptconcept_sFFLHD_RC-figure/unnamed-chunk-3-1.png)

***



```r
s <- sFFLHD.seq(D = 2, L = 3)
l <- 9
plot(NULL, xlim=0:1, ylim=0:1)
abline(h=(0:l)/l, v=(0:l)/l)
for(i in 1:27) points(s$get.batch(), col=i, pch=19)
```

![plot of chunk unnamed-chunk-4](adaptconcept_sFFLHD_RC-figure/unnamed-chunk-4-1.png)

```
[1] "Going one deeper"
```



Example
================



Actual function: Gaussian


```r
contourfilled.func(gaussian1)
```

![plot of chunk unnamed-chunk-6](adaptconcept_sFFLHD_RC-figure/unnamed-chunk-6-1.png)


Adaptive sFFLHD
===============
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="1280px" height="800px" />

a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="1280px" height="800px" />


a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="1280px" height="800px" />

a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="1280px" height="800px" />

a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" width="1280px" height="800px" />


a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="1280px" height="800px" />

a
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="1280px" height="800px" />
