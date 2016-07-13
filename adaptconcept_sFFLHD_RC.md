Adaptive sFFLHD sampling concept
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
           [,1]       [,2]       [,3]
[1,] 0.80704681 0.76735849 0.02498429
[2,] 0.07001103 0.54852614 0.76341669
[3,] 0.61340117 0.82163175 0.58389773
[4,] 0.37208543 0.32099826 0.25862690
[5,] 0.59966032 0.05542548 0.85744879
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



Example: Gaussian
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

a2
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="1280px" height="800px" />


a3
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="1280px" height="800px" />

a4
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="1280px" height="800px" />

a5
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" width="1280px" height="800px" />


a6
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="1280px" height="800px" />

a7
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="1280px" height="800px" />

a8
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" width="1280px" height="800px" />

a9
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="1280px" height="800px" />

a10
===========
title: FALSE


```
[1] "Going one deeper"
```

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" width="1280px" height="800px" />


















Example: Sinumoid
================



Actual function: Sinusoid with a plateau


```r
contourfilled.func(sinumoid)
```

![plot of chunk unnamed-chunk-18](adaptconcept_sFFLHD_RC-figure/unnamed-chunk-18-1.png)


Adaptive sFFLHD
===============
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-19-1.png" title="plot of chunk unnamed-chunk-19" alt="plot of chunk unnamed-chunk-19" width="1280px" height="800px" />

a2
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-20-1.png" title="plot of chunk unnamed-chunk-20" alt="plot of chunk unnamed-chunk-20" width="1280px" height="800px" />


a3
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-21-1.png" title="plot of chunk unnamed-chunk-21" alt="plot of chunk unnamed-chunk-21" width="1280px" height="800px" />

a4
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-22-1.png" title="plot of chunk unnamed-chunk-22" alt="plot of chunk unnamed-chunk-22" width="1280px" height="800px" />

a5
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-23-1.png" title="plot of chunk unnamed-chunk-23" alt="plot of chunk unnamed-chunk-23" width="1280px" height="800px" />


a6
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-24-1.png" title="plot of chunk unnamed-chunk-24" alt="plot of chunk unnamed-chunk-24" width="1280px" height="800px" />

a7
===========
title: FALSE


```
[1] "Going one deeper"
```

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-25-1.png" title="plot of chunk unnamed-chunk-25" alt="plot of chunk unnamed-chunk-25" width="1280px" height="800px" />

a8
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-26-1.png" title="plot of chunk unnamed-chunk-26" alt="plot of chunk unnamed-chunk-26" width="1280px" height="800px" />

a9
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-27-1.png" title="plot of chunk unnamed-chunk-27" alt="plot of chunk unnamed-chunk-27" width="1280px" height="800px" />

a10
===========
title: FALSE

<img src="adaptconcept_sFFLHD_RC-figure/unnamed-chunk-28-1.png" title="plot of chunk unnamed-chunk-28" alt="plot of chunk unnamed-chunk-28" width="1280px" height="800px" />
