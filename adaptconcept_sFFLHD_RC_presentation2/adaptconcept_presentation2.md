Adaptive sampling discussion
========================================================
author: Collin Erickson
date: 7/28/16
autosize: true



What I was doing before seems overcomplicated
==================
- Boxes are arbitrary
- Points taken from boxes are just next ones from design
- E.g. problem: focused on a box because of the corner
- Box levels limit point selection unnecessarily


Goal
=============
- Adaptive design that is space-filling
- Balance between gooal design and focus on interesting areas
- Good design => Just use sFFLHD and take slices
- Adaptive => Pick optimum points based on criterion
- We need to balance these
- Give n options taken from design, let it pick best


Give n options taken from design, let it pick best
=============
- Idea: keep set of n points taken from design, pick "best" to run
- n = 1 => Normal design
- n = $\infty$ => Adaptive with no restrictions


Problem: quickly all n points will uninteresting
============
- Will just use the next point selected
- Will essentially become the design
- Solution: Let n increase over time
- How should n increase?
- Queue will become very long
- Lots of recomputing


Problem: No guarantee on acheiving design stages
===========
- Points could stay in queue forever
- Solution: occasionally take points that have in queue longest
- E.g. every tenth batch (or 10% each batch) is oldest in queue



Problem: Lots to recalculate
===============
- Only recalculate the ones that were nearly selected 
- Or ones that haven't been recalculated recently

Problem: points selected will be close to each other
==================
- Enforce slice/LH in some way?
- Mix in oldest in queue with best
- Active learning
- Mutual information?

PROBLEM!!!
==================
- GP predictive variance doesn't depend on Y values
- Only on design X and parameters (implicit Y)
- So space-filling does best at minimizing pvar (and MSE?)
- Reducing pvar not best use of adaptive
- Not a problem if searching for max or other response-based criterion

Working idea
======
- Take 3 batches of size L at each iteration
- Select L best
- 20% of batches are just oldest points in queue


================








<img src="adaptconcept_presentation2-figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="1280px" height="800px" />


=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="1280px" height="800px" />






=============

<img src="adaptconcept_presentation2-figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="1280px" height="800px" />


=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-9-1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="1280px" height="800px" />


=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-11-1.png" title="plot of chunk unnamed-chunk-11" alt="plot of chunk unnamed-chunk-11" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="1280px" height="800px" />

=============
<img src="adaptconcept_presentation2-figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="1280px" height="800px" />
