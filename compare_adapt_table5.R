# CBE 8/11/19
# This is for the second submission, first to ONR.
# Adding functions that should make us look better.

# Get tables of data for paper,
#  this just runs the func from compare_adapt_table3
source('~/GitHub/DOE-code/compare_adapt_table3.R')

# Need to load the objects listed below from adapt_paper_runs4_addtruegrad.R
#  and adapt_paper_runs3

selfs1 <- list(bran1, banana1, lim1, gramacy2Dexp1)
nams1 <- c("Branin", "Roshan", "Lim", "Gramacy (2D)")
selfs2 <- list(otl1, gramacy6D1, bh1, wing1)
nams2 <- c("OTL", "Gramacy (6D)", "Borehole", "Wing weight")

# 10D
get_xtable2(selfs=selfs1, 
            nams=nams1,
            caption="Comparison of methods on the first four functions after evaluating $10D$ points.",
            label="datatable1_10D", no_candidates=TRUE, digits=2, batches='10D')
get_xtable2(selfs=selfs2, 
            nams=nams2,
            caption="Comparison of methods on the latter four functions after evaluating $10D$ points.",
            label="datatable2_10D", no_candidates=TRUE, digits=2, batches='10D')
# 20D
get_xtable2(selfs=selfs1, 
            nams=nams1,
            caption="Comparison of methods on the first four functions after evaluating $20D$ points.",
            label="datatable1_20D", no_candidates=TRUE, digits=2, batches='20D')
get_xtable2(selfs=selfs2, 
            nams=nams2,
            caption="Comparison of methods on the latter four functions after evaluating $20D$ points.",
            label="datatable2_20D", no_candidates=TRUE, digits=2, batches='20D')
# 40D
get_xtable2(selfs=selfs1, 
            nams=nams1,
            caption="Comparison of methods on the first four functions after evaluating $40D$ points.",
            label="datatable1_40D", no_candidates=TRUE, digits=2, batches='40D')
get_xtable2(selfs=selfs2, 
            nams=nams2,
            caption="Comparison of methods on the latter four functions after evaluating $40D$ points.",
            label="datatable2_40D", no_candidates=TRUE, digits=2, batches='40D')


# 
# # 10D
# get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
#             nams=c('Branin', 'Franke', 'Lim', 'Beam'),
#             caption="Comparison of methods on the first four functions after evaluating $10D$ points.",
#             label="datatable1_10D", no_candidates=TRUE, digits=2, batches='10D')
# get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
#             nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
#             caption="Comparison of methods on the latter four functions after evaluating $10D$ points.",
#             label="datatable2_10D", no_candidates=TRUE, digits=2, batches='10D')
# # 20D
# get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
#             nams=c('Branin', 'Franke', 'Lim', 'Beam'),
#             caption="Comparison of methods on the first four functions after evaluating $20D$ points.",
#             label="datatable1_20D", no_candidates=TRUE, digits=2, batches='20D')
# get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
#             nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
#             caption="Comparison of methods on the latter four functions after evaluating $20D$ points.",
#             label="datatable2_20D", no_candidates=TRUE, digits=2, batches='20D')
# # 40D
# get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
#             nams=c('Branin', 'Franke', 'Lim', 'Beam'),
#             caption="Comparison of methods on the first four functions after evaluating $40D$ points.",
#             label="datatable1_40D", no_candidates=TRUE, digits=2, batches='40D')
# get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
#             nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
#             caption="Comparison of methods on the latter four functions after evaluating $40D$ points.",
#             label="datatable2_40D", no_candidates=TRUE, digits=2, batches='40D')
