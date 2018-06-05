# Get tables of data for paper,
#  this just runs the func from compare_adapt_table3
source('~/GitHub/DOE-code/compare_adapt_table3.R')

# 10D
get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
            nams=c('Branin', 'Franke', 'Lim', 'Beam'),
            caption="Comparison of methods on the first four functions after evaluating $10D$ points.",
            label="datatable1_10D", no_candidates=TRUE, digits=2, batches='10D')
get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
            nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
            caption="Comparison of methods on the latter four functions after evaluating $10D$ points.",
            label="datatable2_10D", no_candidates=TRUE, digits=2, batches='10D')
# 20D
get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
            nams=c('Branin', 'Franke', 'Lim', 'Beam'),
            caption="Comparison of methods on the first four functions after evaluating $20D$ points.",
            label="datatable1_20D", no_candidates=TRUE, digits=2, batches='20D')
get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
            nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
            caption="Comparison of methods on the latter four functions after evaluating $20D$ points.",
            label="datatable2_20D", no_candidates=TRUE, digits=2, batches='20D')
# 40D
get_xtable2(selfs=list(bran1, franke1, lim1, beam1), 
            nams=c('Branin', 'Franke', 'Lim', 'Beam'),
            caption="Comparison of methods on the first four functions after evaluating $40D$ points.",
            label="datatable1_40D", no_candidates=TRUE, digits=2, batches='40D')
get_xtable2(selfs=list(otl1, piston1, bh1, wing1), 
            nams=c('OTL','Piston', 'Borehole', 'Wing weight'),
            caption="Comparison of methods on the latter four functions after evaluating $40D$ points.",
            label="datatable2_40D", no_candidates=TRUE, digits=2, batches='40D')
