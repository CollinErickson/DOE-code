# test data.tree
library('data.tree')
acme <- Node$new('acme inc')
acme$AddChild('bb')
acme$AddChild('ccc')
acme$AddChild('dddd')
e <- acme$AddChild('e')
e$AddChild('ea')
acme$AddChild(49)
e$AddChild(paste(3,2,4))
bo <- e$AddChild('bo')
acme$AddChildNode(bo)
plot(acme)

# test dagR
library('dagR')
dag.draw(demo.dag1())
dag.draw(dag.init(covs=rep(1,4),arcs=c(1,2)))

# test gRbase
library('gRbase')
ug11 <- dag( ~a*b + b*c*d)
plot(ug11)
plot(dag(c(1,2,4,5),c(3,4),c(2,3)))