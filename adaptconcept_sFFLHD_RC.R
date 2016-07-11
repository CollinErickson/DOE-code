setOldClass("sFFLHD.seq")
setOldClass("UGP")
adapt.concept.sFFLHD.RC <- setRefClass("sFFLHD.seq",
  fields = list(
    func = "function", D = "integer", L = "integer", g = "integer", level = "integer",
    lims = "matrix", lims.second = "list",
    X = "matrix", Z = "numeric", Xnotrun = "matrix",
    s = "sFFLHD.seq", mod = "UGP"
  ),
  methods = list()
)