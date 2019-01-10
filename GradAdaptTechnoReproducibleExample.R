# Technometrics requires a reproducible example from the paper.
# This makes the plots for the 2D limnonpoly comparison between our method and standard IMSE
# I pulled all of my functions out of their packages into here.
# Only DoE.base is not included, it needs to be loaded with `library`.

# install.packages("DoE.base")
library(DoE.base)

# Test desirability functions
# A des func for finding large gradient
#' @param return_se whether the se prediction should be returned along with
#'   the des, all will be returned in data.frame, this will save
#'   time if calculating the werror function since it is faster
#'   to predict both at once instead of separately
des_func_grad_norm2_mean <- function(mod, XX, return_se=F) {
  if ("IGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod$grad_norm2_dist(XX=XX)$mean
    des <- mod$mod$grad_norm2_mean(XX=XX)
  } else if ("IGP_laGP_GauPro_kernel" %in% class(mod)) {
    # des <- mod$mod.extra$GauPro$mod$grad_norm2_dist(XX=XX)$mean
    des <- mod$mod.extra$GauPro$mod$grad_norm2_mean(XX=XX)
  } else if ("IGP_LOOEC_GauPro_kernel" %in% class(mod)) {
    des <- mod$mod$grad_norm2_mean(XX=XX)
  } else {
    stop("des_func_grad_norm2_mean only works with GauPro_kernel_model or laGP_GauPro_kernel")
  }
  des
}


msfunc <- function(func1,lims,pow=1L,batch=F, n=1e3) {
  # Find mean square of function over limits using grid sample
  #X1 <- simple.grid(10,nrow(lims),scaledto=lims) # Too big in high dimensions, switching to just random points
  d <- nrow(lims)
  X1 <- simple.random(n=n, d=d, scaledto=lims)
  if(batch) {return(mean(func1(X1)^pow))}
  mean(apply(X1,1,func1)^pow)
}

simple.random <- function(n,d,scaledto=NULL) {
  m <- matrix(runif(n*d), ncol=d, nrow=n)
  if(!is.null(scaledto)) {
    m <- m * matrix(scaledto[,2]-scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T) + matrix(scaledto[,1],nrow=nrow(m),ncol=ncol(m),byrow=T)
  } 
  m
}

msecalc <- function(truefunc, guessfunc,lims, n=500) {
  #X1 <- simple.grid(20,nrow(lims),scaledto=lims)
  #X1 <- lhs::maximinLHS(n, nrow(lims))
  d <- nrow(lims)
  X1 <- matrix(runif(n*d), n, d)
  mean((apply(X1,1,function(xx){truefunc(xx) - guessfunc(xx)}))^2)
}

simple.LHS <- function(n,d,scaled=TRUE,centered=FALSE) {
  m <- matrix(rep(1:n,d),n,d)
  m <- apply(m,2,function(xx){sample(xx)})
  if(scaled) m <- (m - runif(n*d) ) / n
  if(centered) m <- m - ifelse(scaled,.5,n/2+.5)
  m
}

#' Close all open screens
#' 
#' Closes the screens open, which happens
#' when plotting with `split.screen` is interrupted.
#' It often happens when there is a error while plotting.
#' When you try to plot
#' the next thing it gives an error.
#' Running this function will reset the plot screen.
#' It just does `close.screen(all.screens=TRUE)` but is faster to type.
#' 
#' @param silent Should the output of `close.screen` not be returned?
#' 
#' @examples
#' # Split screen into fourths
#' split.screen(c(2,2))
#' hist(rnorm(100))
#' screen(2)
#' hist(runif(100))
#' # Use csa() to go back to normal plotting
#' csa()
#' hist(rexp(100))
#' @export
csa <- function(silent=FALSE) {
  if (silent) { # Suppresses "FALSE" if already closed
    invisible(close.screen(all.screens=TRUE))
  } else {
    close.screen(all.screens=TRUE)
  }
}

#' Makes filled contour plot from function 
#' 
#' A contour plot of the given function without sidebar by default.
#' It calls the function `cf_grid` to make the actual plot.
#' 
#' @param fn0  function to plot, first argument must be two-dimensional
#' @param n  number of points in each dimension
#' @param xlim  x limits for the contour plot
#' @param ylim  y limits for the contour plot
#' @param xylim  x and y limits for the contour plot, use when both are same
#' #@param mainminmax  whether the min and max values should be shown in the title of plot
#' @param batchmax  number of datapoints that can be computed at a time
#' @param out.col.name  if a column needs to be selected from the function, specify it
#' @param out.name Selects with a $ the name from output to be used, for lists and data frames
#' #@param pretitle Text to be preappended to end of plot title
#' #@param posttitle Text to be appended to end of plot title
#' #@param title Title for the plot
#' #@param mainminmax_minmax Whether [min,max]= should be shown in title or just the numbers
#' @param pts Points to plot on top of contour
#' @param use_lines If the contour should be made with lines. Otherwise is made
#' using colors. Defaults to colors.
#' @param ...  Passed to cf_grid
#' @examples 
#' cf_func(function(x){x[1]*x[2]})
#' cf_func(function(x)(exp(-(x[1]-.5)^2-5*(x[2]-.5)^2)))
#' cf_func(function(xx){exp(-sum((xx-.5)^2/.1))}, bar=TRUE)
#' cf_func(function(xx){exp(-sum((xx-.5)^2/.1))}, bar=TRUE, mainminmax=TRUE)
#' cf_func(function(x)(exp(-(x[1]-.5)^2-5*(x[2]-.5)^2)), use_lines=TRUE)
#' @references
#' [1] filled.contour R function, copied function but removed part for sidebar
#' @references
#' [2] http://stackoverflow.com/questions/16774928/removing-part-of-a-graphic-in-r, answer by P Lapointe
#' @export
cf_func <- function(fn0, n=100,
                    xlim=c(0,1), ylim=c(0,1), xylim=NULL,
                    batchmax=1, out.col.name=NULL,
                    out.name=NULL,
                    pts=NULL,
                    use_lines=FALSE,
                    ...) {
  if(!is.null(out.col.name)) {
    fn <- function(xx){fn0(xx)[,out.col.name]}
  } else if (!is.null(out.name)) {
    fn <- function(xx){fn0(xx)[[out.name]]}
  } else {
    fn <- fn0
  }
  if (!is.null(xylim)) {xlim <- ylim <- xylim}
  x <- seq(xlim[1],xlim[2],length.out = n)
  y <- seq(ylim[1],ylim[2],length.out = n)
  z <- eval_over_grid_with_batch(x, y, fn, batchmax)
  if (use_lines) {
    contour(x, y, z, ...)
    points(pts, pch=19)
  } else
    cf_grid(x,y,z, pts=pts, ...)
}

#' Evaluate function over grid of points
#' 
#' `batchmax` gives how many can be evaluated at a time.
#' If more than 1, then the input is given to the function
#' as rows of a matrix.
#' 
#'
#' @param x Vector of x values to evaluate
#' @param y Vector of y values to evaluate
#' @param fn Function that takes in a length two vector if `batchmax` is 1
#' or a matrix with two columns if greater than 1.
#' @param batchmax Number of points that can evaluated simultaneously.
#' If 1, points are passed to `fn` as a vector of length two.
#' If greater than 1, points are passed to `fn` as rows of a matrix.
#'
#' @return Matrix of size `length(x)` by `length(y)`
#' @export
#'
#' @examples
#' eval_over_grid_with_batch(c(0,.5,1), c(10,20,30), function(a)a[1]+a[2], batchmax=1)
#' eval_over_grid_with_batch(c(0,.5,1), c(10,20,30), function(a)a[,1]+a[,2], batchmax=Inf)
eval_over_grid_with_batch <- function(x, y, fn, batchmax) {
  nx <- length(x)
  ny <- length(y)
  if(batchmax<=1) { # calculate single Z value at a time
    #for(xi in 1:n) for(yi in 1:n) {z[xi,yi] <- fn(c(x[xi],y[yi]))}
    fn_outer <- Vectorize(function(xi, yi) {fn(c(x[xi], y[yi]))})
    z <- outer(1:nx, 1:ny, fn_outer)
  } else {
    z <- matrix(NA,nx,ny)
    inbatch <- 0
    for(xi in 1:nx) {
      for(yi in 1:ny) {
        if(inbatch==0) XYi <- matrix(c(xi,yi),ncol=2)
        else XYi <- rbind(XYi,matrix(c(xi,yi),ncol=2))
        inbatch <- inbatch + 1
        if(inbatch == batchmax | (xi==nx & yi==ny)) {
          Zbatch <- fn(matrix(c(x[XYi[,1]],y[XYi[,2]]),ncol=2,byrow=F))
          for(rowbatch in 1:length(Zbatch)) {
            z[XYi[rowbatch,1],XYi[rowbatch,2]] <- Zbatch[rowbatch]
          }
          inbatch <- 0
          rm(XYi)
        }
      }
    }
  }
  z
}

#' Create a contour plot from a grid of data
#' 
#' Makes filled contour plot with an optional sidebar, essentially filled.contour function.
#' This version uses the split.screen() function to add the sidebar if bar is TRUE.
#' By default it won't show the bar but will show the min and max values in the plot title
#' along with their colors.
#' Using this function will make other functions such as points() called afterwards not put points
#' where you expect. Pass anything you want added to the plot area to afterplotfunc
#' as a function to get it to work properly.
#' 
#' @param x  x values, must form grid with y. If not given, it is assumed to be from 0 to 1.
#' @param y  y values, must form grid with x. If not given, it is assumed to be from 0 to 1.
#' @param z  z values at grid locations
#' @param xlim  x limits for the plot.
#' @param ylim  y limits for the plot.
#' @param zlim  z limits for the plot.
#' @param levels  a set of levels which are used to partition the range of z. Must be strictly increasing (and finite). Areas with z values between consecutive levels are painted with the same color.
#' @param nlevels  if levels is not specified, the range of z, values is divided into approximately this many levels.
#' @param color.palette  a color palette function to be used to assign colors
#' in the plot. Defaults to cm.colors. Other options include rainbow,
#' heat.colors, terrain.colors, topo.colors, and function(x) {gray((1:x)/x)}.
#' @param col  an explicit set of colors to be used in the plot. This argument overrides any palette function specification. There should be one less color than levels
#' @param plot.title  statements which add titles to the main plot.
#' @param plot.axes  statements which draw axes (and a box) on the main plot. This overrides the default axes.
#' @param key.title  statements which add titles for the plot key.
#' @param key.axes  statements which draw axes on the plot key. This overrides the default axis.
#' @param asp  the y/x aspect ratio, see plot.window.
#' @param xaxs  the x axis style. The default is to use internal labeling.
#' @param yaxs  the y axis style. The default is to use internal labeling.
#' @param las  the style of labeling to be used. The default is to use horizontal labeling.
#' @param axes  logical indicating if axes should be drawn, as in plot.default.
#' @param frame.plot  logical indicating if a box should be drawn, as in plot.default.
#' @param bar Should a bar showing the output range and colors be shown on the right?
#' @param pts Points to plot on top of contour
#' @param reset.par Should the graphical parameters be reset before exiting? Usually should be
#' unless you need to add something to the plot afterwards and bar is TRUE.
#' @param pretitle Text to be preappended to end of plot title
#' @param posttitle Text to be appended to end of plot title
#' @param main Title for the plot
#' @param mainminmax  whether the min and max values should be shown in the title of plot
#' @param mainminmax_minmax Whether [min,max]= should be shown in title or just the numbers
#' @param afterplotfunc Function to call after plotting, such as adding points or lines.
#' @param cex.main The size of the main title. 1.2 is default.
#' @param par.list List of options to pass to par
#' @param xaxis Should x axis be added?
#' @param yaxis Should y axis be added?
#' @param ...  additional graphical parameters, currently only passed to title().
#' @importFrom grDevices cm.colors
#' @importFrom graphics .filled.contour
#' @importFrom graphics Axis 
#' @importFrom graphics box
#' @importFrom graphics plot.new 
#' @importFrom graphics plot.window 
#' @importFrom graphics points 
#' @importFrom graphics title
#' @importFrom graphics par
#' @importFrom graphics axis layout lcm rect
#' @importFrom graphics split.screen screen close.screen
#' @examples 
#' x <- y <- seq(-4*pi, 4*pi, len = 27)
#' r <- sqrt(outer(x^2, y^2, "+"))
#' cf_grid(cos(r^2)*exp(-r/(2*pi)))
#' cf_grid(r, color.palette=heat.colors, bar=TRUE)
#' cf_grid(r, color.palette=function(x) {gray((1:x)/x)}, bar=TRUE)
#' @references
#' [1] filled.contour R function, copied function but removed part for sidebar
#' @references
#' [2] http://stackoverflow.com/questions/16774928/removing-part-of-a-graphic-in-r, answer by P Lapointe
#' @export
cf_grid <-
  function (x = seq(0, 1, length.out = nrow(z)), 
            y = seq(0, 1,length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes, bar=F, pts=NULL, reset.par=TRUE,
            pretitle="", posttitle="", main=NULL,
            mainminmax=!bar, mainminmax_minmax=TRUE,
            afterplotfunc=NULL,
            cex.main=par()$cex.main,
            par.list=NULL,
            xaxis=TRUE, yaxis=TRUE,
            ...)
  {
    # filled.contour gives unnecessary legend, this function removes it
    # Used P Lapointe's solution from here: http://stackoverflow.com/questions/16774928/removing-part-of-a-graphic-in-r
    #   also had to changed .Internal(fillcontour) to .filled.contour
    #   and change layout to layout(matrix(c(1, 1), ncol = 1L), widths = c(1, lcm(w)))
    # Created 3/28/16 by Collin Erickson
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    #if (reset.par) {on.exit({par(par.orig);close.screen(1)})}#all=TRUE)})}
    
    # This allows user to pass in par.list and it will return it after plotting
    par.names.to.save <- c("mar", "las", "mfrow", names(par.list))
    mar.orig <- (par.orig <- par(par.names.to.save))$mar
    if (!is.null(par.list)) {
      par(par.list)
    }
    
    #on.exit(close.screen(all=TRUE))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(if(bar) 2 else 1, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    start.screen.number <- screen()
    if (bar) {
      #split.screen(c(1,2))
      screen.numbers <- split.screen(matrix(c(0,.85,0,1,.85,1,0,1), ncol=4, byrow=T))
      screen1 <- screen.numbers[1]
      screen2 <- screen.numbers[2]
      screen(screen2)
      mar <- mar.orig
      mar[4L] <- 2.5#mar[2L] # right
      mar[1] <- 2.2 # bottom
      mar[3] <- if (mainminmax | !is.null(main)) 1.3 else .3 #1.3#1.3 # top
      mar[2L] <- 0#1 # left
      par(mar = mar)
      plot.new()
      plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                  yaxs = "i")
      rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
      if (missing(key.axes)) {
        if (axes) 
          axis(4)
      }
      else key.axes
      box()
      if (!missing(key.title))
        key.title
      mar <- mar.orig
      mar[4L] <- 1 # right # Why is this here?
      close.screen(screen2)
      screen(screen1)
      mar[1L] <- 2.2 # bottom
      mar[2L] <- 2.5 # left
      mar[3L] <- if (mainminmax | !is.null(main)) 1.3 else .3 #1.3# 1.3 # top
    }
    if (!bar) {
      # Using screen even with 1 screen to avoid error. Adding points after didn't show up properly
      screen.numbers <- split.screen(c(1,1)) 
      screen1 <- screen.numbers[1]
      screen(screen1)
      # Changing the margin to get bigger and square
      mar <- mar.orig #<- par()$mar
      mar[1] <- 2.2 # bottom
      mar[2] <- 2.5 # left
      mar[3] <- if (mainminmax | !is.null(main)) 1.3 else .3 # top
      mar[4] <- 1 # right
      
      # if (!missing(plot.axes)) {
      #   # TODO I shouldn't use plot.axes like this, FIX THIS
      #   mar[1] <- .3 # bottom
      #   mar[2] <- .3 # left
      #   mar[4] <- .3 # 1 # right
      # }
      if (!xaxis && !yaxis) {
        mar[1] <- .3 # bottom
        mar[2] <- .3 # left
        mar[3] <- if (mainminmax | !is.null(main)) 1.3 else .3 # top
        mar[4] <- .3 # right
      } else if (!xaxis) {
        mar[1] <- 1-.7 # bottom
        mar[3] <- if (mainminmax | !is.null(main)) 1.3 else .3 # top
        mar[4] <- .3 # right
      } else if (!yaxis) {
        mar[2] <- 1-.7 # left
        mar[3] <- if (mainminmax | !is.null(main)) 1.3 else .3 # top
        mar[4] <- .3 # right
      }
    }
    par(mar = mar)
    # par(cex.axis = 2)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    #.Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels),
    #                        col = col))
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                    col = col)
    # Something like this will remove axis numbers and ticks
    # Axis(x, side=1, labels=F, tick=F)
    if (missing(plot.axes)) {
      if (axes) {
        # title(main = "", xlab = "", ylab = "")
        if (xaxis) Axis(x, side = 1)
        if (yaxis) Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot)
      box()
    #if (missing(plot.title))
    #  title(...)
    #else plot.title
    
    if (mainminmax | !is.null(main)) {
      make.multicolor.title(main=main, z=z, pretitle=pretitle, posttitle=posttitle, mainminmax_minmax=mainminmax_minmax, cex.main=cex.main)
    }
    
    if (!is.null(pts)) {
      if (!is.matrix(pts)) { # if not a matrix, make it a matrix by row
        if (is.numeric(pts) && (length(pts)%%2==0)) {
          pts <- matrix(pts, ncol=2, byrow = T)
        }
      }
      points(pts, pch=19)
    }
    
    if (!is.null(afterplotfunc)) {
      afterplotfunc()
    }
    
    reset.par.func <- function() {
      if (T) {close.screen(screen1)}
      if (start.screen.number != FALSE) {screen(start.screen.number, new=FALSE)}
      par(par.orig)
    }
    if (reset.par) {# Either reset parameters
      reset.par.func()
      invisible()
    } else { # or return it so user can do it later
      return(reset.par.func)
    }
  }

make.multicolor.title <- function(main, z, pretitle, posttitle, mainminmax_minmax, cex.main=par()$cex.main) {
  if(is.null(main)) {
    title_text <- c(pretitle)
    title_color <- c(1)
    if (mainminmax_minmax) {
      title_text  <- c(title_text, '[','min',      ', ','max',      '] = ')
      title_color <- c(title_color,1,  "#80FFFFFF",1,   "#FF80FFFF",1)
    }
    title_text  <- c(title_text, "[",signif(min(z),3),', ',signif(max(z),3),']',posttitle)
    title_color <- c(title_color,1,  1,               1,   1,               1,  1)
    multicolor.title(title_text,title_color, cex.main=cex.main)
  } else {
    multicolor.title(main, 1, cex.main=cex.main)
  }
}

#' Makes plot title using specified colors for the text
#' @param main  Text to put in main title of plot
#' @param col.main  Colors to use for the text
#' @param collapse  What to put between elements of main, defaults to "" but " " might be appropriate
#' @param cex.main The size of the main title. 1.2 is default.
#' @examples 
#' plot(1:4)
#' multicolor.title(c('Black, ','red, ','green'),c(1,2,3))
#' @export
multicolor.title <- function(main,col.main, collapse='', cex.main=par()$cex.main) {
  if (length(main) != length(col.main)) {stop('main and col must have same length')}
  n <- length(main)
  if(n==1) {
    title(bquote(.(main[1])),col.main=col.main[1], cex.main=cex.main)
  } else {
    # print first
    title(bquote(.(main[1]) * phantom(.(paste0(main[2:n],collapse=collapse)))),col.main=col.main[1], cex.main=cex.main)
    
    # print middle
    if(n > 2) {
      for(i in 2:(n-1)) {
        title(bquote(phantom(.(paste0(main[1:(i-1)],collapse=collapse))) * .(main[i]) * phantom(.(paste0(main[(i+1):n],collapse=collapse)))),col.main=col.main[i], cex.main=cex.main) 
      }
    }
    
    # print last
    title(bquote(phantom(.(paste0(main[1:(n-1)],collapse=collapse))) * .(main[n])),col.main=col.main[n], cex.main=cex.main)
  }
}

#' Simpler function for making contours with cf package.
#' Won't give argument completion, so all must be specified
#'
#' @param ... Arguments to be passed to cf_func or cf_data based on 
#' data type of first argument. If D is given as argument, then it
#' is passed to cf_highdim.
#'
#' @return Whatever is returned from other function, probably nothing
#' @export
#'
#' @examples
#' cf(function(x){x[1]^2 - x[2]})
#' x <- runif(20)
#' y <- runif(20)
#' z <- exp(-(x-.5)^2-5*(y-.5)^2)# + rnorm(20,0,.05)
#' cf(x,y,z)
#' cf(function(x){x[1]^2 - x[2]}, D=3)
cf <- function(...) {
  dots <- list(...)
  if (is.function(dots[[1]])) {
    if ("D" %in% names(dots)) {
      cf_highdim(...)
    } else {
      cf_func(...)
    }
  } else if (is.numeric(dots[[1]])) {
    cf_data(...)
  } else {
    stop("Data not recognized. Use cf_func for function or 
         cf_data for data or cf_grid for full grid of data.")
  }
}


#' limnonpoly: 2 dimensional function.
#' Equation 28 from Lim et al 2002.
#'
#' @references Lim, Yong B., Jerome Sacks, W. J. Studden, and William J. Welch.
#' "Design and analysis of computer experiments when the output is highly
#' correlated over the input space."
#' Canadian Journal of Statistics 30, no. 1 (2002): 109-126.
#' @export
#' @rdname test_func_apply
#' @examples
#' limnonpoly(runif(2))
limnonpoly <- function(x, scale_it=F, scale_low = c(0,0), scale_high = c(1,1), noise=0) {
  if (is.matrix(x)) apply(x, 1, TF_limnonpoly)
  else TF_limnonpoly(x)
}
#' TF_limnonpoly: A function taking in a single vector.
#' 2 dimensional function.
#' See corresponding function with "TF_" for more details.
#' @export
#' @rdname TF_OTL_Circuit
#' @examples
#' TF_limnonpoly(runif(2))
TF_limnonpoly <- function(x) {
  ((30+5*x[1]*sin(5*x[1]))*(4+exp(-5*x[2])) - 100) / 6
}


#' IGP general function
#'
#' @param package Package to use
#' @param X Design matrix
#' @param Z Response matrix or vector
#' @param ... Arguments passed on to IGP_<package>
#'
#' @return IGP model
#' @export
#'
#' @examples
#' x <- seq(0,1,l=10)
#' y <- abs(sin(2*pi*x))
#' IGP(x,y,'DiceKriging')
IGP <- function(X=NULL, Z=NULL, package=NULL, ...) {
  if (length(package)==0) {
    stop("No package specified Error # 5792324572")
  } else if (package %in% c("laGP", "laGP")) {
    u <- IGP_laGP$new(X=X, Z=Z, ...)
  } else if (package %in% c("GauPro_kernel", "gaupro_kernel")) {
    u <- IGP_GauPro_kernel$new(X=X, Z=Z, ...)
  } else if (tolower(package) %in% c("lagp_gaupro_kernel")) {
    u <- IGP_laGP_GauPro_kernel$new(X=X, Z=Z, ...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}


#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP(package='laGP',X=X1,Z=Z1, corr="gauss")
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL, corr="gauss",
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_base <- R6::R6Class(classname = "IGP",
                        public = list(
                          X = NULL, #"matrix",
                          Z = NULL, #"numeric",
                          package = NULL, #"character",
                          .init = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          .update = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          .predict = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          .predict.se = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          .predict.var = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          #.grad = function(...){stop("This function must be overwritten by subclass")},
                          .delete = function(...){self$mod <- NULL}, #"function",
                          #.theta = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          #.nugget = function(...){stop("This function must be overwritten by subclass")}, #"function",
                          #.mean = function(...){stop("This function must be overwritten by subclass")}, # function that gives mean
                          mod = NULL, #"list", # First element is model
                          mod.extra = list(), #"list", # list to store additional data needed for model
                          n.at.last.update = NULL, #"numeric", # count how many in update, must be at end of X
                          corr = NULL, #"numeric",
                          estimate.nugget = NULL, #"logical", Should the nugget be estimated?
                          nugget0 = NULL, #"numeric" # What value should the nugget be set to? NOT logical. If estimate.nugget==TRUE, then it's the starting value
                          
                          initialize = function(X=NULL, Z=NULL, package=NULL, corr="gauss", estimate.nugget=TRUE, nugget0=1e-8, ...) {#browser()
                            if (!is.null(X)) {self$X <- if (is.matrix(X)) X else matrix(X, ncol=1)} # Add else for 1D data passed as vector
                            if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
                            self$package <- package
                            self$n.at.last.update <- 0
                            self$corr <- corr
                            self$estimate.nugget <- estimate.nugget
                            self$nugget0 <- nugget0
                            
                            #if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
                            if(length(self$X) != 0 & length(self$Z) != 0) {
                              self$init(...)
                            }
                          }, # end initialize
                          init = function(X=NULL, Z=NULL, ...) {#browser()
                            if (!is.null(X)) {self$X <- X}
                            if (!is.null(Z)) {self$Z <- Z}
                            if (length(self$X) == 0 | length(self$Z) == 0) {stop("X or Z not set")}
                            self$n.at.last.update <- nrow(self$X)
                            if (max(self$Z) - min(self$Z) < 1e-8) {warning("Z values are too close, adding noise"); self$Z <- self$Z + rnorm(length(self$Z), 0, 1e-6)}
                            
                            self$.init(...)
                          }, # end init
                          update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
                            if (self$n.at.last.update == 0) {
                              #self$init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
                              x <- if(!is.null(Xall)) Xall else Xnew
                              z <- if (!is.null(Zall)) Zall else Znew
                              self$init(X = x, Z = z)
                            } else {
                              if (!is.null(Xall)) {self$X <- Xall} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew)}
                              if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- c(self$Z, Znew)}
                              self$.update(...)
                            }
                            self$n.at.last.update <- nrow(self$X)
                          }, # end update
                          predict = function(XX, se.fit = FALSE, ...) {#browser()
                            if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                            self$.predict(XX, se.fit=se.fit, ...)
                          },
                          predict.se = function(XX, ...) {
                            if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                            self$.predict.se(XX, ...=...)
                          },
                          predict.var = function(XX, ...) {
                            if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
                            self$.predict.var(XX, ...=...)
                          },
                          grad = function (XX, num=FALSE) {#browser() # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
                            if (!is.matrix(XX)) {
                              if (ncol(self$X) == 1) XX <- matrix(XX, ncol=1)
                              else if (length(XX) == ncol(self$X)) XX <- matrix(XX, nrow=1)
                              else stop('Predict input should be matrix')
                            } else {
                              if (ncol(XX) != ncol(self$X)) {stop("Wrong dimension input")}
                            }
                            if (is.null(self$.grad) | num) { # if no method, use numerical
                              #print('using num')
                              self$grad_num(XX)
                            } else {#print('using package')
                              self$.grad(XX)
                            }
                          },
                          grad_num = function (XX) {
                            grad.func <- function(xx) self$predict(xx)
                            grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
                            grad1 <- apply(XX, 1, grad.apply.func)
                            if (ncol(self$X) == 1) return(grad1)
                            t(grad1)
                          },
                          grad_from_theta = function(XX, theta) {
                            if (missing(theta)) {
                              theta <- self$theta()
                              if (is.null(theta)) {
                                stop("Need theta for grad_from_theta")
                              }
                            }
                            mu <- self$mean()
                            D <- ncol(self$X)
                            N <- nrow(self$X)
                            if (!is.matrix(XX)) {
                              if (D == 1) XX <- matrix(XX, ncol=1)
                              else if (length(XX) == D) XX <- matrix(XX, nrow=1)
                              else stop('Predict input should be matrix')
                            } else {
                              if (ncol(XX) != D) {stop("Wrong dimension input")}
                            }
                            # kx.xx <- self$corr_func(self$X, XX, theta=self$theta)
                            kx.xx <- GauPro::corr_gauss_matrix(self$X, XX, theta)
                            Kx <- GauPro::corr_gauss_matrix_symC(self$X, theta)
                            Kx_nug <- Kx + diag(self$nugget(), nrow(Kx))
                            Kinv_Z_minus_mu <- solve(Kx_nug, self$Z - mu)
                            
                            grad1 <-   vapply(1:nrow(XX),
                                              Vectorize(
                                                function(k) {
                                                  t(-2 * outer(1:N, 1:D, Vectorize(function(i,j) {theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]}))
                                                  )  %*%Kinv_Z_minus_mu
                                                }
                                              )
                                              , numeric(D)
                            )
                            if (D == 1) return(grad1)
                            t(grad1)
                          },
                          grad_norm = function (XX) {#browser()
                            grad1 <- self$grad(XX)
                            if (!is.matrix(grad1)) return(abs(grad1))
                            apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
                          },
                          sample = function(XX, n=1) {
                            if (length(XX) != ncol(self$X)) {stop("Can only sample one point at a time right now error 23537898")}
                            XX.pred <- self$predict(XX=XX, se.fit=T)
                            rnorm(n=n, mean=XX.pred$fit, sd=XX.pred$se.fit)
                          },
                          theta = function() {
                            self$.theta()
                          },
                          nugget = function() {
                            self$.nugget()
                          },
                          s2 = function() {
                            self$.s2()
                          },
                          mean = function() {
                            if (!is.null(self$.mean)) {
                              self$.mean()
                            } else {
                              self$predict(matrix(rep(max(abs(self$X)) * 10,ncol(self$X)), nrow=1))
                            }
                          },
                          max.var = function() {
                            self$predict.var(matrix(rep(max(abs(self$X)) * 10,ncol(self$X)), nrow=1))
                          },
                          at.max.var = function(X, val=.9) {#browser() #logical if pred var at least 90% of max var
                            maxvar = c(self$max.var())
                            self$predict.var(X) > val * maxvar
                          },
                          prop.at.max.var =function(Xlims = matrix(c(0,1), nrow=ncol(self$X), ncol=2, byrow=T), n = 200, val=.9) {#browser()
                            maxvar = c(self$max.var())
                            X <- apply(Xlims, 1, function(Xlim) {runif(n, Xlim[1], Xlim[2])})
                            sum(self$predict.var(X) > val * maxvar) / n
                          },
                          plot = function() {#browser()
                            minx <- min(self$X)
                            maxx <- max(self$X)
                            minxeval <- minx - .03 * (maxx - minx)
                            maxxeval <- maxx + .03 * (maxx - minx)
                            if (ncol(self$X) == 1) {
                              XX <- matrix(seq(minxeval,maxxeval,length.out = 300), ncol=1)
                              pp <- self$predict(XX=XX, se.fit=TRUE)
                              pm <- pp$fit
                              ps <- pp$se.fit
                              phigh <- pm + 2 * ps
                              plow  <- pm - 2 * ps
                              plot(XX, pm, col='white', ylim=c(min(plow), max(phigh)),
                                   xlab="X", ylab="Z")
                              points(XX, phigh, type='l', col=2, lwd=2)
                              points(XX, plow, type='l', col=2, lwd=2)
                              points(XX, pm, type='l', lwd=3)
                              points(self$X, self$Z, pch=19, cex=2)
                            } else {
                              stop("No plot method for higher than 1D")
                            }
                          },
                          delete = function(...) {
                            self$.delete(...=...)
                          },
                          finalize = function(...) {
                            self$delete() # Mostly for laGP to delete, Python should close connection
                          }
                        )
)


#' IGP R6 object for fitting GauPro model
#'
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_GauPro_kernel$new(X=X1,Z=Z1, parallel=FALSE)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/IGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{update(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...)}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_GauPro_kernel <- R6::R6Class(
  classname = "IGP_GauPro_kernel",
  inherit = IGP_base,
  public = list(
    .init = function(..., kernel=NULL, theta=NULL) {
      if (!is.null(kernel)) {
        # kernel will be passed in
      } else if (any(c("R6ClassGenerator", "GauPro_kernel")%in% class(self$corr))) {
        kernel <- self$corr
      } else if (self$corr[[1]] == "gauss") {
        kernel <- GauPro::Gaussian$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "matern32") {
        kernel <- GauPro::Matern32$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "matern52") {
        kernel <- GauPro::Matern52$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "exponential") {
        kernel <- GauPro::Exponential$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "periodic") {
        kernel <- GauPro::periodic$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "rationalquadratic") {
        kernel <- GauPro::RatQuad$new(D=ncol(self$X))
      } else {
        stop("Corr/kernel not recognized in IGP_GauPro_kernel")
      }
      if (!is.null(theta)) {kernel$beta <- log(theta, 10)}
      m <- GauPro::GauPro_kernel_model$new(X=self$X, Z=self$Z, kernel=kernel, nug.est=self$estimate.nugget, nug=self$nugget0, ...)
      self$mod <- m
    }, #"function to initialize model with data
    .update = function(...) {
      self$mod$update(Xall=self$X, Zall=self$Z, ...)
    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      if (se.fit) {
        preds <- self$mod$pred(XX=XX, se.fit=T)
        list(fit=preds$mean, se.fit=preds$se)
      } else {
        c(self$mod$pred(XX=XX))
      }
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$se}, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$s2}, #"function to predict the variance
    .grad = function(XX) {self$mod$grad(XX=XX)}, # function to calculate the gradient
    .delete = function(...){self$mod <- NULL}, #"function to delete model beyond simple deletion
    .theta = function() {10 ^ self$mod$kernel$beta}, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = function() {self$mod$nug}, #"function to get nugget
    .s2 = function() {self$mod$s2_hat},
    .mean = function() {self$mod$trend$m} # function that gives mean (constant, other functions not implemented)
    
  )
)


#' IGP R6 object for fitting laGP model
#'
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_laGP$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/IGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{update(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...)}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_laGP <- R6::R6Class(classname = "IGP_laGP", inherit = IGP_base,
                        public = list(
                          .init = function(..., d=NULL, g=NULL, theta=NULL, nugget=NULL, no_update=FALSE) {
                            if (self$corr[[1]] != "gauss") {
                              stop("laGP only uses Gaussian correlation")
                            }
                            
                            if (is.null(d) & !is.null(theta)) {d <- 1/theta}
                            if (is.null(g) && is.null(nugget) && !is.null(self$nugget0)) {g <- self$nugget0}
                            if (is.null(g) & !is.null(nugget)) {g <- nugget}
                            
                            
                            if (no_update) {
                              if (is.null(d)) {stop("d or theta must be given when using no_update")}
                              if (is.null(g)) {stop("g or nugget must be given when using no_update")}
                              da_start <- d
                              ga_start <- g
                            } else { # Estimating params
                              da <- laGP::darg(list(mle=TRUE), X=self$X)
                              ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = T)
                              if (inherits(ga.try, "try-error")) {
                                # warning("Adding noise to ga in laGP"); # Not too important a warning
                                # Sometimes first try doesn't work, so looping with bigger eps
                                eps_ga <- 1e-2
                                while (TRUE) {
                                  ga <- try(laGP::garg(list(mle=TRUE),
                                                       y=self$Z+rnorm(length(self$Z),0,eps_ga)),
                                            silent=T)
                                  if (!inherits(ga, "try-error")) {break()}
                                  eps_ga <- 2 * eps_ga
                                }
                              }
                              # Follow recommendations for small samples, otherwise use bigger range
                              drange <- if (nrow(self$X)<20) c(da$min, da$max) else c(1e-3,1e4) #c(da$min, da$max), # Don't like these small ranges
                              grange <- c(ga$min, ga$max)
                              # da_start <- if (!is.null(d)) d else da$start
                              # ga_start <- if (!is.null(g)) g else ga$start
                              # Need to make sure starting values are in ranges
                              if (!is.null(d)) {
                                da_start <- d
                                drange <- c(min(drange[1], d), max(drange[2], d))
                              } else {
                                da_start <- da$start
                              }
                              if (!is.null(g)) {
                                ga_start <- g
                                grange <- c(min(grange[1], g), max(grange[2], g))
                              } else {
                                ga_start <- ga$start
                              }
                            }
                            mod1 <- laGP::newGPsep(X=self$X, Z=self$Z, d=da_start, g=ga_start, dK = TRUE)
                            #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
                            if (no_update) { #using d and g given
                              self$mod.extra$theta = as.numeric(1 / d) # store theta params
                              self$mod.extra$nugget = as.numeric(g) # store nugget
                            } else if (!no_update && !self$estimate.nugget) { # Only estimate d/theta
                              mle.out <- laGP::mleGPsep(gpsepi = mod1,
                                                        param="d",
                                                        tmin=drange[1], tmax=drange[2],
                                                        verb=0, maxit=1000)
                              self$mod.extra$theta <- as.numeric(1 / mle.out$d) # store theta params
                              self$mod.extra$nugget <- g
                              # Leave nugget as is
                            } else if (!no_update) { # Update both
                              mle.out <- laGP::jmleGPsep(gpsepi = mod1,
                                                         drange=drange,
                                                         grange=grange,
                                                         #dab=da$ab, gab=ga$ab, # Will use MLE without these
                                                         verb=0, maxit=1000)
                              self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
                              self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
                            } else {stop("Shouldn't be here IGP_laGP #32097555")}
                            self$mod <- mod1
                          }, #"function to initialize model with data
                          .update = function(..., no_update=FALSE) {
                            if (no_update) { # just add data and return
                              laGP::updateGPsep(gpsepi=self$mod,
                                                X=self$X[-(1:self$n.at.last.update), , drop=FALSE],
                                                Z=self$Z[-(1:self$n.at.last.update)])
                              return()
                            }
                            
                            # Start over if not many points, had problems getting stuck in bad spots early
                            if (self$n.at.last.update < 20) {
                              self$.delete()
                              self$.init(...)
                              return()
                            }
                            
                            da <- laGP::darg(list(mle=TRUE), X=self$X)
                            # Had problem with garg that didn't make sense, trying to fix by adding noise to Z
                            ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = TRUE)
                            if (inherits(ga.try, "try-error")) {
                              print("Adding noise to Z so garg doesn't give error")
                              noise.sd <- 1e-8
                              while (TRUE) {
                                ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z + rnorm(length(self$Z), 0, noise.sd)), silent = TRUE)
                                if (!inherits(ga.try, "try-error")) {break}
                                noise.sd <- 10 * noise.sd
                              }
                            }
                            n.since.last.update <- nrow(self$X) - self$n.at.last.update
                            if (n.since.last.update < 1) {
                              warning("Can't update, no new X rows, but can optimize again")
                            } else {
                              if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
                                # start over if too many
                                self$.delete(...=...)
                                self$.init(...=...)
                              } else {
                                lagpupdate.try <- try(
                                  laGP::updateGPsep(gpsepi=self$mod,
                                                    X=self$X[-(1:self$n.at.last.update), , drop=FALSE],
                                                    Z=self$Z[-(1:self$n.at.last.update)])
                                )
                                if (inherits(lagpupdate.try, "try-error")) {
                                  stop("Error in lagpupdate.try #8257")
                                }
                              }
                            }
                            drange <- c(1e-3,1e4)
                            grange <- c(min(sqrt(.Machine$double.eps),self$mod.extra$nugget), max(1,self$mod.extra$nugget))
                            if (no_update) {
                              stop("This is covered above at beginning of .update")
                            } else if (!no_update && !self$estimate.nugget) { # update d/theta but not nugget/g
                              mle.try <- try(mle.out <- laGP::mleGPsep(gpsepi = self$mod,
                                                                       param = "d",
                                                                       tmin=drange[1], tmax=drange[2],
                                                                       verb=0, maxit=1000))
                              if (inherits(mle.try, "try-error")) {
                                # Sometimes gives error: L-BFGS-B needs finite values of 'fn'
                                
                                warning('Restarting laGP model after mle error #40297')
                                self$delete()
                                self$init(..., no_update=no_update)
                                return()
                              }
                              # Update stored parameters for when user calls $theta() or $nugget()
                              self$mod.extra$theta = as.numeric(1 / mle.out$d) # store theta params
                              # leave nugget as it was
                            } else if (!no_update && self$estimate.nugget) {
                              mle.try <- try(mle.out <- laGP::jmleGPsep(gpsepi = self$mod,
                                                                        #drange=c(da$min, da$max), # Getting rid of these here too
                                                                        #grange=c(ga$min, ga$max),
                                                                        drange=drange,
                                                                        grange=grange, # Had error of nugget starting outside bound
                                                                        #dab=da$ab, gab=ga$ab,
                                                                        verb=0, maxit=1000))
                              if (inherits(mle.try, "try-error")) {
                                # Sometimes gives error: L-BFGS-B needs finite values of 'fn'
                                
                                warning('Restarting laGP model after jmle error #19378')
                                self$delete()
                                self$init(...)
                                return()
                              }
                              # Update stored parameters for when user calls $theta() or $nugget()
                              self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
                              self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
                            }
                          }, #"function to add data to model or reestimate params
                          .predict = function(XX, se.fit, ...){
                            if (se.fit) {
                              preds <- laGP::predGPsep(self$mod, XX, lite=TRUE)
                              # Sometimes preds$s2 is negative
                              numneg <- sum(preds$s2<0)
                              if (numneg > 0) {
                                if (nrow(XX) - numneg >= 5) {
                                  newmin <- min(preds$s2[preds$s2 > 0])
                                  preds$s2 <- pmax(newmin, preds$s2)
                                  warning(paste0("laGP gave ",numneg," s2 preds < 0, setting them to min pos s2 pred of ", newmin))
                                } else {
                                  warning(paste0("laGP gave ",numneg," s2 preds < 0, setting them to have pos s2 of ", 1e-16))
                                  preds$s2 <- pmax(1e-16, preds$s2)
                                }
                              }
                              list(fit=preds$mean, se.fit=sqrt(preds$s2))
                            } else {
                              laGP::predGPsep(self$mod, XX, lite=TRUE)$mean
                            }
                          }, #"function to predict at new values
                          .predict.se = function(XX, ...) {
                            # sqrt(laGP::predGPsep(self$mod, XX, lite=TRUE)$s2)
                            
                            # laGP can give neg s2 values, problem with sqrt, so check for it
                            sqrt(self$.predict.var(XX=XX, ...))
                          }, #"function predict the standard error/dev
                          .predict.var = function(XX, ...) {
                            # laGP::predGPsep(self$mod, XX, lite=TRUE)$s2
                            
                            # laGP can give neg s2 values, problem with sqrt, so check for it
                            s2 <- laGP::predGPsep(self$mod, XX, lite=TRUE)$s2
                            numneg <- sum(s2<0)
                            if (numneg > 0) {#print("Fixing laGP var")
                              if (nrow(XX) - numneg >= 5) {
                                newmin <- min(s2[s2 > 0])
                                s2 <- pmax(newmin, s2)
                                warning(paste0("laGP gave ",numneg," s2 preds < 0, setting them to min pos s2 pred of ", newmin))
                              } else {
                                warning(paste0("laGP gave ",numneg," s2 preds < 0, setting them to have pos s2 of ", 1e-16))
                                s2 <- pmax(1e-16, s2)
                              }
                            }
                            s2
                            
                          }, #"function to predict the variance
                          .grad = NULL, # function to calculate the gradient
                          .delete = function(...) {
                            if (!is.null(self$mod)) {
                              laGP::deleteGPsep(self$mod)
                              self$mod <- NULL
                            }
                          }, #"function to delete model beyond simple deletion
                          .theta = function() {self$mod.extra$theta}, #"function to get theta, exp(-theta*(x-x)^2)
                          .nugget = function() {self$mod.extra$nugget}, #"function to get nugget
                          .s2 = function() {self$predict.var(rep(1e4, ncol(self$X)))},
                          .mean = function() {0} # function that gives mean (constant, other functions not implemented)
                          
                        )
)




#' IGP R6 object for fitting laGP_GauPro_kernel model
#'
#' Class providing object with methods for fitting a GP model.
#' This mixes laGP and GauPro with a kernel. It fits the model using laGP,
#' then copies the parameters to a GauPro model for prediction.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_laGP_GauPro_kernel$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' c(u$mod.extra$laGP$theta(), u$mod.extra$laGP$nugget())
#' c(u$mod.extra$GauPro$theta(), u$mod.extra$GauPro$nugget())
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/IGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{update(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...)}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_laGP_GauPro_kernel <- R6::R6Class(
  classname = "IGP_laGP_GauPro_kernel",
  inherit = IGP_base,
  public = list(
    .init = function(...) {
      # Fit model to data with laGP
      self$mod.extra$laGP <- IGP(X=self$X, Z=self$Z, package="laGP",
                                 corr=self$corr,
                                 estimate.nugget=self$estimate.nugget,
                                 nugget0=self$nugget0, ...)
      #self$mod.extra$laGP$init(X=self$X, Z=self$Z, ...)
      
      # Copy params to GauPro, don't fit, use this for predicting
      kern <- GauPro::Gaussian$new(D=ncol(self$X), beta=log(self$mod.extra$laGP$theta(),10),
                                   s2=self$mod.extra$laGP$s2())
      self$mod.extra$GauPro <- IGP(X=self$X, Z=self$Z, package="GauPro_kernel",
                                   corr=self$corr,
                                   estimate.nugget=FALSE,
                                   nugget0=self$mod.extra$laGP$nugget(),
                                   kernel=kern,
                                   # theta=self$mod.extra$laGP$theta(),
                                   #nug=self$mod.extra$laGP$nug,
                                   param.est=FALSE)
      # laGPs2 <- self$mod.extra$laGP$s2()
      # self$mod.extra$GauPro$mod$kernel$s2 <- laGPs2
      # self$mod.extra$GauPro$mod$kernel$logs2 <- log(laGPs2, 10)
      # self$mod.extra$GauPro$mod$s2_hat <- laGPs2
      #self$mod.extra$GauPro$init(X=self$X, Z=self$Z,
      #                          theta=self$mod.extra$laGP$theta,
      #                           nug=self$mod.extra$laGP$nug)
      self$mod <- "laGP model is mod.extra$laGP, GauPro model is mod.extra$GauPro. This fits with laGP but predicts with GauPro"
    }, #"function to initialize model with data
    .update = function(..., no_update=FALSE) {
      # Update model in laGP
      if (!no_update) { # won't have it update data even not updating params since I don't like when it gives issues, and we pass Xall anyways
        self$mod.extra$laGP$update(Xall=self$X, Zall=self$Z,
                                   no_update=no_update,
                                   ...)
      }
      # Pass GauPro new theta and nugget if it was updated
      if (!no_update) {
        self$mod.extra$GauPro$mod$kernel$beta <- log(self$mod.extra$laGP$theta(), 10)
        self$mod.extra$GauPro$mod$nug <- self$mod.extra$laGP$nugget()
        laGPs2 <- self$mod.extra$laGP$s2()
        self$mod.extra$GauPro$mod$kernel$s2 <- laGPs2
        self$mod.extra$GauPro$mod$kernel$logs2 <- log(laGPs2, 10)
        self$mod.extra$GauPro$mod$s2_hat <- laGPs2
      }
      self$mod.extra$GauPro$update(Xall=self$X, Zall=self$Z,
                                   no_update=TRUE)
      
    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      self$mod.extra$GauPro$.predict(XX=XX, se.fit=se.fit, ...)
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {
      self$mod.extra$GauPro$.predict.se(XX=XX, ...)
    }, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {
      self$mod.extra$GauPro$.predict.var(XX=XX, ...)
    }, #"function to predict the variance
    .grad = function(XX) {self$mod.extra$GauPro$grad(XX=XX)}, # function to calculate the gradient
    .delete = function(...){
      if (!is.null(self$mod.extra)) {
        self$mod.extra$laGP$delete()
        self$mod.extra$GauPro$delete()
        self$mod.extra <- NULL
      }
      self$mod <- NULL
    }, #"function to delete model beyond simple deletion
    .theta = function() {self$mod.extra$GauPro$theta()}, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = function() {self$mod.extra$GauPro$nugget()}, #"function to get nugget
    .s2 = function() {self$mod.extra$GauPro$s2()},
    .mean = function() {self$mod.extra$GauPro$trend$m} # function that gives mean (constant, other functions not implemented)
    
  )
)






#' Class providing object with methods for adapt.concept2.sFFLHD.R6
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @importFrom stats optim
#' @keywords data, experiments, adaptive, sequential, simulation,
#' Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for running an 
#' adaptive experiment.
#' @format \code{\link{R6Class}} object.
#' @examples
#' a <- adapt.concept2.sFFLHD.R6$new(D=2,L=3,func=gaussian1,obj="desirability",
#'     des_func=des_func14, n0=12, take_until_maxpvar_below=.9, 
#'     package="GauPro", design='sFFLHD', selection_method="max_des_red")
#' a$run(5)
#' @field X Design matrix
#' @field Z Responses
#' @field b batch size
#' @field func Actual function to get experiment values from
#' @field nb Number of batches, if you know before starting
#' @field D Dimension of data
#' @field Xopts Available points
#' @field X0 Initial design
#' @field package Which GP package to use in IGP
#' @field stats List of tracked stats
#' @field iteration Which iteration
#' @field mod The GP model from 
#' @field func_run_together Whether points should be passed to func separately
#'          as vectors or all together as a matrix whose rows are the points.
#' @field func_fast If the function is fast. If TRUE then full plots are made.
#'          In practice this is alway FALSE.
#' @field new_batches_per_batch How many batches of candidate points are added
#'          for each batch taken.
#' @field X_tracker data.frame tracking the points of X, such as when they were
#'          selected.
#' @field X0 An initial matrix of points to be used.
#' @field Xopts A matrix of candidate (option) points.
#' @field Xopts_tracker A data.frame tracking the points of Xopts.
#' @field batch.tracker Tracks when points were added to Xopts.
#' @field Xopts_removed A matrix of points removed from Xopts.
#' @field s The design object for generating candidate points.
#' @field design A string saying which design object should be used.
#' @field stats A data.frame giving stats for each iteration.
#' @field iteration The current iteration.
#' @field obj A string saying what the objective is.
#' @field obj_func A function for the objective.
#' @field n0 The initial number of points to be selected.
#' @field take_until_maxpvar_below A number, if the proportion of points near 
#'          the maximum variance of the GP model, then it will take 
#'          space-filling points.
#' @field package Which GP package should be used by IGP.
#' @field force_old A number saying how often the oldest candidate points
#'          should be forced into the design.
#' @field force_pvar A number saying how often the points with the highest
#'          predictive variance should be forced into the design.
#' @field des_func The desirability function.
#' @field des_func_fast Whether the des_func is fast for candidate points.
#' @field alpha_des The alpha constant for the weight function.
#' @field actual_des_func The true des func used to evaluate the model,
#' not known in practice
#' @field weight_const The weight constant in the weight function, usually 1.
#' @field selection_method What the selection method should be.
#' @field parallel Should new values be calculated in parallel?
#' @field verbose How much detail should be printed to the console. 0 is
#'          minimal, 1 is medium, 2 is a lot.
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to
#'   https://github.com/CollinErickson/DOE-Code}
#'   \item{\code{new(X, Z, corr="Gauss", verbose=0, separable=T, 
#'   useC=F,useGrad=T,
#'          parallel=T, nug.est=T, ...)}}{This method is used to create object
#'          of this class with \code{X} and \code{Z} as the data.}
#'
#'   \item{\code{update(Xnew=NULL, Znew=NULL, Xall=NULL, Zall=NULL,
#' restarts = 5,
#' param_update = T, nug.update = self$nug.est)}}{This method updates the 
#' model, adding new data if given, then running optimization again.}
#'   }
adapt.concept2.sFFLHD.R6 <- R6::R6Class(classname = "adapt.concept2.sFFLHD.seq",
                                        public = list(
                                          func = NULL, # "function", 
                                          func_run_together = NULL, # Should the matrix of values to be run be passed
                                          #  to func as a matrix or by row?, useful if you
                                          #  parallelize your own function or call another
                                          #  program to get actual values
                                          func_fast = NULL, # Is the func super fast so actual MSE can be calculated?
                                          D = NULL, # "numeric", 
                                          L = NULL, # sFFLHD batch size, probably the same as b, or number of points
                                          #   design gives when taking a single batch
                                          b = NULL, # batch size to add each iteration, probably the same as L
                                          new_batches_per_batch = NULL,
                                          g = NULL, # "numeric", # g not used but I'll leave it for now
                                          X = NULL, # "matrix", Z = "numeric", Xopts = "matrix",
                                          X_tracker = NULL, # tracks points that are in X
                                          X0 = NULL,
                                          Xopts = NULL,
                                          Xopts_tracker = NULL, # Keep track of data about candidate points
                                          batch.tracker = NULL, # tracks when Xopts were added
                                          Xopts_removed = NULL,
                                          Z = NULL,
                                          s = NULL, # "sFFLHD" an object with $get.batch to get batch of points
                                          design = NULL,
                                          stats = NULL, # "list", 
                                          iteration = NULL, # "numeric",
                                          obj = NULL, # "character", 
                                          obj_func = NULL, # "function",
                                          obj_nu = NULL,
                                          n0 = NULL, # "numeric"
                                          stage1batches = NULL, # numeric, number of stage 1 batches
                                          take_until_maxpvar_below = NULL, 
                                          package = NULL, # "character",
                                          force_old = NULL, # "numeric", 
                                          force_pvar = NULL, # "numeric",
                                          useSMEDtheta = NULL, # "logical"
                                          mod = NULL,
                                          #desirability_func = NULL, # args are mod and XX, this was the full 
                                          # weighted error function, poorly named
                                          #actual_desirability_func = NULL, # 
                                          des_func = NULL, # desirability function: args are mod and XX, should be
                                          #   the delta desirability function, output from 0 to 1
                                          des_func_fast = NULL, # If des func is slow (using model, not true), then
                                          #it won't be plotted for other points, eg contour plot
                                          alpha_des = NULL,
                                          actual_des_func = NULL,
                                          actual_werror_func = NULL,
                                          #weight_func = NULL, # weight function: 1 + alpha_des * des_func()
                                          weight_const = NULL,
                                          #werror_func = NULL, # weighted error function: 
                                          #    sigmahat * (1+alpha_des*des_func())
                                          error_power = NULL, # 
                                          selection_method = NULL, # string
                                          nconsider = NULL,
                                          nconsider_random = NULL,
                                          
                                          parallel = NULL, # Should the new values be calculated in parallel?
                                          #   Not for the model, for getting actual new Z values
                                          parallel_cores = NULL, # Number of cores used for parallel
                                          parallel_cluster = NULL, # The object for the cluster currently running
                                          
                                          options = NULL, # A list for holding other things that aren't worth giving 
                                          #   own variable
                                          verbose = NULL, # 0 prints only essential, 2 prints a lot
                                          
                                          initialize = function(D,L,b=NULL, package=NULL, obj=NULL,
                                                                n0=0, stage1batches=NULL,
                                                                force_old=0, force_pvar=0,
                                                                useSMEDtheta=F, 
                                                                func, func_run_together=FALSE, func_fast=TRUE,
                                                                take_until_maxpvar_below=NULL,
                                                                design="sFFLHD",
                                                                selection_method, X0=NULL, Xopts=NULL,
                                                                des_func, des_func_fast=TRUE, alpha_des=1,
                                                                new_batches_per_batch=5,
                                                                parallel=FALSE, parallel_cores="detect",
                                                                nugget=1e-6, estimate.nugget = TRUE,
                                                                verbose = 1,
                                                                design_seed=numeric(0),
                                                                weight_const=0,
                                                                error_power=1,
                                                                nconsider=Inf, nconsider_random=0,
                                                                ...) {
                                            self$iteration <- 1
                                            self$D <- D
                                            self$L <- L
                                            self$b <- if (is.null(b)) L else b
                                            self$new_batches_per_batch <- new_batches_per_batch
                                            self$func <- if (is.function(func)) func 
                                            else if (is.character(func)) get(func)
                                            else stop("Bad func #6293")
                                            self$func_run_together <- func_run_together
                                            self$func_fast <- func_fast
                                            self$force_old <- force_old
                                            self$force_pvar <- force_pvar
                                            self$take_until_maxpvar_below <- take_until_maxpvar_below
                                            self$selection_method <- selection_method
                                            self$des_func_fast <- des_func_fast
                                            self$weight_const <- weight_const
                                            self$error_power <- error_power
                                            if (is.null(error_power) || !(error_power %in% c(1,2))) {
                                              stop("error_power must be 1 or 2")
                                            }
                                            self$verbose <- verbose
                                            self$nconsider <- nconsider
                                            self$nconsider_random <- nconsider_random
                                            
                                            if (any(length(D)==0, length(L)==0)) {
                                              message("D and L must be specified")
                                            }
                                            
                                            self$design <- design
                                            if (self$design == "sFFLHD") {
                                              self$s <- sFFLHD::sFFLHD(D=D, L=L, maximin=T, seed=design_seed)
                                            } else if (self$design == "sFFLHD_Lflex") {
                                              self$s <- sFFLHD::sFFLHD_Lflex$new(D=D, L=L, maximin=T,
                                                                                 seed=design_seed)
                                            } else if (self$design == "random") {
                                              self$s <- random_design$new(D=D, L=L, use_lhs=FALSE, seed=design_seed)
                                            } else if (self$design == "lhd") {
                                              self$s <- random_design$new(D=D, L=L, use_lhs=TRUE, seed=design_seed)
                                            } else if (self$design == "sobol") {
                                              self$s <- sobol_design$new(D=D, L=L, seed=design_seed)
                                            } else if (self$design == "given") { # This means Xopts is given in and
                                              # no new points will be added to design
                                              self$s <- NULL
                                            } else {
                                              stop(paste("Design <", self$design,"> isn't recognized #3285729"))
                                            }
                                            self$X0 <- X0
                                            self$X <- matrix(NA,0,D)
                                            if (is.null(Xopts)) {
                                              self$Xopts <- matrix(NA,0,D)
                                            } else { # Option to give in Xopts
                                              self$Xopts <- Xopts
                                              self$Xopts_tracker_add(Xopts)
                                            }
                                            self$Xopts_removed <- matrix(NA,0,D)
                                            
                                            if(is.null(package)) {self$package <- "laGP"}
                                            else {self$package <- package}
                                            self$mod <- IGP(package = self$package, estimate.nugget=estimate.nugget, 
                                                            nugget=nugget)
                                            self$stats <- list(iteration=c(),n=c(),pvar=c(),mse=c(), ppu=c(), 
                                                               minbatch=c(), pamv=c(), actual_intwerror=c(),
                                                               actual_intwvar=c(), intwerror=c(), intwvar=c(),
                                                               intwerror01=c())
                                            self$obj_nu <- NaN
                                            
                                            # set objective function according to obj
                                            self$obj <- obj
                                            if (is.null(self$obj)) {stop("Must give in obj")}
                                            if (self$obj == "mse") { # The default
                                              self$obj_func <- mod$predict.var
                                            } else if (self$obj == "maxerr") {
                                              self$obj_func <- function(lims) {
                                                maxgridfunc(self$mod$predict.var, lims=lims, batch=T)
                                              }
                                            } else if (self$obj == "grad") {
                                              self$obj_func <- self$mod$grad_norm
                                            } else if (self$obj == "func") {
                                              self$obj_func <- function(xx) pmax(1e-16, self$mod$predict(xx))
                                            } else if (self$obj == "pvar") {
                                              self$obj_func <- function(xx) pmax(1e-16, self$mod$predict.var(xx))
                                            } else if (self$obj == "gradpvarnu") {
                                              self$obj_func <- function(xx) {
                                                if (is.nan(self$obj_nu)) { # if not defined yet, set obj_nu so the two
                                                  #   are balanced
                                                  XXX <- matrix(runif(1e3*self$D), ncol=self$D)
                                                  gn_max  <- max(self$mod$grad_norm(XXX))
                                                  pse_max <- max(self$mod$predict.se(XXX))
                                                  self$obj_nu <- gn_max / pse_max
                                                }
                                                1           *      self$mod$grad_norm(xx) + 
                                                  self$obj_nu *      pmax(1e-16, self$mod$predict.se(xx))
                                              }
                                            } else if (self$obj == "nonadapt") {
                                              # use next batch only
                                            } else if (self$obj %in% c("desirability", "des")) {
                                              self$obj <- "desirability"
                                              if (missing(des_func)) {
                                                stop("Must give in des_func when using desirability")
                                              }
                                              self$des_func <- des_func
                                              # obj_func is used by SMED, give it the weight function
                                              self$obj_func <- function(XX) {self$weight_func(mod=self$mod, XX=XX)}
                                              if (missing(alpha_des)) {stop("alpha_des must be given in")}
                                              self$alpha_des <- alpha_des
                                              if (is.character(self$des_func)) {
                                                if (grepl(pattern="\\(", x=self$des_func)) { # If parentheses, then
                                                  self$des_func <- eval(parse(text=self$des_func))
                                                } else { # Just a function name in quote, eg "des_func_relmax"
                                                  self$des_func <- get(self$des_func)
                                                }
                                              }
                                            }
                                            
                                            # This can be used even when not using desirability in order to
                                            #   make comparisons
                                            if ('actual_des_func' %in% names(list(...))) {
                                              self$actual_des_func <- list(...)$actual_des_func
                                              if (is.character(self$actual_des_func)) {
                                                if (grepl(pattern="\\(", x=self$actual_des_func)) {
                                                  # If parentheses, then
                                                  self$actual_des_func <- eval(parse(text=self$actual_des_func))
                                                } else { # Just a function name in quote, eg "des_func_relmax"
                                                  self$actual_des_func <- get(self$actual_des_func)
                                                }
                                              }
                                            }
                                            if ('actual_intwerror_func' %in% names(list(...))) {
                                              stop("Don't do it like this")
                                            }
                                            if (is.null(self$alpha_des) && !missing(alpha_des)){
                                              self$alpha_des <- alpha_des
                                            }
                                            
                                            self$n0 <- n0
                                            if (F && !is.null(self$X0)) {
                                              self$X <- self$X0
                                              self$Z <- c(self$Z, apply(self$X,1,self$func))
                                              self$mod$update(Xall=self$X, Zall=self$Z)
                                            }
                                            self$stage1batches <- stage1batches
                                            
                                            self$useSMEDtheta <- if (length(useSMEDtheta)==0) {FALSE} 
                                            else {useSMEDtheta}
                                            
                                            # Set up parallel stuff
                                            self$parallel <- parallel
                                            if (self$parallel) {
                                              # Use a list to store info about parallel, such as num nodes, cluster
                                              if (parallel_cores == "detect") {
                                                self$parallel_cores <- parallel::detectCores()
                                              } else {
                                                self$parallel_cores <- parallel_cores
                                              }
                                              # For now assume using parallel package
                                              self$parallel_cluster <- parallel::makeCluster(
                                                spec = self$parallel_cores, type = "SOCK")
                                            }
                                          },
                                          run = function(maxit, plotlastonly=F, noplot=F) {
                                            # Run multiple iterations
                                            i <- 1
                                            while(i <= maxit) {
                                              if (self$verbose >= 1) {
                                                cat(paste('Starting iteration', self$iteration, "at", Sys.time(), "\n"))
                                              }
                                              iplotit <- ((i == maxit) | !plotlastonly) & !noplot
                                              self$run1(plotit=iplotit)
                                              i <- i + 1
                                            }
                                            invisible(self)
                                          },
                                          run1 = function(plotit=TRUE) {
                                            # Run single iteration
                                            if (is.null(self$s)) { # If no design s, then we can only add points when
                                              #  we have enough left, so check to make sure there are at least b left
                                              if (nrow(self$Xopts) + nrow(self$Xopts_removed) < self$b) {
                                                stop("Not enough points left to get a batch #82389, 
                                                     initial design not big enough, b reached")
                                              }
                                              }
                                            self$add_data()
                                            self$update_mod()
                                            self$update_stats()
                                            if (plotit) {
                                              self$plot1()
                                            }
                                            #set_params()
                                            self$iteration <- self$iteration + 1
                                            invisible(self)
                                            },
                                          add_data = function() {
                                            # newL will be the L points selected from Xopts
                                            #   to add to the design
                                            newL <- NULL # Indices of Xopts rows to add
                                            reason <- NA # Reason for selecting those points
                                            
                                            # First check to see if X hasn't been initialized yet
                                            if (nrow(self$X) == 0 ) {
                                              if (!is.null(self$X0)) { # If X0, use it
                                                self$add_newL_points_to_design(newL=NULL, use_X0=TRUE,
                                                                               reason="X0 given")
                                                return()
                                              } else if (!is.null(self$n0) && self$n0 > 0) {
                                                # Take first batches up to n0 and use it
                                                self$add_new_batches_to_Xopts(
                                                  num_batches_to_take = ceiling(self$n0/self$L))
                                                newL <- 1:self$n0
                                                reason <- "Taking first n0 from Xopts since X is empty"
                                              } else { # no X0 or n0, so take first L
                                                self$add_new_batches_to_Xopts(num_batches_to_take = 1)
                                                newL <- 1:self$b
                                                reason <- "Taking first b from Xopts since X is empty"
                                              }
                                              
                                              #
                                            } else if (!is.null(self$stage1batches) && 
                                                       self$iteration <= self$stage1batches) {
                                              cat("stage1batch adding\n")
                                              self$add_new_batches_to_Xopts(num_batches_to_take = 1)
                                              newL <- 1:self$b
                                              reason <- "Taking next b since still stage 1"
                                              # If nonadaptive, just take first L from design
                                            } else if (self$selection_method %in% c("nonadapt", "noadapt")) {
                                              self$add_new_batches_to_Xopts(num_batches_to_take = 1)
                                              newL <- 1:self$b
                                              reason <- "Taking next b since nonadapt"
                                              # If variance is too high across surface, take points
                                            } else if (!is.null(self$take_until_maxpvar_below) && 
                                                       self$mod$prop.at.max.var(val=self$take_until_maxpvar_below) > 0.1) {
                                              #print(paste("Taking until pvar lower: ", 
                                              #      self$mod$prop.at.max.var(val=self$take_until_maxpvar_below)))
                                              if (self$package == 'GauPro') {
                                                cat(paste("Taking until pvar lower: ", 
                                                          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                                                          "   avg t_LOO: ",
                                                          mean(abs(self$mod$mod$pred_LOO(se.fit=T)$t)),
                                                          '\n'))
                                              } else if (self$package == 'laGP_GauPro') {
                                                cat(paste("Taking until pvar lower: ", 
                                                          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                                                          "   avg t_LOO: ",
                                                          mean(abs(
                                                            self$mod$mod.extra$GauPro$mod$pred_LOO(se.fit=T)$t)),
                                                          '\n'))
                                              } else {
                                                cat(paste("Taking until pvar lower: ", 
                                                          self$mod$prop.at.max.var(val=self$take_until_maxpvar_below),
                                                          '\n'))
                                              }
                                              if (FALSE) { # Take oldest points, not doing it
                                                Xnew <- self$s$get.batch()
                                                Znew <- apply(Xnew, 1, self$func)
                                                self$X <- rbind(self$X, Xnew)
                                                self$Z <- c(self$Z, Znew)
                                                return()
                                              } else { # Instead of taking old trying to take space filling
                                                self$add_new_batches_to_Xopts(1)
                                                Xdesign <- self$X
                                                newL <- c()
                                                for (ell in 1:self$b) {
                                                  mindistsq <- apply(self$Xopts, 1,
                                                                     function(xvec) {
                                                                       min(rowSums(sweep(Xdesign, 2, xvec)^2))
                                                                     })
                                                  whichmaxmin <- which.max(mindistsq)
                                                  newL <- c(newL, whichmaxmin)
                                                  Xdesign <- rbind(Xdesign, self$Xopts[whichmaxmin,])
                                                }
                                                reason <- "pvar high so taking maximin dist"
                                              }
                                            }
                                            
                                            # Add new points
                                            
                                            # Add new batches if newL haven't already been selected
                                            if (is.null(newL)) {
                                              self$add_new_batches_to_Xopts()
                                            }
                                            
                                            # Check if forcing old or pvar
                                            # Returns NULL if not selecting, otherwise the L indices
                                            if (is.null(newL)) {
                                              newL <- self$select_new_points_from_old_or_pvar()
                                              if (!is.null(newL)) {reason <- "Taking from old or pvar"}
                                            }
                                            
                                            # If nothing forced, use selection method to get points
                                            if (is.null(newL)) {
                                              if (self$selection_method %in% c("SMED","SMED_true")) {
                                                # standard min energy
                                                newL <- self$select_new_points_from_SMED()
                                                reason <- "SMED"
                                              } else if (self$selection_method %in% c("max_des", "max_des_all",
                                                                                      "max_des_all_best", "ALM",
                                                                                      "ALM_all", "ALM_all_best")) {
                                                # take point with max desirability, update model, requires using se
                                                #   or pvar so adding a point goes to zero
                                                # newL <- self$select_new_points_from_max_des()
                                                # Moved this into des_red even though it isn't a reduction
                                                newL <- self$select_new_points_from_max_des_red()
                                                reason <- "max_des"
                                              } else if (self$selection_method %in% 
                                                         c("max_des_red", "max_des_red_all", "max_des_red_all_best",
                                                           "ALC", "ALC_all", "ALC_all_best")
                                              ) {
                                                # take maximum reduction, update model, requires using se or pvar 
                                                #   so adding a point goes to zero
                                                newL <- self$select_new_points_from_max_des_red()
                                                reason <- "max_des_red or _all or _all_best"
                                              }
                                            }
                                            
                                            self$add_newL_points_to_design(newL = newL, reason=reason)
                                          },
                                          update_obj_nu = function(Xnew, Znew) {
                                            if (is.null(self$mod$X)) {return(rep(NA, nrow(Xnew)))}
                                            if (is.nan(self$obj_nu)) return()
                                            if (is.nan(self$obj_nu)) { # Initialize it intelligently
                                              browser()
                                              self$obj_nu <- .5
                                            }
                                            Zlist <- self$mod$predict(Xnew, se.fit=T)
                                            Zmean <- Zlist$fit
                                            Zse   <- Zlist$se
                                            abs.scores <- abs(Znew - Zmean) / Zse
                                            for (score in abs.scores) {
                                              if (score < 3 && score > .001) { # If score is too close to zero than
                                                # something is wrong? Maybe not, but don't want to reward models
                                                # that just have huge error everywhere
                                                self$obj_nu <- .5 * self$obj_nu
                                              } else {
                                                self$obj_nu <- 2  * self$obj_nu
                                              }
                                            }
                                            print(paste('alpha changed to ', self$obj_nu))
                                          },
                                          update_mod = function() {
                                            # Update GP model for data
                                            self$mod$update(Xall=self$X, Zall=self$Z)
                                          },
                                          set_params = function() {
                                          },
                                          update_stats = function() {
                                            # Keep stats of progress over course of experiment
                                            # self$stats$ <- c(self$stats$, )
                                            self$stats$iteration <- c(self$stats$iteration, self$iteration)
                                            self$stats$n <- c(self$stats$n, nrow(self$X))
                                            self$stats$pvar <- c(self$stats$pvar,
                                                                 msfunc(self$mod$predict.var,
                                                                        cbind(rep(0,self$D),rep(1,self$D))))
                                            self$stats$mse <- c(self$stats$mse, self$mse_func()) 
                                            #msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D))))
                                            self$stats$ppu <- c(self$stats$ppu,
                                                                nrow(self$X) / (nrow(self$X) + nrow(self$Xopts)))
                                            self$stats$minbatch <- c(self$stats$minbatch,
                                                                     if (length(self$batch.tracker>0)) 
                                                                       min(self$batch.tracker) 
                                                                     else 0)
                                            self$stats$pamv <- c(self$stats$pamv, self$mod$prop.at.max.var())
                                            # self$stats$actual_intwerror <- c(self$stats$actual_intwerror,
                                            #                                  self$actual_intwerror_func())
                                            aiwef <- self$actual_intwerror_func(error_power=c(1,2))
                                            self$stats$actual_intwerror <- c(self$stats$actual_intwerror, aiwef[[1]])
                                            self$stats$actual_intwvar <- c(self$stats$actual_intwvar, aiwef[[2]])
                                            if (!is.null(self$des_func)) {
                                              # self$stats$intwerror <- c(self$stats$intwerror, self$intwerror_func())
                                              iwf <- self$intwerror_func(error_power=c(1,2))
                                              self$stats$intwerror <- c(self$stats$intwerror, iwf[[1]])
                                              self$stats$intwvar <- c(self$stats$intwvar, iwf[[2]])
                                              self$stats$intwerror01 <- c(self$stats$intwerror01, NaN) 
                                              # Not using now, should be sped up anyways 
                                              #  by doing at the same time as intwerror, 
                                              #  bad to do it this way 
                                              #  self$intwerror_func(weight_const=0,alpha=1))
                                            } else {
                                              self$stats$intwerror <- c(self$stats$intwerror, NaN)
                                              self$stats$intwvar <- c(self$stats$intwvar, NaN)
                                              self$stats$intwerror01 <- c(self$stats$intwerror01, NaN)
                                            }
                                          },
                                          mse_func = function() {
                                            if (self$func_fast) {
                                              msecalc(self$func,self$mod$predict,cbind(rep(0,self$D),rep(1,self$D)))
                                            } else {
                                              NaN
                                            }
                                          },
                                          plot_mean = function(cex=1, plot.axes=TRUE) {
                                            cf_func(self$mod$predict,batchmax=500, pretitle="Predicted Mean ",
                                                    cex=cex, plot.axes=plot.axes,
                                                    afterplotfunc=function(){
                                                      points(self$X,pch=19)
                                                      if (self$iteration > 1) { # Add points just chosen with yellow
                                                        points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],
                                                               col='yellow',pch=19, cex=.5)} # plot last L separately
                                                    }
                                            )
                                          },
                                          plot_se = function(cex=1, plot.axes=TRUE) {
                                            cf_func(self$mod$predict.se,batchmax=500, pretitle="Predicted SE ",
                                                    cex=cex, plot.axes=plot.axes,
                                                    afterplotfunc=function(){
                                                      points(self$X,pch=19)
                                                      if (self$iteration > 1) { # Plot last L separately
                                                        points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),],
                                                               col='yellow',pch=19, cex=.5)
                                                      }
                                                      points(self$Xopts, col=2); # add points not selected
                                                    }
                                            )
                                          },
                                          plot_abserr = function(cex=1, plot.axes=TRUE) {
                                            cf_func(function(xx){sqrt((
                                              self$mod$predict(xx) - apply(xx, 1, self$func))^2)},
                                              n = 20, mainminmax_minmax = F, pretitle="AbsErr ", batchmax=Inf,
                                              cex=cex, plot.axes=plot.axes)
                                          },
                                          plot_mse = function(statsdf, cex=1) { # Plot MSE and PVar over iterations
                                            par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                            if (missing(statsdf)) {
                                              print("missing statsdf in plot_mse")
                                              statsdf <- as.data.frame(self$stats)
                                            }
                                            plot(rep(statsdf$iter,2), c(statsdf$mse,statsdf$pvar), 
                                                 type='o', log="y", col="white",
                                                 xlab="Iteration", ylab=""
                                            )
                                            legend("topright",legend=c("MSE","PVar"),fill=c(1,2), cex=cex)
                                            points(statsdf$iter, statsdf$mse, type='o', pch=19)
                                            points(statsdf$iter, statsdf$pvar, type='o', pch = 19, col=2)
                                          },
                                          plot_iwe = function(statsdf, cex=1) {
                                            par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                            if (missing(statsdf)) {
                                              print("missing statsdf in plot_iwe")
                                              statsdf <- as.data.frame(self$stats)
                                            }
                                            plot(rep(statsdf$iter,2), c(statsdf$actual_intwerror,statsdf$intwerror), 
                                                 type='o', log="y", col="white",
                                                 xlab="Iteration", ylab=""
                                            )
                                            legend("topright",legend=c("IWE","PIWE"),fill=c(1,2), cex=cex)
                                            points(statsdf$iter, statsdf$actual_intwerror, type='o', pch=19)
                                            points(statsdf$iter, statsdf$intwerror, type='o', pch = 19, col=2)
                                          },
                                          plot_ppu = function(statsdf, cex) {
                                            # Plot percentage of points used over iteration
                                            if (missing(statsdf)) {
                                              print("missing statsdf in plot_ppu")
                                              statsdf <- as.data.frame(self$stats)
                                            }
                                            par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                            plot(statsdf$iter, statsdf$ppu, type='o', pch=19,
                                                 xlab="Iteration")
                                            legend('bottomleft',legend="% pts",fill=1, cex=cex)
                                          },
                                          plot_des_v_acc = function(cex, cex.axis) {
                                            Xplot <- matrix(runif(self$D*100), ncol=self$D)
                                            Xplot_des <- self$des_func(XX=Xplot, mod=self$mod)
                                            Xplot_se <- self$mod$predict.se(Xplot)
                                            
                                            # If func_fast plot des vs se and des vs abserror
                                            if (self$func_fast) {
                                              Xplot_abserror <- abs(self$mod$predict(Xplot) - 
                                                                      apply(Xplot, 1, self$func))
                                              plot(NULL, xlim=c(min(Xplot_des), max(Xplot_des)), 
                                                   ylim=c(min(Xplot_abserror, Xplot_se),
                                                          max(Xplot_abserror, Xplot_se)),
                                                   pch=19, xlab='SE', ylab='Des', cex.axis=cex.axis)#, log='xy')
                                              legend(x = 'topright', legend=c('SE', 'AbsErr'), fill=c(1,2), cex=cex)
                                              points(Xplot_des, Xplot_se, pch=19, col=1)
                                              points(Xplot_des, Xplot_abserror, pch=19, col=2)
                                            } else { # Only plot des vs se
                                              plot(Xplot_des, Xplot_se, pch=19, xlab='SE', ylab='Grad', 
                                                   cex.axis=cex.axis)#, log='xy')
                                            }
                                          },
                                          plot_y_acc = function(residual=FALSE) {
                                            # Plot predicted vs actual with error bars
                                            
                                            # If func_fast, then get for other random points
                                            if (self$func_fast) { # Only do these if fast
                                              Xplot <- matrix(runif(self$D*50), ncol=self$D)
                                              Zplot.pred.all <- self$mod$predict(Xplot, se.fit=TRUE)
                                              Zplot.pred <- Zplot.pred.all$fit
                                              Zplot.se <- Zplot.pred.all$se
                                              Zplot.act <- apply(Xplot,1, self$func)
                                            } else {
                                              Zplot.pred <- c()
                                              Zplot.act <- c()
                                            }
                                            
                                            # Get predictions for points in design
                                            Zused.pred.all <- self$mod$predict(self$X, se.fit=TRUE)
                                            Zused.pred <- Zused.pred.all$fit
                                            Zused.se <- Zused.pred.all$se
                                            
                                            # Make it residuals if residual
                                            if (residual) {
                                              Zused.pred <- Zused.pred - self$Z
                                              if (self$func_fast) {Zplot.pred <- Zplot.pred - Zplot.act}
                                            }
                                            
                                            # Blank plot with right values
                                            if (self$func_fast) {
                                              plot(NULL, xlim=c(min(self$Z, Zplot.act), max(self$Z, Zplot.act)), 
                                                   ylim=c(min(Zused.pred-2*Zused.se, Zplot.pred-2*Zplot.se), 
                                                          max(Zused.pred+2*Zused.se, Zplot.pred+2*Zplot.se)),
                                                   xlab="Z actual", ylab="Z predicted 95%")
                                              legend(x = 'topleft', legend=c("Z", "ZZ"), col = c(2,1), pch=19)
                                            } else {
                                              plot(NULL, xlim=c(min(self$Z), max(self$Z)), 
                                                   ylim=c(min(Zused.pred-2*Zused.se), 
                                                          max(Zused.pred+2*Zused.se)),
                                                   xlab="Z actual", ylab="Z predicted 95%")
                                              legend(x = 'topleft', legend=c("Z", "ZZ"), col = c(2,1), pch=19)
                                            }
                                            
                                            # If fast, then plot values for random points
                                            if (self$func_fast) {
                                              for (i in 1:length(Zplot.se)) {
                                                lines(c(Zplot.act[i],Zplot.act[i]),
                                                      Zplot.pred[i] + 2 * Zplot.se[i] * c(1, -1), col=3)
                                              }
                                            }
                                            for (i in 1:length(Zused.se)) {
                                              lines(c(self$Z[i],self$Z[i]),
                                                    Zused.pred[i] + 2 * Zused.se[i] * c(1, -1), col=4)
                                            }
                                            if (residual) {abline(a=0, b=0)} else {abline(a = 0, b = 1)}
                                            if (self$func_fast) {points(Zplot.act, Zplot.pred, xlab="Z",
                                                                        ylab="Predicted", pch=19)}
                                            points(self$Z, Zused.pred, col=2, pch=19)
                                          },
                                          plot_1D = function() {
                                            x <- matrix(seq(0,1,l=300), ncol=1)
                                            preds <- self$mod$predict(x, se=T)
                                            ylim <- c(min(preds$fit-2*preds$se), max(preds$fit+2*preds$se))
                                            plot(x, preds$fit, ylim=ylim, type='l', col='white', ylab="Z")
                                            points(self$X, rep(ylim[1], length(self$X)))
                                            multicolor.title(c("Actual ","Pred ", "95% ", "Weight ", "Weight*sd"),
                                                             c(3,1,2,6,"cyan"))
                                            if (self$des_func_fast) {
                                              xdes <- self$des_func(mod=self$mod, XX=x)
                                              xdes2 <- ((xdes-min(xdes))/(max(xdes)-min(xdes))) * 
                                                (ylim[2]-ylim[1])*.2 + ylim[1] - .04*diff(ylim)
                                              xwd <- xdes * preds$se
                                              xwd2 <- ((xwd-min(xwd))/(max(xwd)-min(xwd))) * (ylim[2]-ylim[1])*.2 + 
                                                ylim[1] - .04*diff(ylim)
                                              points(x, xdes2, type='l', col=6, lwd=.5)
                                              points(x, xwd2, type='l', col="cyan", lwd=.5)
                                              if (nrow(self$X) > 1) {
                                                dens <- density(self$X)
                                                ydens <- dens$y
                                                ydens2 <- ((ydens-min(ydens))/(max(ydens)-min(ydens))) * 
                                                  (ylim[2]-ylim[1])*.2 + ylim[1] - .04*diff(ylim)
                                                points(dens$x, ydens2, type='l', col="orange", lwd=.5)
                                              }
                                            }
                                            points(x, preds$fit-2*preds$se, type='l', col=2, lwd=2)
                                            points(x, preds$fit+2*preds$se, type='l', col=2, lwd=2)
                                            points(x, preds$fit, type='l', col=1, lwd=3)
                                            if (self$func_fast) {
                                              Zopts <- apply(self$Xopts, 1, self$func)
                                              points(self$Xopts, Zopts, col=4, pch=19)
                                              points(x, apply(x, 1, self$func), type='l', col=3, lwd=3)
                                            }
                                            points(self$X, self$Z, pch=19, cex=2)
                                            if (self$iteration > 1) { # Plot last L separately
                                              points(self$X[(nrow(self$X)-self$b+1):nrow(self$X),], 
                                                     self$Z[(nrow(self$X)-self$b+1):nrow(self$X)],
                                                     col='yellow',pch=19, cex=.5)
                                            }
                                          },
                                          plot_2D = function(twoplot = FALSE, cex=1) {
                                            cex_small = .55 * cex
                                            # twoplot only plots mean and se
                                            if (twoplot) { # Only plot pred surface and pred error
                                              split.screen(matrix(c(0,.5,0,1,.5,1,0,1),byrow=T, ncol=4))
                                              screen(1)
                                              self$plot_mean(cex=cex)
                                              screen(2)
                                              self$plot_se(cex=cex)
                                              close.screen(all=TRUE)
                                              return()
                                            }
                                            
                                            #par(mfrow=c(2,1))
                                            ln <- 5 # number of lower plots
                                            split.screen(matrix(
                                              c(0,.5,.25,1,  .5,1,.25,1,
                                                0,1/ln,0,.25, 1/ln,2/ln,
                                                0,.25, 2/ln,3/ln,
                                                0,.25, 3/ln,4/ln,
                                                0,.25, 4/ln,1,0,.25),
                                              ncol=4,byrow=T))
                                            
                                            # Plot mean
                                            screen(1)
                                            self$plot_mean(cex=cex)
                                            
                                            # Plot se predictions
                                            screen(2)
                                            self$plot_se(cex=cex)
                                            
                                            # Only plot true func if func_fast
                                            if (self$func_fast) {
                                              screen(3) # Actual func
                                              par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                              cf_func(self$func, n = 20, mainminmax_minmax = F, pretitle="Actual ",
                                                      cex=cex_small, plot.axes=FALSE)
                                            }
                                            
                                            # Plot MSE if past iteration 1
                                            if (self$iteration >= 2) {
                                              statsdf <- as.data.frame(self$stats)
                                              screen(4) # MSE plot
                                              self$plot_mse(statsdf=statsdf, cex=cex_small)
                                            }
                                            
                                            if (self$des_func_fast) {
                                              screen(5) # plot des
                                              if (self$des_func_fast && !is.null(self$des_func)) {
                                                # Option to not plot if it is slow
                                                cf_func(function(XX) {self$des_func(XX=XX, mod=self$mod)}, 
                                                        n=20, mainminmax_minmax = F, pretitle="Des ", 
                                                        cex=cex_small, plot.axes=FALSE, batchmax=Inf)
                                              }
                                            }
                                            
                                            if (self$iteration >= 2) {
                                              screen(6) # % of pts used plot 
                                              par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                              if (self$des_func_fast && self$obj == "desirability") {
                                                self$plot_des_v_acc(cex=cex_small, cex.axis = cex_small)
                                              } else {
                                                self$plot_ppu(statsdf=statsdf, cex=cex_small)
                                              }
                                              
                                            }
                                            if (self$func_fast) {
                                              screen(7) # actual squared error plot
                                              # par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                              self$plot_abserr(cex=.55*cex, plot.axes=FALSE)
                                            }       
                                            close.screen(all = TRUE)
                                          },
                                          plot1 = function(twoplot=FALSE, cex=1) {
                                            if (self$D == 1) {
                                              self$plot_1D()
                                            } else if (self$D == 2) {
                                              self$plot_2D(twoplot=twoplot, cex=cex)
                                            } else { # D != 2 
                                              oparmar <- par('mar')
                                              oparmfrow <- par('mfrow')
                                              par(mfrow=c(2,2))
                                              par(mar=c(2,2,0,0.5)) # 5.1 4.1 4.1 2.1 BLTR
                                              statsdf <- as.data.frame(self$stats)
                                              #print(ggplot(statsdf, aes(x=iteration, y=mse, col=level)) + geom_line())
                                              #print(ggplot() + 
                                              #        geom_line(data=statsdf, aes(x=iteration, y=mse, col="red")) + 
                                              #        geom_line(data=statsdf, aes(x=iteration, y=pvar, col="blue"))
                                              #)
                                              if (self$iteration >= 2) {
                                                # 1 mse plot
                                                self$plot_iwe(statsdf=statsdf, cex=cex)
                                                
                                                # 2 yhat vs y plot
                                                self$plot_y_acc(residual=TRUE)
                                                
                                                # 3 % pts used plot
                                                self$plot_ppu(statsdf=statsdf, cex=cex)
                                                
                                                # 4 grad vs pvar
                                                if (!is.null(self$des_func)) {
                                                  self$plot_des_v_acc(cex=cex, cex.axis=cex)
                                                }
                                              }
                                              par(mar=oparmar)
                                              par(mfrow=oparmfrow)
                                            }
                                          },
                                          add_new_batches_to_Xopts = function(
                                            num_batches_to_take=self$new_batches_per_batch) {
                                            if (is.null(self$s)) { # If all options are given by user,
                                              #  don't add new points
                                              return()
                                            }
                                            for (iii in 1:num_batches_to_take) {
                                              Xnew <- self$s$get.batch()
                                              self$Xopts <- rbind(self$Xopts, Xnew)
                                              self$batch.tracker <- c(self$batch.tracker, rep(self$s$b, nrow(Xnew)))
                                              self$Xopts_tracker_add(Xnew) 
                                              #self$Xopts_tracker <- rbind(self$Xopts_tracker, 
                                              #                            self$Xopts_tracker_add(Xnew))
                                            }
                                          },
                                          Xopts_tracker_add = function(Xnew) {
                                            n <- nrow(Xnew)
                                            Xnewdf <- data.frame(iteration_added=rep(self$iteration, n),
                                                                 time_added = rep(Sys.time(), n))
                                            # if (self$obj %in% c("desirability","des")) {
                                            #   if (self$selection_method == "max_des_red") {
                                            #     
                                            #   }
                                            # }
                                            self$Xopts_tracker <- rbind(self$Xopts_tracker, Xnewdf)
                                          },
                                          Xopts_tracker_remove = function(newL) {
                                            # newL is index of pt to remove from Xopts
                                            # Remove from Xopts_tracker, add to X_tracker
                                            removed_rows <- self$Xopts_tracker[newL,, drop=FALSE]
                                            #Zp <- self$predict(newX, se.fit=TRUE)
                                            # self$X_tracker <- rbind(self$X_tracker, newX)
                                            self$Xopts_tracker <- self$Xopts_tracker[-newL,, drop=FALSE]
                                            removed_rows
                                          },
                                          select_new_points_from_old_or_pvar = function() {
                                            newL <- NULL
                                            # Check if forcing old or pvar
                                            if (self$force_old > 0 & self$force_pvar > 0) {
                                              stop("No can force_old and force_pvar")
                                            } else if (self$force_old > 0 & self$force_old <= 1) {
                                              rand1 <- runif(1)
                                              if (rand1 < self$force_old) {newL <- 1:self$b} 
                                            } else if (self$force_old > 1) {
                                              if ((iteration %% as.integer(self$force_old)) == 0) {
                                                newL <- 1:self$b
                                              }
                                            } else if (self$force_pvar > 0 & self$force_pvar <= 1) {
                                              rand1 <- runif(1)
                                              if (rand1 < self$force_pvar) {
                                                newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]
                                              } 
                                            } else if (self$force_pvar > 1) {
                                              if ((iteration %% as.integer(self$force_pvar)) == 0) {
                                                newL <- order(self$mod$predict.var(self$Xopts), decreasing=T)[1:self$b]
                                              }
                                            }
                                            newL
                                          },
                                          select_new_points_from_SMED = function() {
                                            #bestL <- SMED_selectC(f=self$obj_func, n=self$b, 
                                            #                      X0=self$X, Xopt=self$Xopts, 
                                            #                      theta=if (self$useSMEDtheta) {self$mod$theta()}
                                            #                      else {rep(1,2)})
                                            if (self$selection_method == "SMED") {
                                              Yall.try <- try(Yall <- self$obj_func(rbind(self$X, self$Xopts)))
                                            } else if (self$selection_method == "SMED_true") {
                                              Yall.try <- try(Yall <- apply(rbind(self$X, self$Xopts), 1, self$func))
                                            } else {
                                              stop("no SMED #35230")
                                            }
                                            if (inherits(Yall.try, "try-error")) {
                                              browser()
                                              Yall <- self$obj_func(rbind(self$X, self$Xopts))
                                            }
                                            Y0 <- Yall[1:nrow(self$X)]
                                            Yopt <- Yall[(nrow(self$X)+1):length(Yall)]
                                            bestL <- SMED_selectYC(n=self$b, X0=self$X, Xopt=self$Xopts, Y0=Y0,
                                                                   Yopt=Yopt,
                                                                   theta=if (self$useSMEDtheta) {self$mod$theta()}
                                                                   else {rep(1,ncol(self$X))})
                                            newL <- bestL
                                            newL
                                          },
                                          select_new_points_from_max_des = function() {
                                            # take point with max desirability, update model, requires using se or 
                                            #   pvar so adding a point goes to zero
                                            # ALM is active learning Cohn, just picks highest pvar (equiv to se)
                                            gpc <- self$mod$clone(deep=TRUE)
                                            bestL <- c()
                                            for (ell in 1:self$b) {
                                              if (self$selection_method == "ALM") {
                                                objall <- self$mod$predict.se(mod=gpc, XX=rbind(self$X, self$Xopts))
                                              } else {
                                                objall <- self$werror_func(mod=gpc, XX=rbind(self$X, self$Xopts))
                                              }
                                              objopt <- objall[(nrow(self$X)+1):length(objall)]
                                              objopt[bestL] <- -Inf # ignore the ones just selected
                                              bestopt <- which.max(objopt)
                                              bestL <- c(bestL, bestopt)
                                              if (ell < self$b) {
                                                Xnewone <- self$Xopts[bestopt, , drop=FALSE]
                                                Znewone = gpc$predict(Xnewone)
                                                if (self$verbose >= 2) {
                                                  print(Xnewone);print(Znewone);
                                                }
                                                #cf(function(xx) self$desirability_func(gpc, xx),
                                                #  batchmax=1e3, pts=self$Xopts)
                                                gpc$update(Xnew=Xnewone, Znew=Znewone, restarts=0, no_update=TRUE)
                                              }
                                            }
                                            newL <- bestL
                                            rm(gpc, objall, objopt, bestopt, bestL)
                                            newL
                                          },
                                          select_new_points_from_max_des_red = function() {
                                            # Use max weighted error reduction to select batch of points from self$Xopts
                                            # Returns indices of points to use from Xopts
                                            
                                            # gpc is a temp version of the model to add points to
                                            #  now always a GauPro_kernel since it has fastest methods
                                            if (self$package == 'laGP_GauPro_kernel') {
                                              gpc <- self$mod$mod.extra$GauPro$clone(deep=TRUE)
                                            } else if (self$package == "GauPro_kernel") {
                                              gpc <- self$mod$clone(deep=TRUE)
                                            } else {
                                              gpc <- IGP::IGP_GauPro_kernel$new(
                                                X=self$X, Z=self$Z,
                                                kernel=GauPro::Gaussian$new(D=self$D,
                                                                            s2=self$mod$s2(), s2_est=FALSE,
                                                                            beta=log(self$mod$theta(),10),
                                                                            beta_est=F),
                                                trend=GauPro::trend_c$new(D=self$D,
                                                                          m=self$mod$mean(),
                                                                          m_est=FALSE),
                                                no_update=TRUE, nugget=self$mod$nugget(), estimate.nugget=FALSE)
                                            }
                                            # Get indices of points to consider, take most recent
                                            # Xopts_to_consider <- 1:nrow(self$Xopts)
                                            if (self$nconsider[1] < nrow(self$Xopts)) {
                                              Xopts_to_consider <- 1:self$nconsider[1] + nrow(self$Xopts) - 
                                                self$nconsider[1]
                                            } else {
                                              Xopts_to_consider <- 1:nrow(self$Xopts)
                                            }
                                            # Add back in some older points randomly
                                            numrandtoadd <- self$nconsider_random[1]
                                            if (numrandtoadd > 0 &&
                                                length(setdiff(1:nrow(self$Xopts), Xopts_to_consider)) > numrandtoadd) {
                                              Xopts_to_consider <- c(Xopts_to_consider, 
                                                                     sample(setdiff(1:nrow(self$Xopts), 
                                                                                    c(Xopts_to_consider, bestL)), 
                                                                            numrandtoadd, F))
                                            }
                                            # Plot contour function of weighted error function
                                            if (self$D == 2 && self$verbose > 1) {
                                              dontplotfunc <- TRUE
                                              if (dontplotfunc) {
                                                split.screen(matrix(
                                                  c(0,1/2,0,1, 1/2,1,0,1),
                                                  ncol=4,byrow=T))
                                                screen(1)
                                                cf(function(X) {self$werror_func(mod=gpc, XX=X)}, 
                                                   batchmax=Inf, 
                                                   afterplotfunc=function(){
                                                     points(self$X, col=3, pch=2);
                                                     points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
                                                     text(self$Xopts[Xopts_to_consider,])
                                                   },
                                                   main=expression(omega(x)*hat(delta)(x) * "  before")
                                                )
                                              } else {
                                                split.screen(matrix(
                                                  c(0,1/3,0,1, 1/3,2/3,0,1, 2/3,1,0,1),
                                                  ncol=4,byrow=T))
                                                screen(1)
                                                cf(self$mod$predict, batchmax=Inf, pts=self$X)
                                                screen(2)
                                                cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf, 
                                                   afterplotfunc=function(){points(self$X, col=3, pch=2);
                                                     points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
                                                     text(self$Xopts[Xopts_to_consider,])})
                                              }
                                            }
                                            
                                            # Can start with none and select one at time, 
                                            #   or start with random and replace
                                            # TODO make variable for none, all, and best
                                            if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
                                              bestL <- c() # Start with none
                                            } else if (self$selection_method %in% c("max_des_red_all", "ALC_all",
                                                                                    "max_des_all", "ALM_all")) {
                                              # Start with random L and replace
                                              bestL <- sample(Xopts_to_consider, size = self$b, replace = FALSE)
                                            } else if (self$selection_method %in% 
                                                       c("max_des_red_all_best", "ALC_all_best",
                                                         "max_des_all_best", "ALM_all_best")) {
                                              # Start with best L and replace
                                              if (self$selection_method %in% c("ALC_all_best", "ALM_all_best")) {
                                                print("using ALC_all_best")
                                                Xotc_werrors <- self$werror_func(
                                                  XX = self$Xopts[Xopts_to_consider, , drop=FALSE], 
                                                  des_func=function(XX, mod){rep(0, nrow(XX))}, 
                                                  alpha=0, weight_const=1) #, weight_func=self$weight_func)
                                              } else if (self$selection_method %in% c("max_des_red_all_best",
                                                                                      "max_des_all_best")) {
                                                Xotc_werrors <- self$werror_func(
                                                  XX = self$Xopts[Xopts_to_consider, , drop=FALSE])
                                              } else {stop("#92847")}
                                              bestL <- order(Xotc_werrors, decreasing = TRUE)[self$b:1] 
                                              # Make the biggest last so it is least likely to be replaced
                                            } else {
                                              browser("Selection method doesn't match up #92352583")
                                            }
                                            
                                            uses_ALM <- self$selection_method %in% c("ALM", "ALM_all", "ALM_all_best",
                                                                                     "max_des", "max_des_all",
                                                                                     "max_des_all_best")
                                            # TODO put line above earlier, and make sure selection method is valid
                                            if (uses_ALM) { # max_des or ALM
                                              browser()
                                              add_points_weights <- self$weight_func(XX=self$Xopts)
                                              # TODO rename int_werrors_red_func to obj_to_max
                                              int_werrors_red_func <- function(add_points_indices) {
                                                browser()
                                                add_points <- self$Xopts[add_points_indices, ]
                                                if (self$error_power==2) {
                                                  pp <- gpc$predict.var(XX=add_points)
                                                } else {
                                                  pp <- gpc$predict.se(XX=add_points)
                                                }
                                                if (substr(self$selection_method, 1, 3) != "ALM") {
                                                  pp <- pp * add_points_weights[add_points_indices]
                                                }
                                                pp
                                              }
                                            } else {
                                              # Random integration points from simple LHS
                                              int_points <- simple.LHS(1e4, self$D)
                                              
                                              # Make separate int_werror_func for ALC
                                              # if (substr(self$selection_method, 1, 3) == "ALC") {print("Using ALC")
                                              # 
                                              # int_werror_func <- function() {
                                              #   mean(
                                              #     self$werror_func(XX=int_points, mod=gpc, 
                                              #                      des_func=function(XX, mod){rep(0, nrow(XX))},
                                              #                      alpha=0, weight_const=1)
                                              #   )
                                              # }
                                              # } else { # Not ALC, so max_des_red
                                              # There can be alot of variability in calculating the desirability
                                              #   when it involves sampling stuff, so the intwerror values will
                                              #   fluctuate if you recalculate each time. And that is slower.
                                              if (substr(self$selection_method, 1, 3) == "ALC") {print("Using ALC")
                                                int_points_numdes <- rep(1, nrow(int_points))
                                              } else {
                                                int_points_numdes <- self$des_func(XX=int_points, mod=gpc)
                                              }
                                              
                                              # Set function to calculate int_werrors_reduction
                                              if (self$error_power == 1) {
                                                int_werrors_red_func <- function(add_points_indices) {
                                                  # New faster, same results, version
                                                  add_points <- self$Xopts[add_points_indices, , drop=FALSE]
                                                  # Calculate pred vars after adding points
                                                  pvaaps <- gpc$mod$pred_var_after_adding_points_sep(
                                                    add_points=add_points, pred_points=int_points)
                                                  # Some will be a little less than 0, gives NaN for sqrt
                                                  sum_neg <- sum(c(pvaaps)<0)
                                                  # print(eigen(self$mod$mod.extra$GauPro$mod$K, symmetric = T,
                                                  #             only.values = T
                                                  #             )$val %>% {c(min(.), max(.), max(.)/min(.))})
                                                  if (sum_neg > 0) {
                                                    cat("    pvaaps: ",sum_neg,"/",length(pvaaps),
                                                        "are negative, setting to zero", "\n")
                                                  }
                                                  pvaaps <- pmax(pvaaps, 0)
                                                  # Need negative since it isn't reduction, it is total value
                                                  -colMeans(sweep(sqrt(pvaaps), 1, 
                                                                  (self$weight_const+self$alpha_des*int_points_numdes), `*`))
                                                }
                                              } else if (self$error_power == 2) {
                                                int_werrors_red_func <- function(add_points_indices) {
                                                  # New faster, same results, version
                                                  add_points <- self$Xopts[add_points_indices, , drop=FALSE]
                                                  # mean((self$weight_const+self$alpha_des*int_points_numdes)*
                                                  #  gpc$mod$pred_var_reductions(add_points=add_points, 
                                                  #                          pred_points=int_points))
                                                  
                                                  colMeans(sweep(gpc$mod$pred_var_reductions(
                                                    add_points=add_points, pred_points=int_points), 1, 
                                                    (self$weight_const+self$alpha_des*int_points_numdes), `*`))
                                                }
                                              } else {stop("Error power must be 1 or 2 #282362")}
                                            }
                                            # }
                                            # X_with_bestL <- self$X
                                            # Z_with_bestL <- self$Z
                                            Znotrun_preds <- gpc$predict(self$Xopts) # Need to use the predictions
                                            #  before each is added
                                            for (ell in 1:self$b) {
                                              if (self$verbose >= 2) {
                                                cat(paste0('starting iter ', ell,'/',self$b, ', considering ',
                                                           length(unique(Xopts_to_consider,bestL)), "/", 
                                                           nrow(self$Xopts), ', bestL is ', 
                                                           paste0(bestL, collapse = ' '), '\n'))
                                              }
                                              
                                              # The surrogate values
                                              if (exists("use_true_for_surrogates") && use_true_for_surrogates) {
                                                print("cheating")
                                                Znotrun_preds <- apply(self$Xopts, 1, self$func)
                                              }
                                              
                                              int_werror_vals <- rep(Inf, nrow(self$Xopts))
                                              # if (self$selection_method %in% c("max_des_red", "ALC")) {
                                              # Don't have current value, so don't start with anything
                                              # r_star <- NA # Track best index
                                              # int_werror_vals_star <- Inf # Track best value
                                              # } else { # Start with ell and replace
                                              # r_star <- bestL[ell]
                                              # if (ell == 1) { # First time need to calculate current integrated des
                                              # X_with_bestL <- rbind(X_with_bestL, self$Xopts[bestL, , drop=F])
                                              # Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[bestL])
                                              # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                                              #            no_update=TRUE)
                                              # int_werror_vals_star <- int_werror_func()
                                              # } else { # After that it just stays as star value
                                              # essentially int_des_weight_star <- int_des_weight_star
                                              # }
                                              # Need to remove ell of bestL from X since it is being replaced
                                              # Store current selection in IDWs,
                                              #   but not actually using it for anything
                                              # int_werror_vals[bestL[ell]] <- int_werror_vals_star 
                                              # X_with_bestL <- rbind(self$X, self$Xopts[bestL[-ell], , drop=F])
                                              # Z_with_bestL <- c(self$Z, Znotrun_preds[bestL[-ell]])
                                              # }
                                              
                                              # for (r in setdiff(Xopts_to_consider, bestL)) {
                                              #   # gpc$update(Xall = rbind(X_with_bestL, self$Xopts[r, ,drop=F]),
                                              #   #            Zall=c(Z_with_bestL, Znotrun_preds[r]), restarts=0,
                                              #   #            no_update=TRUE)
                                              # 
                                              #   int_werror_vals_r <- int_werror_func()
                                              #   if (int_werror_vals_r < int_werror_vals_star) {
                                              #     int_werror_vals_star <- int_werror_vals_r
                                              #     r_star <- r
                                              #   }
                                              #   int_werror_vals[r] <- int_werror_vals_r
                                              # }
                                              
                                              if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
                                                # if starting with none and adding one
                                                gpc$update(Xall=rbind(self$X, self$Xopts[bestL,,drop=F]),
                                                           Zall=c(self$Z, Znotrun_preds[bestL]),
                                                           no_update=TRUE)
                                                # Consider all except bestL
                                                Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL)
                                              } else { # if starting with L and replacing as go
                                                # Remove one under consideration for replacement
                                                gpc$update(Xall=rbind(self$X, self$Xopts[bestL[-ell],,drop=F]),
                                                           Zall=c(self$Z, Znotrun_preds[bestL[-ell]]),
                                                           no_update=TRUE)
                                                # Consider all except bestL, add back in one that might be replaced
                                                Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL[-ell])
                                              }
                                              # browser()
                                              # Xopts_inds_to_consider <- setdiff(Xopts_to_consider, bestL)
                                              int_werror_red_vals <- int_werrors_red_func(
                                                add_points_indices = Xopts_inds_to_consider)
                                              
                                              r_star <- Xopts_inds_to_consider[which.max(int_werror_red_vals)]
                                              
                                              # print("Here are int_werror_vals")
                                              if (self$verbose >= 2) {
                                                # print(cbind(1:length(int_werror_red_vals), int_werror_red_vals))
                                                print(cbind(Xopts_inds_to_consider, int_werror_red_vals))
                                              }
                                              
                                              # Reduce the number to consider if large
                                              if (ell < self$b) {
                                                numtokeep <- self$nconsider[min(length(self$nconsider), ell+1)] + 
                                                  1 - self$b # b-1 selected that aren't in consideration
                                                Xopts_to_consider <- order(int_werror_vals,
                                                                           decreasing = F)[1:min(length(
                                                                             int_werror_vals), numtokeep)]
                                              }
                                              
                                              # Add back in some random ones
                                              numrandtoadd <- self$nconsider_random[
                                                min(length(self$nconsider_random), ell+1)]
                                              if (numrandtoadd > 0 &&
                                                  length(setdiff(1:nrow(self$Xopts),
                                                                 Xopts_to_consider)) > numrandtoadd) {
                                                Xopts_to_consider <- c(Xopts_to_consider,
                                                                       sample(setdiff(1:nrow(self$Xopts),
                                                                                      c(Xopts_to_consider, bestL)),
                                                                              numrandtoadd, F))
                                              }
                                              
                                              if (self$selection_method %in% c("max_des_red", "ALC", "max_des", "ALM")) {
                                                # if starting with none and adding one
                                                bestL <- c(bestL, r_star)
                                              } else { # if starting with L and replacing as go
                                                if (length(r_star) != 1) {stop("Error in choosing r_star #95882")}
                                                bestL[ell] <- r_star
                                              }
                                              
                                              if (ell < self$b || TRUE) { # REMOVE THIS FOR SPEED
                                                Xnewone <- self$Xopts[r_star, , drop=FALSE]
                                                Znewone <- Znotrun_preds[r_star] #gpc$predict(Xnewone)
                                                # if (self$selection_method %in% c("max_des_red", "ALC")) {
                                                # X_with_bestL <- rbind(X_with_bestL, Xnewone)
                                                # Z_with_bestL <- c(Z_with_bestL, Znewone)
                                                # if (T) { # REMOVE THIS FOR SPEED
                                                #   gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                                                #              no_update=TRUE)
                                                # }
                                                # } else {
                                                # X_with_bestL <- rbind(X_with_bestL, self$Xopts[r_star, ])
                                                # Z_with_bestL <- c(Z_with_bestL, Znotrun_preds[r_star])
                                                # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, restarts=0,
                                                # no_update=TRUE)
                                                # }
                                                if (self$verbose >= 2) {
                                                  cat("\tSelected", r_star, Xnewone, Znewone, "\n")
                                                }
                                              }
                                            }
                                            if (self$verbose >= 2) {
                                              cat("Selected:", bestL, "\n")
                                            }
                                            # If verbose, plot
                                            if (self$D == 2 && self$verbose >= 2) {
                                              # gpc$update(Xall=X_with_bestL, Zall=Z_with_bestL, no_update=TRUE)
                                              gpc$update(Xall=rbind(self$X, self$Xopts[bestL,,drop=F]),
                                                         Zall=c(self$Z, Znotrun_preds[bestL]), no_update=TRUE)
                                              if (dontplotfunc) {
                                                screen(2)
                                                cf(function(X) {self$werror_func(mod=gpc, XX=X)}, 
                                                   batchmax=Inf, 
                                                   afterplotfunc=function(){
                                                     points(self$X, col=3, pch=2);
                                                     points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
                                                     points(self$Xopts[Xopts_to_consider,]);
                                                     points(self$Xopts[bestL,], col=1,pch=19, cex=2);
                                                     text(self$Xopts[bestL,], col=2,pch=19, cex=2)
                                                   },
                                                   main=expression(omega(x)*hat(delta)(x) * "  after")
                                                )
                                              } else {
                                                screen(3)
                                                cf(function(X)self$werror_func(mod=gpc, XX=X), batchmax=Inf,
                                                   afterplotfunc=function(){points(self$X, col=3, pch=2);
                                                     points(self$Xopts[-Xopts_to_consider,], col=4,pch=3);
                                                     points(self$Xopts[Xopts_to_consider,]);
                                                     points(self$Xopts[bestL,], col=1,pch=19, cex=2);
                                                     text(self$Xopts[bestL,], col=2,pch=19, cex=2)})
                                              }
                                              close.screen(all=TRUE)
                                            }
                                            
                                            # Return bestL
                                            bestL
                                          },
                                          int_werror_after_adding = function(Xnew=NULL, Znew=NULL, Xnew_Xoptsrow=NULL, 
                                                                             n=1e4, int_points=NULL, 
                                                                             seed=NULL,
                                                                             ...
                                          ) {
                                            if (!is.null(seed)) {set.seed(seed)}
                                            if (is.null(int_points)) {
                                              int_points = simple.LHS(n=n, d=self$D)
                                            }
                                            if (is.null(Xnew)) {
                                              if (!is.null(Xnew_Xoptsrow)) {
                                                Xnew <- self$Xopts[Xnew_Xoptsrow, ]
                                              } else {
                                                return("Need to give Xnew or Xnew_Xoptsrow")
                                              }
                                            }
                                            if (is.null(Znew)) {
                                              Znew = self$mod$predict(Xnew)
                                            }
                                            mod <- IGP::IGP(X=rbind(self$X, Xnew), Z=c(self$Z, Znew), 
                                                            package = "GauPro",
                                                            theta=self$mod$theta(), param.est=FALSE,
                                                            set.nugget=self$mod$nugget(), estimate.nugget=FALSE)
                                            mn <- mean(self$werror_func(XX=int_points, mod=mod, ...))
                                            mn
                                          },
                                          add_newL_points_to_design = function(newL=NULL, use_X0=FALSE, reason) {
                                            if (use_X0) { # If X0 given and first iter, add them
                                              Xnew <- self$X0
                                              removed_tracker_rows <- data.frame(iteration_added=0,
                                                                                 time_added=Sys.time())
                                            } else { # Else newL must be given
                                              if (length(newL) != self$b) { 
                                                if (length(newL) != self$n0  || nrow(self$X)!=0) {
                                                  browser()
                                                  stop("Selected newL not of length L #84274")
                                                }
                                              }
                                              removed_tracker_rows <- self$Xopts_tracker_remove(newL=newL)
                                              Xnew <- self$Xopts[newL, , drop=FALSE]
                                              self$Xopts <- self$Xopts[-newL, , drop=FALSE]
                                              self$batch.tracker <- self$batch.tracker[-newL]
                                            }
                                            Znew <- self$calculate_Z(Xnew)
                                            if (any(duplicated(rbind(self$X,Xnew)))) {browser()}
                                            self$X <- rbind(self$X,Xnew)
                                            self$Z <- c(self$Z,Znew)
                                            # Track points added
                                            pred <- if (nrow(self$X) == length(newL) || use_X0) { # Model not fit yet
                                              fakelen <- if (use_X0) {nrow(Xnew)} else {length(newL)}
                                              data.frame(fit=rep(NA, fakelen), se.fit=rep(NA, fakelen))
                                              # data.frame(fit=rep(NA, length(newL)), 
                                              #               se.fit=rep(NA, length(newL)))
                                            } else{
                                              self$mod$predict(Xnew, se.fit=TRUE)
                                            }
                                            tracker_rows <- data.frame(
                                              iteration_added_to_opts=removed_tracker_rows$iteration_added, 
                                              time_added_to_opts=removed_tracker_rows$time_added,
                                              iteration_added = self$iteration,
                                              time_added = Sys.time(),
                                              Z=Znew, Zpred=pred$fit, sepred=pred$se.fit,
                                              t=(Znew-pred$fit)/pred$se.fit, reason=reason)
                                            self$X_tracker <- rbind(self$X_tracker, tracker_rows)
                                            
                                            self$update_obj_nu(Xnew=Xnew, Znew=Znew)
                                          },
                                          calculate_Z = function(X) {
                                            # Used to just be apply(self$X, 1, self$func)
                                            if (self$parallel && inherits(self$parallel_cluster, "cluster")) {
                                              # parallel::clusterApply(cl = self$parallal_cluster, x = 1:nrow(X))
                                              parallel::parRapply(cl = self$parallel_cluster, x = X, self$func)
                                            } else if (self$func_run_together) {
                                              self$func(X)
                                            } else {
                                              apply(X, 1, self$func)
                                            }
                                          },
                                          # The weight function 1 + alpha * delta()
                                          weight_func = function(..., XX, mod=self$mod, des_func=self$des_func,
                                                                 alpha=self$alpha_des, weight_const=self$weight_const) {
                                            if (is.function(des_func)) {
                                              weight_const + alpha * des_func(XX=XX, mod=mod)
                                            } else if (is.numeric(des_func)) { 
                                              weight_const + alpha * des_func
                                            } else {
                                              browser("Shouldn't be here error #132817585")
                                            }
                                          },
                                          # The weighted error function sigmahat * (1 + alpha * delta())
                                          werror_func = function(..., XX, mod=self$mod, des_func=self$des_func,
                                                                 alpha=self$alpha_des, weight_const=self$weight_const,
                                                                 weight_func=self$weight_func,
                                                                 error_power=self$error_power) {
                                            err <- mod$predict.se(XX)
                                            if (exists("use_true_for_error") && use_true_for_error) {
                                              if (runif(1) < .01) print("Using true error #9258332")
                                              err <- abs(mod$predict(XX) - apply(XX, 1, self$func))
                                            }
                                            weight_func_out <- weight_func(XX=XX, mod=mod, des_func=des_func,
                                                                           alpha=alpha,weight_const=weight_const)
                                            if (length(error_power) == 1 && error_power == 1) {
                                              err * weight_func_out
                                            } else if (length(error_power) == 1 && error_power == 2) {
                                              err^2 * weight_func_out
                                            } else if (length(error_power) == 2 && error_power == c(1,2)) {
                                              list(err * weight_func_out,
                                                   err^2 * weight_func_out)
                                            } else {stop("error_power not recognized in werror_func #825376")}
                                          },
                                          intwerror_func = function(..., XX=NULL, N=1e4, mod=self$mod,
                                                                    des_func=self$des_func, alpha=self$alpha_des,
                                                                    weight_const=self$weight_const,
                                                                    weight_func=self$weight_func,
                                                                    error_power=self$error_power){
                                            # use self$func instead of self$mod to get actual value
                                            if (is.null(XX)) {
                                              XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
                                            }
                                            werror_func_out <- self$werror_func(XX=XX, mod=mod, des_func=des_func,
                                                                                alpha=alpha,weight_const=weight_const,
                                                                                error_power=error_power)
                                            if (is.list(werror_func_out)) {
                                              sapply(werror_func_out, mean)
                                            } else {
                                              mean(werror_func_out)
                                            }
                                          },
                                          int_pvar_red_for_opts = function(..., Xopts, XX=NULL, N=1e4, mod=self$mod,
                                                                           des_func=self$des_func,
                                                                           alpha=self$alpha_des,
                                                                           weight_const=self$weight_const,
                                                                           weight_func=self$weight_func,
                                                                           delta_pvar_func=mean){
                                            # use self$func instead of self$mod to get actual value
                                            if (is.null(XX)) {
                                              XX <- simple.LHS(N, self$D) #matrix(runif(n*self$D), ncol=self$D)
                                            }
                                            K_X_XX <- mod$mod$corr_func(self$X, XX)
                                            to <- apply(X=Xopts, MARGIN=1, FUN=self$int_pvar_red_for_one, X_=self$X, 
                                                        XX=XX, corr_func=mod$mod$corr_func, Kinv=self$mod$mod$Kinv, 
                                                        s2=self$mod$mod$s2_hat, K_X_XX=K_X_XX,
                                                        delta_pvar_func=delta_pvar_func)
                                            to
                                          },
                                          int_pvar_red_for_one = function(v, X_, XX, corr_func, Kinv, s2, K_X_XX,
                                                                          delta_pvar_func=mean) {
                                            X <- X_ # can't pass X through apply since it matches first arg
                                            vmatrix <- matrix(v, nrow=1)
                                            Kxv <- as.numeric(corr_func(X, vmatrix))
                                            Kvv <- as.numeric(corr_func(vmatrix, vmatrix))
                                            Kxinv_Kxv <- as.numeric(Kinv %*% Kxv) # convert to vector to be faster
                                            s2_over_bottom <- as.numeric(s2/ (Kvv - t(Kxv) %*% Kxinv_Kxv))
                                            reds <- sapply(1:nrow(XX), function(i) {
                                              zmatrix <- XX[i, , drop=F]
                                              # zmatrix <- matrix(z, nrow=1)
                                              # Kxz <- corr_func(X, zmatrix)
                                              Kxz <- K_X_XX[, i]
                                              Kvz <- corr_func(vmatrix, zmatrix)
                                              t1 <- s2_over_bottom * (sum(Kxz * Kxinv_Kxv) - Kvz)^2
                                              if (is.na(t1)) {browser()}
                                              t1
                                            })
                                            # Before was just taking mean
                                            # mean(reds)
                                            # Now letting you pass in func, can weight them, or sqrt * weight
                                            delta_pvar_func(reds)
                                          },
                                          actual_intwerror_func = function(..., N=2e3, mod=self$mod, f=self$func,
                                                                           error_power=self$error_power) {
                                            if (is.null(self$actual_des_func)) {
                                              # Return NaN if user doesn't give actual_des_func
                                              return(rep(NaN, length(error_power)))
                                            }
                                            XX <- simple.LHS(n = N,d = self$D)
                                            ZZ <- mod$predict(XX)
                                            ZZ.actual <- apply(XX, 1, f)
                                            abserr <- abs(ZZ - ZZ.actual)
                                            # TODO LATER Have actual_des_func return ZZ to save time
                                            weight <- self$weight_const +
                                              self$alpha_des * self$actual_des_func(XX=XX, mod=mod)
                                            if (length(error_power) == 1 && error_power == 1) {
                                              mean(weight * abserr)
                                            } else if (length(error_power) == 1 && error_power == 2) {
                                              mean(weight * abserr^2)
                                            } else if (length(error_power) == 2 && error_power == c(1,2)) {
                                              list(mean(weight * abserr),
                                                   mean(weight * abserr^2))
                                            } else {stop("error_power not recognized in actual_intwerror_func #20497")}
                                          },
                                          print_results = function() {
                                            best_index <- which.max(self$Z)
                                            bestZ <- self$Z[best_index]
                                            bestX <- self$X[best_index, ]
                                            cat("Best design point is ", signif(bestX, 3),
                                                " with objective value ", bestZ, '\n')
                                            
                                          },
                                          delete = function() {
                                            self$mod$delete()
                                            if (self$parallel) {
                                              print("Deleting cluster")
                                              parallel::stopCluster(cl = self$parallel_cluster)
                                            }
                                          }
                                        )
)



# Plots should be in color but also work if printed in black in white
# Set color.palette = heat.colors is okay, but a lot of bright yellow/red
# Can remove some of dark red with color.palette=function(x) {heat.colors((x+6))[-(1:6)]},
#  vary 6 to change range.
# Use rev() to reverse high/low
# Could use color.palette=function(x) {gray((1:x+10)/(x+10))} to do all in gray

color.palette=function(x) {(heat.colors((x+6))[-(1:6)])}

# Save images with width 620 and height 600.


# Using our adaptive method

# limnonpoly, grad_norm2_mean, laGP_GauPro_kernel
set.seed(1); csa(); a_IMVSE <- adapt.concept2.sFFLHD.R6$new(
  D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = T,
  obj="desirability", des_func=des_func_grad_norm2_mean,
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=2, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='sFFLHD',
  error_power=2,
  selection_method="max_des_red_all_best"
  # selection_method="ALC_all_best"
); a_IMVSE$run(15, noplot = T)
# Show predicted mean and batch when pts selected
# cf(a$mod$predict, batchmax=Inf,
#    afterplotfunc=function() {
#      text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
#    xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=T)
# Show true function and batch when pts selected
cf(limnonpoly, batchmax=Inf,
   color.palette=color.palette,
   afterplotfunc=function() {
     text(a_IMVSE$X[,1], a_IMVSE$X[,2], a_IMVSE$X_tracker$iteration_added)},
   xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)


# Using the standard IMSE method

# limnonpoly, IMSE
set.seed(10); csa(); a_IMSE <- adapt.concept2.sFFLHD.R6$new(
  D=2,L=3,func=limnonpoly, nugget = 1e-7,estimate.nugget = T,
  obj="desirability", des_func=des_func_grad_norm2_mean,
  actual_des_func=NULL,#get_num_actual_des_func_grad_norm2_mean(),
  stage1batches=2, alpha_des=1, weight_const=0,
  package="laGP_GauPro_kernel", design='sFFLHD',
  error_power=2,
  # selection_method="max_des_red_all_best"
  selection_method="ALC_all_best"
); a_IMSE$run(15, noplot = T)
# Show predicted mean and batch when pts selected
# cf(a$mod$predict, batchmax=Inf,
#    afterplotfunc=function() {
#      text(a$X[,1], a$X[,2], a$X_tracker$iteration_added)},
#    xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=T)
# Show true function and batch when pts selected
cf(limnonpoly, batchmax=Inf,
   color.palette=color.palette,
   afterplotfunc=function() {
     text(a_IMSE$X[,1], a_IMSE$X[,2], a_IMSE$X_tracker$iteration_added)},
   xlim=c(-.02,1.02), ylim=c(-.02,1.02), bar=F, mainminmax=F)
