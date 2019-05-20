
###################################################################################################
#' Surface plot of a function
#'
#' A (filled) contour plot or perspective / surface plot of a function.
#'
#' @param f function to be plotted. The function should either be able to take two vectors or one matrix specifying sample locations. i.e. \code{z=f(X)} or \code{z=f(x2,x1)} where Z is a two column matrix containing the sample locations \code{x1} and \code{x2}.
#' @param lower boundary for x1 and x2 (defaults to \code{c(0,0)}).
#' @param upper boundary (defaults to \code{c(1,1)}).
#' @param type string describing the type of the plot:  \code{"filled.contour"} (default), \code{"contour"}, 
#' \code{"persp"} (perspective), or \code{"persp3d"} plot.
#' Note that "persp3d" is based on the plotly package and will work in RStudio, but not in the standard RGui.
#' @param s number of samples along each dimension. e.g. \code{f} will be evaluated \code{s^2} times.
#' @param xlab lable of first axis
#' @param ylab lable of second axis
#' @param zlab lable of third axis
#' @param color.palette colors used, default is \code{terrain.color}
#' @param title of the plot 
#' @param levels number of levels for the plotted function value. Will be set automatically with default NULL.. (contour plots  only)
#' @param points1 can be omitted, but if given the points in this matrix are added to the plot in form of dots. Contour plots and persp3d only. Contour plots expect matrix with two columns for coordinates. 3Dperspective expects matrix with three columns, third column giving the corresponding observed value of the plotted function.
#' @param points2 can be omitted, but if given the points in this matrix are added to the plot in form of crosses. Contour plots and persp3d only.  Contour plots expect matrix with two columns for coordinates. 3Dperspective expects matrix with three columns, third column giving the corresponding observed value of the plotted function.
#' @param pch1 pch (symbol) setting for points1 (default: 20). (contour plots only)
#' @param pch2 pch (symbol) setting for points2 (default: 8). (contour plots only)
#' @param lwd1 line width for points1 (default: 1). (contour plots only)
#' @param lwd2 line width for points2 (default: 1). (contour plots only)
#' @param cex1 cex for points1 (default: 1). (contour plots only)
#' @param cex2 cex for points2 (default: 1). (contour plots only)
#' @param col1 color for points1 (default: "black"). (contour plots only)
#' @param col2 color for points2 (default: "black"). (contour plots only)
#' @param theta angle defining the viewing direction. theta gives the azimuthal direction and phi the colatitude. (persp plot only) 
#' @param phi angle defining the viewing direction. theta gives the colatitude. (persp plot only) 
#' @param ... additional parameters passed to \code{contour} or \code{filled.contour}
#'
#' @examples
#' plotFunction(function(x){rowSums(x^2)},c(-5,0),c(10,15))
#' plotFunction(function(x){rowSums(x^2)},c(-5,0),c(10,15),type="contour")
#' plotFunction(function(x){rowSums(x^2)},c(-5,0),c(10,15),type="persp")
#'
#' @seealso \code{\link{plotData}}, \code{\link{plotModel}}
#'
#' @export
#' @importFrom plotly plot_ly
#' @importFrom plotly add_trace
#' @importFrom plotly %>%
###################################################################################################
plotFunction <- function(f=function(x){rowSums(x^2)}, 
                                lower=c(0,0) , upper=c(1,1) , 
																type="filled.contour",
																s=100, 
                                xlab="x1",ylab="x2", zlab="y",
                                color.palette = terrain.colors, 
                                title=" ",  levels=NULL, 
                                points1, points2, pch1=20, pch2=8, lwd1=1, lwd2=1, cex1=1, cex2=1, col1="red", col2="black",
																theta=-40,phi=40,
																...){
  x <- seq(lower[1], upper[1], length = s)  
  y <- seq(lower[2], upper[2], length = s) 
  if(length(formals(f))==1){
    fn <- function(a,b){f(cbind(a,b))}	
    z <- outer(x, y, fn)
  }else if(length(formals(f))==2){
    z <- outer(x, y, f)	
  }		
  
  if(is.null(levels))
    levels=pretty(range(z[!is.na(z)]),20)
  
  if(type=="filled.contour"){
    if(missing(points1)&missing(points2)){
      filled.contour(x, y, z, color.palette=color.palette, levels=levels,
                     key.title = title(main = zlab),
                     plot.title=title(title,
                                      xlab=xlab,
                                      ylab=ylab),
										 ...
										 )
    }else if(missing(points1)&!missing(points2)){
      filled.contour(x, y, z, color.palette=color.palette, levels=levels,
                     key.title = title(main = zlab),
                     plot.title=title(title,
                                      xlab=xlab,
                                      ylab=ylab),
                     plot.axes = { points(points2,pch=pch2,col=col2,cex=cex2,lwd=lwd2); axis(1); axis(2);	},
										 ...
										 )
    }else if(!missing(points1)&missing(points2)){
      filled.contour(x, y, z, color.palette=color.palette, levels=levels,
                     key.title = title(main = zlab),
                     plot.title=title(title,
                                      xlab=xlab,
                                      ylab=ylab),
                     plot.axes = { points(points1,pch=pch1,cex=cex1,lwd=lwd1,col=col1); axis(1); axis(2);	 },
										 ...
										 )
    }else{
      filled.contour(x, y, z, color.palette=color.palette, levels=levels,
                     key.title = title(main = zlab),
                     plot.title=title(title,
                                      xlab=xlab,
                                      ylab=ylab),
                     plot.axes = { points(points1,pch=pch1,cex=cex1,lwd=lwd1,col=col1); points(points2,pch=pch2,col=col2,cex=cex2,lwd=lwd2);axis(1); axis(2); },
										 ...
										 )
    }
  }else if(type=="contour"){ #not filled
    contour(x,y,z,
            xlab = xlab,
            ylab = ylab,
            main = "",
            key.title = title(main = zlab),...)
    if(!missing(points1))
       points(points1,pch=pch1,cex=cex1,lwd=lwd1,col=col1)
    if(!missing(points2))
      points(points2,pch=pch2,col=col2,cex=cex2,lwd=lwd2)
  }else if(type=="persp"){ #perspective
		# Color palette terrain, 100 different colors.
		colors<-terrain.colors(100)
		# height of facets (for coloring each facet)
		z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
		# Range of the facet center on a 100-scale (number of colors)
		z.facet.range<-cut(z.facet.center, 100)
		# Plot
		persp(x=x, y=y, z=z, 
					xlab=xlab,
					ylab=ylab,
					zlab=zlab,
					main=title,
					col=colors[z.facet.range],
					theta=theta,phi=phi,...)	
  }else if(type=="persp3d"){ #perspective plot with plotly
		p <- plot_ly(z = ~t(z), x = x, y = y,type = "surface")# %>% add_surface()
    if(!missing(points1))
      p <- p %>% add_trace(data=points1,x=points1[,1],z=points1[,3],y=points1[,2], mode = "markers", type = "scatter3d", 
            marker = list(size = 5, color = col1, symbol = 200))
    if(!missing(points2))
      p <- p %>% add_trace(data=points2,x=points2[,1],z=points2[,3],y=points2[,2], mode = "markers", type = "scatter3d", 
            marker = list(size = 5, color = col2, symbol = 102))
		p
  }
}
