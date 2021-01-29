#' @title AdaTiSS
#' @description To obtain tissue specificity scores.
#' @author Meng Wang
#' \email{mengw1@stanford.edu}
#' @param X The expression matrix (prefered in log scale).
#' @param tiss.abd The expression summarized in tissue level (default: NULL)
#' @return The score matrix in sample level (ada.s) and in tissue level (ada.z) and population fitting info (pop.fit.mx). Note: to take another care on the genes with 'pi0.hat' <= 0.5 in the pop.fit.mx.
#' @export

AdaTiSS = function(X, tiss.abd=NULL) {

	pop.fit.mx = matrix(NA, nrow(X), 6)
	rownames(pop.fit.mx) = rownames(X)
	colnames(pop.fit.mx) = c("n.observed", "gam.sel", "mu0.hat", "sd0.hat" , "pi0.hat" ,"crt")
	id.ls.1 = rownames(X)[rowSums(!is.na(X)) >= 20]
	length(id.ls.1)
	for (i in 1:length(id.ls.1)) { 
	  if(i %% 500 == 0) print(i)
	  id = id.ls.1[i]
	  x.0 = X[id, ] 
	  gam.limit = ifelse(sum(!is.na(x.0)) <= 100, 1, 3)
	  result.x = adapt.gam.rob.fit.fn(x.0, gam.seq=seq(0,gam.limit,by=0.1), bin.num=round(length(x.0)/10))
	  pop.fit.mx[id,] = result.x[["est.hat"]]
	}
	
	id.ls.2 = setdiff(rownames(X), id.ls.1)
	length(id.ls.2)
	pop.info.2 = apply(X[id.ls.2, ], 1, function(x) c(sum(!is.na(x)), median(x, na.rm=TRUE), mad(x, na.rm=TRUE), sum(abs(x-median(x, na.rm=TRUE)) <= 2*mad(x, na.rm=TRUE), na.rm=TRUE )/sum(!is.na(x)) ))
    pop.info.2 = t(pop.info.2)
    pop.fit.mx[id.ls.2, c("n.observed", "mu0.hat", "sd0.hat" , "pi0.hat")] = pop.info.2
    
    pop.fit.mx[, 'sd0.hat'] = pmax(pop.fit.mx[, 'sd0.hat'], 0.01)
    ada.s = (X - outer(pop.fit.mx[rownames(X), "mu0.hat"], rep(1, ncol(X))))/outer(pop.fit.mx[rownames(X), "sd0.hat"], rep(1, ncol(X)))
    if (!is.null(tiss.abd)) {
    	ada.z = (tiss.abd - outer(pop.fit.mx[rownames(tiss.abd), "mu0.hat"], rep(1, ncol(tiss.abd))))/outer(pop.fit.mx[rownames(tiss.abd), "sd0.hat"], rep(1, ncol(tiss.abd)))
    } else {
    	ada.z = NULL
    }


    return(list(ada.s = ada.s, ada.z=ada.z, pop.fit.mx=pop.fit.mx))
}


# ------------------------------------------------------------------------------
# data-adaptive selection procedure
adapt.gam.rob.fit.fn = function (x.00, gam.seq, step=50, mu.fix=NULL, var.fix=NULL, bin.num=NULL) {

	x.0 = x.00[!is.na(x.00)]
	nm = c("mu0.hat", "sd0.hat", "pi0.hat", 'crt')
	par.hat = matrix(NA, length(gam.seq), length(nm))
	rownames(par.hat) = gam.seq
	colnames(par.hat) = nm
	
	for (i in 1:length(gam.seq)) {
		  	gam = gam.seq[i]
			x.mu = ifelse(is.null(mu.fix), mean(x.0), mu.fix)
		  	x.var = ifelse(is.null(var.fix), var(x.0), var.fix)
			result = est.fn(x.0, x.mu, x.var, gam, fix.mu=!is.null(mu.fix), fix.var=!is.null(var.fix), step=step)		
			mu.hat = result$mu.est
			var.hat = result$var.est
			if (!is.na(var.hat)) {
				est.result = efdr.0.fn(x.0, mu.hat, var.hat, gam, bin.num)
				par.hat[i, ] = est.result[colnames(par.hat)]
			}

	}
	crt.hat.0 = abs(pmin(par.hat[,'crt'], 10) - 1)
	ind.comp = !is.na(crt.hat.0)
	gam.comp = gam.seq[ind.comp]
	if (length(gam.comp)  == 0 ) {
		 est.hat=NA
		 x.n=NA
		 w=NA
		 gam.comp=NA
		 crt.hat.0=NA
		 para.hat.mx=NA
	} else {
		crt.hat.0 = crt.hat.0[ind.comp]
		par.hat = matrix(par.hat[ind.comp, ], ncol=length(nm))
		colnames(par.hat) = nm
		rownames(par.hat) = gam.comp
		
		crt.est.0 = pmin(par.hat[,'crt'], 10)
		gam.sel = gam.comp[which.min(crt.hat.0)]
		
		gam.sel.char = as.character(gam.sel)
		est.hat = c(length(x.0), gam.sel, par.hat[gam.sel.char, c("mu0.hat", "sd0.hat", "pi0.hat", "crt")])
		names(est.hat)[1:2] = c("n.observed", "gam.sel")
		
		w.nu = dnorm(x.00, est.hat["mu0.hat"], est.hat["sd0.hat"])^est.hat["gam.sel"]
		w.nu[is.na(x.00)] = NA
		w = w.nu/sum(w.nu, na.rm=TRUE)
	}
                     
	return(list(est.hat=est.hat, x.w=w, gam.comp=gam.comp, crt.hat.0=crt.hat.0,  para.hat.mx=par.hat )) 
	
}

# ------------------------------------------------------------------------------
# expected of fdr criterion
efdr.0.fn = function (x, mu.hat, var.hat, gam, bin.num=NULL) {
	    x = x[!is.na(x)]
	    den.fit = dnorm(x, mu.hat, sqrt(var.hat))
		frac.hat= mean(den.fit^gam)*sqrt(2*pi*var.hat)^gam * sqrt(1 + gam)
			
		my.hist = bk.cnt.fn(x, bin.num)
		bin.bk = my.hist[[1]]
		cnt.bk = my.hist[[2]]

        p0.hat.bin  = numeric(length(bin.bk)-1)
        p0.hat.bin[1] = pnorm(bin.bk[2], mu.hat, sqrt(var.hat))
        p0.hat.bin[length(bin.bk)-1] = 1 - pnorm(bin.bk[length(bin.bk)-1], mu.hat, sqrt(var.hat))
	    for ( j in 3:(length(bin.bk)-1)) {
	      	   p0.hat.bin[j-1] = pnorm(bin.bk[j], mu.hat, sqrt(var.hat)) - pnorm(bin.bk[j-1], mu.hat, sqrt(var.hat))
	     }
	    p.hat.bin =  cnt.bk/sum(cnt.bk)
	    null.efdr.hat = min(1,frac.hat) * sum(p0.hat.bin^2/p.hat.bin)	
	    est.sum = c(gam, mu.hat, sqrt(var.hat),  frac.hat, null.efdr.hat)
	    names(est.sum) = c("gamma", "mu0.hat", "sd0.hat", "pi0.hat",  'crt')
        return(est.sum)
}



# ------------------------------------------------------------------------------
# estimation under a fixed gamma
est.fn = function(x, mu.0, var.0, gam,  tol=10^(-4), step=step, fix.mu=FALSE, fix.var=FALSE) {
	          x = x[!is.na(x)]
	          n = length(x)
	   	      dum = dnorm(x, mu.0, sqrt(var.0))^gam
	          w.0 = dum/sum(dum)
	         
	       	  int = 1
	          flag = FALSE
	    	  diff.par.int = c()
	          mu.int = mu.0
	          var.int = var.0
	          while ( flag == FALSE) {
	          		     mu.1 = ifelse(fix.mu, mu.0, sum(w.0 * x))
	          	         var.1 =  ifelse(fix.var, var.0, (1+gam) * sum(w.0 * (x-mu.1)^2)) 
	          	         if (var.1 < 10^(-4)) {
							  mu.0 = NA
							  var.0 = NA
							  diff.par = NA
							  diff.par.int = c( diff.par.int, diff.par)
   							  break;
						 }
     	          	     	 
	          	         diff.par = abs(mu.1 - mu.0) + abs(sqrt(var.1) - sqrt(var.0))
	          	         diff.par.int = c( diff.par.int, diff.par)
	          	         if (diff.par < tol | int > step ){  
	          	         	flag = TRUE
	          	         	break;
	          	         } else {
	          	       	 	mu.0 = mu.1
	          	         	var.0 = var.1
	          	         	dum = dnorm(x, mu.0, sqrt(var.0))^gam
	                        w.0 = dum/sum(dum)
	                        mu.int = c(mu.int, mu.0)
	                        var.int = c(var.int, var.0)
	                        int = int + 1
	          	         }
	          }
	          return(list(mu.est=mu.0, var.est=var.0, w=w.0, diff.par.est=diff.par.int))
}


# ------------------------------------------------------------------------------
# to merge intervals s.t. each bin has positive number of data points
bk.cnt.fn = function (x, bin.num=NULL) {
	if (is.null(bin.num)) {
		if (length(x) > 1000) {
			bk.num = 20
		} 
		if (length(x) <= 1000 & length(x) > 500) {
			bk.num = 10
		} 
		if (length(x) <= 500){
			bk.num= 5
		}
	} else {
		bk.num = bin.num
	}
	h = hist(x, breaks=bk.num, plot=FALSE)
	# h$counts
	# h$breaks
	ind.zero = (1:length(h$counts))[h$counts == 0]
	if (length(ind.zero) != 0) {
	   bk.start = h$breaks[ind.zero]
		bk.end = h$breaks[ind.zero + 2]
		bk.update = c(bk.start[1])
		if (length(bk.start) > 1) {
				for (i in 1: (length(bk.start)-1) ) {
					  if (bk.end[i] < bk.start[i+1]) bk.update = c(bk.update, bk.end[i])					      
				}
	    	}
	   bk.update = c(bk.update,  bk.end[length(bk.end)] )
		bk.pts = sort(unique(c(h$breaks[-c(ind.zero, ind.zero+1)], bk.update)))
		cnt = unname(table(cut(x, breaks=bk.pts)))
	} else {
		bk.pts = h$breaks
		cnt = h$counts
	}
    return( list( bk.pts, cnt))
}

