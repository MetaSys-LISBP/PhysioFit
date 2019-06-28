# Copyright 2008-2013, INRA, France
# Author: Serguei SOKOL, sokol@insa-toulouse.fr
# License: GPL v2 http://www.gnu.org/licenses/gpl-2.0.html

# This an eclectic collection of R functions useful in variuous
# contexts.

# map values on color table for plot min->col[1], max->col[n]
val2col = function(val, col=rainbow(16L,start=0,end=4/6)) {
# creates vectors of length=NROW(val) translating val to colors
# min_val is mapped to col[1] and
# max_val is mapped to col[maxcol]
# default 16 rainbow colors min in red, max in blue
   val.min = min(val)
   val.max = max(val)
   nval = NROW(val)
   ncol = NROW(col)
   res = rep("", nval)
   for (i in 1L:nval) {
      res[i] = col[1L+floor(ncol*(val[i]-val.min)/(1.001*(val.max-val.min)))]
   }
   res
}
norm2<-function(v) crossprod(as.numeric(v))[1L]
norm<-function(v) {
   sqrt(norm2(v))
}
cnorm2=function(v) {
   Re(sum(Adj(v)%*%(v)))
}
house<-function(v,x,st=1L) {
   # apply housholder transform defined by v to x or
   # to a part of x starting from st
   n=NROW(x)
   x=as.matrix(x,n)
   m=NCOL(x)
   if (st>n) {
      # hmm, do nothing
      return(x)
   }
   v=c(v,rep(0,n-st+1L-NROW(v)))
   for (j in 1L:m) {
      x[st:n,j]=x[st:n,j]-2*v*(v%*%x[st:n,j])
   }
   x
}
#pos=function(x) 0.5*(x+abs(x))
pos=function(x) {x[x<0]=0.; return(x)}
neg=function(x) {x[x>0]=0.; return(-x)}
spick_det=function(sp, npick, lag=5) {
   # singulet pick detection in a nmr spectrum
   # returns a param array suitable for deconvolution
   n=length(sp)
   cs=exp_smooth(cumsum(sp))
   di=diff(cs, 3, lag=lag)
   o=order(di)+30
   param=c(0., 0., 0.0001); # sigma and lambda detected later, gauss fraction set very small
   for (i in 1L:npick) {
      # next greatest pick
      param=c(param, o[1L]); # add position of the pick
      # detect half height pick support
      maxp=sp[o[1L]]
      il=o[1L]:max(o[1L]-10*lag,1)
      lpos=which(sp[il]<=0.5*maxp)[1L]
      ir=o[1L]:min(o[1L]+10*lag,n)
      rpos=which(sp[ir]<=0.5*maxp)[1L]
      hw=0.5*(rpos+lpos)
      amp=pi*hw*maxp
      param=c(param, amp); # add amplitude of the pick
      # cancel 4*hw from o
      il=max(o[1L]-2*hw,1)
      ir=min(o[1L]+2*hw,n)
      o=o[!(o %in% il:ir)]
   }
   param[1L] = param[2L] = hw
}
exp_smooth=function(y, lambda=1./3.) {
   # smooth by convolution with normfactor*exp(-lambda*|x|)
   # lambda is in [0;1]
   # 0 => very smoothed (constant=leftmost value, not normalized)
   # 1 => no smooth at all
   n=NROW(y)
   ys=y
   l1=(1.-lambda)
   # forward march
   for (i in 2:n) {
      ys[i]=lambda*ys[i]+l1*ys[i-1L]
   }
   # backward march
   for (i in (n-1L):1L) {
      ys[i]=lambda*ys[i]+l1*ys[i+1L]
   }
   return(ys)
}
proc2par=function(dirn) {
   # in dirname dirn get procs and proc2s files, extract
   # "xdim", "si"  and some other params, then transform them in
   # a list with members: dims, ppm, axename
   # dims[1] is moving fastest
   fn=c("procs", "proc2s")
   
   # place holders
   xdim=1L:2L
   si=1:2
   ppm=cbind(start=c(0, 0), end=c(0, 0))
   sf=c(0, 0)
   axename=c("", "")
   
   for (i in 1L:2L) {
      f=file.path(dirn, fn[i])
      data=as.matrix(read.table(f, comment.char="", fill=TRUE))
      
      ir=which(data[,1L]=="##$XDIM=")
      xdim[i]=as.integer(data[ir,2L])
      
      ir=which(data[,1L]=="##$SI=")
      si[i]=as.integer(data[ir,2L])
      
      ir=which(data[,1L]=="##$AXNAME=")
      axename[i]=data[ir,2L]
      
      ir=which(data[,1L]=="##$F1P=")
      ppm[i,]=as.double(c(data[ir+1L,2L], data[ir,2L]))
      
      ir=which(data[,1L]=="##$SF=")
      sf[i]=as.double(c(data[ir, 2]))
   }
   rownames(ppm)=axename
   return(list(dims=c(xdim, floor(si/xdim)), ppm=ppm, sf=sf, axename=axename))
}
bf2mtx=function(file_name, dims=proc2par(dirname(file_name))$dims, ty="int") {
   # read a block structured file (like 2rr of topspin)
   # cnct is a connection
   # dims is a vector of dimensions of size 4
   # dims[1] is for index moving fastest
   d=array(readBin(file_name, ty, prod(dims)), dims)
   dc=matrix(d, dims[1L]*dims[3L], dims[2L]*dims[4L])
   # repack in contigous way
   if (dims[3L]==1 && dims[4L]==1) {
      return(dc)
   }
   for (i4 in 1:dims[4L]) {
      for (i2 in 1:dims[2L]) {
         dc[,i2+(i4-1)*dims[2L]]=d[,i2,,i4]
      }
   }
   return(dc)
}
read.2Dbruk=function(file_name, dims=proc2par(dirname(file_name))$dims, ty="int") {
   # read bruker formated file of 2D spectrum (like 2rr, 2ri and so on)
   # In the same directory it must be files procs and proc2s or
   # parameter dims provided explicitly.
   # Return a 4-dimenesionnal array, first two dimensions for
   # a block size, and the last two for block numbers in F1 and F2
   # directions.
   # use matrix(aperm(result, c(1,3,2,4)), nrow=dim1*dim3) to transofrm
   # a resulting array into contigous matrix
   array(readBin(file_name, ty, prod(dims)), dims)
}
read.1Dbruk=function(dir_name, ty="int") {
   # read bruker formated file of 1D spectrum real=1r, imaginary=1i
   # In the same directory it must be files 1r and 1i
   # Return a complex vector.
   fr=file.path(dir_name, "1r")
   fi=file.path(dir_name, "1i")
   nb_pnt=file.info(fr)$size/4
   complex(
      re=readBin(fr, ty, nb_pnt),
      im=readBin(fi, ty, nb_pnt)
   )
}
terre3d=function(z3, ...) {
   library(rgl)
   x=1:nrow(z3)
   y=1:ncol(z3)
   zlim <- range(z3)
   fact=max(x,y)/(zlim[2L]-zlim[1L])/2
   zlen=16
   colorlut=terrain.colors(zlen)
   col=colorlut[round((z3-zlim[1L])/(zlim[2L]-zlim[1L])*(zlen-1)) + 1 ]
   surface3d(x, y, fact*z3, color=col, back="lines", ...)
}
qcut=function(x, qntl=0.5, val=0) {
   # quantile cut to val of x columns
   # v[v<=quantile(v,qntl)]=val
   x=as.matrix(x, nrow=NROW(x))
   return(apply(x,2,function(v){v[v<=quantile(v,qntl)]=val; return(v)}))
}
cmpcut=function(x, lag=round(0.5*NROW(x)), val=0, app=2) {
   # compact cut to val of x columns (app=2) or rows (app=1)
   # start=which.max(diff(cumsum(v),lag=lag))
   # v[1:start,(start+lag):end]=val
   # NB: the support of non cut values is equal lag+1
   x=as.matrix(x, nrow=NROW(x))
   x=apply(x,app,function(v){n=NROW(v); if (lag>=n) { s=1; } else { s=which.max(diff(cumsum(v),lag=lag)); }; if (s > 1) { v[1:(s-1)]=val; }; if ((s+lag) < n) { v[(s+lag+1):n]=val; }; return(v)})
   if (app==1) {
      x=t(x)
   }
   x[x<val]=val
   return(x)
}
smooth025=function(v, repet=1) {
   # smooth by convolution with 0.25;0.5;0.25
   # NB: the returned vector is two points shorter than initial one
   # (one point out from each end)
   v=as.vector(v)
   n=NROW(v)
#cat("n=", n, "\n", sep="")
   if (n<3 || repet < 1) {
      # the vector is too short or no smoothing is asked
      # send it back as is
      return(v)
   }
   v=v[1:(n-1)]+v[2:n]
   v=0.25*(v[1:(n-2)]+v[2:(n-1)])
   if (repet<=1) {
      return(v)
   } else {
      return(smooth025(v, repet-1))
   }
}
smooth725=function(v, repet=1) {
   # smooth by convolution with 0.25;0.5;0.25
   # NB: the returned vector have the same length as original one.
   # Extream points have weights 0.75, 0.25 on left and 0.25;0.75 on right
   v=as.vector(v)
   n=NROW(v)
#cat("n=", n, "\n", sep="")
   if (n<2 || repet < 1) {
      # the vector is too short or no smoothing is asked
      # send it back as is
      return(v)
   }
   # extend v beyond interval and NA gaps
   # indexes of valid value having NA on its left
   nal=c(0, which(is.na(head(v,-1)) & !is.na(v[-1L])))+1
   # indexes of valid value having NA on its right
   nar=c(which(is.na(v[-1L]) & !is.na(head(v,-1))), n)
   ve=c(NA, v, NA)
   ve[nal]=v[nal]
   ve[nar+2L]=v[nar]
   ve=head(ve,-1)+ve[-1L]
   ve=0.25*(head(ve,-1)+ve[-1L])
   if (repet<=1) {
      return(ve)
   } else {
      return(smooth725(ve, repet-1))
   }
}
smdec=function(v) {
   # smooth by convolution with 0.25;0.5;0.25 and decimate
   # by taking every second point
   n=NROW(v)
#cat("n=", n, "\n", sep="")
   if (n<3) {
      # the vector is too short
      # send it back as is
      return(v)
   }
#cat("n2=", n, "\n", sep="")
   return(smooth025(v)[seq(1,n-2,2)])
}

smdec2d=function(z) {
   z1=apply(z,2,smdec)
   return(t(apply(z1,1,smdec)))
}

coarse_n=function(v, nsc=ceiling(log2(length(v))), f=smdec) {
   # return a list of scaled down vectors
   # first list element is original vector
   # second one is scaled down once and so on.
   # If all scales down are possible
   # then the list will have nsc+1 items
   # otherwise the last item is the vector
   # which could not be scaled down by lack of elements
   # f is a function used to passe from fine
   # to coarse level (smdec by default)
   res=list(v)
   if (nsc <= 0) {
      return(res)
   }
   for (i in 1:nsc) {
      sc=f(res[[i]])
      if (NROW(sc)==NROW(res[[i]]) || (NCOL(res[[i]]) > 1 && (NCOL(sc)==NCOL(res[[i]])))) {
         # no more possible scale down
         return(res)
      }
      res=append(res, list(sc))
   }
}

refine=function(vc,nf) {
   # refine the coarse vector vc to the size nf.
   # This is inverse of the coarsening (smooth+decimation) operation
   # refinement is done by linear interpolation.
   nc=NROW(vc)
   if (nf < 2*nc+1) {
      # bad coherence between vector lengths
      # do nothing
      return(vc)
   }
   vf=rep(0., nf)
   dvc=0.5*diff(c(0., vc, 0.))
   i2=seq(2,nf-1,by=2)
   # inject coarse values on even nodes
   vf[i2]=vc
   # interpolate on odd nodes
   vf[i2-1L]=vc-dvc[1:nc]
   # interpolate at the right end
   vf[i2[nc]+1L]=vc[nc]+dvc[nc+1L]
   if (i2[nc]!=nf-1) {
      # very far right end is just put to zero
      vf[nf]=0.
   }
   return(vf)
}

refine2x1d=function(zc,nf1,nf2) {
   # return refined zc to n1 x n2 field
   zf=matrix(0, nf1, nf2)
   nc1=nrow(zc)
   nc2=ncol(zc)
   for (i2 in 1:nc2) {
      # refine along the columns
      zf[,i2]=refine(zc[,i2], nf1)
   }
   for (i1 in 1:nf1) {
      # refine along the rows
      zf[i1,]=refine(zf[i1,1:nc2], nf2)
   }
   return(zf)
}

ipeaks=function(y, decreasing=TRUE) {
   # return a vector of indexes correponding to peaks in a vector y
   # By default, the peaks are sorted in decreasing order.
   # The peaks are searched only in the form /\ (mountain) and not in /-\ (plato)
   d=sign(diff(y))
   # make the same length as y (linear extrapolation)
   d=c(1,d,-1); # doute benefice on the ends
   d=diff(d)
   ip=which(d==-2)
   o=order(y[ip], decreasing=decreasing)
   return(ip[o])
}
peak_tv=function(y, ip=ipeaks(y)) {
   # y is data vector, ip are indexes of peaks
   # total variation of a peak=2*y-yleft-yright where yleft
   # and yright are closest valley's bottoms on left and right from peak
   # return tv vector of the same length as ip and in the same order
   # than peaks
   # cf. ivalleys() and ihhw_peak()
   o=order(ip)
   n=length(y)
   ipb=c(1, ip[o], n)
   np=length(ip)
   yl=apply(t(1:np), 2, function(i) {max(y[ipb[i+1L]]-y[ipb[i]:ipb[i+1L]], na.rm=T)})
   yr=apply(t(1+(1:np)), 2, function(i) {max(y[ipb[i]]-y[ipb[i]:ipb[i+1L]], na.rm=T)})
   tv=1:np
   tv[o]=yr+yl
   return(tv)
}
ivalleys=function(y, ip=ipeaks(y)) {
   # y is data vector, ip are indexes of peaks
   # return indexes of deepest points between peaks and index limits
   o=order(ip)
   n=length(y)
   ipb=c(1, ip[o], n)
   np=length(ip)
   il=ir=numeric(np)
   il[o]=ipb[1:np]-1+apply(t(1:np), 2, function(i) {which.max(y[ipb[i+1L]]-y[ipb[i]:ipb[i+1L]])})
   ir[o]=ipb[1+(1:np)]-1+apply(t(1+(1:np)), 2, function(i) {which.max(y[ipb[i]]-y[ipb[i]:ipb[i+1L]])})
   return(cbind(left=il, center=ip, right=ir))
}
ihhw_peak=function(y, ip=ipeaks(y)) {
   # y is data vector, ip are indexes of peaks
   # return indexes of half-height width or deepest points between
   # peaks and index limits whichever comes first.
   o=order(ip)
   n=length(y)
   ipb=c(1, ip[o], n)
   np=length(ip)
   il=ir=numeric(np)
   il[o]=ipb[1:np]-1+apply(t(1:np), 2, function(i) {v=y[ipb[i]:ipb[i+1L]]; which.last(y[ipb[i+1L]]/2.>=v | v==min(v, na.rm=T))})
   ir[o]=ipb[1+(1:np)]-1+apply(t(1+(1:np)), 2, function(i) {v=y[ipb[i]:ipb[i+1L]]; which.first(y[ipb[i]]/2.>=v | v==min(v, na.rm=T))})
   return(cbind(left=il, center=ip, right=ir))
}
ipeaks2d=function(y, decreasing=TRUE) {
   # return a matrix n x 2 who's rows are composed of
   # indexe couples correponding to peaks in a matrix y
   # By default, the peaks are sorted in decreasing order.
   # The peaks are searched only in the form /\ (2d mountain in x,y direction)
   # and not in /-\ (plato) form.
   # The matrix y is prolongated by 0-order approximation
   n=dim(y)
   # make the same length as y (const extrapolation)
   d1=diff(sign(diff(rbind(y[1,],y,y[n[1L],]))))
   d2=t(diff(sign(diff(rbind(y[,1L],t(y),y[,n[2L]])))))
   i=d1*d2==4
   o=order(y[i], decreasing=decreasing)
   return(matrix(which(i, arr.ind=T)[o,], length(o)))
}
ipeak_close=function(d, iref, thre=0.5, repet=0, ip=ipeaks(smooth025(d-min(d),repet))) {
   # among all peaks in d chose one whose index
   # is closest to iref and whose height is above
   # thre*highest peak
   n=length(d)
   # shift d to make it positive then smooth
   d=smooth025(d-min(d),repet)
   # keep only peaks who's tv is above threshold*highest_tv
   tv=peak_tv(d, ip)
   ip=ip[tv>=thre*tv[1L]]+repet
   irt=ip[which.min(abs(ip-iref))]
   return(irt)
}
lorentz=function(w,w0=0,lambda=1) {
   l=1./lambda
   w=l*(w-w0)
   res=l/(pi*(w**2+1.))
   return(res)
   # (w-w0) at half height = lambda
}

discr_lorentz=function(w,w0=0,lambda=1) {
   # discreet approximation of lorentz function
   # w must be a vector of midpoints on which the
   # mean lorentz values are calculated
   l=1./lambda
   w=l*(w-w0)
   n=NROW(w)
   if (n<2) {
      stop("discr_lorentz: vector dimension is too low")
   }
   w=0.5*c(3*w[1L]-w[2L],w[1:(n-1)]+w[2:n],3*w[n]-w[n-1L])
   return(diff(atan(w))/diff(w)/pi)
}

solve_toeplitz=function(Toe,y) {
   # solves toeplitz system Toe%*%x=y
   # Toe is a list with fileds d, col, row
   # defining a toeplitz matrix (n+1)x(n+1) where n=length(col)
   # The Levison-Durbin algorithm is used.
   # x is returned
   col=Toe$col
   ns=length(col)
   d=Toe$d
   row=Toe$row
   x=f=b=rep(0., ns+1)
   f[1L]=b[1L]=1./d
   x[1L]=y[1L]*b[1L]
   for (i in 1:ns) {
      ii=1:i
      ii1=c(ii,i+1)
      eb=row[ii]%*%b[ii]
      ef=rev(col[ii])%*%f[ii]
      ex=rev(col[ii])%*%x[ii]
      deti=1./(1.-ef*eb)
      ft=deti*f[ii1]
      bt=deti*c(0.,b[ii])
      f[ii1]=ft-ef*bt
      b[ii1]=bt-eb*ft
      x[ii1]=x[ii1]+(y[i+1L]-ex)*b[ii1]
   }
   return(x)
}

cqr_gs=function(a, tol=1.e-10) {
   # QR decomposition of a (complex) matrix a by Gram-Schmidt algorithm
   # without pivoting
   # NB. For some matrices it can be numerically unstable.
   # Q is normal and orthogonal Q* %*% Q=diag(1)
   # dim(a) = c(n,p), n>=p
   # dim(Q) = c(n,p1) with p1 < p if a is rank deficient
   # dim(R) = c(p1,p)
   # rank(a) = p1
   a=as.matrix(a)
   n=dim(a)[1L]
   p=dim(a)[2L]
   p1=p
   Q=matrix(0, n, p)
   R=matrix(0, p, p)
   aux=rep(0, p)
   for (ic in 1:p) {
      u=Adj(a[,ic])
      if (ic > 1) {
         icm1=min(ic-1,p1)
         #r=Adj(u%*%Q[,1:icm1])
         # reorthogonalisation
         r=rep(0,icm1)
         for (icc in 1:icm1) {
            r[icc]=Adj(u%*%Q[,icc])
            u=u-Adj(r[icc]*Q[,icc])
         }
         #prj=Q[,1:icm1]%*%r
         #u=u-Adj(prj)
         R[1:icm1,ic]=r
      }
      nrm=as.vector(sqrt(Mod(u%*%Adj(u))))
      if (p1==p && Mod(nrm)<=tol) {
         # rank deficient
         p1=ic-1
      }
      if (Mod(nrm)>tol) {
         u=u/nrm
         R[ic,ic]=u%*%a[,ic]
         Q[,ic]=Adj(u)
      }
      aux[ic]=nrm
   }
   return(list(Q=Q[,1:p1], R=R[1:p1,], rank=p1, nrm=aux))
}
cqr_gs.solve=function(cqr, b) {
   # solve cqr*x=b where cqr is a qr decomposition of some rectangular
   # complex matrix
   
   # QR*x=b => R*x=Adj(Q)*b => x=inv(R)*Adj(Q)*b
   if (!is.list(cqr)) {
      cqr=cqr_gs(cqr)
   }
   b=as.matrix(b)
   n=dim(b)[1L]
   nb=dim(b)[2L]
   b=Adj(Adj(b)%*%cqr$Q)
   p1=cqr$rank
   R=as.matrix(cqr$R)
   p=dim(R)[2L]
   x=matrix(0., p, nb)
   for (ir in p1:1) {
      if (ir < p) {
         b[1:ir,]=b[1:ir,]-R[1:ir,ir+1L]%o%x[ir+1,]
      }
      x[ir,]=b[ir,]/R[ir,ir]
   }
   return(x)
}
vander=function(x, nc=NROW(x)) {
   # construct vandermonde matrix
   outer(x, 0:(nc-1), `^`)
}
# <-- 2013-05-02 from http://rwiki.sciviews.org/doku.php?id=misc:r_accuracy:high_precision_arithmetic
# the following function is written to work both
# with double precision and with mpfr numbers (cf Rmpfr package)
inv.vandermonde <- function(a)
{
    n <- length(a)
    B <- array(a - a, dim=c(n, n)) # zeros with the same accuracy
    B[1, ] <- 1
    num <- B[1, ]
    for (i in 1:n)
    {
        C <-  - a[i] * B
        C[2:n, ] <- C[2:n, ] + B[1:(n-1), ]
        B[, - i] <- C[, - i]
        p <- a - a[i]
        p[i] <- 1
        num <- num * p
    }
    B / array(rep(num, each=n), dim=c(n, n))
}
# -->
Adj=function(m) {
   # adjoint matrix (Conj + transp)
   return(Conj(t(m)))
}
solve_vander_t=function(v, b) {
   # solve transpose of square vandermonde
   # system (ones are on the first row) t(V(v))*x=b
   # defined by a vector v and rhs b.
   # b might be a matrix of multiple rhs.
   # algorithm of Bjork
   n=length(v)
   if (n==1) {
      return(b)
   }
   nm1=n-1
   b=as.matrix(b)
   for (k in 1:nm1) {
      for (i in n:(k+1)) {
         b[i,]=b[i,]-v[k]*b[i-1,]
      }
   }
   for (k in nm1:1) {
      for (i in (k+1):n) {
         b[i,]=b[i,]/(v[i]-v[i-k])
      }
      for (i in k:nm1) {
         b[i,]=b[i,]-b[i+1,]
      }
   }
   return(b)
}
gsum=function(z,n) {
   # calculate the sum of geometric series
   # 1+z+z^2+...+z^(n-1) (n terms)
   if (z==1) {
      return(n)
   }
   return((1-z^n)/(1-z))
}
eval_formula=function(f, l) {
   # evaluate an expression of rhs for formula f with values from a list l
   f=as.formula(f)
   return(with(l, eval(f[[length(f)]])))
}
# robust estimators
# weighted estimators are used
# weights are calculated with biquadratique function 1/(1+x^2)^2
biquadw=function(v, center=median(v), scale=mad(v,center=center,constant=3)){
   scale=max(scale, max(abs(v-center))/length(v), 1.e-14) # avoid div by 0
   v=(v-center)/scale
   v=1+v*v
   1./(v*v)
}
sigmaw=function(v, center=median(v), scale=mad(v,center)) {
   # weight function in sigmoid form. The sigmoid is obtained with pnorm()
   if (scale==0)
      scale=1
   v=abs(v-center)/scale
   pnorm(2*(v-2.), lower.tail=F)
}
rcov=function(m,...){
   # robust estimation of covariance matrix (covariance of columns)
   m=matrix(m,NROW(m))
   wt=apply(m,2,biquadw);	# weight per element
   wt=apply(wt,1,prod)**(1/NCOL(m));	# weight per row (individum)
   wt=wt/sum(wt)
   cov.wt(m,wt)
}
rz.pval.mono=function(x,...){
   # consider x as gausian sample, estimate its mean and sdev in robust way
   # and return p-values mono-tail.
   wt=biquadw(x,...)
   wt=wt/sum(wt)
   par=cov.wt(matrix(x, NROW(x)),wt)
   pnorm((x-par$center)/sqrt(par$cov), lower.tail=FALSE)
}
rz.pval.bi=function(x, fweight=sigmaw, ...){
   # consider x as gausian sample, estimate its mean and sdev in robust way
   # and return p-values bi-tail (NA is returned for NA in x)
   # treat only valid entries (skip NA)
   iva=!is.na(x)
   res=rep(NA, length(x))
   if (sum(iva) == 0) {
      return(res)
   }
   x=x[iva]
   wt=fweight(x,...)
   wt=wt/sum(wt)
   par=cov.wt(matrix(x, NROW(x)),wt)
   sd=sqrt(par$cov)
   if (sd==0.)
      sd=1.
   res[iva]=pnorm(abs(x-c(par$center))/c(sd), lower.tail=FALSE)*2
   attr(res, "robust_mean") = par$center
   attr(res, "robust_sd") = sd
   attr(res, "wt") = par$wt
   return(res)
}
z.pval.mono=function(x,...){
   # consider x as gausian sample, estimate its mean and sdev
   # and return p-values mono-tail.
   pnorm((x-mean(x))/sd(x), lower.tail=FALSE)
}
rmean_sd=function(x,...) {
   # estimate mean and sd in robust way
   wt=biquadw(x,...)
   wt=wt/sum(wt)
   p=cov.wt(matrix(x, NROW(x)),wt)
   return(list(m=p$center, sd=sqrt(p$cov)))
}
join=function(sep, v, ...) {
   # convert elements of vector v (and all following arguments)
   # in strings and join them using sep as separator.
   return(paste(c(v,...), sep="", collapse=sep))
}
qr_ceq=function(a, b, ceq=NULL, deq=NULL) {
   # solve least square linear system ax=b with equality constrainsts cx=d
   # method: null space base is used to reduce the dimension
   # of original system x=n*y+x_p where x_p is particular solution
   # for cx_p=d
   if (is.null(ceq) || is.null(deq)) {
      # plain qr problem
      return(qr.solve(a, b))
   }
   #library(MASS)
   ns=Null(t(ceq))
   # particular solution of constraints
   x_p=qr.solve(ceq, deq)

   # reduced least square problem
   ar=a%*%ns
   br=b-a%*%x_p
   y=qr.solve(ar, br)
   xr=ns%*%y+x_p
   return(xr)
}
int2bit=function(i, len=31, onech="1", zch="0") {
   # return a string with a bit mask of an integer i
   mask=bitwShiftL(1,(len-1):0)
   paste(ifelse(bitwAnd(i, mask)>0, onech, zch), sep="", collapse="")
}
str2ind=function(qry, ref) {
   # find indexes of strings from qry in ref string vector
   # It can be usefull for retriving row/column indexes
   # from row/column names for a given set of names
   # return integer vector of indexes (NA if an entry is not in ref)
   return(apply(as.matrix(qry), 1, function(s)which(ref==s)[1L]))
}
which.first=function(lv) {
   # return the first index of logical lv vector
   # which is true
   return(head(which(lv),1))
}
which.last=function(lv) {
   # return the last index of logical lv vector
   # which is true
   return(tail(which(lv),1))
}
which.contrib=function(v, thre=0.95) {
   # return the indexes of v who's cumulated sum of squares is
   # over norm2(v)*thre
   v2=v*v
   o=order(v2, decreasing=T)
   cs=cumsum(v2[o])
   i=o[which(cs<=tail(cs,1)*thre)]
   if (length(i)==0) {
      return(head(o,1))
   } else {
      return(i)
   }
}

# custom operators
`%dot%`<-function(a,b)tcrossprod(a,t(b)) # a%dot%b===a%*%b sometimes more efficient than plain %*%
# you might want to set
#`%*%`<-`%dot%` # use rm(`%*%`) to turn back to the classic %*%
# test it for performance before any use in production

`%tmm%`<-function(a,b)crossprod(a,b) # a%tmm%b===t(a)%*%b
`%mmt%`<-function(a,b)tcrossprod(a,b) # a%mmt%b===a%*%t(b)
`%mrv%`<-function(m, v) {
   # multiply each row of m by a vector v term by term
   # length(v) must be == ncol(m)
   if (inherits(m, "simple_triplet_matrix")) {
      m$v=m$v*v[m$j]
      return(m)
   } else {
      return(m*rep(v, rep(nrow(m), length(v))))
   }
}
`%rep%`<-function(v, n) {
   # repeat vector v n times
   rep.int(v, n)
}
`%s+%`<-function(s1, s2) {
   # concatenate two strings by paste()
   paste(s1, s2, sep="")
}
`%s*%`<-function(s, n) {
   # repeat string s n times by paste() with collapse
   paste0(rep(s, n), collapse="")
}

Heaviside=function(x) {
   return(0.5+0.5*sign(x))
}

Nulla=function (M, rcond=1.e10) {
    # use Lapack for null space basis
    # derived from MASS::Null
    tmp <- qr(as.matrix(M), LAPACK=TRUE)
    d=abs(diag(tmp$qr))
    n=length(d)
#browser()
    if (d[1L]==0.) {
       tmp$rank = 0
    } else {
       tmp$rank = sum(d/d[1L] > 1./rcond)
    }
    set <- if (tmp$rank == 0L)
        1L:n
    else
        -(1L:tmp$rank)
    bn=qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
    attr(bn, "qr")=tmp
    return(bn)
}
jordan=function(a, rcond=1.e10) {
   # similar to eigen()
   # additional field is $jordan
   # which is a 2 column integer matrix of eigenvalue indexes relative to $values (first column)
   # and their jordan order (second column, 1 means no jordan block)
   # potential for optimization:
   #  - projections to eigen spaces
   
   # 2012-04-20 sokol: first version
   # 2013-09-12 sokol: jordan info to 2 column matrix
   
   a=as.matrix(a)
   n=nrow(a)
   stopifnot(n==ncol(a))
   if (n==1) {
      return(list(values=c(a), vectors=as.matrix(1.), jordan=t(c(1L, 1L))))
   }
   tol=1./rcond
   i=1:n
   ea=eigen(a)
   u=ea$vectors
   if (typeof(u)=="complex") {
      cu=Conj(u)%tmm%u
   } else {
      cu=crossprod(u)
   }
   # start with no jordan blocks, later if found, extra rows are deleted
   jordan=cbind(i, jorder=1L)
   
   # overdiagonal correlation==1 indicates jordan bloc of order > 1
   dve=which(abs(abs(cu[cbind(i[-n], i[-n]+1)])-1.)<=tol)
   if (length(dve)==0) {
   } else if (length(dve)==1) {
      jordan[dve,2L]=2
      jordan=jordan[-dve-1L,,drop=F] # delete extra row
      i=dve:(dve+1L)
      ea$values[i]=mean(ea$values[i])
   } else {
      # find left and right end of each linear interval in dve
      dd=diff(dve)
      le=c(dve[1L], dve[1L+which(dd!=1L)]) # startin rows
      ri=c(dve[which(dd!=1L)],tail(dve,1L)) # ending rows of the overdiag
      jordan[le,2L]=ri-le+2
      iextra=c()
      for (i in 1L:length(le)) {
         ii=le[i]:ri[i]
         iextra=c(iextra, ii)
         ii=c(ii, ri[i]+1)
         ea$values[ii]=mean(ea$values[ii]) # replace by mean close values
      }
      iextra=1L+iextra
      jordan=jordan[-iextra,,drop=F]
   }
   ea$jordan=jordan
   # get basis vectors for each block
   for (i in 1L:nrow(jordan)) {
      v=jordan[i,]
      if (v[2L] == 1) {
         # not jordan, skip it
         next
      }
      qa=qr(as.matrix(a)-diag(ea$values[v[1L]], n), LAPACK=T)
      # calculate u_j as solution of alam%*%u_j=u_(j-1)
      for (ii in v[1L]+(1L:(v[2L]-1L))) {
         #cat("i=", i, "\n")
         ev=ls_ln(qa, ea$vectors[,ii-1L])
         ea$vectors[,ii] = ev #/sqrt(norm2(ev))
      }
   }
   return(ea)
}
bjord=function(values, jordan) {
   # prepare jordan form with vector values on diagonal and 1 on over diagonal
   # where appropriate (indicated by last column in jordan, cf jordan())
   # For test: aj=jordan(a); a must be close to
   # aj$vectors%*%bjord(aj$values, aj$jordan)%*%solve(aj$vectors)
   res=diag(values, length(values))
   if (length(jordan) == 0) {
      return(res)
   }
   iover=which(jordan[,2L]>1)
   if (length(iover) == 0L) {
      return(res)
   }
   ile=jordan[iover,1L]
   iri=ile+jordan[iover,2L]-2L
   i=unlist(lapply(1L:length(iover), function(i) ile[i]:iri[i]))
   res[cbind(i, i+1L)]=1
   return(res)
}
expjord=function(values, jordan, t=1) {
   # calculate exponential of jordan form exp(J*t)
   # given by diagonal vector values and
   # a jordan structure (cf jordan())
   # time moments are in the vector t.
   
   # return a three dimensional array, first two dimensions
   # are for blocks and the third for time t.
   # 2012-04-21 sokol
   # 2013-09-12 sokol: jordan is passed to a 2 column matrix
   
   n=length(values)
   nt=length(t)
   # get maximal size of jordan block
   maxj=max(jordan[,2L])
   mj1=maxj-1L
   tv=vander(t, maxj)%mrv%(1./cumprod(c(1., iseq(mj1))))
   expvt=exp(values%o%t)
   res=vapply(1L:nt, function(iti) diag(expvt[,iti], n), numeric(n*n))
   if (maxj==1L) {
      # plain diagonal matrices
      return(array(res, dim=c(n, n, nt)))
   }
   # fill overdiagonals by t^i/i!: i=0 on diagonal, i=1 on first overdiagonal and so on
   # make index row-col for overdiag terms, one per ioff
   irc=list()
   for (ioff in iseq(mj1)) {
      ib=which(jordan[,2L]-ioff > 0)
      ile=jordan[ib,1L]
      iri=ile+(jordan[ib,2L]-ioff-1L)
      irc=append(irc, list(unlist(lapply(1L:length(ib), function(i) {
         (ile[i]:iri[i])*(n+1L)+(ioff-1L)*n
      }))))
   }
   for (ioff in iseq(mj1)) {
      for (iti in 1L:nt) {
         i=irc[[ioff]]
         res[i,iti]=tv[iti, ioff+1]
         # multiply by exp
         ie=(i-1)%%n+1
         res[i,iti]=res[i,iti]*expvt[ie,iti]
      }
   }
   if (nt==1) {
      dim(res)=c(n, n)
   } else {
      dim(res)=c(n, n, nt)
   }
   return(res)
}
expmj=function(a, t=1., x=NULL) {
   # calculate a series of exp(a*t) for multiple t
   # by using jordan canonical form
   # Optionally exp(a*t) is multiplied by a vector x.
   # Return a 3d array if x==NULL or a matrix
   # where each column correspond to a time moment t: exp(a*t)*x
   # 2012/04/22 sokol
   
   nt=length(t)
   aj=jordan(a)
   u=aj$vectors
   ej=expjord(aj$values, aj$jordan, t)
   if (is.null(x)) {
      uinv=solve(u)
      res=vapply(1:nt, function(iti) {
         c(u%*%ej[,,iti]%*%uinv)
      }, c(ej[,,1L]))
   } else {
      ux=solve(u, x)
      res=vapply(1:nt, function(iti) {
         c(u%*%(ej[,,iti]%*%ux))
      }, x)
   }
   if (nt==1) {
      res=matrix(res, nrow=nrow(a))
   } else {
      res=array(res, dim=dim(ej))
   }
   return(res)
}
cubic=function(x, x1, x2, y1, y2, k1, k2) {
   # return a matrix of cubic curves, time runs along columns
   t=(x-x1)/(x2-x1)
   a=k1*(x2-x1)-(y2-y1)
   b=-k2*(x2-x1)+(y2-y1)
   q=(1-t)%o%y1 + t%o%y2 + (t*(1.-t))*((1.-t)%o%a + t%o%b)
   return(q)
}
powmx=function(a, pow, x=NULL) {
   # power square matrix a**pow optionnaly multiplyed by a vector x
   # pow is rounded to an integer.
   # due to a binary decomposition of n, it can be less multiplications
   # than n.
    if (!is.null(x)) {
       ans=x
    } else {
       ans="identity"
    }
    pow=round(pow)
    if (pow==0L) {
       if (ans[1L]=="identity") {
          return(diag(nrow(a)))
       }
       return(ans)
    } else if (pow < 0L) {
       a=solve(a)
       pow=-pow
    }
    mmcount=0
    while (pow > 0L) {
        if (pow%%2L) {
           if (ans[1L]=="identity") {
              ans=a
           } else {
              ans=a%*%ans
              mmcount=mmcount+1L
           }
        }
        pow=pow%/%2L
        if (pow > 0L) {
           a=a%*%a
           mmcount=mmcount+1L
        }
    }
    #cat("mmcount=", mmcount, "\n", sep="")
    return(ans)
}
is.diff=function(v1, v2, tol=1.e-7) {
   # test if two numerical vectors of the same length (not checked here)
   # are different. NA are silently omitted.
   # return TRUE if max(abs(v1-v2)) > tol
   # If any of v1 or v2 is NULL, TRUE is returned
   if (is.null(v1) || is.null(v2)) {
      return(TRUE)
   }
   return(max(abs(v1-v2), na.rm=T) > tol)
}
iseq=function(n) {
   # positive sequence of integer numbers 1:n
   # if n=0 then integer(0) is returned
   seq.int(from=1L, to=n, length=pos(n))
}
erf=function(x) 2.*pnorm(x*sqrt(2.))-1.
erfc=function(x) 1.-erf(x)

`%stw%`=function(a, b) {
   # check if the string(s) a starts with the string(s) b
   substr(a, 1L, nchar(b))==b
}
memuse=function(...) {
   # Return a vector of used memory in bytes for each object returned by ls(...)
   # in the calling frame
   # The vector is named by object names. The "..." is passed through to ls()
   env=parent.frame()
   sapply(ls(env=env, ...), function(nam) object.size(get(nam, env=env)))
}
conf=function(x, f) {
   # convolves a signal x with a filter f
   # return a vector of the same length as x
   convolve(x, rev(f), t="o")[length(f)/2+iseq(length(x))]
}
funique=function(x, tol=1.e-7) {
   # like unique() but on float numbers which are considered the same
   # if their abs diff is less or equal to tol
   n=length(x)
   inx=iseq(n)
   # find same numbers (only lower tirangle is accounted)
   #a=abs(outer(x, x, `-`))<=tol
   a=abs(rep(x, each=n) - x)<=tol
   b=outer(inx, inx, `>=`)
   dim(a)=dim(b)
   same=a & b
   # a unique index
   iu=unique(apply(same, 2, function(v) max(which(v))))
   # how many times repeated
   #rep=rowSums(same)
   return(x[iu])
}
adot=function(a1, a2, d1=NULL, d2=NULL) {
   # dot product for two arrays a1 and a2 over
   # dimensions d1 and d2 respectively.
   # If not given, d1 and d2 are take to be
   # the last and first indexes respectively.
   # For example, a classical matrix dot product
   # would result in a call adot(m1, m2, 2, 1) or
   # shorter adot(m1, m2)
   
   # Method: a1 and a2 are permuted to make d1 the last and d2 the
   # first index in the array, then arrays are reshaped to matrices
   # to use a classical dot product then the result is reshaped
   # back to a due array.
   
   dim1=dim(a1)
   nd1=length(dim1)
   dim2=dim(a2)
   if (is.null(dim2)) {
      dim2=c(length(a2), 1L)
   }
   nd2=length(dim2)
   if (is.null(d1)) {
      d1=nd1
   }
   if (is.null(d2)) {
      d2=1L
   }
   stopifnot(dim1[d1] == dim2[d2])
   # permute a1 to make d1 the last index
   if (d1 != nd1 && d1 != 1L) {
      a1=aperm(a1, c((1L:nd1)[-d1],d1))
   }
   a1t=d1==1L
   a2t=F
   if (! a1t) {
      a2t=d2==nd2
   }
   # permute a2 to make d2 the first index
   if (d2 != 1L && !a2t) {
      a2=aperm(a2, c(d2, (1L:nd2)[-d2]))
   }
   # reshape
   if (a1t) {
      # keep running dimension at first place
      dim(a1)=c(dim1[d1], prod(dim1[-d1]))
   } else {
      dim(a1)=c(prod(dim1[-d1]), dim1[d1])
   }
   if (a2t) {
      # keep running dimension at last place
      dim(a2)=c(prod(dim2[-d2]), dim2[d2])
   } else {
      dim(a2)=c(dim2[d2], prod(dim2[-d2]))
   }
   # dot product
   if (a1t) {
      res=a1%tmm%a2
   } else if (a2t) {
      res=a1%mmt%a2
   } else {
      if (dim1[2L]==1L) {
         res=c(a1)*rep(a2, each=dim1[1L])
      } else {
         res=a1%*%a2
      }
   }
   # back-reshape
   #res=as.array(res)
   dim(res)=c(dim1[-d1], dim2[-d2])
   return(res)
}
fmatch=function(x, table, tol=1.e-7) {
   # like match() but works on floating numbers considered as matching
   # if their abs diff is less or equal to tol
   # unmatched values have NA in their corresponding positions
   di=abs(outer(x, table, "-")) <= tol
   dim(di)=c(length(x), length(table))
   res=vapply(1:nrow(di), function(i) which(di[i,])[1L], 1L)
   return(res)
}
extrem_parab=function(x, y) {
   # find the extremum point of a quadratic parabola
   # defined by three distinct points (length(x)==length(y)==3)
   # Return a three component vector named x, y and d2
   # The "d2" component is set to the second derivative,
   # which helps to know if it's a minimum or maximum.
   # In the case of a straight line, d2 is set to 0 and x=y=NA
   # To get the value of this parabola in any point xa from the retrived
   # vector res, one can do
   # res["y"]+0.5*res["d2"]*(xa-res["x"])**2
   # test: extrem_parab(1:3, (1:3)**2) must return c(0, 0, 2)
   
   stopifnot(length(x)==length(y))
   stopifnot(length(x)==3L)
   dx=diff(x)
   dy=diff(y)
   d=dy/dx
   if (d[1L]==d[2L]) {
      # it's a straight line :(
      res=c(x=NA, y=NA, d2=0.)
      return(res)
   }
   
   # diff of middles of the intervals
   dxmid=mean(dx)
   
   # second finite difference
   d2=diff(d)/dxmid
   type=if (d2 > 0.) "min" else "max" # d2==0 is already treated
   xmid=mean(x[1L:2L])
   xe=xmid-d[1L]/d2 # this is the "x"
   # find the first derivative at x1
   dy1=d[1]-d2*dx[1]*0.5
   dxe=xe-x[1]
   ye=y[1]+dy1*dxe+0.5*d2*dxe*dxe
   return(c(x=xe, y=ye, d2=d2))
}
vars2list=function(nm_vars) {
   # Put variables whose names are in the vector nm_vars as items of
   # the returned list.
   stopifnot(is.character(nm_vars))
   res=list()
   for (nm in nm_vars) {
      res[[nm]]=get(nm, env=parent.frame())
   }
   return(res)
}
list2vars=function(li) {
   # Set variables in the calling frame from the items of named list li
   # Return NULL
   stopifnot(is.list(li))
   for (nm in names(li)) {
      assign(nm, li[[nm]], env=parent.frame())
   }
   return(NULL)
}
stop_mes=function(mes="", file=stderr()) {
   # print the error message mes to the error file (by default stderr)
   # and exit with status 1 without environement saving
   # If the session is intercative and the file is not stderr,
   # the message is duplicated on the screen via stop(mes)
   cat(mes, "\n", sep="", file=file)
   if (isatty(stdin()) && file!=stderr()) {
      stop(mes, call.=F)
   } else {
      q("no", status=1)
   }
}
lusolve=function(lua, b, perm=NULL) {
   # solve a normal system a*x=b when lu(a) is already available.
   # b may be a matrix.
   # Return the solution x.
   db=dim(b)
   b=as.matrix(b)
   class(b)="double"
   clua=class(lua)
   if (clua=="sparseLU") {
      # sparse LU is detected
      x=backsolve(lua@U, forwardsolve(lua@L, b[lua@p+1L,,drop=F]))
      x[lua@q+1L,]=x
   } else if (clua=="denseLU") {
      # dense LU here, use lapack
      if (!is.loaded("dgetrs")) {
         lapack.path <- file.path(R.home(), ifelse(.Platform$OS.type == "windows",
            file.path("bin", .Platform$r_arch, "Rlapack"), file.path("lib", "libRlapack")))
         dyn.load(paste(lapack.path,.Platform$dynlib.ext, sep=""))
      }
      n=nrow(b)
      info=-1L
      trans="N"
      res=.Fortran("dgetrs", trans, n, ncol(b), lua@x, n, lua@perm, b, n, info)
      x=res[[7L]]
   } else if (clua=="magmaLU") {
      #cat("magmaLU\n")
      x=magma::solve(lua, b)
   } else if (clua=="matrix" || clua=="numeric") {
      # dense plain matrix
      # pivot is supposed to be given. if null => plain 1,2,3,...
      lua=as.matrix(lua)
      if (is.null(perm)) {
         perm=seq_len(ncol(lua))
      }
      if (!is.loaded("dgetrs")) {
         lapack.path <- file.path(R.home(), ifelse(.Platform$OS.type == "windows",
            file.path("bin", .Platform$r_arch, "Rlapack"), file.path("lib", "libRlapack")))
         dyn.load(paste(lapack.path,.Platform$dynlib.ext, sep=""))
      }
      n=nrow(b)
      info=-1L
      trans="N"
      res=.Fortran("dgetrs", trans, n, ncol(b), lua, n, perm, b, n, info)
      x=res[[7L]]
   } else {
      stop(sprintf("Expected an object of class 'sparseLU', 'denseLU' or 'magmaLU'. Instead, got a '%s' class.", class(lua)))
   }
   dim(x)=db
   return(x)
}
imat2iv=function(imat, dim) {
   # translate matrix indeces to a flat vector ones
   # e.g. if imat is a two column matrix, its first column is supposed
   # to be row numbers and the second one the column numbers.
   # the results will be a vector of indeces pointing to the same
   # matrix elements as imat but treating mat as a flat vector
   # (column wise stored as usual in R)
   # length(dim) must be equal to ncol(imat)
   # dim is a vector of dimension in a matrix to which imat is supposed
   # point to.
   # Return a vector iv with length(iv)=nrow(imat)
   # NB If dim(imat) is not set, i.e. it is a plain vector,
   # it is returned as is.
   # If dim(imat) is set, it's length must equal 2, i.e. imat must be
   # a matrix, not nD array. if the matrix imat has n columns, its elements
   # are pointing to a nD array.
   
   stopifnot(ncol(imat)==length(dim))
   if (is.null(dim(imat))) {
      return(imat)
   }
   stopifnot(length(dim(imat))==2L)
   imat[,-1L]=imat[,-1L]-1L
   iv=rowSums(imat%mrv%c(1L, cumprod(dim[-length(dim)])))
   return(iv)
}
dgemm=function(a, b, transa=F, transb=F, alpha=1., beta=0., c=matrix(0., ifelse(transa, ncol(a), nrow(a)), ifelse(transb, nrow(b), ncol(b)))) {
   if (!is.loaded("dgemm")) {
      lapack.path <- file.path(R.home(), ifelse(.Platform$OS.type == "windows",
         file.path("bin", "Rlapack"), file.path("lib", "libRlapack")))
      #dyn.load(paste(lapack.path,.Platform$dynlib.ext, sep=""))
      #dyn.load("/opt/intel/mkl/lib/intel64/libmkl_gf_ilp64.so")
      #dyn.load("/opt/intel/mkl/lib/intel64/libmkl_gnu_thread.so")
      #dyn.load("/opt/intel/mkl/lib/intel64/libmkl_core.so")
      #dyn.load("/usr/local/src/xianyi-OpenBLAS-9c51cdf/libopenblas.so")
      #dyn.load("/opt/OpenBLAS/lib/libopenblas.so")
      #dyn.load("/usr/lib64/atlas/libf77blas.so")
      #dyn.load("/usr/local/atlas/lib64/libsatlas.so")
      #dyn.load("/usr/local/src/ATLAS-3.11.22/build/lib/libsatlas.so")
      #dyn.load("/usr/local/atlas/lib64/libf77blas.so")
      #dyn.load("/usr/local/src/ATLAS-3.10.2/build/lib/libsatlas.so")
      dyn.load("/usr/local/src/ATLAS-3.11.39/build/lib/libsatlas.so")
   }
   m=ifelse(transa, ncol(a), nrow(a))
   n=ifelse(transb, nrow(b), ncol(b))
   k=ifelse(transa, nrow(a), ncol(a))
   res=.Fortran("dgemm", ifelse(transa, "T", "N"), ifelse(transb, "T", "N"), m, n, k, as.double(alpha), as.double(a), nrow(a), as.double(b), nrow(b), as.double(beta), as.double(c), nrow(c))
   c=res[[12L]]
   dim(c)=c(m, n)
   #dimnames(c)=list(if (transa) colnames(a) else rownames(a), if (transb) rownames(b) else colnames(b))
   return(c)
}
ci.plot=function(x, y=NULL, ylow=NULL, yup=NULL, length=0.1, ...) {
   # error bars
   if (missing(yup)) {
      yup=ylow
      ylow=y
      y=x
      x=seq_len(length(y))
   }
   if (!is.null(yup))
      arrows(x, y, x, yup, angle=90, code=2, length=length, ...)
   if (!is.null(ylow))
      arrows(x, y, x, ylow, angle=90, code=2, length=length, ...)
}
