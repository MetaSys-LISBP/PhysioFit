# define cBind(), rBind() if Matrix is not loaded
if (length(find("cBind")) == 0) {
   #cat("cBind not found\n")
   cBind=cbind
   rBind=rbind
}
nlsic=function(par, r, u=NULL, co=NULL, control=list(), e=NULL, eco=NULL, flsi=lsi,...) {
   # solve non linear least square problem min ||r(par,...)$res||
   # with inequality constraints u*par>=co
   # and optional equality constraints e*par=eco
   # Parameters:
   # - par is initial values for parameter vector (can be in non feasible domain)
   # - r is a function calculating residual vector
   # by a call to r(par, cjac=T|F, ...)
   # par is a current paramenter vector,
   # ... are params passed through to r()
   # and cjac is set to TRUE if jacobian is required
   # The call to r() must return a list having a field "res" containing the residual vector and optionnaly a field "jacobian"
   # when cjac is set to TRUE.
   # - u linear constrant matrix in u*par>=co
   # - co constraint vector
   # controls (=by default values):
   # - errx=1.e-7 error on l2 norm of the iteration step sqrt(pt*p).
   # - maxit=100 maximum of newton iterations
   # - btfrac=0.5 (0;1) by this value we diminish the step till tight up
   # to the quadratic model of norm reduction in backtrack (bt) iterations
   # - btdesc=0.5 (0;1) how good we have to tight up to the quadratic model
   # 0-we are very relaxe, 1 we are very close (may need many bt iterations)
   # - btmaxit=15 maximum of backtrack iterations
   # - rcond=1.e10 maximal condition number in QR decomposition
   #   starting from which a matrix is considered as numerically rank
   #   deficient.
   # - ci = list of options relative to confidence interval estimation
   # - history=TRUE report or not the convergence history
   # - adaptbt=FALSE use or not adaptive backtracking
   # Output:
   # - list(par, <some convergence information>)
   # Method:
   # sequential lsi globalized by backtracking technique.
   # If e, eco are not NULL, reduce jacobian before lsi() call.
   # Notes:
   # If function r() return a list having field "jacobian" it is
   # supposed to be equal to jaconian dr_dpar.
   # If not, numerical derivation numDeriv::jacobian()
   # is used for its construction.
   # If equality constraints are present, MASS library is necessary
   # for Null() function.
   
   n=length(par)
   m=NROW(co)
   nm_par=names(par)
   
   
   # predefined controls, then overwritten by actual values
   con=list(errx=1.e-7, maxit=100, btstart=1., btfrac=0.5, btdesc=0.5,
      btmaxit=15, btkmin=1.e-7, trace=0, rcond=1.e10, ci=list(p=0.95, report=T),
      history=F, adaptbt=F)
   nmsC=names(con)
   nmsci=names(con$ci)
   # predefined confidence interval parameters, then overwritten
   ci=con$ci
   if (!is.null(control$ci)) {
      ci[(namci <- names(control$ci))] <- control$ci
   } else {
      namci <- names(con$ci)
   }
   con[(namc <- names(control))] <- control
   con$ci=ci
   if (length(noNms <- namc[!namc %in% nmsC])) {
      warning("nlsic: unknown names in control: ", paste(noNms, collapse = ", "))
   }
   if (length(noNms <- namci[!namci %in% nmsci])) {
      warning("nlsic: unknown names in control$ci: ", paste(noNms, collapse = ", "))
   }
   
   # history matrices
   hpar=hp=hstep=hres=c()
   if (con$history) {
      hist=list()
   } else {
      hist=NULL
   }
   
   if (!is.null(e)) {
      e=matrix(e, ncol=n)
      ne=nrow(e)
   } else {
      ne=0
   }
   econstr=!is.null(e) && !is.null(eco) && nrow(e)>0 && ncol(e)==length(par)

   # make u matrix if u=NULL
   if (is.null(u)) {
      u=matrix(0, nrow=0, ncol=n)
      co=c()
   }
   if (econstr) {
      # affine transform epar -> par s.t. e*par=eco
      # particular solution
      parp=qr.solve(e, eco)
      # par=parp+Null(t(e))%*%epar, where epar is new parameter vector
      nte=Null(t(e))
      ue=u%*%nte
      epar=qr.solve(nte, par-parp)
      par[]=nte%*%epar+parp
   }
   ine=co-u%*%par
   if (any(ine > 1./con$rcond)) {
      # project par into feasible domain
      # in case if residual function is defined only on feasible domain
      if (econstr) {
         x=ldp(ue, ine, con$rcond)
         if (is.null(x)) {
            return(list(
               par=NULL,
               error=1,
               mes="nlsic: Unfeasible constraints at starting point"
            ))
         }
         laststep=nte%*%x
         par[]=par+laststep
      } else {
         x=ldp(u, ine, con$rcond)
         if (is.null(x)) {
            mes="nlsic: Unfeasible constraints at starting point"
            if (!is.null(rownames(u))) {
                  mes=join("\n", c(mes, rownames(u)))
            }
            return(list(
               par=NULL,
               error=1,
               mes=mes
            ))
         }
         laststep=x
         par[]=par+laststep
      }
   }
   
   # newton globalized iterations
   it=0
   btit=0
   converged=F
   while (!converged && it < con$maxit) {
      mes=""
      if (it == 0) {
         # residual vector
         lres=r(par, cjac=T, ...)
         if ((!is.null(lres$err) && lres$err) || any(is.na(lres$res))) {
            return(list(
               par=c(par),
               retres=lres,
               it=it,
               btit=btit,
               error=1,
               mes=paste("nlsic: Problem in first residual calculation",
                  lres$mes, sep="\n"))
            )
         }
         res=lres$res
         b=-res
         norb=sqrt(norm2(b))
         #cat("norb0=", norb, "\n", sep="")
         norres=norb
         klast=0
         if (con$trace) {
            cat("it=", it, "\tres=", norb, "\n", sep="")
         }
      } else {
         b=-res
         norb=norres
      }
      # jacobian
      if ("jacobian" %in% names(lres)) {
         a=r(par, cjac=TRUE, ...)$jacobian
      } else {
         # needs something like library(numDeriv)
         a=jacobian(function(x)r(x, cjac=F, ...)$res, par)
      }
      # solve linear ls
      if (econstr) {
         p=flsi(a%*%nte, b, ue, co-u%*%par, rcond=con$rcond); # newton direction
         p=nte%*%p
      } else {
         p=flsi(a, b, u, co-u%*%par, rcond=con$rcond); # newton direction
      }
      names(p)=nm_par
#browser()
      if (con$history) {
         hpar=cBind(hpar, par)
         hp=cBind(hp, c(p))
         hres=cBind(hres, res)
         hist=list(par=hpar, p=hp, res=hres, step=hstep)
      }
      if (any(is.na(p))) {
         names(par)=nm_par
         return(list(
            par=c(par),
            it=it,
            btit=btit,
            error=1,
            history=hist,
            jacobian=a,
            mes=paste("nlsic: Problem in solving linear problem with constraint",
            attr(p, "mes"), sep="\n")))
      } else if (! is.null(attr(p, "mes"))) {
         mes=attr(p, "mes")
      }
      
      normp=sqrt(norm2(p))
      converged=(normp <= con$errx)
      if (converged) {
         # no need for backtracking at this stage
         laststep=p
         par[]=par+laststep
         it=it+1
         btit=0
         resdecr=NULL
         res=r(par, cjac=FALSE, ...)$res
         norres=sqrt(norm2(res))
         if (con$history) {
            hstep=cBind(hstep, c(laststep))
            hist=list(par=hpar, p=hp, res=hres, step=hstep)
         }
         if (con$trace && (it < 10 || !it%%10)) {
            cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
         }
         break
      }
      # check if the direction p is descending
      ap=as.double(a%*%p)
      n2ap=norm2(ap)
      resap=crossprod(res, ap)
      if (resap > 0.) {
         laststep=p*NA
         mes=paste("nlsic: LSI returned not descending direction",
            attr(p, "mes"), sep="\n")
         break
      }
      k=con$btstart; # fraction of p
      # backtrack iterations
      btit=0
      descending=F
      while (!descending && btit < con$btmaxit && k >= con$btkmin) {
         laststep=k*p
         lres=r(par+c(laststep), cjac=FALSE, ...)
         if ((!is.null(lres$err) && lres$err) || any(is.na(lres$res))) {
            return(list(
               par=c(par),
               laststep=c(laststep),
               retres=lres,
               it=it,
               btit=btit,
               error=1,
               mes=paste("nlsic: Problem in residual calculation",
                  lres$mes, sep="\n"))
            )
         }
         res=lres$res
         norres=sqrt(norm2(res))
         #cat("norres=", norres, "\n", sep="")
         scaprod=crossprod(res, ap)-resap
         #resdecr=sqrt(scaprod); # residual decreasing
         descending=(
            scaprod >=0. &&
            sqrt(scaprod) >= con$btdesc*sqrt(n2ap*k) &&
            norres < norb
         )
         klast=k
         if (con$adaptbt) {
            h=max(min(2., 1./(2.-(norres-norb)*(norres+norb)/(resap*k))), con$btfrac)
            k=k*h
            #cat("h=", h, "klast=", klast, "k=", k, "\n")
         } else {
            k=k*con$btfrac
         }
         k=max(k, con$btkmin)
         btit=btit+1
      }
      par[]=par+laststep
      it=it+1
      if (con$trace && (it < 10 || !it%%10)) {
         cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
      }
      if (con$history) {
         hstep=cBind(hstep, laststep)
         hist=list(par=hpar, p=hp, res=hres, step=hstep)
      }
   }
   if (con$trace && !(it < 10 || !it%%10)) {
      cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
   }
   if (it >= con$maxit) {
      mes=paste(mes, "nlsic: Maximal non linear iteration number is achieved", sep="\n")
   }
   if (btit >= con$btmaxit) {
      mes=paste(mes, "nlsic: Maximal backtrack iteration number is achieved", sep="\n")
   }
   # restore names
   names(par)=names(laststep)=names(p)=nm_par
   # confidence interval
   #browser()
   nr=length(res)
   fdeg=max(1., nr-n+ne) # freedom degrees
   if (con$ci$report) {
      sd_res=sqrt(sum(res**2)/fdeg)
      if (econstr) {
         jac=a%*%nte
      } else {
         jac=a
      }
      d=try(sqrt(diag(solve(crossprod(jac)))), silent=T)
      if (inherits(d, "try-error")) {
         sv=svd(jac)
         d=sqrt(diag(sv$v%*%diag(1./pmax(sv$d, 1./con$rcond)**2, length(sv$d))%*%t(sv$v)))
      }
      if (econstr) {
         d=abs(as.numeric(nte%*%d))
      }
      hci=(-qt((1.-con$ci$p)/2., fdeg)*sd_res)*d
      names(hci)=nm_par
   } else {
      sd_res=hci=NULL
   }
   return(list(
      par=c(par),
      lastp=c(p),
      hci=hci,
      ci_p=con$ci$p,
      sd_res=sd_res,
      laststep=c(laststep),
      normp=normp,
      res=c(res),
      prevres=-c(b),
      jacobian=a,
      retres=lres,
      it=it,
      btit=btit,
      history=hist,
      error=0,
      mes=mes))
}
lsi=function(a, b, u=NULL, co=NULL, rcond=1.e10) {
   # solve linear least square problem (min ||Ax-b||)
   # with inequality constraints ux>=co
   # Method:
   # 1. reduce the pb to ldp (min(xat*xa) => least distance programming)
   # 2. solve ldp
   # 3. change back to x
   
   if (class(a)!="qr") {
      n=ncol(a)
      aq=qr(a, LAPACK=T)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=diag(aq$qr)
   aq$rank=sum(abs(d)>abs(d[1])/rcond)
   if (aq$rank < n) {
      x=rep(NA, n)
      attr(x, "mes")=
         paste("lsi: Rank deficient matrix in least squares\n",
         n-aq$rank,
         " unsolvable variable(s):\n",
         paste(dimnames(a)[[2]][aq$pivot[-(1:aq$rank)]],
         aq$pivot[-(1:aq$rank)], sep="\t", collapse="\n"), "\n",
         sep="")
      return(x)
   }
   x0=qr.solve(aq, b)
   if (!is.null(u) && nrow(u)>0) {
      # we do have inequalities
      # prepare variable change
      ut=t(backsolve(aq$qr, t(u[, aq$pivot, drop=F]), trans=T))
      xa=ldp(ut, co-u%*%x0, rcond)
      if (is.null(xa)) {
         x=rep(NA, n)
         attr(x, "mes")="lsi: Infeasible constraints detected by ldp()\n"
         return(x)
      }
      x=xa
      xa=backsolve(aq$qr, xa)
      x[aq$pivot]=xa
      x=x+x0
      cou=co-u%*%x
      if (any(cou > 1./rcond)) {
         # round errors made fail inequality constraints
         # force them as much as possible even by degrading the residual
         x=x+ldp(u, cou, rcond)
      }
   } else {
      # plain QR
      x=x0
   }
   return(x)
}
lsi_ln=function(a, b, u=NULL, co=NULL, rcond=1.e10) {
   # solve linear least square problem (min ||Ax-b||)
   # with inequality constraints ux>=co
   # If A is rank deficient, least norm solution is used
   # to reduce the pb to lsi with not rank deficient matrix
   
   mes=""
   if (class(a)!="qr") {
      n=ncol(a)
      aq=qr(a, LAPACK=T)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=diag(aq$qr)
   aq$rank=sum(abs(d)>abs(d[1])/rcond)
   rdefic=aq$rank < n
   if (!rdefic) {
      # just passthrough params to lsi()
      return(lsi(aq, b, u, co, rcond))
   }
   # prepare free variable substitution
   ndefi=n-aq$rank
   i1=1:aq$rank
   i2=(aq$rank+1):n
   r=qr.R(aq)[i1,,drop=F]
   r1=r[,i1,drop=F]
   r2=r[,i2,drop=F]
   qrr=qr(t(r))
   # null space and its complement bases
   ba=qr.Q(qrr, complete=T)
   bal=ba[,i1,drop=F]
   ban=ba[,i2,drop=F]
   # least norm without constraints
   bt=qr.qty(aq, b)
   #x0=bal%*%backsolve(qr.R(qrr), bt, trans=T)
   x0=backsolve(r1, bt) # implicitly, free variables are set to zero
   x=c(x0, rep(0., ndefi)) # we unpivot just before return
   #browser()
   if (!is.null(u) && nrow(u)>0) {
      # solve ldp with null part scaled by sqrt(rcond)
      up=u[,aq$pivot,drop=F]
      cou=co-up%*%x
      if (any(cou>0.)) {
         sca=sqrt(rcond)
         uz=(up%*%ban)*sca
         ut=cBind(t(backsolve(r1, t(u[, aq$pivot[i1], drop=F]), trans=T)),uz)
         ya=ldp(ut, cou, rcond)
         if (!is.null(ya)) {
            x[i1]=x[i1]+backsolve(r1, ya[i1])
            x=x+(ban%*%ya[i2])*sca
         }
         cou=co-up%*%x
      }
      # now minimize the norm of x while keeping the norm of ya
      #z0=crossprod(ban, x)
      #x_ln=x-ban%*%z0
      #cot=co-up%*%x_ln
      #if (any(cot>0.)) {
      #   z=ldp(uz, cot, rcond)
      #   if (!is.null(z)) {
      #      x=x_ln+ban%*%z
      #   }
      #} else {
      #   x=x_ln
      #}
      if (any(cou > 1./rcond)) {
         # round errors might make fail inequality constraints
         # force them as much as possible even by degrading the residual
         dx=ldp(up, cou, rcond)
         if (!is.null(dx)) {
            x=x+dx
         } else {
            mes="lsi_ln: Unstable ldp() solution"
         }
      }
   }
   x[aq$pivot]=x
   attr(x, "mes")=
   paste("lsi_ln: Rank deficient matrix in least squares\n",
      ndefi,
      " free variable(s):\n",
      paste(dimnames(a)[[2]][aq$pivot[i2]],
      aq$pivot[-(1:aq$rank)], sep="\t", collapse="\n"),
      "\nLeast L2-norm solution is provided.\n",
      mes, "\n",
      sep="")
   return(x)
}
ldp=function(u, co, rcond=1.e10) {
   # solve least distance programing: find x satisfying u*x>=co and s.t. min(||x||)
   # by passing to nnls
   # (non negative least square) problem
   # return x or Null (in case of non-feasible inequalities)
   m=NROW(u)
   n=ncol(u)
   maxu=max(abs(u))
   maxco=max(abs(co))
   rcond=abs(rcond)
   if (m==0 || (maxu==0. && maxco==0.) || max(co)<=0.) {
      # trivial case
      x=rep(0., n)
      return(x)
   }
   e=rBind(t(u), t(co))
   f=c(rep(0., n), 1.)
   resnnls=nnls(e, f)
   feasible=sqrt(resnnls$deviance) > 1./rcond && resnnls$residuals[n+1] != 0.
   if (feasible) {
      x=resnnls$residuals[1:n]/(-resnnls$residuals[n+1])
      # check for numerical stability problems
      ux=u%*%x
      cou=co-ux
      if (any(cou > 1./rcond)) {
         # second trial
         e[n+1,]=cou
         rn=nnls(e, f)
         if (rn$residuals[n+1]!=0.) {
            x=x-rn$residuals[1:n]/rn$residuals[n+1]
         } else {
            x=NULL
         }
      } else if (all(cou<0) && !all(x==0.)) {
         # round error pushed the solution inside of feasible domain
         # shorten it till the most close border
         alpha=ux/co
         alpha=min(alpha[alpha>=1.], na.rm=T)
         x=x/alpha
      }
   } else {
      x=NULL
   }
   return(x)
}
if (length(find("norm2"))==0) {
   norm2=function(v)crossprod(as.numeric(v))[1]
}
ls_ln=function(a, b, rcond=1.e10) {
   # least square with least norm
   if (class(a)!="qr") {
      a=qr(a, LAPACK=T)
   }
   n=ncol(a$qr)
   d=diag(a$qr)
   a$rank=sum(abs(d)>abs(d[1])/rcond)
   rdefic=a$rank < n
   i1=1:a$rank
   i2=(a$rank+1):n
   r=qr.R(a)[i1,,drop=F]
   if (!rdefic) {
      # plain ls
      x=backsolve(r, qr.qty(a, b))
      x[a$pivot]=x
      return(x)
   }
   # prepare free variable substitution
   r1=r[,i1,drop=F]
   r2=r[,i2,drop=F]
   qrr=qr(t(r))
   # null space and its complement bases
   ba=qr.Q(qrr, complete=T)
   bal=ba[,i1,drop=F]
   # least norm
   x=bal%*%backsolve(qr.R(qrr), qr.qty(a, b), trans=T)
   x[a$pivot]=x
   return(x)
}
lsi_reg=function(a, b, u, co, rcond=1.e10) {
   # regularized lsi problem
   # regularization is made by concatenating diag(lam,n) to a
   # and calling lsi
   if (class(a)!="qr") {
      n=ncol(a)
      aq=qr(a, LAPACK=T)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=diag(aq$qr)
   aq$rank=sum(abs(d)>abs(d[1])/rcond)
   rdefic=aq$rank < n
   if (!rdefic) {
      # just passthrough params to lsi()
      return(lsi(aq, b, u, co, rcond))
   }
   
   # back to a
   if (class(a)=="qr") {
      p=a$pivot
      a=qr.Q(a)%*%qr.R(a)
      a[,p]=a
   }
   lam=1./sqrt(rcond)
   a=rBind(a, diag(lam*d[1], n))
   return(lsi(a, c(b, rep(0., n)), u, co, rcond))
}
