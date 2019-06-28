# Copyright 2011-2014, INRA, France
# Author: Serguei SOKOL, sokol@insa-toulouse.fr
# License: GPL v2 http://www.gnu.org/licenses/gpl-2.0.html

# define cBind(), rBind() if Matrix is not loaded
if (!exists("cBind", mode="function")) {
   #cat("cBind not found\n")
   cBind=cbind
   rBind=rbind
}
if (!exists("join", mode="function")) {
   #cat("cBind not found\n")
   join=function(sep, v, ...) {
      # convert elements of vector v (and all following arguments)
      # in strings and join them using sep as separator.
      return(paste(c(v,...), sep="", collapse=sep))
   }
}

nlsic=function(par, r, u=NULL, co=NULL, control=list(), e=NULL, eco=NULL, flsi=lsi,...) {
   # solve non linear least square problem min ||r(par,...)$res||
   # with optional inequality constraints u*par>=co
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
   # - btstart=1 staring value for backtracking coefficient coefficient
   # - btfrac=0.5 (0;1) by this value we diminish the step till tight up
   # to the quadratic model of norm reduction in backtrack (bt) iterations
   # - btdesc=0.1 (0;1) how good we have to tight up to the quadratic model
   # 0-we are very relaxe, 1 we are very tight (may need many bt iterations)
   # - btmaxit=15 maximum of backtrack iterations
   # - btkmin=1.e-7 a floor value for backtracking fractioning
   # - trace=0 no information is printed during iterations (1 - print is active)
   # - rcond=1.e10 maximal condition number in QR decomposition
   #   starting from which a matrix is considered as numerically rank
   #   deficient. The inverse of this number is also used as a measure of very
   #   small number.
   # - ci = list of options relative to confidence interval estimation
   #  + p=0.95 p-value of confidence interval. If you need a plain sd value,
   #    set p-value to 0.6826895
   #  + report=T report (or not and hence calculate or not) confidence
   #       interval information (in the field hci, as 'half confidence interval' width)
   # - history=TRUE report or not the convergence history
   # - adaptbt=FALSE use or not adaptive backtracking
   # - mnorm=NULL a norm matrix for a case sln==TRUE (cf next parameter)
   # - sln=FALSE use or not (default) least norm of the solution
   #   (instead of least norm of the increment)
   # - maxstep=NULL a maximal norm for increment step (if not NULL), must be positive
   # - monotone=FALSE should or not the cost decrease be monotone. If TRUE, then at
   #    first non decrease of the cost, the iterations are stopped with a warning message.
   # Output:
   # - list(par, <some convergence information>)
   # Method:
   # sequential lsi globalized by backtracking technique.
   # If e, eco are not NULL, reduce jacobian before lsi() call.
   # Notes:
   # If the function r() returns a list having a field "jacobian" it is
   # supposed to be equal to the jacobian dr/dpar.
   # If not, numerical derivation numDeriv::jacobian()
   # is used for its construction.
   
   n=length(par)
   m=NROW(co)
   nm_par=names(par)
   
   # predefined controls, then overwritten by actual values
   con=list(errx=1.e-7, maxit=100, btstart=1., btfrac=0.5, btdesc=0.1,
      btmaxit=15, btkmin=1.e-7, trace=0, rcond=1.e10, ci=list(p=0.95, report=T),
      history=F, adaptbt=F, mnorm=NULL, sln=FALSE, maxstep=NULL, monotone=FALSE)
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
   mes=NULL
   if (length(noNms <- namc[!namc %in% nmsC])) {
      mes=join("", "nlsic: unknown name(s) in control: ", join(", ", noNms))
      return(list(par=NULL, error=1, mes=mes))
   }
   if (length(noNms <- namci[!namci %in% nmsci])) {
      mes=join("", "nlsic: unknown names in control$ci: ", join(", ", noNms))
      return(list(par=NULL, error=1, mes=mes))
   }
   if (!is.null(con$maxstep) && con$maxstep <= 0.) {
      mes=sprintf("nlsic: maxstep parameter in control list must be positive, got %g instead", con$maxstep)
      return(list(par=NULL, error=1, mes=mes))
   }
   if (!is.null(con$maxstep) && con$maxstep <= con$errx) {
      mes=sprintf("nlsic: maxstep parameter in control list is less than errx (maxstep=%g, errx=%g)", con$maxstep, con$errx)
      return(list(par=NULL, error=1, mes=mes))
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
      nr=nrow(e)
      nc=ncol(e)
      if (nr > nc) {
         return(list(par=NULL, error=1, mes="nlsic: matrix e is over determined."))
      }
      # affine transform epar -> par s.t. e*par=eco
      # par=parp+Null(t(e))%*%epar, where epar is new parameter vector
      nte=Nulla(t(e))
      qe=attr(nte, "qr")
      if (qe$rank < nr) {
         return(list(par=NULL, error=1, mes="nlsic: matrix e is rank deficient."))
      }
      # particular solution
      parp=qr.qy(qe, c(backsolve(qe$qr, eco[qe$pivot], trans=T), double(n-nr)))
      ue=u%*%nte
      epar=crossprod(nte, c(par)-parp)
      par[]=nte%*%epar+parp
   }
   ine=if (nrow(u) > 0) co-u%*%c(par) else NULL
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
         par[]=c(par)+laststep
      } else {
         x=ldp(u, ine, con$rcond)
         if (is.null(x)) {
            mes="nlsic: Unfeasible constraints at starting point (u%*%par-co>0):"
            if (!is.null(rownames(u))) {
               mes=join("\n", c(mes, paste(rownames(u), " (", -ine, ")", sep="")))
            }
            return(list(
               par=NULL,
               error=1,
               mes=mes
            ))
         }
         laststep=x
         par[]=c(par)+laststep
      }
   }
   
   # newton globalized iterations
   it=0
   btit=0
   converged=F
   laststep=par
   laststep[]=0.
   p=laststep
   while (!converged && it < con$maxit) {
      mes=NULL
      if (it == 0) {
         # residual vector
         lres=r(par, cjac=T, ...)
         if (!is.null(lres$err) && lres$err) {
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
         res=as.numeric(lres$res)
         iva=which(!is.na(res)) # valid entries in res
         if (length(iva)==0) {
            # we check it only at it=0 hoping that NA-structure remains the same
            # all the long ietrations
            return(list(
               par=NULL,
               error=1,
               mes="nlsic: no valid residual value."
            ))
         }
         vres=res[iva]
         b=-res
         norb=sqrt(norm2(vres))
         #cat("norb0=", norb, "\n", sep="")
         norres=norb
         klast=0
         if (con$trace) {
            cat("it=", it, "\tres=", norb, "\n", sep="")
         }
         if (norb < 1.e-20) {
            # we are already at solution point
            converged=TRUE
            normp=0
            a=NULL
            next
         }
      } else {
         b=-res
         norb=norres
      }
      # jacobian
      a=NULL
      if ("jacobian" %in% names(lres) && !is.null(lres$jacobian)) {
         a=lres$jacobian[iva,,drop=FALSE]
         colnames(a)=nm_par
      } else {
         # may be we are here because the last call was with cjac=FALSE
         if (it > 0) {
            lres=r(par, cjac=TRUE, ...)
            res=as.numeric(lres$res)
            iva=!is.na(res)
            vres=res[iva]
            b=-res
            norres=sqrt(norm2(vres))
            norb=norres
            if ("jacobian" %in% names(lres) && !is.null(lres$jacobian)) {
               a=lres$jacobian[iva,,drop=FALSE]
               colnames(a)=nm_par
            }
         }
      }
      if (is.null(a)) {
         # needs something like library(numDeriv) before calling nlsic()
         a=jacobian(function(x) as.numeric(r(x, cjac=F, ...)$res), par)[iva,,drop=FALSE]
         colnames(a)=nm_par
      }
      # solve linear ls
      
      ine=if (nrow(u) > 0) co-u%*%c(par) else NULL
      # newton direction
      if (econstr) {
         if (con$sln) {
            p=flsi(a%*%nte, -vres, ue, ine, rcond=con$rcond, con$mnorm, -c(par))
         } else {
            p=flsi(a%*%nte, -vres, ue, ine, rcond=con$rcond)
         }
         p=nte%*%p
      } else {
         if (con$sln) {
            p=flsi(a, -vres, u, ine, rcond=con$rcond, con$mnorm, -c(par))
         } else {
            p=flsi(a, -vres, u, ine, rcond=con$rcond)
         }
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
            retres=lres,
            mes=paste("nlsic: Problem in solving linear problem with constraint",
               attr(p, "mes"), sep="\n")))
      } else if (! is.null(attr(p, "mes"))) {
         mes=join("\n", mes, attr(p, "mes"))
      }
      
      normp=sqrt(norm2(p))
      converged=(normp <= con$errx)
      if (normp == 0) {
         ############### no need for backtracking at this stage
         laststep=p
         par[]=c(par)+laststep
         it=it+1
         btit=0
         resdecr=NULL
         res=as.numeric(r(par, cjac=FALSE, ...)$res)
         iva=!is.na(res)
         vres=res[iva]
         norres=sqrt(norm2(vres))
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
      resap=crossprod(vres, ap)
      if (resap > 0.) {
         laststep=p*NA
         mes=paste("nlsic: LSI returned not descending direction",
            attr(p, "mes"), sep="\n")
         break
      }
      k=con$btstart; # fraction of p
      if (!is.null(con$maxstep)) {
         kms=con$maxstep/normp
         if (kms < k) {
            k=kms
         }
      }
      # backtrack iterations
      btit=0
      descending=F
      while (!descending && btit < con$btmaxit && k >= con$btkmin) {
         laststep=k*p
         lres=r(par+c(laststep), cjac=FALSE, ...)
         if (!is.null(lres$err) && lres$err) {
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
         res=as.numeric(lres$res)
         iva=!is.na(res) # valid entries in res
         vres=res[iva]
         norres=sqrt(norm2(vres))
         #cat("norres=", norres, "\n", sep="")
         #scaprod=crossprod(vres, ap)-resap
         #resdecr=sqrt(scaprod); # residual decreasing
         realdecr=crossprod(b[iva]-vres, vres+b[iva]) # normally it is positive
         lindecr=-k*(2*resap+k*n2ap)
         #descending=(
         #   scaprod >=0. &&
         #   sqrt(scaprod) >= con$btdesc*sqrt(n2ap*k) &&
         #   norres < norb
         #)
         descending=(realdecr>=con$btdesc*lindecr)
         klast=k
         if (con$adaptbt) {
            h=max(min(2., 1./(2.-(norres-norb)*(norres+norb)/(resap*k))), con$btfrac)
            # check that p*k*h won't violate inequalities (for h > 1 only)
            if (h > 1. && nrow(u) > 0) {
               if (any(u%*%(c(par)+(k*h)*p)-co <= -1.e-10)) {
                  h=con$btfrac
               }
            }
            k=k*h
            #cat("h=", h, "klast=", klast, "k=", k, "\n")
         } else {
            k=k*con$btfrac
         }
         k=max(k, con$btkmin)
         btit=btit+1
      }
      par[]=c(par)+laststep
      it=it+1
      if (con$monotone) {
         if (norres > norb) {
            converged=T
            mes=join("\n", mes, "nlsic: Non monotone descrease of cost function has occured.")
         }
      }
      if (con$trace && ((it < 10 || !it%%10) || converged)) {
         cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
      }
      if (con$history) {
         hstep=cBind(hstep, laststep)
         hist=list(par=hpar, p=hp, res=hres, step=hstep)
      }
   }
   if (con$trace && !(it < 10 || !it%%10) && !converged) {
      cat("it=", it, "\tres=", norres, "\tnormstep=", normp, "\tbtk=", klast, "\n", sep="")
   }
   if (it >= con$maxit) {
      mes=join("\n", mes, "nlsic: Maximal non linear iteration number is achieved")
   }
   if (btit >= con$btmaxit) {
      mes=join("\n", mes, "nlsic: Maximal backtrack iteration number is achieved")
   }
   # restore names
   names(par)=names(laststep)=names(p)=nm_par
   # confidence interval
   #browser()
   nr=length(vres)
   fdeg=max(1., nr-n+ne) # freedom degrees
   if (con$ci$report) {
      sd_res=sqrt(sum(vres**2)/fdeg)
      if (is.null(a)) {
         a=r(par, TRUE, ...)$jacobian
         if (is.null(a)) {
            # needs something like library(numDeriv) before calling nlsic()
            a=jacobian(function(x) as.numeric(r(x, cjac=F, ...)$res), par)
            colnames(a)=nm_par
         }
      }
      if (econstr) {
         jac=a%*%nte
      } else {
         jac=a
      }
      covpar=try(solve(crossprod(jac)), silent=T)
      if (inherits(covpar, "try-error")) {
         sv=svd(jac)
         covpar=tcrossprod(sv$v*rep(1./pmax(sv$d, 1./con$rcond), rep(nrow(sv$v), ncol(sv$v))))
      }
      if (econstr) {
         covpar=nte%*%tcrossprod(covpar, nte)
      }
      dimnames(covpar)=list(nm_par, nm_par)
      d=sqrt(diag(covpar))
      hci=(-qt((1.-con$ci$p)/2., fdeg)*sd_res)*d
      names(hci)=nm_par
   } else {
      sd_res=hci=covpar=ci_fdeg=NULL
   }
   return(list(
      par=c(par),
      lastp=c(p),
      hci=hci,
      ci_p=con$ci$p,
      ci_fdeg=fdeg,
      sd_res=sd_res,
      covpar=covpar,
      laststep=c(laststep),
      normp=normp,
      res=c(res),
      prevres=c(vres),
      jacobian=a,
      retres=lres,
      it=it,
      btit=btit,
      history=hist,
      error=0,
      mes=mes)
   )
}
lsi=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   # solve linear least square problem (min ||Ax-b||)
   # with inequality constraints ux>=co
   # Method:
   # 1. reduce the pb to ldp (min(xat*xa) => least distance programming)
   # 2. solve ldp
   # 3. change back to x
   # mnrom, and x0 are dummy parameters which are here to make lsi()
   # compatible with lsi_ln() argument list
   
   if (! is.qr(a)) {
      n=ncol(a)
      aq=qr(a, LAPACK=T)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=abs(diag(aq$qr))
#cat("d=", d, "\n")
   aq$rank=sum(d>d[1]/rcond)
   if (is.na(aq$rank)) {
      x=rep(NA, n)
      attr(x, "mes")="lsi: Rank could not be estimated in least squares\n"
      return(x)
   }
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
      ut=t(backsolve(aq$qr, t(as.matrix(u[, aq$pivot, drop=FALSE])), trans=T))
      xa=ldp(ut, co-u%*%x0, rcond)
      if (is.null(xa)) {
         # Infeasible constraints detected by ldp
         xa=rep(NA, n)
         attr(xa, "mes")="lsi: ldp() revealed unfeasible constrants"
         return(xa)
      }
      x=xa
      xa=backsolve(aq$qr, xa)
      x[aq$pivot]=xa
      x=x+x0
      cou=co-u%*%x
      if (any(cou > 1.e-10)) {
         # round errors made fail inequality constraints
         # force them as much as possible even by degrading the residual
         dx=ldp(u, cou, rcond)
         if (!is.null(dx)) {
            x=x+dx
         } # else leave x as is
      }
   } else {
      # plain QR
      x=x0
   }
   return(x)
}
lsi_ln=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   # solve linear least square problem (min ||A*x-b||)
   # with inequality constraints ux>=co
   # If A is rank deficient, least norm solution ||mnorm*(x-x0)|| is used.
   # If the matrix mnorm is NULL, it is supposed to be an identity matrix.
   # If the vector x0 is NULL, it is treated as 0 vector.
   
   mes=""
   if (! is.qr(a)) {
      a=as.matrix(a)
      n=ncol(a)
      aq=qr(a, LAPACK=T)
   } else {
      n=ncol(a$qr)
      aq=a
   }
   d=diag(aq$qr)
   aq$rank=sum(abs(d)>abs(d[1])/rcond)
   if (aq$rank == 0) {
      # matrix a is zero => no x can diminish the residual norm
      if (!is.null(u) && nrow(u) > 0) {
         # there are some inequalities to deal with
         if (any(co > 1.e-10)) {
            # at least we satisfy inequalities
            x=ldp(u, co)
            if (is.null(x)) {
               x=rep(NA, n)
               attr(x, "mes")="lsi_ln: matrix a is zero rank and ldp() revealed unfeasible constrants"
            }
            return(x)
         } else {
            # they are already satisfied => return 0 vector
            return(double(n))
         }
      } else {
         # no inequalities
         return(double(n))
      }
   }
   # here after the rank is > 0
   rdefic=aq$rank < n
   if (!rdefic) {
      # just passthrough params to lsi()
      return(lsi(aq, b, u, co, rcond))
   }
   # prepare free variable substitution
   ndefic=n-aq$rank
   i1=seq(len=aq$rank)
   ndefic=n-aq$rank
   i2=seq(aq$rank+1, n, len=ndefic)
   r=qr.R(aq)[i1,,drop=FALSE]
   r1=r[,i1,drop=FALSE]
   r2=r[,i2,drop=FALSE]
   qrr=qr(t(r), LAPACK=T)
   # null space and its complement bases
   ba=qr.Q(qrr, complete=T)
   bal=ba[,i1,drop=FALSE]
   ban=ba[,i2,drop=FALSE]
   # least norm without constraints
   x=bal%*%backsolve(qrr$qr, qr.qty(aq, b)[qrr$pivot], trans=T)
   #x0=backsolve(r1, bt) # implicitly, free variables are set to zero
   #x=c(x0, rep(0., ndefic)) # we unpivot just before return
   #browser()
   if (!is.null(u) && nrow(u)>0) {
      up=u[,aq$pivot,drop=FALSE]
      cou=co-up%*%x
      if (any(cou > 1.e-10)) {
         uz=up%*%ban
         # minimize ||y;z||
         ut=t(backsolve(qrr$qr, t(as.matrix(up%*%bal)), trans=F))
         ut=cBind(ut, uz)
         ya=ldp(ut, cou, rcond)
         if (!is.null(ya)) {
            # move along the borders to diminish ||a*x-b||
            # the solution is searched for in the form
            # y=ya+bau*yz
            # ia are indexes of border equations to add as equality constraints (active set)
            cou=as.numeric(cou-ut%*%ya)
            ia=which(cou > -1.e-10) # must not be empty if we are here
            if (length(ia) == 0L) {
               # but it can happen for badly conditioned ut
               # let take the inequality closest to 0
               ia=which.max(cou)
            }
            bau=Nulla(t(ut[ia,,drop=FALSE]))
            yz=lsi_ln(bau[i1,,drop=FALSE], -ya[i1], ut[-ia,,drop=FALSE]%*%bau, cou[-ia], rcond=rcond)
            if (!any(is.na(yz))) {
               # got a solution
               ya=ya+bau%*%yz
               dx=bal%*%backsolve(qrr$qr, ya[i1], trans=T)+ban%*%ya[i2]
               x=x+dx
               attr(x, "mes")=attr(yz, "mes")
            } else {
               x=rep(NA, n)
               attr(x, "mes")=join("\n", attr(x, "mes"), "lsi_ln: failed to glide along borders.")
               return(x)
            }
         } else {
            # unfeasible constraints
            x=rep(NA, n)
            attr(x, "mes")="lsi_ln: unfeasible constraints"
            return(x)
         }
         # solve lsi_ln within just Null space (without modifying ||y||)
         cou=co-up%*%x
         if(is.null(mnorm)) {
            mban=ban
            if (is.null(x0)) {
               x0=-x
            } else {
               x0=x0-x
            }
         } else {
            mban=mnorm%*%ban
            if (is.null(x0)) {
               x0=-mnorm%*%x
            } else {
               x0=mnorm%*%(x0-x)
            }
         }
         z=lsi_ln(mban, x0, uz, cou)
         if (!any(is.na(z))) {
            x=x+ban%*%z
            attr(x, "mes")=join("\n", attr(x, "mes"), attr(z, "mes"))
         } else {
            #print(list(a=a, b=b, u=u, co=co))
#browser()
            attr(x, "mes")=join("\n", attr(x, "mes"), "lsi_ln: Unfeasible constraints in null space. (It must not be. Round off errors struck again.)")
         }
      }
   }
   x[aq$pivot]=x
   if (!is.null(u)) {
      cou=co-u%*%x
      if (any(cou > 1.e-10)) {
         # round errors made fail inequality constraints
         # force them as much as possible even by degrading the residual
         dx=ldp(u, cou, rcond)
         if (!is.null(dx)) {
            x=x+dx
         } # else leave x as is
      }
   }
   attr(x, "mes")=join("\n", attr(x, "mes"), paste(
      "lsi_ln: Rank deficient matrix in least squares\n",
      ndefic, " free variable(s):\n",
      paste(dimnames(a)[[2]][aq$pivot[i2]],
      aq$pivot[-i1], sep="\t", collapse="\n"),
      "\nLeast L2-norm solution is provided.",
      sep=""))
   return(x)
}
ldp=function(u, co, rcond=1.e10) {
   # solve least distance programing: find x satisfying u*x>=co and s.t. min(||x||)
   # by passing to nnls
   # (non negative least square) problem
   # return x or Null (in case of non-feasible inequalities)
   m=NROW(u)
   n=ncol(u)
   if (m == 0) {
      # no inequality to satisfy => trivial case
      return(double(n))
   }
   rcond=abs(rcond)
   # eliminate 0 rows from u
   maxu=apply(u, 1, function(v) max(abs(v)))
   i0=which(maxu < 1.e-13)
   if (length(i0) > 0) {
      if (all(co[i0] < 1.e-10)) {
         u=u[-i0,,drop=FALSE]
         co=co[-i0]
         m=m-length(i0)
         if (m == 0) {
            # no inequality to satisfy => trivial case
            return(double(n))
         }
      } else {
         # non feasible 0 constraints
         return(NULL)
      }
   }
   
   maxco=max(co)
   if (maxco<=1.e-10) {
      # all rhs are < 0 => trivial case
      return(double(n))
   }
   e=rBind(t(u), t(co))
   f=double(n+1L)
   f[n+1]=1.
   resnnls=nnls::nnls(e, f)
   feasible=sqrt(resnnls$deviance) > 1.e-10 && resnnls$residuals[n+1] != 0.
   if (feasible) {
      x=resnnls$residuals[1:n]/(-resnnls$residuals[n+1])
      # check for numerical stability problems
      ux=u%*%x
      cou=co-ux
      if (F && any(cou > 1.e-10)) { # F because it can worsen the situation
         # second trial
         e[n+1,]=cou
         rn=nnls::nnls(e, f)
         if (rn$residuals[n+1]!=0.) {
            x=x-rn$residuals[1:n]/rn$residuals[n+1]
         } else {
            x=NULL
         }
      } else if (all(cou<0) && !all(x==0.)) {
         # round off error pushed the solution inside of feasible domain
         # shorten it till the closest border
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
ls_ln_old=function(a, b, rcond=1.e10) {
   # least squares with least norm
   # no LAPACK in the second qr()
   if (! is.qr(a)) {
      a=qr(a, LAPACK=T)
   }
   n=ncol(a$qr)
   d=diag(a$qr)
   a$rank=sum(abs(d)>abs(d[1])/rcond)
   rdefic=a$rank < n
   i1=1:a$rank
   i2=(a$rank+1):n
   r=qr.R(a)[i1,,drop=FALSE]
   if (!rdefic) {
      # plain ls
      x=backsolve(r, qr.qty(a, b))
      x[a$pivot]=x
      return(x)
   }
   # prepare free variable substitution
   r1=r[,i1,drop=FALSE]
   r2=r[,i2,drop=FALSE]
   qrr=qr(t(r), LAPACK=F)
   # null space and its complement bases
   ba=qr.Q(qrr, complete=T)
   bal=ba[,i1,drop=FALSE]
   # least norm
   x=bal%*%backsolve(qr.R(qrr), qr.qty(a, b), trans=T)
   x[a$pivot]=x
   return(x)
}
ls_ln=function(a, b, rcond=1.e10) {
   # least squares with least norm
   # LAPACK is called in the second qr() too.
   if (! is.qr(a)) {
      a=qr(a, LAPACK=T)
   }
   n=ncol(a$qr)
   d=diag(a$qr)
   a$rank=sum(abs(d)>abs(d[1])/rcond)
   if (a$rank == 0) {
      return(double(n)) # return plain 0 vector
   }
   rdefic=a$rank < n
   i1=seq(len=a$rank)
   r=qr.R(a)[i1,,drop=FALSE]
   if (!rdefic) {
      # plain ls
      x=backsolve(r, qr.qty(a, b))
      x[a$pivot]=x
      return(x)
   }
   ndefic=n-a$rank
   i2=seq(a$rank+1, n, len=ndefic)
   # prepare free variable substitution
   qrr=qr(t(r), LAPACK=T)
   # least norm
   x=qr.qy(qrr, c(backsolve(qrr$qr, qr.qty(a, b)[qrr$pivot], trans=T), rep.int(0., ndefic)))
   x[a$pivot]=x
   return(x)
}
lsie_ln=function(a, b, u=NULL, co=NULL, e=NULL, ce=NULL, rcond=1.e10) {
   # solve linear least square problem (min ||A*x-b||)
   # with inequality constraints u*x>=co and equality constraints e*x=ce
   # Method:
   # reduce the pb to lsi_ln on the null-space of e
   if (is.qr(a)) {
      n=ncol(a$qr)
      aqr=T
   } else {
      n=ncol(a)
      aqr=F
   }
   if (!is.null(e) && nrow(e) > 0) {
      # x is searched in a form xp+bn*z
      # where xp is a particular solution of e*x=ce
      # bn is a basis of e null space and z is a
      # new unknown vector.
      nr=nrow(e)
      nc=ncol(e)
      if (nr > nc) {
         x=rep(NA, n)
         attr(x, "mes")="lsie_ln: matrix e is over determined."
         return(x)
      }
      bn=Nulla(t(e))
      qe=attr(bn, "qr")
      if (qe$rank < nr) {
         x=rep(NA, n)
         attr(x, "mes")="lsie_ln: matrix e is rank deficient."
         return(x)
      }
      # particular solution
      xp=qr.qy(qe, c(backsolve(qe$qr, ce[qe$pivot], trans=T), double(n-nr)))
      if (aqr) {
         R=qr.R(a, complete=T)
         b=b-qr.qy(a, R%*%xp[a$pivot])
         a=qr.qy(a, R%*%bn[a$pivot,,drop=FALSE])
      } else {
         b=b-a%*%xp
         a=a%*%bn
      }
      if (!is.null(u)) {
         co=co-u%*%xp
         u=u%*%bn
      }
      z=lsi_ln(a, b, u, co, rcond)
      x=xp+bn%*%z
      if (!is.null(attr(z, "mes"))) {
         attr(x, "mes")=attr(z, "mes")
      }
      return(x)
   }
}
ls_ln_svd=function(a, b, rcond=1.e10) {
   # least squares of least norm by svd(a)
   sa=svd(a)
   d=sa$d
   rank=sum(d > d[1L]/rcond)
   if (rank == 0L) {
      return(rep(0., ncol(a)))
   }
   i=seq_len(rank)
   x=sa$v[,i]%*%(crossprod(sa$u[,i], b)/d[i])
   return(x)
}
tls=function(a, b) {
   # total least square by svd
   sab=svd(cbind(a, b))
   n=ncol(a)
   v=sab$v[,length(sab$d)]
   return(v[seq_len(n)]/(-v[n+1L]))
}
lsi_reg=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   # solve linear least square problem (min_x ||a*x-b||)
   # with inequality constraints ux>=co
   # If a is rank deficient, regularization term lambda^2*||mnorm*(x-x0)||^2
   # is added to ||a*x-b||^2. The rank is estimated as number of singular values
   # above d[1]*1.e-10 where d[1] is the highest singular value.
   # The scalar lambda is an optional positive parameter
   # and is treated as relative to the highest singular value of a.
   # Its value is calculated from rcond parameter (which is preserved
   # for compatibility with others lsi_...() functions
   # At return, lambda can be found in attributes of the returned vector x.
   # NB. lambda is set to NA
   # if rank(a)==0 or a is of full rank or if there is no inequality.
   # If the matrix mnorm is NULL, it is supposed to be an identity matrix.
   # If the vector x0 is NULL, it is treated as 0 vector.
   # Return: the solution vector x with attribute "lambda"
   
   mes=""
   lambda=1./sqrt(rcond)
   if (is.qr(a)) {
      # get back the matrix from qr
      a=qr.X(a)
   }
   n=ncol(a)
   m=nrow(a)
   if (m < n) {
      # make nrow(a) >= ncol(a) by repeating a and b as many times as needed
      nrep=ceiling(n/m)
      a=aperm(array(a, dim=c(m, n, nrep)), c(1,3,2))
      dim(a)=c(m*nrep, n)
      b=rep(b, nrep)
   }
   sva=svd(a)
   d=sva$d
   arank=sum(d > d[1]*1.e-10)
   x0=if (is.null(x0)) double(n) else x0
   if (arank == 0) {
      # matrix a is zero => no x can diminish the residual norm
      if (!is.null(u) && nrow(u) > 0) {
         # there are some inequalities to deal with
         co0=co-u%*%x0
         if (any(co0 > 1.e-10)) {
            # at least we satisfy inequalities
            if (is.null(mnorm)) {
               dx=ldp(u, co0)
               if (is.null(dx)) {
                  x=rep(NA, n)
                  attr(x, "mes")="lsi_reg: matrix a is zero rank and ldp() revealed unfeasible constrants"
               }
            } else {
               x=lsi_reg(mnorm, mnorm%*%x0, u, co, rcond, NULL, NULL)
            }
         } else {
            # inequalities are already satisfied
            x=x0
         }
      } else {
         # no inequalities
         x=x0
      }
      attr(x, "lambda")=NA
      return(x)
   }
   # hereafter, the arank is > 0
   i=seq_len(arank)
   if (NROW(u) == 0) {
      # no inequality
      # plain ls_ln
      if (arank == n) {
         x=sva$v[,i]%*%(crossprod(sva$u[,i,drop=FALSE], b)/d[i])
         attr(x, "lambda")=NA
      } else {
         mes=sprintf("lsi_reg: Rank deficient matrix in least squares
There is (are) %d free variable(s).
Regularized L2-norm solution is provided.", n-arank)
         if (is.null(mnorm)) {
            mnorm=diag(1., n)
         }
         lambda=d[1]*lambda
         lam2=lambda**2
         d2=diag(d**2, n)
         mv=mnorm%*%sva$v
         mv2=crossprod(mv)
         bx0=crossprod(mv, mnorm%*%x0)
         bt=d*crossprod(sva$u, b)+lam2*x0
         x=sva$v%*%solve(d2+lam2*mv2, d*bt)
         attr(x, "lambda")=lambda
         attr(x, "mes")=mes
      }
   } else {
      # there are inequalities
      if (arank == n) {
         # a is of full rank => no lambda to use, plain lsi
         invd=1./d
         bt=crossprod(sva$u, b)
         ut=(u%*%sva$v)*rep(invd, rep(n, n))
         y=ldp(ut, co-ut%*%bt)
         x=sva$v%*%(invd*(y+bt))
         attr(x, "lambda")=NA
      } else {
         mes=sprintf("lsi_reg: Rank deficient matrix in least squares
There is (are) %d free variable(s).
Regularized L2-norm solution is provided.", n-arank)
         if (is.null(mnorm)) {
            mnorm=diag(1., n)
         }
         # transformed rhs
         mv=mnorm%*%sva$v
         bx0=crossprod(mv, mnorm%*%x0)
         bt=d*crossprod(sva$u, b)
         # semi-transformed u (it lacks a factor 1/(d**2+(lambda*mv)**2) )
         uv=u%*%sva$v
         mv2=crossprod(mv)
         # ldp to find approximate solution
         lambda=d[1]*lambda
         lam2=lambda**2
         d2=diag(d**2, n)+lam2*mv2
         ut=t(solve(t(d2), t(uv)))
         btt=bt+lam2*bx0
         y=ldp(ut, co-ut%*%btt)
         x=sva$v%*%solve(d2, btt+y)
         attr(x, "lambda")=lambda
         attr(x, "mes")=mes
      }
   }
   return(x)
}
uplo2uco=function(param, upper=NULL, lower=NULL, linear=NULL) {
   # Return a list with a matrix u and a vector co such that u%*%param-co>=0
   # translates the inequalities param <= upper and param >= lower
   # Names in upper and lower must be present in names of param
   nm_par=names(param)
   nlo=length(lower)
   nup=length(upper)
   nli=length(linear)
   u=matrix(0., nrow=nup+nlo+nli, ncol=length(param))
   co=numeric(nrow(u))
   colnames(u)=nm_par
   rownames(u)=c(
      if (nup) paste(names(upper), " <= ", upper, sep="") else NULL,
      if (nlo) paste(names(lower), " >= ", lower, sep="") else NULL,
      linear
   )
   names(co)=rownames(u)
   # fill u and co
   if (nup) {
      u[seq_len(nup), names(upper)]=diag(-1, nup)
      co[seq_len(nup)]=-upper
   }
   u[nup+seq_len(nlo), names(lower)]=diag(1, nlo)
   co[nup+seq_len(nlo)]=lower

   # parse inequalities
   if (nli) {
      vema=try(equa2vecmat(nm_par, linear, sep=">="))
      if (inherits(vema, "try-error")) {
         stop("Error in parsing inequalities")
      }
      i=nup+nlo+seq_len(nli)
      u[i, ]=vema[,-1L]
      co[i]=vema[,1L]
   }
   return(list(u=u, co=co))
}
equa2vecmat=function(nm_par, linear, sep="=") {
   # parse a text vector of linear equations and produce a corresponding
   # matrix and right hand side vector
   # Input:
   #  - nm_par a text vector of variable names. It will be used in the symbolic
   # derivation.
   #  - linear is a text vector of linear equations like "a+2*c+3*b = 0"
   #  - sep is the separator of two parts of equations. Use for example
   #    ">=" for linear inequalities
   # Return: an augmented matrix. Its first column is the rhs vector.
   # Other columns are named by nm_par.
   # If the vector linear is NULL or its content is empty a NULL is returned
   
   # Sanity check
   stopifnot(length(sep)==1 && nchar(sep) > 0)
   stopifnot(length(nm_par)>=1 && all(nchar(nm_par)) > 0)
   
   vlin=sapply(linear, function(it) strsplit(it, sep)[[1]])
   if (length(vlin) && !is.matrix(vlin)) {
      stop(sprintf("Linear (in)equalities are expected to have '%s' sign in them", sep))
   }
   if (length(vlin) == 0) {
      # we are set
      return(NULL)
   }
   # each column in vlin is an (in)equality
   # first row is left part, the second row is the right part
   # derive left and right parts to get the matrix
   ze=double(length(nm_par))
   names(ze)=nm_par
   ze=as.list(ze)
   de=apply(vlin, 2L, function(ineq) {
      le=with(ze, eval(deriv(parse(text=ineq[1]), nm_par)))
      ri=with(ze, eval(deriv(parse(text=ineq[2]), nm_par)))
      v=c(ri-le, attr(le, "gradient")-attr(ri, "gradient"))
      return(v)
   })
   rownames(de)=c("rhs", nm_par)
   return(t(de))
}
lsi_lim=function(a, b, u=NULL, co=NULL, rcond=1.e10, mnorm=NULL, x0=NULL) {
   suppressWarnings(limSolve::lsei(A=a, B=b, G=u, H=co)$X)
}
