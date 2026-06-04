subroutine dran_ini(iseed0)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(np=14)
  parameter(nbit=31)
  parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
  integer ix(ip)
  dimension g(0:m)

  data c0,c1,c2/2.515517,0.802853,0.010328/
  data d1,d2,d3/1.432788,0.189269,0.001308/
  data pi/3.141592653589793d0/

  common /ixx/ ix
  common /icc/ ic     
  common /gg/ g

  dseed=iseed0
  do i=1,ip
     ix(i)=0
     do j=0,nbit-1
        if(rand_xx(dseed).lt.0.5) ix(i)=ibset(ix(i),j)
     enddo
  enddo
  ic=0

  do i=m/2,m
     p=1.0-dble(i+1)/(m+2)
     t=sqrt(-2.0*log(p))
     x=t-(c0+t*(c1+c2*t))/(1.0+t*(d1+t*(d2+t*d3)))
     g(i)=x
     g(m-i)=-x
  enddo

  u2th=1.0-dble(m+2)/m*sqrt(2.0/pi)*g(m)*exp(-g(m)*g(m)/2)
  u2th=nn1*sqrt(u2th)
  g=g/u2th

  return
end subroutine dran_ini

subroutine dran_read(iunit)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(np=14)
  parameter(m=2**np)
  integer ix(ip)
  dimension g(0:m)
  common /ixx/ ix
  common /icc/ ic
  common /gg/ g
  read (iunit,*) ic
  read (iunit,*) (ix(i),i=1,ip)
  read (iunit,*) (g(i),i=0,m)
  return
end subroutine dran_read

subroutine dran_write(iunit)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(np=14)
  parameter(m=2**np)
  integer ix(ip)
  dimension g(0:m)
  common /ixx/ ix
  common /icc/ ic      
  common /gg/ g
  write (iunit,*) ic
  write (iunit,*) (ix(i),i=1,ip)
  write (iunit,*) (g(i),i=0,m)
  return
end subroutine dran_write

function i_dran(n)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(is=ip-iq)
  integer ix(ip)
  common /ixx/ ix
  common /icc/ ic
  ic=ic+1
  if (ic > ip) ic=1
  if (ic > iq) then
     ix(ic)=ieor(ix(ic),ix(ic-iq))
  else
     ix(ic)=ieor(ix(ic),ix(ic+is))
  endif
  i_dran=ix(ic)
  if (n > 0) i_dran=mod(i_dran,n)+1
  return
end function i_dran

function dran_u()
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(is=ip-iq)
  parameter (rmax=2147483648.0d0)
  integer ix(ip)
  common /ixx/ ix
  common /icc/ ic
  ic=ic+1
  if (ic > ip) ic=1
  if (ic > iq) then
     ix(ic)=ieor(ix(ic),ix(ic-iq))
  else
     ix(ic)=ieor(ix(ic),ix(ic+is))
  endif
  dran_u=(dble(ix(ic))+0.5d0)/rmax
  return
end function dran_u


function dran_g()
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(np=14)
  parameter(nbit=31)
  parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
  parameter(is=ip-iq)

  integer ix(ip)
  dimension g(0:m)

  common /ixx/ ix
  common /icc/ ic
  common /gg/ g

  ic=ic+1
  if(ic > ip) ic=1
  if(ic > iq) then
     ix(ic)=ieor(ix(ic),ix(ic-iq))
  else
     ix(ic)=ieor(ix(ic),ix(ic+is))
  endif
  i=ishft(ix(ic),-np1)
  i2=iand(ix(ic),nn)
  dran_g=i2*g(i+1)+(nn1-i2)*g(i)
  return
end function dran_g


function dran_gbmw()
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(is=ip-iq)
  parameter (rmax=2147483648.0d0)
  integer ix(ip)
  integer, save :: icount=1
  double precision, save :: u,v
  common /ixx/ ix
  common /icc/ ic
  data pi2 /6.283185307179586d0/
  if (icount.eq.1) then
     ic=ic+1
     if (ic > ip) ic=1
     if (ic > iq) then
	ix(ic)=ieor(ix(ic),ix(ic-iq))
     else
	ix(ic)=ieor(ix(ic),ix(ic+is))
     endif
     u=pi2*dble(ix(ic)+0.5d0)/rmax
     ic=ic+1
     if(ic > ip) ic=1
     if(ic > iq) then
	ix(ic)=ieor(ix(ic),ix(ic-iq))
     else
	ix(ic)=ieor(ix(ic),ix(ic+is))
     endif
     v=(dble(ix(ic))+0.5d0)/rmax
     v=dsqrt(-2.0d0*log(v))
     dran_gbmw=dcos(u)*v
     icount=2
  else
     dran_gbmw=dsin(u)*v
     icount=1
  endif
  return
end function dran_gbmw


subroutine dran_gv(u,n)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(np=14)
  parameter(nbit=31)
  parameter(m=2**np,np1=nbit-np,nn=2**np1-1,nn1=nn+1)
  parameter(is=ip-iq)
  dimension g(0:m)
  dimension u(n)
  dimension ix(ip)
  common /gg/ g
  common /ixx/ ix
  common /icc/ic

  n1=0
  do while (n1 < n)
     if (ic.lt.iq) then
        kmax=min(n-n1,iq-ic)
        do k=1,kmax
           ic=ic+1
           ix(ic)=ieor(ix(ic),ix(ic+is))
           i=ishft(ix(ic),-np1)
           i2=iand(ix(ic),nn)
           u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)
        enddo
     else
        kmax=min(n-n1,ip-ic)
        do k=1,kmax
           ic=ic+1
           ix(ic)=ieor(ix(ic),ix(ic-iq))
           i=ishft(ix(ic),-np1)
           i2=iand(ix(ic),nn)
           u(n1+k)=i2*g(i+1)+(nn1-i2)*g(i)
        enddo
     endif
     if(ic.ge.ip) ic=0
     n1=n1+kmax
  enddo

  return
end subroutine dran_gv

subroutine dran_uv(u,n)
  implicit double precision(a-h,o-z)
  parameter(ip=1279)
  parameter(iq=418)
  parameter(is=ip-iq)
  parameter (rmax=2147483648.0d0)
  dimension u(n)
  dimension ix(ip)
  common /ixx/ ix
  common /icc/ic

  n1=0
  do while (n1 < n)
     if (ic.lt.iq) then
        kmax=min(n-n1,iq-ic)
        do k=1,kmax
           ic=ic+1
           ix(ic)=ieor(ix(ic),ix(ic+is))
           u(n1+k)=(dble(ix(ic))+0.5d0)/rmax
        enddo
     else
        kmax=min(n-n1,ip-ic)
        do k=1,kmax
           ic=ic+1
           ix(ic)=ieor(ix(ic),ix(ic-iq))
           u(n1+k)=(dble(ix(ic))+0.5d0)/rmax
        enddo
     endif
     if(ic.ge.ip) ic=0
     n1=n1+kmax
  enddo

  return
end subroutine dran_uv



function rand_xx(dseed)
  double precision a,c,xm,rm,dseed,rand_xx
  parameter (xm=2.d0**32,rm=1.d0/xm,a=69069.d0,c=1.d0)
  dseed=mod(dseed*a+c,xm)
  rand_xx=dseed*rm
  return
end function rand_xx


function dran_gamma(alpha)
!!$ Added on May 3rd 2016
!!$ Function to generate random numbers according to the gamma distribution
!!$ f(x)=exp[-x]x^(alpha-1)/gamma[alpha], for alpha>0,
!!$ using the rejection method explained in exercise 4.11 of Stochastic Numerical Methods course.
!!$ There are two different formulas whether alpha<1 or alpha>=1
!!$ For alpha>3.5 it is better to use the Numerical Recipes method using
!!$ a Lorentzian proposal

  implicit double precision (a-h,o-z)
  data alpha0/-1.0d0/
  save alpha0,alpha1,alpha2,x0,ca,cb,cd,am,s
  if (alpha.le.0.0d0) stop 'alpha not valid'
  if (alpha.ne.alpha0) then
!!$ New value of alpha. Need to compute parameters again
     x0=1.0d0-alpha
     if (alpha.lt.1.0d0) then
        ca=x0**(1.0d0-1.0d0/alpha)
        cb=x0*dlog(x0)-x0
        alpha1=1.0d0/alpha
     else if (alpha.lt.3.5d0) then
        cd=x0*dlog(alpha)-x0
        alpha2=1.0d0/alpha-1.0d0
     else
        am=alpha-1.0d0
        s=dsqrt(2.0d0*am+1.0d0)
     endif
     alpha0=alpha
  endif

  if (alpha.lt.1.0) then
2    u=dran_u()
     if (u.lt.x0) then
        dran_gamma=ca*u**alpha1
        if (dran_u().lt.dexp(-dran_gamma)) return
     else
        dran_gamma=x0-dlog(alpha1*(1.0d0-u))
        h=cb-x0*dlog(dran_gamma)
        if (dran_u().lt.dexp(h)) return
     endif
     goto 2
  else if (alpha.lt.3.5d0) then
3    u=dran_u()
     dran_gamma=-alpha*dlog(u)
     h=cd+alpha2*dran_gamma-x0*dlog(dran_gamma)
     if(dran_u().lt.dexp(h)) return
     goto 3
  else
1    v1=dran_u()
     v2=2.0d0*dran_u()-1.0d0
     if (v1*v1+v2*v2.gt.1.0d0) goto 1
     y=v2/v1
     dran_gamma=s*y+am
     if(dran_gamma.le.0.0d0) goto 1
     h=(1.0d0+y*y)*dexp(am*dlog(dran_gamma/am)-s*y)
     if (dran_u().lt.h) return
     goto 1
  endif
end function dran_gamma

!!$ Added on October 26th 2021
!!$ Function to generate random numbers according to the Poisson distribution
!!$ p_i=e^(-lambda)lambda^i/i!
!!$ using the method explained in my notes Poisson.pdf located at Stoch/SimMeth/LectureNotes/PoissonDistribution
!!$ There are two different formulas whether:
!!$ lambda<50 (from the exponential distribution of times)
!!$ lambda>=50 (rejection using a trial exponential distribution) 

function iran_poisson(lambda)
  implicit double precision (a-h,o-z)
  double precision lambda
  if (lambda.ge.50.0d0) then
1    u=dran_u()
     if (u.lt.0.5d0) then
        x=lambda+dsqrt(1.5d0*lambda)*dlog(2.0d0*u)
     else
        x=lambda-dsqrt(1.5d0*lambda)*dlog(2.0d0*(1.0d0-u))
     endif
     i=x
     if (i.lt.0) goto 1
     h=1.5186836405255704d0*&
          dexp(-lambda+dsqrt(2./(3*lambda))*abs(x-lambda)+(i+0.5d0)*dlog(lambda)-gammln(dble(i+1)))
     if (dran_u().gt.h) goto 1
     iran_poisson=i
  else
     F=dran_u()
     a=dexp(-lambda)
     iran_poisson=0
     do  while(F.gt.a)
        F=F*dran_u()
        iran_poisson=iran_poisson+1
     enddo
  endif
end function iran_poisson


!!$ Generation of binomial distribution
!!$ Copied from BNLDEV of Numerical Recipes and "translated" to modern fortran of double precision
!!$ Modified to use dran routines
!!$ Added March 20th 2023
function iran_bin(pp,n)
  implicit double precision (a-h,o-z)
  parameter (pi=3.1415926535897932385d0)
  SAVE nold,pold,pc,plog,pclog,en,oldg
  DATA nold /-1/, pold /-1.0d0/
  if (pp.le.0.5d0)then
     p=pp
  else
     p=1.0d0-pp
  endif
  am=n*p
  if (n.lt.25)then
     iran_bin=0
     do  j=1,n
        if(dran_u().lt.p) iran_bin=iran_bin+1
     enddo
  else if (am.lt.1.0d0) then
     g=dexp(-am)
     t=1.0d0
     do j=0,n
        t=t*dran_u()
        if (t.lt.g) go to 1
     enddo
     j=n
1    iran_bin=j
  else
     if (n.ne.nold) then
        en=n
        oldg=gammln(en+1.0d0)
        nold=n
     endif
     if (p.ne.pold) then
        pc=1.0d0-p
        plog=dlog(p)
        pclog=dlog(pc)
        pold=p
     endif
     sq=dsqrt(2.0d0*am*pc)
2    y=dtan(pi*dran_u())
     em=sq*y+am
     if (em.lt.0.0d0.or.em.ge.en+1.0d0) go to 2
     em=int(em)
     t=1.20d0*sq*(1.0d0+y**2)*dexp(oldg-gammln(em+1.0d0)-gammln(en-em+1.0d0)+em*plog+(en-em)*pclog)
     if (dran_u().gt.t) go to 2
     iran_bin=int(em)
  endif
  if (p.ne.pp) iran_bin=n-iran_bin
  return
end function iran_bin

function gammln(xx)  
  double precision gammln,xx  
  integer j  
  double precision ser,stp,tmp,x,y,cof(6)  
  SAVE cof,stp  
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0, &
       -1.231739572450155d0,.1208650973866179d-2, -.5395239384953d-5,2.5066282746310005d0/  
  x=xx  
  y=x  
  tmp=x+5.5d0  
  tmp=(x+0.5d0)*log(tmp)-tmp  
  ser=1.000000000190015d0  
  do j=1,6  
     y=y+1.d0  
     ser=ser+cof(j)/y  
  enddo
  gammln=tmp+log(stp*ser/x)  
end function gammln
