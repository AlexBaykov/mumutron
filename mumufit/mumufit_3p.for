      program mumufit
      implicit none
c --- this file is needed here to get access to some FIT variables
      include 'fitdata1.inc'
c approximation of e+e- -> KSKL pi0 -> 7 gamma with fitting method: 1
c 14 parameters:
c      MsOm
c      MsA0
c      CrsRho1
c      MsRho1
c      GRho1
c      PhRho1
c      CrsRho2
c      MsRho2
c      GRho2
c      PhRho2
c      CrsRho3
c      MsRho3
c      GRho3
c      PhRho3
c initial values are defined in "test.dat"
c
c === fitbegin data ===
c
c --- a number of components
      integer inproc
      parameter ( inproc = 1 )
c --- auxiliary parameters: number of points per component
      integer in1
      parameter ( in1 = 3 )
c --- a period of radcor calculations, usual values:
c     < -10000 derivatives at the first pass only, no radcor
c     -1       derivatives at every step, no radcor
c     0        no radcor, no derivatives
c     1        radcor and derivatives at every step
c     > 10000  radcor and derivatives at the first pass
      integer initer
      parameter ( initer = 0 )
c --- radcor cut-off energy, E_max, should be treated carefully
      real*8 demaxrg(inproc) / 500.0 /
c --- par-component settings: (fitmode)+10*(superidx)+1000*(radtype)
c fitmode, now available:
c    0 - Process is not considered
c    1 - Asymmetric Chi2 (derrnevdn/derrnevup should be supplied)
c    2 - Pseudo-Poisson  (disabled)
c    3 - Poisson         (for small stat)
c    4 - Symmetric Chi2  (for high stat)
c    5 - Auto-selection
c superidx - can be 0..99
c radtype, is one of the following:
c    0 - usual operation, according initer
c    1 - no radcor, only derivatives
c    2 - no radcor, no derivatives
      integer ifitmode(inproc)
      data ifitmode / 3 /
c --- number of points for each component
      integer inener(inproc)
      data inener / in1 /
c --- the length of data array
      integer iarrlen
      parameter ( iarrlen = in1 )
c --- energy in c.m.s., it will be filled from debeam array
      real*8 decm(iarrlen)
      data decm/
     &  213.300,213.306,213.316/
c --- error of decm, here it is considered to be 50 keV
      real*8 derre2(iarrlen)
      data derre2 / iarrlen*0.0 /
c --- spread of decm
      real*8 dsige2(iarrlen)
      data dsige2 / iarrlen*0.0 /
c --- integrated luminosity and its error
      real*8 dlum(iarrlen)
      data dlum /
c     &  288.,288.,288.,288.,288./
     &  288.,288.,288./
      real*8 derrlum(iarrlen)
      data derrlum / iarrlen*0.0 /
c --- live time of measurement
      real*8 dlivtim(iarrlen)
      data dlivtim / iarrlen*1000000.0 /
c --- efficiency and its error (obtained from MC simulation)
      real*8 deff(iarrlen) /iarrlen*1.0 /
      real*8 derreff(iarrlen) /iarrlen*0.0 /
c --- cross section (nb)
      real*8 dnev(iarrlen)
      data dnev /
     &  99704.73, 99851.73, 100096.18/
c --- arrays for user supplied errors of number of events
c     they are actually not used in this program
      real*8 derrnevdn(iarrlen), derrnevup(iarrlen)
      data derrnevdn /
c     &  0.,0.,0.,0.,0./
     &  0.,0.,0./
      data derrnevup /
c     &  0.,0.,0.,0.,0./
     &  0.,0.,0./   
c --- filename of minuit input command file
      character*(*) cminfile
      parameter ( cminfile = 'mumufit_3p' )
c --- filename of output files
      character*(*) chisfile
      parameter ( chisfile = 'mumufit_3p' )
c === end of fitbegin data ===
c
c --- local variables
      integer i
c
c --- call to inifit, initialization of auxiliary library RADCOR
c      call inifit()
c       call rkkstar
      CALL csmumu
      thanxtoandy = .true.
c --- the main call
      call fitbegin(
     &  inproc,
     &  initer,
     &  ifitmode,
     &  inener,
     &  iarrlen,
     &  dnev,
     &  decm,
     &  derre2,
     &  dsige2,
     &  dlum,
     &  derrlum,
     &  dlivtim,
     &  deff,
     &  derreff,
     &  derrnevup,
     &  derrnevdn,
     &  demaxrg,
     &  cminfile,
     &  chisfile)
      END
c---------------------------------------------------------------------
c a subroutine to calculate a theoretical cross section.
c should be written by user.
c input parameters:
c   iproc - component index,
c   dener - energy in c.m.s.
c   parmin() - values of parameters, are transferred through common block.
c output parameters:
c   d_crspr - calculated cross section
c---------------------------------------------------------------------
      subroutine crs_proc(iproc,dener,d_crspr)
      implicit none
      include 'fitdata1.inc'
c
      integer iproc
      real*8 dener,d_crspr
      double precision RELTOL,ABSTOL,ERR
      integer nseg
      double precision de,sigb,e,xmin,xmax,csrx,csfun,cs
      common /cscom/e,sigb
      external csfun

      de = parmin( 1)
      sigb  = parmin( 2)
      cs =parmin( 3)
      e=dener+de
      xmin=e-4.
      xmax=e+4
      RELTOL=1.D-10
      ABSTOL=0.d0
      NSEG=1
      CALL dADAPT(csfun,xmin,xmax,NSEG,RELTOL,ABSTOL,csrx,ERR)
      d_crspr=csrx*cs
      end

      subroutine csmumu
      implicit none
      double precision x,f,e,sigb
      double precision rmass,pi,alpha,mevnb
      parameter (alpha=1/137.036)
      parameter (rmass=105.658)
      parameter (mevnb=0.389379338d12)
      parameter (pi=3.14159265d0)
      double precision  xmin,xmax,csrad,cs,csr(10001)
      double precision wcs,sw,bw,cs0,en(10001),csrx
      common /cscom/e,sigb
      common /funmumu/en,csr
      external wcs
      double precision RELTOL,ABSTOL,ERR
      integer i,nseg

      do i=1,10001
        e=2.*rmass+0.001d0*(i-1)+1.d-7
        en(i)=e
        xmin=1.0d-14
        xmax=1.-4.d0*rmass**2/e**2
        RELTOL=1.D-10
        ABSTOL=0.d0
        NSEG=1
        CALL dADAPT(wcs,xmin,xmax,NSEG,RELTOL,ABSTOL,csrx,ERR)
        csr(i)=csrx+sw(e,xmin)*bw(e)
c        if(i.le.100)print *,e-2.*rmass,csr(i)
      enddo
      end

      double precision function csfun(x)
      implicit none
      double precision x,en(10001),csr(10001),f,sigb,pi,e
c      parameter (sigb=0.400)
      parameter (pi=3.14159265d0)
      common /cscom/e,sigb
      common /funmumu/en,csr
      integer i

      if(x.lt.en(1))then
        f=0.
      else if(x.ge.en(10001))then
        f=csr(10001)
      else
        i=min(int((x-en(1))/0.001d0)+1,10000)
        f=csr(i)+((x-en(1))/0.001d0+1-i)*(csr(i+1)-csr(i))
      endif
      csfun=f/sigb/sqrt(2.d0*pi)*exp(-(x-e)**2/2.d0/sigb**2)
      end

      double precision function wcs(x)
      implicit none
      double precision x,e,bw
      common /cscom/e
      double precision alpha,pi,a2,a3
      parameter (alpha=1.d0/137.036d0)
      parameter (pi=3.14159265d0)
      parameter (a2=1.64493407d0)
      parameter (a3=1.2020569d0)
      double precision s,L,beta,delta2,delta,w

      s=e*e
      L=log(s/(0.511d0)**2)
      beta=2.*alpha/pi*(L-1.d0)
      delta2=(1.125d0-2.d0*a2)*L**2-
     &  (2.8125d0-5.5d0*a2-3.d0*a3)*L-
     &  1.2d0*a2**2-4.5d0*a3-6.d0*a2*log(2.d0)+
     &  0.375d0*a2+4.75d0
      delta=1.d0+alpha/pi*(1.5d0*L+pi**2/3.d0-2.d0)+
     &  (alpha/pi)**2*delta2
      w=delta*beta*x**(beta-1.d0)-beta/2.d0*(2.d0-x)+
     &  beta**2/8.d0*((2.d0-x)*(3.d0*log(1.d0-x)-4.d0*log(x))-
     &  4.d0*log(1.d0-x)/x-6.d0+x)
      wcs=w*bw(e*sqrt(1.d0-x))
      end

      double precision function sw(e,x)
      implicit none
      double precision e,s,x
      double precision alpha,pi,a2,a3
      parameter (alpha=1.d0/137.036d0)
      parameter (pi=3.14159265d0)
      parameter (a2=1.64493407d0)
      parameter (a3=1.2020569d0)
      double precision L,beta,delta2,delta,w

      s=e*e
      L=log(s/(0.511d0)**2)
      beta=2.d0*alpha/pi*(L-1.d0)
      delta2=(1.125d0-2.d0*a2)*L**2-
     &  (2.8125d0-5.5d0*a2-3.d0*a3)*L-
     &  1.2d0*a2**2-4.5d0*a3-6.d0*a2*log(2.d0)+
     &  0.375d0*a2+4.75d0
      delta=1.d0+alpha/pi*(1.5d0*L+pi**2/3.d0-2.d0)+
     &  (alpha/pi)**2*delta2
      sw=delta*x**beta
     &  -beta*x
     &  -beta**2*(x*log(x)-0.75d0*x)

      end

      double precision function bw(e)
      implicit none
      double precision rmass,pi,alpha,mevnb
      parameter (pi=3.14159265d0)
      parameter (alpha=1/137.036)
      parameter (rmass=105.658)
      parameter (mevnb=0.389379338d12)
      double precision e,beta,cs0,cs,y,c

      if(e.lt.2.*rmass)then
        bw=0
      else if(e.eq.2.*rmass)then
        bw=2.d0*pi**2*alpha**3/e**2*mevnb
      else
        beta=sqrt(1.d0-4.d0*rmass**2/e**2)
        y=pi*alpha/beta
        c=y/(1.d0-exp(-y))
        cs0=2.d0*pi*alpha**2*beta/e**2*(1.d0-beta**2/3)*mevnb
        bw=cs0*C
      endif
      end

      double precision function cs0(e)
      implicit none
      double precision rmass,pi,alpha,mevnb
      parameter (pi=3.14159265d0)
      parameter (alpha=1/137.036)
      parameter (rmass=105.658)
      parameter (mevnb=0.389379338d12)
      double precision e,beta,cs0,cs,y,c

      if(e.le.2.*rmass)then
        cs0=0
      else
        beta=sqrt(1.d0-4.d0*rmass**2/e**2)
        cs0=2.d0*pi*alpha**2*beta/e**2*(1.d0-beta**2/3)*mevnb
      endif
      end

      subroutine crs_bgrd(iproc,dener,d_crsbgl,d_crsbgt)
      implicit none
      integer iproc
      real*8 dener,d_crsbgl,d_crsbgt
c
c      if ( .not.(iproc.eq.1) ) then
c        write(6,*) 'crs_bgrd unknown process id'
c        stop
c      end if
      d_crsbgl=0.0d0
      d_crsbgt=0.0d0
      END
c      include 'mumu_fit.f'

