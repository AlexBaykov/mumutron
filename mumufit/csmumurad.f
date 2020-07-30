      program csmumu
      implicit none
      double precision rmass,pi,alpha,mevnb
      parameter (alpha=1/137.035999074)
      parameter (rmass=105.6583745)
      parameter (mevnb=0.389379338d12)
      parameter (pi=3.14159265d0)
      double precision  e,xmin,xmax,csrad,cs,csr
      double precision wcs,sw,bw,cs0
      common /cscom/e
      external wcs
      double precision RELTOL,ABSTOL,ERR
      integer i,nseg

      do i=1,2001
        e=2.*rmass+0.001d0*(i-1)+1.d-7
        xmin=1.0d-14
        xmax=1.-4.d0*rmass**2/e**2
        RELTOL=1.D-10
        ABSTOL=0.d0
        NSEG=1
        CALL dADAPT(wcs,xmin,xmax,NSEG,RELTOL,ABSTOL,csr,ERR)
        csr=csr+sw(e,xmin)*bw(e)
        cs=bw(e)
        print *,e-2.*rmass,cs0(e),cs,csr
      enddo

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


      
      
