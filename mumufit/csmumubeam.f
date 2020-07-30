      program csmum
      implicit none
      integer i,nseg,j,npoi
      double precision rmass
      parameter (rmass=105.658)
      parameter(npoi=10001)
      double precision  e,xmin,xmax,st, cs
      double precision csfun, x, f(npoi),dsimps
      external csfun
      common /cscom/e
      
      call csmumu
      
c      do i=1,4001
      do i=1,3
c        e=2.*rmass+0.001d0*(i-1)+1.d-6-2.

        if(i.eq.1)then
          e=209.756
        else if(i.eq.2  )then
          e=213.300
        else
          e=211.316 + 2 
        endif  
        xmin=e-2.
        xmax=e+2
        st=(xmax-xmin)/(npoi-1)
        do j=1,npoi
          x=xmin+st*(j-1)
          f(j)=csfun(x)
        enddo
c        print *,f
        cs=DSIMPS(F,xmin,xmax,npoi-1)
c       print *,e-2.*rmass,cs
        print *,e,cs,cs*0.08*3600.
      enddo
      end


      subroutine csmumu
      implicit none
      double precision x,f,e
      double precision rmass,pi,alpha,mevnb
      parameter (alpha=1/137.036)
      parameter (rmass=105.658)
      parameter (mevnb=0.389379338d12)
      parameter (pi=3.14159265d0)
      double precision  xmin,xmax,csrad,cs,csr(10001)
      double precision wcs,sw,bw,cs0,en(10001),csrx
      common /cscom/e
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
      parameter (sigb=0.400)
      parameter (pi=3.14159265d0)
      common /cscom/e
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


      
      
