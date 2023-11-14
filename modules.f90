!================================================================
!  Modules and subroutines to be used with sfec13d.f90
!  modules.f90
!================================================================
 module constants  
! 
! Some relevant constant 
!#########################
   use geometry
   implicit none

   real(kind=8), parameter :: yr=5.9355d7,dvfac=0.3d0
   real(kind=8), parameter :: pi=3.14159265358979d0,pi4=4d0*pi,deg2rad=pi/180d0

  contains

!=============================================================================
 subroutine readinput(string,value,string_value)
! Search the input file and assigns a value to a selected string variable
!=============================================================================
   implicit none
   integer :: k
   real (kind=8) :: v,value
   character(len=*) :: string,string_value
   character(len=100) :: s,sv
!=============================================================================

   k=0
   open(1,file="sfec.inp")
   
   do while(.true.)

     read(1,*,end=99) s
     if(s(1:1).ne.'#') then
       backspace(1)
       read(1,*,end=99) s,v,sv
     endif
     
     if (s.eq.string.and.sv.eq.'no') then
       value=v
       k=1
     elseif (s.eq.string.and.int(v).eq.999) then
       string_value=sv
       k=1
     endif
     
   end do
99 continue

   if(k.ne.1) then
     print *, 'Variable ',string,' missing in sfec.inp'
     stop
   endif
   
   close(1)
 
 end subroutine readinput
!===================================================================


 end module constants   
!#####################


!################################################################
 module grid
! 
! Module used to build the discrete radial and angular structure 
!################################################################
   implicit none
   integer :: nt,nlayer,ndima
   real(kind=8) :: rcr
   real(kind=8), dimension(:), allocatable   :: rad,eta,roots,wghts

 contains

!==============================================================
	SUBROUTINE gauleg(x1,x2,x,w)
!	
! The subroutine generates the n-dimensional vector x of roots
! of Legendre polynomials and w of the corresponding weights
!
!  see NRF90 Ch. B4, pag. 1059
!==============================================================	
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: xl,xm
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauleg')
	m=(n+1)/2
	xm=0.5_dp*(x2+x1)
	xl=0.5_dp*(x2-x1)
	z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			end where
		end do
		where (unfinished)
			pp=n*(z*p1-p2)/(z*z-1.0_dp)
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	x(1:m)=xm-xl*z
	x(n:n-m+1:-1)=xm+xl*z
	w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	w(n:n-m+1:-1)=w(1:m)
	END SUBROUTINE gauleg
!======================================================
 end module grid 
!#################


!###########################################################
 module legendre
! 
! This module is used to perform the angular integration 
!###########################################################
   implicit none
   real(kind=8), dimension(:,:), allocatable :: plg
 contains   
!==========================================================================
 subroutine dpnm(x,n,pol)
! 
! This subroutine returns the array of the fully normalized
! associated legendre functions
!
! see: Z.Martinec, "Program to calculate the spectral harmonic expansion 
!                   coefficients of the two scalar fields product"
!      Comp.Phys.Comm. (1989), 54, 177-182
! N.B. Originally in F77, contains some obsolescent instruction
!==========================================================================
   implicit none   
   integer :: n,iflag,ms,np1,m,i,j,i1,i2,mp2
   real(kind=8), dimension (:)  ::  pol
   real(kind=8) :: sth,x,pi4,sou,zn,f1,f2,f3,sthm
!      
      sth=dsqrt(1d00-x*x)
      pi4=16d00*datan(1d00)
      pol(1)=1d00/dsqrt(pi4)
      if(n.le.0) return
      sou=1d00
      zn=-1d00
      np1=n+1
      iflag=0
      do 3 ms=1,np1
        m=ms-1
        f1=dfloat(m+m)
        i=(m+1)*(m+2)/2
        if(m.eq.0) goto 1
        sou=sou*(f1+1d00)/f1
        if(iflag-1.lt.0) goto 10
   10   sthm=sth**m
        if(sthm.lt.1.d-55) iflag=1
        goto 25
   20   sthm=0d00
   25   continue
        pol(i)=zn*dsqrt(sou/pi4)*sthm
        if(m.eq.n) goto 3
        zn=-zn
    1   i1=i+m+1
        pol(i1)=dsqrt(f1+3d00)*x*pol(i)
        mp2=m+2
        if(mp2.gt.n) goto 3
        do 2 j=mp2,n
          f1=dfloat((j-m)*(j+m))
          f2=dfloat((j+j-1)*(j+j+1))
          f3=dfloat((j+j+1)*(j-m-1)*(j+m-1))/dfloat(j+j-3)
          f2=dsqrt(f2/f1)*x
          f3=dsqrt(f3/f1)
          i=j*(j+1)/2+m+1
          i1=i-j
          i2=i1-j+1
          pol(i)=f2*pol(i1)-f3*pol(i2)
    2   continue
    3 continue
 end subroutine dpnm	  
!====================================================================
 subroutine cyjm(j,m,leg,ph,yjm)      
!
! Complex spherical harmonic Y_{jm}(th,phi) for a given j and m
! (-j.le.m.le.j)
!====================================================================
   implicit none
   integer :: j,m,jm
   real(kind=8) :: zn,ph,leg
   complex(kind=8) :: cunit,csupl,yjm
!
   cunit=(0d0,1d0)
   jm=(j*(j+1))/2+iabs(m)+1
   zn=1d0
   if(m.lt.0) zn=(-1d0)**m
   csupl=cunit*dfloat(m)*ph
   yjm=zn*leg*cdexp(csupl)
 end subroutine cyjm
!===================================================================

 end module legendre
!###################### 


!############################################################## 
 module matrix
! 
! This module is used to setup the matrix of the system 
!##############################################################
   implicit none
   integer :: iccs
   integer, dimension(:), allocatable   :: row_ind,col_ptr
   real(kind=8), dimension(:), allocatable   :: val,q1,q2,dq1,dq2
   real(kind=8), dimension(:,:), allocatable :: u,v,pr,du,dv,dpr
   complex (kind=8), dimension(:), allocatable :: rhs
   complex (kind=8), dimension(:,:), allocatable  :: crho
 contains   
!==============================================================
 subroutine msfe(jmin,jmax)
! 
! Matrix of the System of Finite Elements
!==============================================================
   use grid
   
   implicit none    
   integer :: j,m,k,jm,jmin,jmax,jmmax,ip
   
   jmmax=jmax*(jmax+1)/2+jmax+1
   allocate(q1(jmmax),q2(jmmax),u(jmmax,nlayer+1),v(jmmax,nlayer+1),pr(jmmax,nlayer+1))
   allocate(dq1(jmmax),dq2(jmmax),du(jmmax,nlayer+1),dv(jmmax,nlayer+1),dpr(jmmax,nlayer+1))
   
   q1=0d0; q2=0d0; u=0d0; v=0d0; pr=0d0
   
   iccs=0
   ip=0
   do j=jmin,jmax
     write(*,*) 'Current j in msfe:',j
     do m=0,j
       jm=(j*(j+1))/2+m+1
       do k=1,nlayer+1
! block relative to the Lagrange multiplier "q1" (cmb vertical stress)
         if(k.ne.1) goto 2
         ip=ip+1
         q1(jm)=1d0
         call vef(k,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
         q1(jm)=0d0
         call ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)	 
! block relative to vertical velocity "u"
 2       ip=ip+1
         u(jm,k)=1d0
         call vef(k,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
         u(jm,k)=0d0
         call ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)	 
! block relative to horizontal velocity "v"  
         ip=ip+1
         v(jm,k)=1d0
         call vef(k,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
         v(jm,k)=0d0
         call ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)	 
! block relative to the Lagrange multiplier "q2" (top vertical stress)
         if(k.ne.nlayer+1) goto 3
         ip=ip+1
         q2(jm)=1d0
         call vef(k,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
         q2(jm)=0d0
         call ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)	 
! block relative to the Lagrange multiplier "p" (pressure)      
 3       if(k.eq.nlayer+1) exit
         ip=ip+1
         pr(jm,k)=1d0
         call vef(k,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
         pr(jm,k)=0d0
         call ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)	 
       end do
     end do
   end do
   col_ptr(ip+1)=iccs+1

end subroutine msfe

!========================================================================
 subroutine vef(kref,jm,j,q1,q2,pr,u,v,dq1,dq2,dpr,du,dv)
! 
!  Variation of the Energy Functional (1-D case) with respect to
!    - L. multiplier representing cmb radial stress (dq1) 
!    - L. multiplier representing top radial stress (dq2)
!    - L. multiplier representing pressure          (dpr)
!    - vertical velocity                            (du)
!    - horizontal velocity                          (dv)
!========================================================================
   use constants
   use grid
   use legendre
   
   implicit none
   integer :: k,kref,j,jm
   real (kind=8) :: jj,p1n,p2n,p3n
   real (kind=8) :: a,b,I1,I6,I3a,I3b,K2a,K2b,pom1,pom2
   real (kind=8), dimension(:)   :: q1,q2,dq1,dq2
   real (kind=8), dimension(:,:) :: u,v,pr,du,dv,dpr
      
   dq1=0d0; dq2=0d0; du=0d0; dv=0d0; dpr=0d0
   
   jj=dble(j*(j+1))
!----------------------------------------
!  Variation of the functional
!----------------------------------------
   du(jm,1)=q1(jm)*radc*radc     
   dv(jm,1)=0d0
   do k=1,nlayer ! k-loop over finite elements
     if(k.eq.kref-1.or.k.eq.kref.or.k.eq.kref+1) then      	      

! Results of the ortogonality relations
       p1n=eta(k)                
       p2n=jj*eta(k)
       p3n=(j-1)*jj*(j+2)*eta(k)
	 
! Integrals over finite elements (notation according to Martinec(2001), GJI)
       b=rad(k+1)
       a=rad(k)                          
       I1=(b*b+a*b+a*a)/(b-a)/3d0     !  I1(k,k)=-I1(k,k+1)=I1(k+1,k+1)=-K1(k)=K1(k+1)
       I6=(b-a)/3d0                   !  I6(k,k)=2*I6(k,k+1)=I6(k+1,k+1)
       I3a=-(b+a+a)/6d0               !  I3(k,k)=-I3(k,k+1)
       I3b=(b+b+a)/6d0                !  I3(k+1,k+1)=-I3(k+1,k)
       K2a=(b-a)*(b+a+a)/6d0          !  K2(k)
       K2b=(b-a)*(b+b+a)/6d0          !  K2(k+1)
	 
!-------------------------------------------------------
!  Variation with respect to the two Lagrange multipliers
!-------------------------------------------------------
       dq1(jm)=radc*radc*u(jm,1)  
       dq2(jm)=erad*erad*u(jm,nlayer+1)
	 
!--------------------------------------------------
!  Variation with respect to pressure
!--------------------------------------------------
       dpr(jm,k)=-a*a*u(jm,k)+b*b*u(jm,k+1)-jj*(K2a*v(jm,k)+K2b*v(jm,k+1))
	 
!--------------------------------------------------
!  Variation with respect to vertical velocity
!--------------------------------------------------
       du(jm,k)=du(jm,k)-pr(jm,k)*rad(k)*rad(k)
       du(jm,k+1)=pr(jm,k)*rad(k+1)*rad(k+1)
       pom1= u(jm,k)*I1-u(jm,k+1)*I1
       pom2=-u(jm,k)*I1+u(jm,k+1)*I1
       du(jm,k)=du(jm,k)+2d0*p1n*pom1
       du(jm,k+1)=du(jm,k+1)+2d0*p1n*pom2
       pom1=u(jm,k)*I6+u(jm,k+1)*I6/2d0
       pom2=u(jm,k)*I6/2d0+u(jm,k+1)*I6
       du(jm,k)=du(jm,k)+4d0*p1n*pom1
       du(jm,k+1)=du(jm,k+1)+4d0*p1n*pom2
       pom1=v(jm,k)*I3a-v(jm,k+1)*I3a
       pom2=-v(jm,k)*I3b+v(jm,k+1)*I3b
       du(jm,k)=du(jm,k)+p2n*pom1
       du(jm,k+1)=du(jm,k+1)+p2n*pom2
       pom1=(-v(jm,k)-0.5d0*v(jm,k+1)+u(jm,k)+0.5d0*u(jm,k+1))*I6
       pom2=(-0.5d0*v(jm,k)-v(jm,k+1)+0.5d0*u(jm,k)+u(jm,k+1))*I6
       du(jm,k)=du(jm,k)+p2n*pom1
       du(jm,k+1)=du(jm,k+1)+p2n*pom2
       pom1=(-v(jm,k)-0.5d0*v(jm,k+1))*I6
       pom2=(-0.5d0*v(jm,k)-v(jm,k+1))*I6
       du(jm,k)=du(jm,k)+2d0*p1n*jj*pom1
       du(jm,k+1)=du(jm,k+1)+2d0*p1n*jj*pom2
       du(jm,nlayer+1)=du(jm,nlayer+1)+q2(jm)*erad*erad
	 
!-------------------------------------------------------
!  Variation with respect to horizontal velocity
!-------------------------------------------------------
       dv(jm,k)=dv(jm,k)-pr(jm,k)*K2a*jj
       dv(jm,k+1)=-pr(jm,k)*K2b*jj
       pom1=v(jm,k)*I1-v(jm,k+1)*I1-v(jm,k)*I3a+v(jm,k+1)*I3b+u(jm,k)*I3a-u(jm,k+1)*I3b
       pom2=-v(jm,k)*I1+v(jm,k+1)*I1+v(jm,k)*I3a-v(jm,k+1)*I3b-u(jm,k)*I3a+u(jm,k+1)*I3b
       dv(jm,k)=dv(jm,k)+p2n*pom1
       dv(jm,k+1)=dv(jm,k+1)+p2n*pom2
       pom1=v(jm,k)*I3a-v(jm,k+1)*I3a-v(jm,k)*I6-v(jm,k+1)*I6/2d0+u(jm,k)*I6+u(jm,k+1)*I6/2d0
       pom2=-v(jm,k)*I3b+v(jm,k+1)*I3b-v(jm,k)*I6/2d0-v(jm,k+1)*I6+u(jm,k)*I6/2d0+u(jm,k+1)*I6
       dv(jm,k)=dv(jm,k)-p2n*pom1
       dv(jm,k+1)=dv(jm,k+1)-p2n*pom2
       pom1=(v(jm,k)+v(jm,k+1)/2d0)*jj-(u(jm,k)+u(jm,k+1)/2d0)*2d0
       pom2=(v(jm,k)/2d0+v(jm,k+1))*jj-(u(jm,k)/2d0+u(jm,k+1))*2d0
       dv(jm,k)=dv(jm,k)+p1n*jj*I6*pom1
       dv(jm,k+1)=dv(jm,k+1)+p1n*jj*I6*pom2
       pom1=v(jm,k)+v(jm,k+1)/2d0
       pom2=v(jm,k)/2d0+v(jm,k+1)
       dv(jm,k)=dv(jm,k)+p3n*I6*pom1
       dv(jm,k+1)=dv(jm,k+1)+p3n*I6*pom2

     endif      
   end do 

end subroutine vef                  


!==================================================================
subroutine ccstore(ip,iccs,jmin,jmax,dq1,dq2,dpr,du,dv)
!
! The subroutine stores the matrix compactly into the arrays     
! val,row_ind, col_ptr 
!
! Only non-zero elements are stored according to
! CCS (Compressed Column Storage) method
!
! see Barret et al. "Templates for the solution of linear systems:
!                     building blocks for iterative methods"
!     http://www.netlib.org/templates/Templates.html   
!==================================================================
   use grid
   
   implicit none
   integer :: j,m,jm,k,icol,iq,ip,jmin,jmax,iccs
   real(kind=8), dimension(:)   :: dq1,dq2
   real(kind=8), dimension(:,:) :: du,dv,dpr
      
   icol=0   
   iq=-4
   do j=jmin,jmax      
     do m=0,j
       jm=j*(j+1)/2+m+1
       do k=1,nlayer+1
         iq=iq+4	  
         if(k.eq.1) then     
           if(dq1(jm).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dq1(jm)
             row_ind(iccs)=iq+1
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(du(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=du(jm,k)
             row_ind(iccs)=iq+2
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dv(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dv(jm,k)
             row_ind(iccs)=iq+3
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dpr(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dpr(jm,k)
             row_ind(iccs)=iq+4
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
         elseif(k.eq.nlayer+1) then
           if(du(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=du(jm,k)
             row_ind(iccs)=iq+1
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dv(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dv(jm,k)
             row_ind(iccs)=iq+2
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dq2(jm).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dq2(jm)
             row_ind(iccs)=iq+3
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
         else
           if(du(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=du(jm,k)
             row_ind(iccs)=iq+1
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dv(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dv(jm,k)
             row_ind(iccs)=iq+2
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
           if(dpr(jm,k).ne.0d0) then
             icol=icol+1
             iccs=iccs+1
             val(iccs)=dpr(jm,k)
             row_ind(iccs)=iq+3
             if(icol.eq.1) col_ptr(ip)=iccs
           endif
         endif
         if(k.ne.1) iq=iq-1
       end do
     end do  
   end do       

end subroutine ccstore      


!=============================================================
 subroutine rhstm(nltm,indld,jmin,jmax,jlmin,jlmax,rhs)
! 
!  Right Hand Side Term of the linear system
!  3D load of the tomographic model
!=============================================================
   use constants, only: gref
   use grid
   
   implicit none
   integer :: k,nl,j,m,jm,jmin,jmax,jlmin,jlmax,kld,nltm,i1,i2
   integer, dimension(:) :: indld
   real (kind=8) :: a,b,qk3,qk3p
   complex (kind=8), dimension(:) :: rhs
   
   nl=0
   do j=jmin,jlmax
     do m=0,j
       nl=nl+1
       jm=(j*(j+1))/2+m+1
       do k=1,nltm
         kld=indld(k)
         a=rad(kld)
         b=rad(kld+1)                               ! Integrals:
         qk3=(b-a)*(b*b+2.d0*a*b+3.d0*a*a)/12.d0    !  k3(k)
         qk3p=(b-a)*(3.d0*b*b+2.d0*a*b+a*a)/12.d0   !  K3(k+1)        
         i1=(nl-1)*(3*nlayer+4)+(3*kld-1)
         i2=(nl-1)*(3*nlayer+4)+(3*kld+2)
         rhs(i1)=qk3*gref*crho(k,jm)
         rhs(i2)=qk3p*gref*crho(k,jm)
       end do	 
     end do  
   end do  

 end subroutine rhstm
 

!=============================================================

 end module matrix
!##########################


!#################################################################
 module matrixmakeup
! 
! The module is used to perform some operation on the matrix, on 
! the preconditioner and to solve the linear system 
!#################################################################
  implicit none
  integer :: m1,m2
  integer, dimension (:), allocatable :: indx,indxtr
  real (kind=8), dimension (:), allocatable   :: x
  real (kind=8), dimension (:,:), allocatable :: acomp,al,acomptr,altr

contains 

!=============================================================
	SUBROUTINE bandec(a,m1,m2,al,indx,d)
!	
! LU decomposition of a band diagonal matrix
!
! see NRF90 Ch. B2, pag. 1020
!=============================================================
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,swap,arth
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), INTENT(IN) :: m1,m2
	REAL(DP), DIMENSION(:,:), INTENT(OUT) :: al
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(DP), INTENT(OUT) :: d
	REAL(DP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	REAL(DP) :: dum
	n=assert_eq(size(a,1),size(al,1),size(indx),'bandec: n')
	mm=assert_eq(size(a,2),m1+m2+1,'bandec: mm')
	mdum=assert_eq(size(al,2),m1,'bandec: mdum')
	a(1:m1,:)=eoshift(a(1:m1,:),dim=2,shift=arth(m1,-1,m1))
	d=1.0
	do k=1,n
		l=min(m1+k,n)
		i=imaxloc(abs(a(k:l,1)))+k-1
		dum=a(i,1)
		if (dum == 0.0) a(k,1)=TINY
		indx(k)=i
		if (i /= k) then
			d=-d
			call swap(a(k,1:mm),a(i,1:mm))
		end if
		do i=k+1,l
			dum=a(i,1)/a(k,1)
			al(k,i-k)=dum
			a(i,1:mm-1)=a(i,2:mm)-dum*a(k,2:mm)
			a(i,mm)=0.0
		end do
	end do
	END SUBROUTINE bandec

!=============================================================================
	SUBROUTINE linbcg(b,x,itol,tol,itmax,iter,err)
!	
!  Solution of a linear system by Preconditioned Biconiugate Gradient method
!	
! see NRF90 Ch. B2, pag. 1034
!=============================================================================
        USE matrix
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: b
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
	INTEGER(I4B), INTENT(IN) :: itol,itmax
	REAL(DP), INTENT(IN) :: tol
	INTEGER(I4B), INTENT(OUT) :: iter
	REAL(DP), INTENT(OUT) :: err
	REAL(DP), PARAMETER :: EPS=1.0e-14_dp
	INTEGER(I4B) :: n
	REAL(DP) :: ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm
	REAL(DP), DIMENSION(size(b)) :: p1,pp,r,rr,z,zz
	n=assert_eq(size(b),size(x),'linbcg')
	iter=0
	call atimes(x,r,0)
	r=b-r
	rr=r
!
!	call atimes(r,rr,0)
	select case(itol)
		case(1)
			bnrm=snrm(b,itol)
			call asolve(r,z,0)
		case(2)
			call asolve(b,z,0)
			bnrm=snrm(z,itol)
			call asolve(r,z,0)
		case(3:4)
			call asolve(b,z,0)
			bnrm=snrm(z,itol)
			call asolve(r,z,0)
			znrm=snrm(z,itol)
		case default
			call nrerror('illegal itol in linbcg')
	end select
	do
		if (iter > itmax) exit
		iter=iter+1
		call asolve(rr,zz,1)
		bknum=dot_product(z,rr)
		if (iter == 1) then
			p1=z
			pp=zz
		else
			bk=bknum/bkden
			p1=bk*p1+z
			pp=bk*pp+zz
		end if
		bkden=bknum
		call atimes(p1,z,0)
		akden=dot_product(z,pp)
		ak=bknum/akden
		call atimes(pp,zz,1)
		x=x+ak*p1
		r=r-ak*z
		rr=rr-ak*zz
		call asolve(r,z,0)
		select case(itol)
			case(1)
				err=snrm(r,itol)/bnrm
			case(2)
				err=snrm(z,itol)/bnrm
			case(3:4)
				zm1nrm=znrm
				znrm=snrm(z,itol)
				if (abs(zm1nrm-znrm) > EPS*znrm) then
					dxnrm=abs(ak)*snrm(p1,itol)
					err=znrm/abs(zm1nrm-znrm)*dxnrm
				else
					err=znrm/bnrm
					cycle
				end if
				xnrm=snrm(x,itol)
				if (err <= 0.5_dp*xnrm) then
					err=err/xnrm
				else
					err=znrm/bnrm
					cycle
				end if
		end select
!		write (16,*) ' iter=',iter,' err=',err
		if (err <= tol) exit
	end do
	write (*,*) ' iter=',iter,' err=',err
	END SUBROUTINE linbcg

!================================================	
	FUNCTION snrm(sx,itol)
!	
! Compute one of two norms for a vector sx(1:n), 
! as signaled by itol. Used by linbcg. 
!
! see NRF90 Ch.B2, pag. 1036
!================================================	
	USE nrtype
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: sx
	INTEGER(I4B), INTENT(IN) :: itol
	REAL(DP) :: snrm
	if (itol <= 3) then
		snrm=sqrt(dot_product(sx,sx))
	else
		snrm=maxval(abs(sx))
	end if
	END FUNCTION snrm

!===================================================
 subroutine atimes(x,y,idor)
! 
! Mutiplication y = A * x    if(idor.eq.0)
! Mutiplication y = A^t * x  if(idor.eq.1)
! A is stored in "val,icol_ind,irow_ptr" CCS-way
!===================================================
   use matrix 
   use grid, only: ndima
!   
   implicit none
   integer :: idor,i,j
   real (kind=8), dimension(:) :: x,y
!      
   if(idor.eq.0) then
     y=0d0
     do j=1,ndima
       do i=col_ptr(j),col_ptr(j+1)-1	 
         y(row_ind(i))=y(row_ind(i))+val(i)*x(j)
       end do
     end do
   else
     do i=1,ndima
       y(i)=0d0
       do j=col_ptr(i),col_ptr(i+1)-1
         y(i)=y(i)+val(j)*x(row_ind(j))
       end do	
     end do
   endif
 end subroutine atimes

!=============================================================
 subroutine asolve(b,x,idor)
! 
! Solution of the system M * x = b    if(idor.eq.0)
! Solution of the system M^t * x = b  if(idor.eq.1), 
! where M is the preconditioner, stored in 'acomp' and 'al' and
! M^t its transposed, stored in 'acomptr' and 'altr'
!=============================================================
   use grid
   implicit none
   integer :: idor
   real (kind=8), dimension(:) :: b,x
   real (kind=8), dimension(size(b)) :: baux
!   
   baux=b
   if(idor.eq.0) then 
     call banbks(acomp,m1,m2,al,indx,baux)
   elseif(idor.eq.1) then
     call banbks(acomptr,m1,m2,altr,indxtr,baux)
   endif
   x=baux
 end subroutine asolve

!=======================================================================
	SUBROUTINE banbks(a,m1,m2,al,indx,b)
!	
! Solution of a linear system with band diagonal LU-decomposed matrix
!
! see NRF90 Ch. B2, pag. 1021
!=======================================================================
	USE nrtype; USE nrutil, ONLY : assert_eq,swap
	IMPLICIT NONE
	REAL(DP), DIMENSION(:,:), INTENT(IN) :: a,al
	INTEGER(I4B), INTENT(IN) :: m1,m2
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,k,l,mdum,mm,n
	n=assert_eq(size(a,1),size(al,1),size(b),size(indx),'banbks: n')
	mm=assert_eq(size(a,2),m1+m2+1,'banbks: mm')
	mdum=assert_eq(size(al,2),m1,'banbks: mdum')
	do k=1,n
		l=min(n,m1+k)
		i=indx(k)
		if (i /= k) call swap(b(i),b(k))
		b(k+1:l)=b(k+1:l)-al(k,1:l-k)*b(k)
	end do
	do i=n,1,-1
		l=min(mm,n-i+1)
		b(i)=(b(i)-dot_product(a(i,2:l),b(1+i:i+l-1)))/a(i,1)
	end do
	END SUBROUTINE banbks
!=============================================================

!=========================================================
SUBROUTINE gaussj(a,b)
! Linear system solution via Gauss-Jordan elimination
!=========================================================
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
  INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
  LOGICAL(LGT), DIMENSION(size(a,1)) :: lpiv
  REAL(DP) :: pivinv
  REAL(DP), DIMENSION(size(a,1)) :: dumc
  INTEGER(I4B), TARGET :: irc(2)
  INTEGER(I4B) :: i,l,n
  INTEGER(I4B), POINTER :: irow,icol
  n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
  irow => irc(1)
  icol => irc(2)
  ipiv=0
  do i=1,n
        lpiv = (ipiv == 0)
        irc=maxloc(abs(a),outerand(lpiv,lpiv))
        ipiv(icol)=ipiv(icol)+1
        if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
        if (irow /= icol) then
                call swap(a(irow,:),a(icol,:))
                call swap(b(irow,:),b(icol,:))
        end if
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol) == 0.0) &
                call nrerror('gaussj: singular matrix (2)')
        pivinv=1.0_sp/a(icol,icol)
        a(icol,icol)=1.0
        a(icol,:)=a(icol,:)*pivinv
        b(icol,:)=b(icol,:)*pivinv
        dumc=a(:,icol)
        a(:,icol)=0.0
        a(icol,icol)=pivinv
        a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
        b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
        a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
        b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
  end do
  do l=n,1,-1
        call swap(a(:,indxr(l)),a(:,indxc(l)))
  end do

END SUBROUTINE gaussj


 end module matrixmakeup

!############################## 

