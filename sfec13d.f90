!=====================================================================================================
program sfec13d
!  sfec13d.f90 
!  Spectral Finite Element Convection with 1D viscosity and 3D response
!  Solves the Stokes equation using density anomalies derived from a tomographic
!  model and computes the resulting geoid and dynamic topographies of the
!  surface and CMB
!=====================================================================================================
   use constants
   use grid
   use legendre
   use matrix
   use matrixmakeup
   
   implicit none   
   integer :: i, k, l, j, m, jm, k1, k2, kc
   integer :: jmin, jmax, jlmin, jlmax, jmmax, jmlmax
   integer :: nltm, nd, nlprov, nz, itol, itmax, iter, kph, sg
   integer, dimension(:), allocatable :: indld
   real (kind=8) :: readval
   real (kind=8) :: rerho, imrho, dz, dr, suld, xth, tol, d, dtr, err, fatt, th, ph, leg, zn
   real (kind=8), dimension(:), allocatable :: z, pnm, Rex, Imx
   real (kind=8), dimension(:,:), allocatable :: acomppom, pot, potld, geoid, topotop, topocore, mtrx
   real (kind=8), dimension(4) :: Re_sgX, Im_sgX
   real (kind=8), dimension(4,2) :: sgB
   real (kind=8), dimension(4,4) :: sgM
   complex (kind=8) :: yjm
   complex (kind=8), dimension(:), allocatable :: tauc, taut, topocore_jm, topotop_jm, potsld_jm, potcld_jm, pot_jm
   complex (kind=8), dimension(:,:), allocatable :: vflow, hflow, pres
   character (len=30) :: tomography,readsval
!=====================================================================================================
   
!-------------------------------------------------------      
! Radial and spectral structure of model and load
!-------------------------------------------------------
! Minimum degree for:      
   call readinput("Min_degree",readval,readsval)
   jmin=int(readval) ! response
   jlmin=jmin        ! load
   write(*,*) "Min. harmonic degree:", jmin

!---------------------------------------------------------    
! Read Re & Im coeff. of the tomographic model     
!---------------------------------------------------------    
   call readinput("Tomography_file",readval,readsval)
   tomography=readsval
   write(*,*) "Tomographic model:  ",tomography
   open(1,file=tomography)
   read(1,*) nltm,jlmax        ! number of layers of the tomographic model and its max. degree
   jmlmax=jlmax*(jlmax+1)/2+jlmax+1

   allocate(z(nltm+1))          ! radii of the tomographic model
   allocate(crho(nltm,jmlmax))  ! complex density
   do k=1,nltm
     read(1,*) z(k)
     z(k)=z(k)*1d3
     do j=0,jlmax
       do m=0,j
         jm=j*(j+1)/2+m+1
         read(1,*) rerho,imrho           ! Re and Im part of tomographic densities
         crho(k,jm)=dcmplx(rerho,imrho)   
       end do 
     end do
   end do
   close(1)

! Maximum degree for       
   call readinput("Max_degree",readval,readsval)
   jlmax=int(readval) ! response
   jmax=jlmax         ! load
   write(*,*) "Max. harmonic degree:", jmax

! The complex density is multiplied by the thickness of the tomographic layer
   dz=z(2)-z(1)   ! thickness of one tomographic layer (constant)
   crho=dz*crho  
   
!----------------------------------------------------------
! Setup of the radial grid of finite elements
!----------------------------------------------------------
   write(*,*) 'Number of layers of the tomographic model=',nltm
   call readinput("N_node_layer",readval,readsval)
   nd=int(readval)

   nlprov=4*nltm*nd  ! temporary estimate of the size of rad
   allocate(rad(0:nlprov),indld(nltm+1)) 
   rad(0)=radc
   rad(1)=radc
   suld=1d0  ! support for the piecewise const. func. of the load (must be 1!)
   kc=1
   do k1=1,nltm+1
     if(k1.eq.1) then
       dr=(z(1)-radc)/nd
     elseif(k1.eq.nltm+1) then
       dr=(erad-z(nltm))/nd
     else  
       dr=(z(k1)-z(k1-1))/nd
     endif
     do k=1,nd
       kc=kc+1
       if(k.eq.1) then
         rad(kc)=rad(kc-2)+dr
       else  
         rad(kc)=rad(kc-1)+dr
       endif  
       if(k.eq.nd) rad(kc)=z(k1)
       if(k.eq.nd.and.k1.eq.nltm+1) rad(kc)=erad
       if(k.eq.nd) indld(k1)=kc  ! index associated to the load position
     end do
     if(k1.eq.nltm+1) exit  
     kc=kc+1
     rad(kc)=z(k1)+suld
   end do
   nlayer=kc-1   ! final number of layer

!  CONTROLLARE SE NECESSARIO   
!  per far coincidere rm e rlith con un nodo della griglia   
   do k=1,nlayer
     if(rad(k).le.rm.and.rm.le.rad(k+1)) rad(k)=rm
     if(rad(k).le.rlith.and.rlith.le.rad(k+1)) rad(k)=rlith
   end do
   write(*,*) 'Total number of layers of the model:',nlayer
   
!----------------------------------------------      
! Setup of viscosity on the 1-D radial grid
! (3-layers are used here)
!----------------------------------------------
   allocate(eta(nlayer+1))
   do k=1,nlayer+1

     ! Lower mantle 
     if(rad(k).ge.radc.and.rad(k).lt.rm) then
       call readinput("Visc_LM",readval,readsval)
       eta(k)=readval
     ! Upper mantle mantle 
     elseif(rad(k).ge.rm.and.rad(k).lt.rlith) then
       call readinput("Visc_UM",readval,readsval)
       eta(k)=readval
     else
     ! Lithosphere 
       call readinput("Visc_Lith",readval,readsval)
       eta(k)=readval
     endif

   end do  

!----------------------------------------------------------------------------
! Evaluation of Legendre polynomials on the theta-grid and relative storage
!---------------------------------------------------------------------------- 
   call readinput("Theta_points",readval,readsval)
   nt=int(readval)
   allocate(roots(nt),wghts(nt))
   call gauleg(-1d0,1d0,roots,wghts)      
   roots=-roots
   jmmax=jmax*(jmax+1)/2+jmax+1
   allocate(pnm(jmmax+1),plg(nt,jmmax+1))
   do l=1,nt
     xth=roots(l)
     call dpnm(xth,jmax,pnm)
     do j=0,jmax  
       do m=0,j
         jm=j*(j+1)/2+m+1
         plg(l,jm)=pnm(jm)
       end do
     end do  
   end do

!--------------------------------------
! Setup of the linear system 
!--------------------------------------
! Dimension of the matrix of the system
   ndima=(3*nlayer+4)*((jmax+1)*(jmax+2)/2-1)
   nz=0.01*(ndima*ndima) 
   write(*,*) 'Setup of the matrix of the system'
   write(*,*) ' Matrix dimension:',ndima
   
   allocate(val(nz),row_ind(nz),col_ptr(ndima+1))
   call msfe(jmin,jmax)   

!----------------------------------------------------
! Uncomment the following block if you need to store 
! the full matrix into array mat(:,:)
!----------------------------------------------------
!   allocate(mtrx(ndima,ndima))
!   mtrx=0d0
!   do i=1,ndima
!     do j=col_ptr(i),col_ptr(i+1)-1
!       mtrx(row_ind(j),i)=val(j)
!     end do
!   end do
!   open(25,file='matrix.dat')
!   do i=1,ndima
!     write(25,98) (mtrx(i,j), j=1,ndima)
!   end do
!98 format(10000E14.3)   
!   close(1)
!   stop
!----------------------------------------------------   
   
   write(*,*) 'Total number of matrix elements:', ndima*ndima
   write(*,*) 'Number of non-zero matrix elements:', iccs
   write(*,*) 'Compact storage of preconditioner' 
   
! The principal diagonal of the matrix forms a band that is used to build 
! the preconditioner and is stored compactly into the array acomp
   m1=4         ! number of subdiagonals
   m2=4         ! number of superdiagonals
   allocate(acomp(ndima,m1+m2+1),acomptr(ndima,m1+m2+1),acomppom(m1+m2+1,ndima))
   acomppom=0d0
   acomp=0d0
   
! The three loops that follow recover the elements of the principal
! band of the matrix from its CCS, saving them into acompt
   do i=1,ndima
     k=i-m1-1
     do m=col_ptr(i),col_ptr(i+1)-1
       do j=max(1,1-k),min(m1+m2+1,ndima-k)
           if(k+j.eq.row_ind(m)) then
           acomppom(j,i)=val(m)
         endif
       end do
     end do
   end do

! The following two loops rearrange acompt into acomp that has
! the suitable format to be read from bandec
   do j=1,m1+m2+1
     k=j-m1-1
     do i=max(1,1-k),min(ndima,ndima-k)
       acomp(i+k,m1+m2+2-j)=acomppom(j,i)
     end do
   end do

! The following two loops build the transpose acomptr of acomp in
! the suitable format to be read from bandec
   do i=1,m1+m2+1
     k1=m1+1-i
     k2=ndima+m1-i
     do j=max(1,k1+1),min(ndima,k2+1)
       acomptr(j,i)=acomp(j-m1-1+i,m1+m2+2-i)
     end do
   end do  
   write(*,*) 'LU-decomposition of the preconditioner' 
   
! LU-decomposition of the compactly stored preconditioner acomp and its transpose acomptr
   allocate(indx(ndima),al(ndima,m1))
   call bandec(acomp,m1,m2,al,indx,d)
   allocate(indxtr(ndima),altr(ndima,m1))
   call bandec(acomptr,m1,m2,altr,indxtr,dtr)

! Setup of the complex r.h.s. term    
   allocate(rhs(ndima))
   rhs=0d0
   call rhstm(nltm,indld,jmin,jmax,jlmin,jlmax,rhs)

!--------------------------------
! Solution of the system
!--------------------------------
   call readinput("CG_tol_test",readval,readsval)
   itol=int(readval)  ! tolerance test
   call readinput("CG_tol",readval,readsval)
   tol=readval       ! tolerance 
   call readinput("CG_maxit",readval,readsval)
   itmax=int(readval) ! highest number of allowed iterations

   write(*,*) 'Start of the iterations (solution for real part)'
! Starting guess for the solution Rex (the null vector is used)
   allocate(Rex(ndima))
   Rex=0d0
! Solution for real part
   call linbcg(dble(rhs),Rex,itol,tol,itmax,iter,err) 

   write(*,*) 'Start of the iterations (solution for imaginary part)'
! Starting guess for the solution Imx (the null vector is used)
   allocate(Imx(ndima))
   Imx=0d0
! Solution for imaginary part
   call linbcg(dimag(rhs),Imx,itol,tol,itmax,iter,err) 
   
!------------------------------------------------
! Storage of the solutions of the system
!------------------------------------------------
   allocate(tauc(jmmax),taut(jmmax),topocore_jm(jmmax),topotop_jm(jmmax),     &
            vflow(jmmax,nlayer+1),hflow(jmmax,nlayer+1),pres(jmmax,nlayer+1))
   i=0
   do j=jmin,jmax

     do m=0,j
       jm=(j*(j+1))/2+m+1

       do k=1,nlayer+1

         if(k.eq.1) then
           i=i+1
           tauc(jm)=dcmplx(Rex(i),Imx(i))    ! radial stress at cmb 

         endif

         i=i+1
         vflow(jm,k)=dcmplx(Rex(i),Imx(i))   ! vertical flow 

         i=i+1
         hflow(jm,k)=dcmplx(Rex(i),Imx(i))   ! horizontal flow

         if(k.le.nlayer) then
           i=i+1
           pres(jm,k)=dcmplx(Rex(i),Imx(i))  ! pressure 
         endif

         if(k.eq.nlayer+1) then
           i=i+1
           taut(jm)=dcmplx(Rex(i),Imx(i))    ! radial stress at top
         endif

       end do
     end do
   end do  



!------------------------------
! Geoid computation
!------------------------------
   allocate(potsld_jm(jmmax),potcld_jm(jmmax),pot_jm(jmmax))
   allocate(potld(nt,0:360),pot(nt,0:360),geoid(nt,0:360),topotop(nt,0:360),topocore(nt,0:360))
   potsld_jm=0d0 
   potcld_jm=0d0 
   pot_jm=0d0   
   topotop_jm=0d0
   topocore_jm=0d0
   
! Coefficients of surface potential (potsld_jm) and 
! core potential (potcld_jm) due to the internal density anomalies only
   do j=jlmin,jlmax   
     do m=0,j
       jm=(j*(j+1))/2+m+1

! Radial integration       
       do k=1,nltm
         potsld_jm(jm)=potsld_jm(jm)+z(k)*(z(k)/erad)**(dble(j+1d0))*crho(k,jm)
         potcld_jm(jm)=potcld_jm(jm)+z(k)*(radc/z(k))**dble(j)*crho(k,jm)       

       end do
     end do
   end do

! Self-gravitation
   call readinput("Self_gravitation",readval,readsval)
   sg=int(readval)
   if(sg.eq.1) then
     write(*,*) "Self-gravitating model"

     topotop_jm=0d0
     topocore_jm=0d0
     sgM=0d0; sgB=0d0
     do j=jmin,jmax   
       fatt=pi4*gi/(2d0*dble(j)+1)

       sgM(1,1)=1d0
       sgM(1,2)=0d0     
       sgM(1,3)=-fatt*radc*derhoc 
       sgM(1,4)=-fatt*erad*(radc/erad)**dble(j)*derhot

       sgM(2,1)=0d0      
       sgM(2,2)=1d0
       sgM(2,3)=-fatt*radc*(radc/erad)**(dble(j)+1d0)*derhoc
       sgM(2,4)=-fatt*erad*derhot 

       sgM(3,1)=-1d0/gref      
       sgM(3,2)=0d0
       sgM(3,3)=1d0
       sgM(3,4)=0d0
  
       sgM(4,1)=0d0
       sgM(4,2)=1d0/gref
       sgM(4,3)=0d0
       sgM(4,4)=1d0

       do m=0,j
         jm=(j*(j+1))/2+m+1

! Real part of the RHS
         sgB(1,1) = fatt*dble(potcld_jm(jm))
         sgB(2,1) = fatt*dble(potsld_jm(jm))
         sgB(3,1) = -dble(tauc(jm))/derhoc/gref
         sgB(4,1) = -dble(taut(jm))/derhot/gref       

! Imaginary part of the RHS
         sgB(1,2) = fatt*dimag(potcld_jm(jm))
         sgB(2,2) = fatt*dimag(potsld_jm(jm))
         sgB(3,2) = -dimag(tauc(jm))/derhoc/gref
         sgB(4,2) = -dimag(taut(jm))/derhot/gref

! Linear system solution
         call gaussj(sgM,sgB)
         Re_sgX=sgB(:,1)
         Im_sgX=sgB(:,2)

! Dynamic topographies corrected for self-gravitation
         topocore_jm(jm)=dcmplx(Re_sgX(3),Im_sgX(3))
         topotop_jm(jm)=dcmplx(Re_sgX(4),Im_sgX(4))  

       end do
     end do

! No self-gravitation     
   else

     do j=jmin,jmax   
       fatt=pi4*gi/(2d0*dble(j)+1)

       do m=0,j
         jm=(j*(j+1))/2+m+1

         topotop_jm(jm)=-taut(jm)/derhot/gref
         topocore_jm(jm)=-tauc(jm)/derhoc/gref
	
       end do
     end do  

   endif

! Calculate the coefficients of total surface potential 
! due to internal load, surface and CMB dynamic topographies
   do j=jmin,jmax   
     fatt=pi4*gi/(2d0*dble(j)+1)

     do m=0,j
       jm=(j*(j+1))/2+m+1

       pot_jm(jm)=fatt*(erad*derhot*topotop_jm(jm)              &
                                            +potsld_jm(jm)) 

 !      write(*,*) 'Core Not Taken into account'
 !      write(*,*) fatt*(-radc*(radc/erad)**(dble(j)+1d0)*derhoc*topocore_jm(jm))
 !      write(*,*) pot_jm(jm)
     end do
   end do  
! -radc*(radc/erad)**(dble(j)+1d0)*derhoc*topocore_jm(jm)  &

! Spherical harmonic synthesis of:
!  total geoid (geoid)
!  static geoid (potld)
!  surface topography (topotop)
!  core toporgaphy (topocore)

   pot=0d0; potld=0d0
   geoid=0d0; topotop=0d0; topocore=0d0; 

   write(*,*) 'Writing output files...'
   open(2,file='geoid.dat')

   do l=1,nt
     th=dacos(roots(l))/deg2rad

     do kph=0,360,2
       ph=kph*deg2rad

       do j=jmin,jmax
         fatt=pi4*gi/(2d0*dble(j)+1)

         do m=-j,j
           jm=j*(j+1)/2+iabs(m)+1           
           leg=plg(l,jm)
           call cyjm(j,m,leg,ph,yjm)

           if(m.ge.0) then
             potld(l,kph)=potld(l,kph)+fatt*potsld_jm(jm)*yjm
             pot(l,kph)=pot(l,kph)+pot_jm(jm)*yjm
             topotop(l,kph)=topotop(l,kph)+topotop_jm(jm)*yjm
             topocore(l,kph)=topocore(l,kph)+topocore_jm(jm)*yjm
           else  
             zn=(-1d0)**m      
             potld(l,kph)=potld(l,kph)+zn*fatt*potsld_jm(jm)*conjg(yjm)
             pot(l,kph)=pot(l,kph)+zn*pot_jm(jm)*conjg(yjm)
             topotop(l,kph)=topotop(l,kph)+zn*topotop_jm(jm)*conjg(yjm)
             topocore(l,kph)=topocore(l,kph)+zn*topocore_jm(jm)*conjg(yjm)
           endif

           geoid(l,kph)=pot(l,kph)/gref


         end do  
       end do  

       write(2,90) real(kph),90-th,geoid(l,kph),potld(l,kph)/gref,topotop(l,kph),topocore(l,kph)
     end do  
   end do 
   close(2)

   write(*,*) '...done.'


90 format(6E15.4) 

!----------------------------------------------------------------------
! ZM modification:
! BR modification - including the potsld_jm to output
   write(*,*) 'Writing spherical harmonics of the dynamic geoid, surface topography'
   open(3,file='geoid_coeff.dat')
       do j=jmin,jmax
         fatt=pi4*gi/(2d0*dble(j)+1)
!         do m=0,j
!           jm=j*(j+1)/2+m+1  
	 do m=-j,j
           jm=j*(j+1)/2+iabs(m)+1          
           write(3,91) j,m,pot_jm(jm),topotop_jm(jm),topocore_jm(jm),fatt*potsld_jm(jm)
         end do
       end do
   close(3)
91 format(2i4,2E17.8,2x,2E17.8,2x,2E17.8,2x,2E17.8) 
!----------------------------------------------------------------------
  

end program sfec13d 

















