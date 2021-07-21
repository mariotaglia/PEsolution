!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)

      use mkinsol
      use system
      implicit none

      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)

      common /psize/ neq

      do  i = 1, neq
         vv(i) = vv(i) * pp(i)
      enddo
      ier = 0

      return
      end

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     The routine kpreco is the preconditioner setup routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

      subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

      use mkinsol
      use system
      implicit none
      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)


      common /psize/ neq

      do i = 1, neq
         pp(i) = 0.1  / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
      enddo

!      do i = dimz+1, 2*dimz
!          pp(i) = 1.0
!  !        pp(i)= 0.1 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
!      enddo

!      do i = 2*dimz+1,3*dimz
!         pp(i) = 0.1 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
!      enddo

      ier = 0

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine calls kinsol solver
!

      subroutine call_kinsol(x1, xg1, ier)

      use solver   
      use system
      implicit none

      real*8 x1(2), xg1(2)
      integer*8 iout(15) ! Kinsol additional output information
      real*8 rout(2) ! Kinsol additional out information
      integer*8 msbpre
      real*8 fnormtol, scsteptol
      real*8 uscale(2)
      real*8 fscale(2)
      real*8 constr(2)
      integer*4  globalstrat, maxl, maxlrst
      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations
      integer i
      integer max_niter
      real*8 rerr_func
      real*8 eta_const

      common /psize/ neq ! Kinsol

! INIT KINSOL

      neq = 2
      msbpre  = 5 ! maximum number of iterations without prec. setup (?)
      fnormtol = 1.0d-10 ! Function-norm stopping tolerance
      scsteptol = 1.0d-10 ! Function-norm stopping tolerance

      maxl = 5000 ! maximum Krylov subspace dimesion
      maxlrst = 2 ! maximum number of restarts
      globalstrat = 0
      
      max_niter = 10000

      call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
      if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
        print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
        stop
      endif

      call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
      if (ier .ne. 0) then
         print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
         stop
      endif

      call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      call fkinsetrin('SSTEP_TOL', scsteptol, ier)
      call fkinsetiin('MAX_NITERS', max_niter, ier)

         do i = 1, neq  !constraint vector
         constr(i) = 1.0	!mayor igual que cero
         enddo

         call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector

!       CALL FKINSPTFQMR (MAXL, IER)
if(linearsolver.eq.0)call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system 
if(linearsolver.eq.1)call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG
if(linearsolver.eq.2)call fkindense(neq,ier)
 
      if (ier .ne. 0) then
      print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
         call fkinfree ! free memory
         stop
      endif
!      call fkinspilssetprec(1, ier) ! preconditioner


      do i = 1, neq ! scaling vector
      uscale(i) = 1.0E5
      fscale(i) = 1.0
      enddo

      do i = 1, neq ! Initial guess
      x1(i) = xg1(i)
      enddo


!      print*, 'Entrada', 10**(-x1(1)), 10**(-x1(2)), csal
      call fkinsol(x1, globalstrat, uscale, fscale, ier)         ! call kinsol
!      print*, 'Salida', 10**(-x1(1)), 10**(-x1(2)), norma

      


!      if (ier .lt. 0) then
!      print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
!      print*, 'Salida', 10**(-x1(1)), 10**(-x1(2)), norma
!      print*, 'Linear Solver returned IER = ', iout(9) 

      call fkinfree ! free memory

!      endif
!      stop

      return
      end subroutine



