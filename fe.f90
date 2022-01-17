subroutine fe(x,elib,xmsalt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine calculates the free energy of the system and chemical potential 
! of chains
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use system
!use kai
!use bulk
!use molecules
use const
implicit none

real*16 x(2)
real*16 elib,xmsolv,xmOHmin,xmHplus, xmNaCl, xmsalt
real*16 Free_energy 
integer i, iz, iiz
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solventireal
real*16 xphiA,xphiB,xsolv,fa_A,fa_B,xmphiA,xmphiB,fnc_a,fnc_B,fc_A,fc_B,fNa_a,fCl_b
real*16 frac(8)
real*16 aa,bb,cc

! Calculation of xsolvent
xphiA=x(1)
xphiB=x(2)

aa = expmupos*expmuneg*Ksal*2.
bb = (1.+expmupos+expmuneg+expmuHplus+expmuOHmin)
cc = -(1.0 -xphiA-xphiB)

xsolv = (-bb + sqrt(bb**2 - 4.*aa*cc))/(2.0*aa) ! solvent volume fraction
xmsolv=xsolv/vs ! solvent number density

xmHplus=xmsolv*expmuHplus ! number density
xmOHmin=xmsolv*expmuOHmin ! number density
xmNaplus=xmsolv*expmupos  ! number density
xmClmin=xmsolv*expmuneg   ! number density
xmNaCl = ksal/vs*expmupos*expmuneg ! number density

call fracasos(x,frac)
fa_A=frac(1)
fa_B=frac(2)
fnc_A=frac(3)
fnc_B=frac(4)
fc_A=frac(5)
fc_B=frac(6)
fNa_a=frac(7)
fCl_b=frac(8)

xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)

Free_Energy = 0.0

! 1. solvent entropy

  Free_Energy=Free_Energy-xmphiA-xmphiB-xmsolv-xmHplus-xmOHmin-xmNaplus-xmClmin-xmNaCl

  Free_Energy=Free_Energy +xmphiA*Ma*fa_A+(log(xsolv))/vs

elib=Free_Energy
xmsalt = (xmNaplus+xmClmin)/2.0 + xmNaCl
return 
end subroutine

