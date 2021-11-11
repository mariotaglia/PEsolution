subroutine fracasos(x,frac)
use system
use const
 implicit none
integer i
real*16 frac(8),x(2)
real*16 xphiA,xphiB,xmphiA,xmphiB,xsolv
real*16 M,Kap,Kbp,Kop,xOHmin,xHplus
real*16 KNap,KClp
real*16 alfaA, alfaB, betaA, betaB, gamaA, gamaB, auxB, auxC

xphiA=x(1);
xphiB=x(2);
xmphiA=xphiA/(Ma*vp);
xmphiB=xphiB/(Mb*vp);

xsolv=(1.0 -xphiA-xphiB)/(1.+expmupos+expmuneg+expmuHplus+expmuOHmin)
xHplus=xsolv*expmuHplus
xOHmin=xsolv*expmuOHmin
xNaplus=xsolv*expmupos
xClmin=xsolv*expmuneg
!print* ,x

betaA=K0A*xsolv/xHplus
betaB=K0B*xsolv/xOHmin

KNap=xNaplus*K0A/(xHplus*K0ANa)
KClp=xClmin*K0B/(xOHmin*K0BCl)

  auxC = xphiB/xphiA

  alfaA = K0ANa*xNaplus/xsolv
  alfaB = K0BCl*xClmin/xsolv

  gamaA = 1.0+alfaA+1.0/betaA
  gamaB = 1.0+alfaB+1.0/betaB

  auxB = -1.0 -auxC -gamaA*gamaB/Kbind0/(xphiA/vpol)

  frac(1) = (-auxB - SQRT(auxB**2 - 4.0*auxC))/2.0 ! fA_a
  frac(2) = frac(1)/auxC ! fB_a

  frac(5) = (1.0 - frac(1))/betaA/gamaA ! fAc
  frac(6) = (1.0 - frac(2))/betaB/gamaB ! fBc

  frac(7) = (1.0 - frac(1))*alfaA/gamaA ! fANa
  frac(8) = (1.0 - frac(2))*alfaB/gamaB ! fBCl

  frac(3)=(1.-frac(1)-frac(5)-frac(7)) !fAnc
  frac(4)=(1.-frac(2)-frac(6)-frac(8)) !fBnc

! viejo

!Knew=Kop*Kap*Kbp/((1.0+Kap+KNap)*(1.0+Kbp+KClp))
!frac(2)=0.5*(1.0+etap+etap/Knew )-( (0.5*(1.0+etap+etap/Knew))**2-etap)**0.5;!fB_a
!frac(1)=frac(2)*eta;!fA_a
!frac(3)=(1.-frac(1))/(1.+Kap) ;!fAnc
!frac(4)=(1.-frac(2))/(1.+Kbp) ;!fBnc
!frac(5)=(1.-frac(1))*Kap/(1.+Kap+KNap);!fAc
!frac(6)=(1.-frac(2))*Kbp/(1.+Kbp+KClp);!fBc
!frac(7)=KNap*(1.-frac(1)-frac(5))/(1.+KNap) !fANa
!frac(8)=KClp*(1.-frac(2)-frac(6))/(1.+KClp) !fBCl
!frac(3)=(1.-frac(1)-frac(5)-frac(7)) !fAnc
!frac(4)=(1.-frac(2)-frac(6)-frac(8)) !fBnc

!print* ,xsolv
!print* , 'f',frac
!stop
!do i = 1,6
!  if (frac(i) .lt. 10**(-15))then
!    frac(i)=10**(-15)
! endif
!enddo
                                  
end subroutine
