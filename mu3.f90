
subroutine mu3(x,potquim3)
use const
use system
implicit none
real*16 x(2)
real*16 potquim3,fa_A,fa_B,fc_a,fc_b,fnc_A,fnc_b,fcl_b,fNa_a
real*16 xphiA,xphiB,xmphiA,xmphiB,xsolv,xOHmin
real*16 frac(8)
real*16 aa,bb,cc

xphiA=x(1)
xphiB=x(2)
xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)

aa = expmupos*expmuneg*Ksal*2.
bb = (1.+expmupos+expmuneg+expmuHplus+expmuOHmin)
cc = -(1.0 -xphiA-xphiB)

xsolv = (-bb + sqrt(bb**2 - 4.*aa*cc))/(2.0*aa) ! solvent volume fraction

call fracasos(x,frac)
fa_A=frac(1)
fa_B=frac(2)
fnc_a=frac(3)
fnc_b=frac(4)
fc_A=frac(5)
fc_B=frac(6)
fNa_a=frac(7)
fCl_b=frac(8)
potquim3= log(xmphiB*vs)-(Mb*vp/vs)*log(xsolv)+Mb*(log(fc_B))-chi*Ma*Ma*(xmphiA+xmphiB)
!return (potquim3)
!print* ,potquim3
!stop
end subroutine 
