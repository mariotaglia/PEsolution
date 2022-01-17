
subroutine mu2(x,potquim2)
use const
use system
implicit none
real*16 x(2)
real*16 potquim2,fa_A,fa_B,fnc_B,fnc_a,fc_A,fc_b,fNa_a,fCl_b
real*16 xphiA,xphiB,xmphiA,xmphiB,xsolv,xHplus
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
!frac=fracasos(x,xphiB)
fa_A=frac(1)
fa_B=frac(2)
fnc_A=frac(3)
fnc_B=frac(4)
fc_a=frac(5)
fc_b=frac(6)
fNa_a=frac(7)
fCl_b=frac(8)
potquim2= log(xmphiA*vs)-(Ma*vp/vs)*log(xsolv)+Ma*(log(fc_A))!-log(xHplus)+log(xsolv))

!print* ,potquim2
!stop
end  subroutine 
