subroutine solve(flagcrash)

!use pks 
use system
use results
use solver
use const

implicit none
real*16 x2beta,x2alpha
real*16 x2beta0,x2alpha0
real*16 x3beta,x3alpha
real*16 x3beta0,x3alpha0
real*16 tolerancia,criterio,check_Ka_alpha,check_ka_beta,checkb_Kai,KK0check,KKaAcheckplus,kkaBcheckmin
real*16 check_Ko_alpha,check_Ko_beta,check_Kb_alpha,check_Kb_beta
integer nmax
real*16 xiteri
real*16 xiterii
real*16 xiterj
real*16 xiterjj
real*8 x1(3)
real*8 x1g(3)
integer ier, i,newt, j, ii, jj
integer flagcrash
real*16 potquim2,potquim3,mu2alpha,mu2beta
real*16 xatest(2),xbtest(2),segsolalpha(2),segsolbeta(2),yy(2),yybet(2)
real*16 Mtemp,frac(8),fracbet(8)
real*16 phi0,xphisales,xsolvit
integer flag
integer ngrid
!real*16 arrayalphagrid(2,50000), arraybetagrid(2,50000)
integer gridpoints
real*16 xiter,xsolventalpha,xsolventbeta,xHplusalpha,xHplusbeta,xNaplusalpha,xNaplusbeta
real*16 xClminalpha,xClminbeta
real* 16 testarray(2),testmu2,testmu3,elib1,testfrac(8),xOhminalpha,xOHminbeta
real* 16 testarrayb(2),testmub2,testmub3,elib2,testfracb(8)

integer pos,neg,li
!   phimin=phiminread; !valor minimo del exponente
!   phimax=phimaxread ! valor maximo del exponente
pos=0
neg=0
   nmax=npasos     ! npasos por consola i
   ngrid = npasosgrid
   criterio=1E-8!criterio pra la norma
   tolerancia=1E-3!criterio pra la diferencia de concentraciones relativa
   gridpoints = 0

! #### GRID SEARCH

linearsolver = 2
flaggg=0
! ida

   do i=0,ngrid
   do j=i+1,ngrid
  
      x2alpha  = 10**xiteri ! x2phialpha
      x2beta  = 10**xiterj

      x2alpha = float(i)
      x2beta = float(j)
 

      x1(1)=x2alpha
      x1g(1)=x1(1)
      x1(2)=x2beta     !x2phibeta inicial
      x1g(2)=x1(2)

      call call_kinsol(x1, x1g, ier)

       if ((norma.lt.criterio).and.(abs(x1(1)-x1(2)/(x1(1)+x1(2))).gt.tolerancia))exit

   enddo
       if ((norma.lt.criterio).and.(abs(x1(1)-x1(2)/(x1(1)+x1(2))).gt.tolerancia))exit
   enddo

       if ((norma.lt.criterio).and.(abs(x1(1)-x1(2)/(x1(1)+x1(2))).gt.tolerancia)) then
flaggg=1
         gridpoints = gridpoints + 1
         print*,'Grid Point OK',gridpoints,yes
         print*,'csal,x2alpha,x2beta',csal,min(10**(-x1(1)),10**(-x1(2))),max(10**(-x1(1)),10**(-x1(2)))

         arrayalphagrid(1,gridpoints)=x1(1)
         arrayalphagrid(2,gridpoints)=x1(1)
         arraybetagrid(1,gridpoints)=x1(2)
         arraybetagrid(2,gridpoints)=x1(2)
!arraycsal(1,gridpoints)=csal
!arrayalpha(1,gridpoints)=min(10**(-x1(1)),10**(-x1(2)))
!arraybeta(1,gridpoints)=max(10**(-x1(1)),10**(-x1(2)))
xsolvit=(1.0-2*min(10**(-arrayalphagrid(1,1)),10**(-arraybetagrid(1,1))))
xsolvit=xsolvit/(1.+expmuHplus+expmuOHmin+expmuneg+expmupos)
xphisales=xsolvit*(expmupos+expmuneg)

print*,'xphisales',xphisales

       endif

return
end subroutine

