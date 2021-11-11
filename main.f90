  !###############################################################################!     
  !     Simple brush: Standard Molecular Theory Program 
  !    
  !     Calculates a weak polyelectrolyte brush in poor sv conditions 
  !     Calculates free-energy
  !     Solve BINODAL METHOD
  !     MARCH 2020
  !         
  !###############################################################################
!use pks
use system
use const
use results
  implicit none
  integer i, flagcrash,iii

  real*16 pKw,vsal,KAinput,Kbind,KBinput,KANainput,KBClinput,kte

  real*16 xposbulk,xnegbulk,xHplusbulk,xOHminbulk,Kw,xsalt
  real*16 xphisales,xsolvit,phialphamol,phibetamol
  real*16 saleihalphamol,saleihbeta,testsolvbeta,testsolvalpha 
  real*16 xsolvalpha,xsolvbeta,saleihbetamol
!print*, 'Program Simple Brush'
  !print*, 'GIT Version: ', _VERSION
  yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
  flagcrash=1
  call readinput  ! reads input variables from file
!  call allocation ! allocates memory

  Kbind=10**(-pKeo)
  KAinput=10**(-pKaA)
  KBinput=10**(-pKaB)
  KANainput=10**(-pKaNa)
  KBClinput=10**(-pKaCl)

  pKw=14.0
  kW=10**(-pKw)
  pOHbulk=pKw-pHbulk
  cOHminbulk=10**(-pOHbulk)
  cHplusbulk=10**(-pHbulk)

  vp=vpol!
  vsal=vs/vs! Necesario segun notas !
  
  xposbulk=phi_sal  !NUEVO
  xnegbulk=phi_sal  !NUEVO
  
  xHplusbulk = (cHplusbulk*Na/(1.0d24))*(vs)
  xOHminbulk = (cOHminbulk*Na/(1.0d24))*(vs) 

  do i = 1, ncsal ! loop in csal

  csal = csalini + (csalfin-csalini)/float(ncsal-1)*float(i-1) 


  xsalt=csal*vsal*vs*6.02/10.0

  if(pHbulk.le.7) then  ! pH<= 7
     xposbulk=xsalt
     xnegbulk= xsalt +(xHplusbulk -xOHminbulk)*vsal! NaCl+ HCl  
  else                  ! pH >7 
     xposbulk=xsalt +(xOHminbulk -xHplusbulk)*vsal ! NaCl+ NaOH   
     xnegbulk=xsalt
  endif

  xsolbulk=1.0 -xHplusbulk -xOHminbulk - xnegbulk -xposbulk

  K0A = (KAinput*vs/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
  K0B = (Kw/KBinput*vs/xsolbulk)*(Na/1.0d24)

  K0ANa = (KANainput/vs)/(Na/1.0d24)! Esta definida igual que en el paper JCP
  K0BCl = (KBClinput/vs)/(Na/1.0d24)! Esta definida igual que en el paper JCP

  Kbind0 = (Kbind)*(1.0d24/Na) !!!!!!!!!!!!!!!!!!!!!!!
  

!  cHplusbulk=xHplusbulk
!  cOHminbulk=xOHminbulk
!csalt=xsalt/(vsal*vs*6.02/10)
expmupos=xposbulk/xsolbulk**vsal
expmuneg=xnegbulk/xsolbulk**vsal
expmuHplus=xHplusbulk/xsolbulk ! vHplus=vsol
expmuOHmin=xOHminbulk/xsolbulk ! vOHminus=vsol

  call solve(flagcrash)
!call fe(cc, ccc)         ! calculates and saves free energy to disk
!  call salvando(flagcrash)
!  print*, 'Save OK',yes

if (flaggg==1)then
yes=yes+1


arraycsalmol(1,yes)=csal
arraycsal(1,yes)=xsalt
arrayalpha(1,yes)=min(10**(-arrayalphagrid(1,1)),10**(-arraybetagrid(1,1)))
arraybeta(1,yes)=max(10**(-arrayalphagrid(1,1)),10**(-arraybetagrid(1,1)))

phialphamol=arrayalpha(1,yes)*10/(6.02*vp)
arrayalphamol(1,yes)=phialphamol*2
phibetamol=arraybeta(1,yes)*10/(6.02*vp)
arraybetamol(1,yes)=phibetamol*2

saleihalphamol=arrayalpha(1,yes)*10/(6.02*0.195)  !usando que vp =0.195  de saleih
arraysaleihalpha(1,yes)=saleihalphamol*2

saleihbetamol=arraybeta(1,yes)*10/(6.02*vp)    !usando que vp=0.195 de saleih
arraysaleihbeta(1,yes)=saleihbetamol*2

xsolvalpha=(1.-arrayalpha(1,yes)-arrayalpha(1,yes)-xsalt)/(1.+expmuHplus+expmuohmin)
testsolvalpha=(1.-arrayalpha(1,yes)-arrayalpha(1,yes)-xsalt)
arraysolvalpha(1,yes)= xsolvalpha*10/(6.02*0.03)

xsolvbeta=(1.-arraybeta(1,yes)-arraybeta(1,yes)-xsalt)/(1.+expmuHplus+expmuohmin)
testsolvbeta=(1.-arraybeta(1,yes)-arraybeta(1,yes)-xsalt)
arraysolvbeta(1,yes)= xsolvbeta*10/(6.02*0.03)

print*,'Solvent, testsolvent',xsolvalpha,testsolvalpha
print*,'Solvent, testsolvent',xsolvbeta,testsolvbeta

arrayconstalpha(1,yes)=2*arrayalpha(1,yes)+xsolvalpha+xsalt
arrayconstbeta(1,yes)=2*arraybeta(1,yes)+xsolvbeta+xsalt

porcenpolalpha(1,yes)=100*arraysaleihalpha(1,yes)/(arraysaleihalpha(1,yes)+csal+arraysolvalpha(1,yes))
porcenpolbeta(1,yes)=100*arraysaleihbeta(1,yes)/(arraysaleihbeta(1,yes)+csal+arraysolvbeta(1,yes))
porcensalalpha(1,yes)=100*csal/(arraysaleihalpha(1,yes)+csal+arraysolvalpha(1,yes))
porcensalbeta(1,yes)=100*csal/(arraysaleihbeta(1,yes)+csal+arraysolvbeta(1,yes))

print*, 'xphisales, phialpha, phibeta', arraycsal(1,yes),arrayalpha(1,yes),arraybeta(1,yes)
print*,'csalmol, alphamol,betamol',arraycsal(1,yes),arraybetamol(1,yes),arrayalphamol(1,yes)
endif


enddo ! i

open (unit=1,file='alphaarray.txt',status='replace')

do iii=1,yes
   write (1,*) arrayalpha(1,iii), arraycsal(1,iii)
end do

open (unit=2,file='betaarray.txt',status='replace')

do iii=1,yes
   write (2,*) arraybeta(1,iii), arraycsal(1,iii)
end do


open (unit=3,file='csalmolalpha.txt',status='replace')

do iii=1,yes
   write (3,*) arrayalphamol(1,iii), arraycsalmol(1,iii)
end do

open (unit=4,file='csalmolbeta.txt',status='replace')

do iii=1,yes
   write (4,*) arraybetamol(1,iii), arraycsalmol(1,iii)
end do


open (unit=5,file='saleihalpha.txt',status='replace')

do iii=1,yes
   write (5,*) arraysaleihalpha(1,iii), arraycsal(1,iii)
end do

open (unit=6,file='saleihbeta.txt',status='replace')

do iii=1,yes
   write (6,*) arraysaleihbeta(1,iii), arraycsal(1,iii)
end do


open (unit=7,file='porcalpha.txt',status='replace')

do iii=1,yes
   write (7,*) porcenpolalpha(1,iii), porcensalalpha(1,iii)
end do

open (unit=8,file='porcbeta.txt',status='replace')

do iii=1,yes
   write (8,*) porcenpolbeta(1,iii), porcensalbeta(1,iii)
end do

  call endall     ! clean up and terminate


end 
subroutine endall
stop
end subroutine
