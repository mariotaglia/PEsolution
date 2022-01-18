!  ###############################################################################!     
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

  real*16 rhoNa, rhoCl, rhoNaCl
  real*16 pKw,vsal,KAinput,Kbind,KBinput,KANainput,KBClinput,kte
  real*16 aa, bb, cc, xNaClbulk
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

! Concentration of free anf paired Na+ and Cl- in bulk reference

  rhoNa = xposbulk/vs
  rhoCl = xnegbulk/vs

  aa = 1.
  bb = -1.*(rhoNa+rhoCl+1./Ksal/vs)
  cc = rhoNa*rhoCl

  rhoNaCl = ((-bb - sqrt(bb**2 - 4.*aa*cc))/(2.0*aa))

  rhoNa = rhoNa - rhoNaCl
  rhoCl = rhoCl - rhoNaCl
 
!  print*, rhoNa*rhoCl/rhoNaCl

  xposbulk = rhoNa*vs
  xnegbulk = rhoCl*vs
  xNaClbulk = rhoNaCl*2.*vs

!  print*, '!!!!', xposbulk, xnegbulk, xNaClbulk
!  print*, '$$$$', (xposbulk/vs)*(xnegbulk/vs)/(xNaClbulk/2.0/vs), Ksal

  if(xNaClbulk.lt.0.0)cycle ! if negative, go to next salt concentration

!!!!!!!!!!!!!!!!!!!!!

  xsolbulk=1.0 -xHplusbulk -xOHminbulk - xnegbulk -xposbulk -xNaClbulk

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

if (flaggg==1)then
yes=yes+1

arraycsalmol(1,yes)=xmsaltalpha/Na*1.d24 
arraycsalmol(2,yes)=xmsaltbeta/Na*1.d24 

arraycsalmolbulk(1,yes)=((xnegbulk+xposbulk)/2.0+xNaClbulk)/Na*1.d24/vs 
arraycsalmolbulk(2,yes)=((xnegbulk+xposbulk)/2.0+xNaClbulk)/Na*1.d24/vs

arraycsal(1,yes)=xmsaltalpha
arraycsal(2,yes)=xmsaltbeta

arraypolmol(1,yes)=xpolalpha*1.e24/Na/vp*2  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
arraypolmol(2,yes)=xpolbeta*1.e24/Na/vp*2

endif


enddo ! i

open (unit=3,file='csalmolalpha.txt',status='replace')

do iii=1,yes
   write (3,*) arraypolmol(1,iii), arraycsalmol(1,iii)
end do

open (unit=4,file='csalmolbeta.txt',status='replace')

do iii=1,yes
   write (4,*) arraypolmol(2,iii), arraycsalmol(2,iii)
end do


open (unit=30,file='csalbulkalpha.txt',status='replace')

do iii=1,yes
   write (30,*) arraypolmol(1,iii), arraycsalmolbulk(1,iii)
end do

open (unit=40,file='csalbulkbeta.txt',status='replace')

do iii=1,yes
   write (40,*) arraypolmol(2,iii), arraycsalmolbulk(2,iii)
end do

 
  call endall     ! clean up and terminate


end 
subroutine endall
stop
end subroutine
