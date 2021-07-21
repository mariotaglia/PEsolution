subroutine salvar(flagcrash)
use results
implicit none
integer i,j,flagcrash
real*16 l

open (unit=1,file='alphaarray.txt',status='replace')

do i=1,conteo
   write (1,*) arrayalpha(1,i), arrayalpha(2,i)
end do

open (unit=1,file='betaarray.txt',status='replace')

do j=1,conteo
   write (1,*) arraybeta(1,j), arraybeta(2,j)
end do

open (unit=2,file='alphacarga.txt',status='replace')

do i=1,conteo
   write (2,*) cargaalpha(1,i)
end do

open (unit=2,file='betacarga.txt',status='replace')

do j=1,conteo
   write (2,*) cargabeta(1,j)
end do

open (unit=4,file='arraysolvalpha.txt',status='replace')

do j=1,conteo
   write (4,*) arraysolvalpha(1,j)
end do

open (unit=4,file='arraysolvbeta.txt',status='replace')

do j=1,conteo
   write (4,*) arraysolvbeta(1,j)
end do


open (unit=5,file='arrayNaplusalpha.txt',status='replace')

do j=1,conteo
   write (5,*) arrayalpha(1,j), arrayNaplusalpha(1,j)
end do

open (unit=5,file='arrayNaplusbeta.txt',status='replace')

do j=1,conteo
   write (5,*) arraybeta(1,j), arrayNaplusbeta(1,j)
end do


open (unit=5,file='sumarrayNaClalpha.txt',status='replace')

do j=1,conteo
   l=arrayalpha(1,j)+arrayalpha(2,j)
   write (5,*) l, sumalphaNacl(1,j)
end do

open (unit=5,file='sumarrayNaClbeta.txt',status='replace')

do j=1,conteo
   l=arraybeta(1,j)+arraybeta(2,j)
   write (5,*) l, sumbetaNacl(1,j)
end do


open (unit=6,file='arrayClminalpha.txt',status='replace')

do j=1,conteo
   write (6,*) arrayClminalpha(1,j)
end do

open (unit=6,file='arrayClminbeta.txt',status='replace')

do j=1,conteo
   write (6,*) arrayClminbeta(1,j)
end do

open (unit=7,file='arrayHplusalpha.txt',status='replace')

do j=1,conteo
   write (7,*) arrayHplusalpha(1,j)
end do

open (unit=7,file='arrayHplusbeta.txt',status='replace')

do j=1,conteo
   write (7,*) arrayHplusbeta(1,j)
end do



open (unit=8,file='arrayOHminalpha.txt',status='replace')

do j=1,conteo
   write (8,*) arrayOHminalpha(1,j)
end do

open (unit=8,file='arrayOHminbeta.txt',status='replace')

do j=1,conteo
   write (8,*) arrayOHminbeta(1,j)
end do



open (unit=9,file='array_fchA_alpha.txt',status='replace')

do j=1,conteo
   write (9,*) array_fchA_alpha(1,j)
end do

open (unit=9,file='array_fchA_beta.txt',status='replace')

do j=1,conteo
   write (9,*) array_fchA_beta(1,j)
end do


open (unit=9,file='array_fchB_alpha.txt',status='replace')

do j=1,conteo
   write (9,*) array_fchB_alpha(1,j)
end do

open (unit=9,file='array_fchB_beta.txt',status='replace')

do j=1,conteo
   write (9,*) array_fchB_beta(1,j)
end do


open (unit=10,file='array_fucA_alpha.txt',status='replace')

do j=1,conteo
   write (10,*) array_fucA_alpha(1,j)
end do

open (unit=10,file='array_fucA_beta.txt',status='replace')

do j=1,conteo
   write (10,*) array_fucA_beta(1,j)
end do


open (unit=10,file='array_fucB_alpha.txt',status='replace')

do j=1,conteo
   write (10,*) array_fucB_alpha(1,j)
end do

open (unit=10,file='array_fucB_beta.txt',status='replace')

do j=1,conteo
   write (10,*) array_fucB_beta(1,j)
end do



open (unit=11,file='array_fIPA_alpha.txt',status='replace')

do j=1,conteo
   write (11,*) array_fIPA_alpha(1,j)
end do

open (unit=11,file='array_fIPA_beta.txt',status='replace')

do j=1,conteo
   write (11,*) array_fIPA_beta(1,j)
end do


open (unit=11,file='array_fIPB_alpha.txt',status='replace')

do j=1,conteo
   write (11,*) array_fIPB_alpha(1,j)
end do

open (unit=11,file='array_fIPB_beta.txt',status='replace')

do j=1,conteo
   write (11,*) array_fIPB_beta(1,j)
end do


open (unit=12,file='array_fPPA_alpha.txt',status='replace')

do j=1,conteo
   write (12,*) array_fPPA_alpha(1,j)
end do

open (unit=12,file='array_fPPA_beta.txt',status='replace')

do j=1,conteo
   write (12,*) array_fPPA_beta(1,j)
end do


open (unit=12,file='array_fPPB_alpha.txt',status='replace')

do j=1,conteo
   write (12,*) array_fPPB_alpha(1,j)
end do

open (unit=12,file='array_fPPB_beta.txt',status='replace')

do j=1,conteo
   write (12,*) array_fPPB_beta(1,j)
end do


open (unit=112,file='fe_alpha_array.txt',status='replace')

do i=1,conteo
   write (112,*) arrayalpha(1,i) ,array_fe_alpha(1,i) 

end do

open (unit=112,file='fe_beta_array.txt',status='replace')

do i=1,conteo
   write (112,*) arraybeta(1,i) ,array_fe_beta(1,i)

end do

open (unit=112,file='mu2_alpha_array.txt',status='replace')

do i=1,conteo
   write (112,*) arrayalpha(1,i) ,array_mu2_alpha(1,i)

end do

open (unit=112,file='mu2_beta_array.txt',status='replace')

do i=1,conteo
   write (112,*) arraybeta(1,i) ,array_mu2_beta(1,i)

end do


open (unit=112,file='mu3_alpha_array.txt',status='replace')

do i=1,conteo
   write (112,*) arrayalpha(1,i) ,array_mu3_alpha(1,i)

end do

open (unit=112,file='mu3_beta_array.txt',status='replace')

do i=1,conteo
   write (112,*) arraybeta(1,i) ,array_mu3_beta(1,i)

end do

!open (unit=2,file='betaarrayspi.txt',status='replace')
!
!do j=1,cont
!   write (2,*) arraybetaspi(1,j), arraybetaspi(2,j)
!end do


end subroutine

