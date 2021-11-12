subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use pks
use system
use solver
!use kai

implicit none
integer i
character basura

! read starts here, not that read is performed sequentially! 

read(8,*), basura
read(8,*), Ma    ! Ma  polA

read(8, *), basura
read(8, *), Mb  !   Mb  polB

read(8, *), basura! vp
read(8, *), vpol !

read(8, *), basura! rsal
read(8, *), rsal !

read(8, *), basura!
read(8, *), phi_sal

read(8, *), basura!
read(8, *), phi_pol

read(8, *), basura
read(8, *), npasos  !   Npasos

read(8, *), basura
read(8, *), npasosgrid  !   Npasos

read(8, *), basura
read(8, *), phimin, phimax  !   Npasos


read(8, *), basura
read(8, *), csalini,csalfin,ncsal  ! salt concentration in bulk (Molar)
!
read(8, *), basura
read(8, *), pKaA    ! pKaA of weak polyacid segments
!
read(8, *), basura
read(8, *), pKaB    ! pKaB of weak polyacid segments
!
read(8, *), basura
read(8, *), pHbulk    ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKeo     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKaNa     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKaCl     ! polymer-polymer attraction strenght in kBT

read(8, *), basura
read(8, *), ksal  ! cutoff for porr sv interaction in lattice sites

!read(8, *), basura
!read(8, *), infile ! read input from file?

end subroutine

