
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine quaternion      
use MDPARAM,only : NSITE,C,CNEW
implicit none   
real8 :: kai,eta,xi,zeta,p,F
integer :: k
real8,dimension(3,NSITE) :: AM

CNEW=0.0d0
1 call random_number(p)
  kai=2.0d0*p-1.0d0
  call random_number(p)
  eta=2.0d0*p-1.0d0
  call random_number(p)
  xi=2.0d0*p-1.0d0
  F =kai**2+eta**2+xi**2
if(F>=1.0d0) goto 1
zeta=sqrt(1.0d0-F)
AM(1,1:3)=(/ -xi**2+eta**2-zeta**2+kai**2,2.0d0*(zeta*kai-xi*eta),2.0d0*(eta*zeta+xi*kai) /)
AM(2,1:3)=(/ -2.0d0*(xi*eta+zeta*kai),xi**2-eta**2-zeta**2+kai**2,2.0d0*(eta*kai-xi*zeta) /)
AM(3,1:3)=(/ 2.0d0*(eta*zeta-xi*kai),-2.0d0*(xi*zeta+eta*kai),-xi**2-eta**2+zeta**2+kai**2 /)

CNEW=matmul(AM,C)

return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
