********************************************************
subroutine EulerAngler
use MDPARAM,only : NSITE,C,CNEW
implicit none
real8 :: cosphi,sinphi,cospsi,sinpsi,costhe,sinthe
integer :: i
call random_number(cosphi)
call random_number(costhe)
call random_number(cospsi)
cosphi=2.d0*cosphi-1.d0
costhe=2.d0*costhe-1.d0
cospsi=2.d0*cospsi-1.d0
sinphi=sqrt(1.d0-cosphi**2)
sinthe=sqrt(1.d0-costhe**2)
sinpsi=sqrt(1.d0-cosphi**2)

do i=1,NSITE
  CNEW(1,i)=(cospsi*cosphi-costhe*sinphi*sinpsi)*C(1,i) &
    +(cospsi*sinphi+costhe*cosphi*sinpsi)*C(2,i) &
    +sinpsi*sinthe*C(3,i)
  CNEW(2,i)=(-sinpsi*cosphi-costhe*sinphi*cospsi)*C(1,i) &
    +(-sinpsi*sinphi+costhe*cosphi*cospsi)*C(2,i) &
    +cosphi*sinthe*C(3,i)
  CNEW(3,i)=sinthe*sinphi*C(1,i) &
    -sinthe*cosphi*C(2,i) &
    +costhe*C(3,i)
enddo
return
end
********************************************************
