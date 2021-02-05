!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*
*
*  H2O 
*
*
* long-range calculation
*
* 
* 
*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module EWALD
implicit none
integer :: TOTAL_K
integer,parameter :: MAX_K=1000,KMAX=5,KSQMAX=20
real8,dimension(MAX_K) :: KVEC ! store h-vector
real8 :: alpha
end module EWALD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EWALD_SETUP
use MDPARAM,only : PI,RBOX,IOUT
use EWALD
implicit none
integer :: kx,ky,kz   
integer ::ksq
real8 :: rkx,rky,rkz
real8 :: rksq
real8 :: tmp
real8 :: AL
alpha=5.0d0/RBOX
AL=0.25d0/(alpha**2)
TOTAL_K=0  ! number of k
tmp=2.0d0*PI/RBOX
do kx=-kmax,kmax
  rkx=tmp*dble(kx)
  do ky=-kmax,kmax
    rky=tmp*dble(ky)
    do kz=-kmax,kmax
      rkz=tmp*dble(kz)
      ksq=kx**2+ky**2+kz**2
      if( (ksq<ksqmax) .and. (ksq /= 0)) then 
        total_k=total_k+1
        if(total_k > max_k) stop
        rksq=rkx**2+rky**2+rkz**2 ! |G|^2  G=2*pi/RBOX*k
        kvec(total_k)=exp(-rksq*AL)/rksq
      endif
    enddo
  enddo
enddo

write(IOUT,*) "# TOTAL_K:                ",total_k
end   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine R_EWALD(FR,UR)
use MDPARAM,only : N,NSITE,Q,PI,A,B,RDFFLG,RBOX,RCUT
use VARIABLES,only : R
use EWALD,only : alpha
implicit none
integer :: i,j,k,i2,i3,num
real8,intent(out),dimension(3,NSITE,N) :: FR
real8,intent(out) :: UR
real8 :: r12,ftmp,Q2,RCUT2
real8 :: rij(3)
real8 :: alphar
!LJ
real8 :: rij2inv,Arij12inv,Brij6inv,rr
FR=0.0d0
UR=0.0d0
RCUT2=RCUT**2
do i=1,N-1!triangle loop
  do j=i+1,N
    do i2=1,NSITE
      Q2=Q(i2)
      do i3=1,NSITE
        do k=1,3 !xyz
          rij(k)=r(k,i2,i)-r(k,i3,j)
          rij(k)=rij(k)-anint(rij(k)/RBOX)*RBOX
        enddo
        r12=dot_product(rij,rij)
        if(r12>=RCUT2) cycle
        rij2inv=1.d0/r12
        r12=sqrt(r12)
        alphar=alpha*r12
        UR=UR+Q2*Q(i3)*erfc(alphar)/r12
        ftmp=Q2*Q(i3)*(erfc(alphar)/(r12**3)+2.0d0*alpha/sqrt(PI)&
              *exp(-(alphar)**2)/(r12**2))
        do k=1,3
          FR(k,i2,i)=FR(k,i2,i)+ftmp*rij(k)
          FR(k,i3,j)=FR(k,i3,j)-ftmp*rij(k)
        enddo
       !LJ potension
        if(i2==1 .and. i3==1) then
          Arij12inv = A*rij2inv**6
          Brij6inv  = B*rij2inv**3
          UR = UR + Arij12inv - Brij6inv
          rr = 12.0d0*Arij12inv*rij2inv - 6.0d0*Brij6inv*rij2inv
          do k=1,3
            FR(k,1,i)=FR(k,1,i)+rr*rij(k)
            FR(k,1,j)=FR(k,1,j)-rr*rij(k)
          enddo
        endif
      enddo
    enddo
  enddo
enddo
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine I_EWALD(FI,UI)   
use MDPARAM,only : N,NSITE,Q,PI,RBOX
use VARIABLES,only : R
use EWALD
implicit none
integer :: k,kx,ky,kz,i,i2
real8,intent(out) :: UI
real8,intent(out),dimension(3,NSITE,N) :: FI
real8 :: U3,tmp,tmp2,ksq
complex8,dimension(1:NSITE,1:N,-kmax:kmax) :: eikx,eiky,eikz !Euler.Eq.
complex8,dimension(NSITE,N) :: eik
complex8 :: sum_e
UI=0.0d0
FI=0.0d0
U3=0.0d0
TOTAL_K=0
tmp=2.0d0*PI/RBOX
do i=1,N
  do i2=1,NSITE
    U3=U3-Q(i2)**2
    !G*R set up
    !h=0
    eikx(i2,i,0)=cmplx(1.0d0,0.0d0,kind(0d0))
    eiky(i2,i,0)=cmplx(1.0d0,0.0d0,kind(0d0))
    eikz(i2,i,0)=cmplx(1.0d0,0.0d0,kind(0d0))
    !h=1
    eikx(i2,i,1)=cmplx(cos(tmp*R(1,i2,i)),sin(tmp*R(1,i2,i)),kind(0d0))
    eiky(i2,i,1)=cmplx(cos(tmp*R(2,i2,i)),sin(tmp*R(2,i2,i)),kind(0d0))
    eikz(i2,i,1)=cmplx(cos(tmp*R(3,i2,i)),sin(tmp*R(3,i2,i)),kind(0d0))
    !h=-1
    eikx(i2,i,-1)=conjg(eikx(i2,i,1))
    eiky(i2,i,-1)=conjg(eiky(i2,i,1))
    eikz(i2,i,-1)=conjg(eikz(i2,i,1))
  enddo
enddo

do k=2,kmax
  do i=1,N
    do i2=1,NSITE
      
      !kx
      eikx(i2,i,k)=eikx(i2,i,k-1)*eikx(i2,i,1)
      eikx(i2,i,-k)=conjg(eikx(i2,i,k))
      !ky
      eiky(i2,i,k)=eiky(i2,i,k-1)*eiky(i2,i,1)
      eiky(i2,i,-k)=conjg(eiky(i2,i,k))
      !kz
      eikz(i2,i,k)=eikz(i2,i,k-1)*eikz(i2,i,1)
      eikz(i2,i,-k)=conjg(eikz(i2,i,k))
    enddo
  enddo
enddo

do kx=-kmax,kmax
  do ky=-kmax,kmax
    do kz=-kmax,kmax
      ksq=kx**2+ky**2+kz**2
      if( (ksq<ksqmax) .and. (ksq /= 0)) then !EWALD_INIT
        total_k=total_k+1
        SUM_e=cmplx(0.0d0,0.0d0,kind(0d0))
        do i=1,N
          do i2=1,NSITE
            eik(i2,i)=eikx(i2,i,kx)*eiky(i2,i,ky)*eikz(i2,i,kz)! exp(2*pi*hvec/RBOX)
            sum_e=sum_e+Q(i2)*eik(i2,i)
          enddo
        enddo
        UI=UI+kvec(total_k)*(sum_e*conjg(sum_e))
        do i=1,N
          do i2=1,NSITE
            tmp2=kvec(total_k)*(Q(i2)*aimag(eik(i2,i))*real(sum_e)-Q(i2)*real(eik(i2,i))*aimag(sum_e))
            FI(1,i2,i)=FI(1,i2,i)+tmp2*kx*tmp 
            FI(2,i2,i)=FI(2,i2,i)+tmp2*ky*tmp
            FI(3,i2,i)=FI(3,i2,i)+tmp2*kz*tmp
          enddo
        enddo
      endif
    enddo
  enddo
enddo
UI=UI*tmp/(RBOX**2)
FI=FI*4.0d0*PI/(RBOX**3)
UI=UI+U3*alpha/sqrt(PI)
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
