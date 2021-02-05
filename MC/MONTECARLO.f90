************************************************

 We calculate RDF using  montecarlo

 Moleculer : Ar 

 NVT montecarlo

 Judge metropolis

***********************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module RDFPARAM
implicit none
integer,parameter :: NBIN=1000
real8 :: RBIN  
real8,dimension(NBIN) :: HIST 
real8 :: ave,dis
end module RDFPARAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MCNVT(ISTEP)
use MCPARAM,only : N,drmax
use VARIABLES,only : adopt,reject,R,U
implicit none
integer,intent(in) :: istep
integer :: i,j,k
real8 :: rand(3)
real8 :: deltaU,Ut
real8,dimension(3) :: rnew,rnow
real8 :: Unew,Uold
Ut=0.d0
do i=1,N
  rnow(1:3)=r(1:3,i)
  call potential(rnow,i,Uold)
  call random_number(rand)
  do k=1,3
    rnew(k) = rnow(k) + (2.d0*rand(k)-1.d0)*drmax
  enddo
  call potential(rnew,i,unew)
  deltaU=Unew-Uold
  call random_number(rand)
  if(deltaU<0.d0) then
    do k=1,3
      r(k,i)=rnew(k)
    enddo
    Ut=Ut+deltaU
    adopt=adopt+1
  else if(exp(-beta*deltaU) > rand(1)) then
    do k=1,3
      r(k,i)=rnew(k)
    enddo
    Ut=Ut+deltaU
    adopt=adopt+1
  else
    reject=reject+1
  endif
enddo
U=U+Ut
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine potential(Ri,i,Ui)
use MCPARAM,only : sig,N,RBOX,RCUT,eps
implicit none
integer,intent(in) :: i
real8,dimension(3),intent(in) :: Ri
real8,intent(out) :: Ui

integer :: j,k
real8 :: rcut2,rij2,SS6
real8,dimension(3) :: rij

rcut2=rcut**2
Ui=0.d0
do j=1,N
  if(j==i) cycle
  do k=1,3
    rij(k)=ri(k)-r(k,j)
    rij(k)=rij(k)-anint(rij(k)/RBOX)*RBOX
  enddo
  rij2=dot_product(rij,rij)
  if(rij2 >= rcut2 ) cycle
  SS6=sig**6/rij2**3
  Ui=Ui + 4.d0*EPS*(SS6**2-SS6)
enddo
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RDF(istep)
use MCPARAM,only : N,RBOX
use VARIABLES,only : R
use RDFPARAM
implicit none
integer,intent(in) :: istep
integer :: i,j,k,num
real8,dimension(3) :: rij
real8 :: rij2,rijsq
real8 :: tmp

do i=2,N
  do j=1,i-1
    do k=1,3
      rij(k)=r(k,j)-r(k,i)
      rij(k)=rij(k)-anint(rij(k)/RBOX)*RBOX
    enddo
    rij2=dot_product(rij,rij)
    rijsq=sqrt(rij2)
    num=int(rijsq/RBIN)+1
    if (num>nbin) cycle
    hist(num)=hist(num)+1.d0
    ave=ave+rijsq
    dis=dis+rij2
  enddo
enddo

return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RDFOUTPUT
use MCPARAM,RDFOUT,FRDF,N,BOHR
use VARIABLES
use RDFPARAM
implicit none
integer :: i
real8 :: tmp

  open(RDFOUT,file=FRDF,form='formatted')
  ave=ave/dble((maxstp+1)*n*(n-1)/2)
  dis=dis/dble((maxstp+1)*n*(n-1)/2)-ave**2

  hist=hist/dble(maxstp+1)
   
  tmp=1.d0/(4.d0*PI)*dble(NBIN**3)/(dble(n*(n-1)))
  do i=1,NBIN
    hist(i)=hist(i)*tmp
    write(RDFOUT,*) dble(i)*RBIN*BOHR,HIST(i)/dble(i**2)*8.d0 ! rbox = 2*rcut
  enddo
  write(RDFOUT,*) '#AVG,DIS=',ave,dis
  close(RDFOUT)

return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SETMC 
use MCPARAM,only : N,RBOX,RCUT,rijmin,SIG,EPS
use VARIABLES,only : R,U
implicit none
integer :: i,j,k,num
real8 :: SS
real8,dimension(3) :: r12
real8 :: rij,rr,rcut2,rijmin2
rijmin2=rijmin**2
rcut2=rcut**2
U=0.0d0

do i=2,N
  do j=1,i-1
    do k=1,3
      r12(k)=r(k,i)-r(k,j)
      r12(k)=r12(k)-rbox*anint(r12(k)/rbox)
    enddo
    rij=dot_product(r12,r12)
    if(rij<rijmin2) then
      write(*,*) 'initial configulation is overlapped'
      stop
    else if(rij<RCUT2) then
      rij=sqrt(rij)
      SS=SIG/rij
      U=U+4.0d0*EPS*(SS**12-SS**6)
    endif
  enddo
enddo
return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DRSCALE
use MCPARAM,only : drmax
use VARIABLES,only : adopt,reject
implicit none
real8 :: ans

ans=dble(adopt)/dble(adopt+reject)
if(ans > 0.55d0) then
  drmax=drmax*1.05d0
else if(ans < 0.45d0) then
  drmax=drmax*0.95d0
endif
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
