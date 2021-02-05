! Box-Muller法のGauss型乱数発生
module sub
  implicit none
contains
  subroutine BOX(X)
  real(8) :: x
  integer :: iset
  real(8) :: fac,gset,rsq,v1,v2
  save iset,gset
  data iset/ 0/
  if(iset==0)then
  1 call random_number(v1)
    call random_number(v2)
   v1=v1*2.0d0-1.0d0
   v2=v2*2.0d0-1.0d0
    rsq=v1**2+v2**2
    if(rsq<0.0d0 .or. rsq>1.0d0) goto 1
    fac=sqrt(-2.d00*log(rsq)/rsq)
    gset=v1*fac
    x=v2*fac
    iset=1
  else
    x=gset
    iset=0
  endif
  end subroutine BOX
end module sub

program main
  use sub
  implicit none
  integer,parameter :: n=1000000,Nbin=200
  real(8),parameter :: W=0.02d0
  integer,dimension(-Nbin:Nbin) :: L

  integer,dimension(2) :: ISEED=(/80,1001 /)
  integer :: i,j,fo=10

  real(8) :: x,ave,dis
  open(fo,file='output.d')
  L=0
  ave=0.0d0
  dis=0.0d0
  call random_seed(put=ISEED)
  do i=1,n
    call BOX(X)
    ! j=int(x/w)+1　切り捨て
    j=nint(x/w)　!四捨五入
    ! j=(x+w-1)/w  切り上げ
    if(j>nbin) j=nbin
    if(j<-nbin) j=-nbin
    L(j)=L(j)+1
    ave=ave+x
    dis=dis+x*x
  enddo
  ave=ave/dble(n)
  dis=dis/dble(n)-ave**2

  do  j=-nbin,nbin
    write(fo,*) dble(j)*w,L(J)
  enddo
  write(fo,*) '#ave,dis=',ave,dis
  close(fo)
end program main


