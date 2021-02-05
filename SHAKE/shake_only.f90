!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*
*   calculation H2O dynamics 1: Oxygen 2,3:Hydrogen
*
*   using shake 
*
*   NVE
*
*   
*
*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NEXT
use VARIABLES,only : R,ROLD,RNEW
implicit none
ROLD=R
R=RNEW
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SHAKE  !including verler 
use MDPARAM,only : N,NSITE,dt,W,RBOX,Dij
use VARIABLES,only : R,RNEW,ROLD,V,T
implicit none
integer :: i,j,i2,i3,k
real8,dimension(NSITE) :: rdashij !no constrain
real8,dimension(NSITE) :: rij
real8 :: rdashijSq,tmp_inner,diff_sq
real8 :: mass_i,mass_j !inverse
!! shake
real8 :: GAM
!! iteration judge
logical :: iteflg
logical,dimension(NSITE)  :: moving,moved
integer :: cnt,max_cnt
real8 :: er0=1.d-6  ! tolerance 
T=0.d0
max_cnt=30

do i=1,N 
  !verlet   only interforce
  do i2=1,NSITE
    do k=1,3  
      RNEW(k,i2,i)=2.d0*R(k,i2,i)-ROLD(k,i2,i)+dt**2*F(k,i2,i)/W(i2)
    enddo
  moving(i2)=.false.
  moved(i2)=.true.
  enddo
  !end verlet
  iteflg = .false.
  cnt=0
  !iteration method
  do while(.NOT. iteflg)
    iteflg=.true.
    do i2=1,NSITE     ! corretion position
      !i3=mod(i2,NSITE)+1
      i3=i2+1
      if(i3> nsite) i3=1
      if((moved(i2) .or. moved(i3)) ) then !cycle ! judge which changed position
        do k=1,3 !xyz
          rdashij(k)=RNEW(k,i2,i)-RNEW(k,i3,i) !拘束前の距離
          rdashij(k)=rdashij(k)-anint(rdashij(k)/RBOX)*RBOX
        enddo
        rdashijSq=dot_product(rdashij,rdashij)
        diff_sq=rdashijsq-Dij(i2)**2 ! 拘束条件との差
        if( abs(diff_sq) >  er0*(Dij(i2))**2 )then! cycle   ! judge satisfied constrain
          do k=1,3
            rij(k)=R(k,i2,i)-R(k,i3,i)
            rij(k)=rij(k)-anint(rij(k)/RBOX)*RBOX
          enddo
          tmp_inner=dot_product(rij,rdashij)
          mass_i=1.d0/W(i2)
          mass_j=1.d0/W(i3)
          GAM=diff_sq/(2.d0*(mass_i+mass_j)*tmp_inner)   ! ignore the term of square
          do k=1,3
            RNEW(k,i2,i)=RNEW(k,i2,i)-GAM*mass_i*rij(k)
            RNEW(k,i3,i)=RNEW(k,i3,i)+GAM*mass_j*rij(k)
          enddo
          moving(i2)=.true.  ! because changed position by constrain
          moving(i3)=.true.
          iteflg=.false.
        endif
      endif
    enddo
    ! end corretion position
    do i2=1,NSITE  ! set up  next iteration
      moved(i2)=moving(i2)
      moving(i2)=.false.
    enddo
    cnt=cnt+1
    if(cnt>max_cnt) then
      write(*,*) " cnt = ",cnt," is too many iteration ! I'm tired!!"
      stop
    endif
  enddo ! end loop

  !!velocity cal
  do i2=1,NSITE
    do k=1,3
      V(k,i2,i)=(RNEW(k,i2,i)-ROLD(k,i2,i))/(2.d0*DT)
      T=T+(V(k,i2,i)**2)*W(i2)*0.5d0
    enddo
  enddo

enddo  

 ! end shake
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

