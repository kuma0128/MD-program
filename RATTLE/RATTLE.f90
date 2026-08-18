!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 Calculation H2O dynamics 

 using RATTLE 

 Listen to

 https://youtu.be/8B7xr_EjbzE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rattle_moveA  !including moveA 
use MDPARAM
use VARIABLES
implicit none
integer :: i,j,i2,i3,k
real8,dimension(NSITE) :: rdashij !no constrain
real8,dimension(NSITE) :: rij
real8 :: rdashijSq,tmp_inner,diff_sq
real8 :: mass_i,mass_j !inverse
real8 :: dt2
real8,dimension(3,NSITE) :: Qnew  ! Rnew 
!! shake
real8 :: GAM
!! iteration judge
logical :: iteflg
logical,dimension(NSITE)  :: moving,moved
integer :: cnt,max_cnt
real8 :: er0=1.d-5  ! tolerance

dt2=0.5d0*dt**2
max_cnt=100

do i=1,N ! in  a molecure 
  !moveA   onlu interforce
  do i2=1,NSITE
    do k=1,3  
      Qnew(k,i2)=r(k,i2,i)+v(k,i2,i)*dt+dt2*f(k,i2,i)/W(i2)
      v(k,i2,i)=v(k,i2,i)+0.5d0*dt*F(k,i2,i)/W(i2)
    enddo
    moving(i2)=.false.
    moved(i2)=.true.
  enddo
  !finish moveA
  iteflg = .false.
  cnt=0
  !iteration method
  do while(.NOT. iteflg)
    iteflg=.true.
    do i2=1,NSITE     ! corretion position
      i3=mod(i2,NSITE)+1
      if((moved(i2) .or. moved(i3)) ) then !cycle ! judge which changed position
        do k=1,3 !xyz
          rdashij(k)=QNEW(k,i2)-QNEW(k,i3) 
          rdashij(k)=rdashij(k)-anint(rdashij(k)/RBOX)*RBOX
        enddo
        rdashijSq=dot_product(rdashij,rdashij)
        diff_sq=rdashijsq-Dij(i2)**2 
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
            Qnew(k,i2)=Qnew(k,i2)-GAM*mass_i*rij(k)*0.5d0
            Qnew(k,i3)=Qnew(k,i3)+GAM*mass_j*rij(k)*0.5d0
            V(k,i2,i)=V(k,i2,i)-0.5d0*mass_i*GAM*rij(k)/dt
            V(k,i3,i)=V(k,i3,i)+0.5d0*mass_j*GAM*rij(k)/dt
          enddo
          moving(i2)=.true.  ! because changed position by constrain
          moving(i3)=.true.
          iteflg=.false.
        endif
      endif
    enddo
    ! finish corretion position
    do i2=1,NSITE  ! set up  next iteration
      moved(i2)=moving(i2)
      moving(i2)=.false.
    enddo
    cnt=cnt+1
    if(cnt>max_cnt) then
      write(*,*) " cnt = ",cnt," is too many iteration ! I'm tired!!"
      stop
    endif
  enddo
  ! finish loop
!  write(*,*) "cnt=",cnt 
  !!R update
  do i2=1,NSITE
    do k=1,3
      R(k,i2,i)=Qnew(k,i2)
    enddo
  enddo
  !!

enddo  
!finish rattle_moveA
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rattle_moveB  !including moveB 
use MDPARAM
use VARIABLES
implicit none
integer :: i,j,i2,i3,k
real8,dimension(NSITE) :: rdashij !no constrain
real8,dimension(NSITE) :: rij,vij 
real8 :: rdashijSq,tmp_inner,diff_sq
real8 :: mass_i,mass_j !inverse
real8,dimension(3,NSITE) :: Qnew  ! Rnew 
!! rattle moveB
real8 :: GAM2
!! iteration judge
logical :: iteflg
logical,dimension(NSITE)  :: moving,moved
integer :: cnt,max_cnt
real8 :: er0=1.d-5  !tolerance 

T=0.d0
max_cnt=100

do i=1,N ! in  a molecure 
  !moveB   only constrain force
  do i2=1,NSITE
    do k=1,3  
      v(k,i2,i)=v(k,i2,i)+0.5d0*dt*F(k,i2,i)/W(i2)
    enddo
    moving(i2)=.false.
    moved(i2)=.true.
  enddo
  !end moveB
  iteflg = .false.
  cnt=0
  !iteration method
  do while(.NOT. iteflg)
    iteflg=.true.
    do i2=1,NSITE     ! corretion position
      i3=mod(i2,NSITE)+1
      if((moved(i2) .or. moved(i3)) ) then !cycle ! judge which changed position
        do k=1,3 !xyz
          rdashij(k)=R(k,i2,i)-R(k,i3,i) 
          rdashij(k)=rdashij(k)-anint(rdashij(k)/RBOX)*RBOX
        enddo
        rdashijSq=dot_product(rdashij,rdashij)
        diff_sq=rdashijsq-Dij(i2)**2 
        if( abs(diff_sq) >  er0*(Dij(i2))**2 )then! cycle   ! judge satisfied constrain
          do k=1,3
            vij(k)=V(k,i2,i)-V(k,i3,i)
          enddo
          tmp_inner=dot_product(rij,vij)
          mass_i=1.d0/W(i2)
          mass_j=1.d0/W(i3)
          GAM2=tmp_inner/((mass_i+mass_j)*Dij(i2)**2)   ! ignore the term of square
          do k=1,3
            V(k,i2,i)=V(k,i2,i)-mass_i*GAM2*rij(k)
            V(k,i3,i)=V(k,i3,i)+mass_j*GAM2*rij(k)
          enddo
          moving(i2)=.true.  ! because changed position by constrain
          moving(i3)=.true.
          iteflg=.false.
        endif
      endif
    enddo
    ! finish corretion position
    do i2=1,NSITE  ! set up  next iteration
      moved(i2)=moving(i2)
      moving(i2)=.false.
    enddo
    cnt=cnt+1
    if(cnt>max_cnt) then
      write(*,*) " cnt = ",cnt," is too many iteration ! I'm tired!!"
      stop
    endif
  enddo
  ! finish loop
  !!T cal
  do i2=1,NSITE
    do k=1,3
      T=T+V(k,i2,i)**2*0.5d0*W(i2)
    enddo
  enddo
  !!

enddo  
!finish rattle_moveB
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
