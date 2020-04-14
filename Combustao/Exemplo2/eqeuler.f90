!Domain Divition for different equations 
!______________________________________________________
!             |                         |              |
!   NRBC      |         NRBC            |     NRBC     |   
!  streching  |       streching         |  streching   |
!  corner  1  |       xlayer_top        |   corner  3  | 
!_____________|_________________________|______________|
!             |                         |              |
!   NRBC      |                         |    NRBC      |
! streching   |                         |  streching   |   
! ylayer_left |       NS equations      | ylayer_right |   
!             |                         |              |   
!             |                         |              |
!             |                         |              |
!_____________|_________________________|______________|
!   NRBC      |                         |    NRBC      |
! streching   |         NRBC            | streaching   |   
!   corner 2  |      xlayer_bottom      |   corner 4   |
!_____________|_________________________|______________|

module  eqeuler 

use global
USE DIFF
!USE MPIOWN
contains

subroutine tau_ij()!its sended by the global.
IMPLICIT NONE
INTEGER::I,J
!Claudio REAL(KIND=IP),DIMENSION(1:imax,1:jmax)::diver
!TAU(1,...) = TAU_XX, // TAU(2,...) = TAU_XY // TAU(3,...) = TAU_YY,
!TO CALCULATE THE TENSOR TAU(IJ)
! !$omp parallel do private(i)

! !$acc data copyin(mi(1:imax,1:jmax),dudm(2:3,1:imax,1:jmax),dudn(2:3,1:imax,1:jmax)) copyout(tau(1:3,1:imax,1:jmax))
!$acc parallel loop collapse(2)
do j = 1,jmax
  do i = 1,imax
            !Claudio diver(i,j) =  - 2.0D0/3.0D0*(MI(I,J))*(DUDM(2,I,J)+DUDN(3,I,J))      
            !TAU(1,I,J) =    diver(i,j)+2.0d0*(mi(i,j))*(dudm(2,i,j))
            !TAU(2,I,J) =   (MI(I,J))*(DUDN(2,I,J)+DUDM(3,I,J))                 
            !TAU(3,I,J) =    diver(i,j)+2.0d0*(mi(i,j))*(dudn(3,i,j))             
            TAU(1,I,J) =   (- 2.0D0/3.0D0*(MI(I,J))*(DUDM(2,I,J)+DUDN(3,I,J)))+2.0d0*(mi(i,j))*(dudm(2,i,j))
            TAU(2,I,J) =   (MI(I,J))*(DUDN(2,I,J)+DUDM(3,I,J))                 
            TAU(3,I,J) =   (- 2.0D0/3.0D0*(MI(I,J))*(DUDM(2,I,J)+DUDN(3,I,J)))+2.0d0*(mi(i,j))*(dudn(3,i,j))             
  end do
end do
! !$omp end parallel do

!!!! Stress tensor derivates x direction  
  call deropm(DTAUDM,TAU,1,imax,1,jmax,size(tau,1))
!!!! Stress tensor derivates y direction 
  call deropn(DTAUDN,TAU,1,imax,1,jmax,size(tau,1))

end subroutine

!#####################
!#####################
!#####################

SUBROUTINE EULERNONLINEAR(DUDT,U,NK)

IMPLICIT NONE                     
INTEGER::K,L,I,J,NK
REAL(KIND=IP),DIMENSION(1:nk,1:imax,1:jmax)::U,DUDT

!LOOP TO CALCULATE THE RHS OF THE EULER EQUATIONS (D/DT=RHS) WITH AND WITHOUT PML 

    ! !$omp parallel do private(I)
!    !$acc data copyin(DUDM,DUDN,A,B,Re,U,DTAUDM,DTAUDN,s_acou,gamma,PR,D2TEMPDM,D2TEMPDN) copyout(DUDT)
    !$acc parallel loop collapse(2)
    DO J=1,jmax
      DO I=1,imax
               DUDT(1,I,J) = - (  A(1,1,I,J)*(DUDM(1,I,J))          &
                                + A(1,2,I,J)*(DUDM(2,I,J))          &
                                + B(1,1,I,J)*(DUDN(1,I,J))          &
                                + B(1,3,I,J)*(DUDN(3,I,J))          )  

               DUDT(2,I,J) = - (  A(2,2,I,J)*(DUDM(2,I,J))          &
                                + A(2,4,I,J)*(DUDM(4,I,J))          &
                                + B(2,2,I,J)*(DUDN(2,I,J))          &
                                - 1.d0/Re/(U(1,I,J))*               &
                                 (DTAUDM(1,I,J) + DTAUDN(2,I,J))    )

               DUDT(3,I,J) = - (  A(3,3,I,J)*(DUDM(3,I,J))          &
                                + B(3,3,I,J)*(DUDN(3,I,J))          &
                                + B(3,4,I,J)*(DUDN(4,I,J))          &
                                - 1.d0/Re/(U(1,I,J))*               &
                                 (DTAUDM(2,I,J) + DTAUDN(3,I,J))    )            

               DUDT(4,I,J) = - (  A(4,2,I,J)*(DUDM(2,I,J))          &
                                + A(4,4,I,J)*(DUDM(4,I,J))          &
                                + B(4,3,I,J)*(DUDN(3,I,J))          &
                                + B(4,4,I,J)*(DUDN(4,I,J))          &
                                - s_acou(i,j)                       &
                                - (gamma-1.0D0)/Re * (              &
                                + (DUDM(2,I,J))*TAU(1,I,J)          &
                                + (DUDM(3,I,J))*TAU(2,I,J)          &
                                + (DUDN(2,I,J))*TAU(2,I,J)          &
                                + (DUDN(3,I,J))*TAU(3,I,J)   )      & 
				- 1.0D0/RE/PR * (                   &
                                  D2TEMPDM(1,i,j) +                   &
                                  D2TEMPDN(1,i,j)             )       )    

               DUDT(5,I,J) = - (U(2,i,j)*DUDM(5,i,j)+U(3,i,j)*DUDN(5,i,j)&
			     -1.0d0/Pe/U(1,i,j)*(d2psiFdm(1,i,j)+d2psiFdn(1,i,j)) )

       end do
    end do
   ! !$omp end parallel do
end subroutine 

end module
