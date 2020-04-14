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
USE MPIOWN
contains

subroutine tau_ij()!its sended by the global.
IMPLICIT NONE
INTEGER::I,J
!TAU(1,...) = TAU_XX, // TAU(2,...) = TAU_XY // TAU(3,...) = TAU_YY,
!TO CALCULATE THE TENSOR TAU(IJ)
do j = s_s, e_s
  do i = i_s, d_s
      IF( (I>=D+1).AND.(I<=IMAX-D))THEN
        IF((J>=D+1).AND.(J<=JMAX-D))THEN
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         * ((DUDN(2,I,J)+(DUDNB(2,I,J))                         &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))
        END IF
      END IF
    !***********************ylayer_left***********************************
       if((i>=D+1).and.(i<=imax-D))then
          if((j>=1).and.(j<=D))then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))

            TAU(2,I,J) =  (MI(I,J)+MIb(I,J))                                    &
                         *((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                  &
                         + (DUDM(3,I,J)+DUDMB(3,I,J)))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
          end if
       end if

    !***********************ylayer_right***********************************
         if((i>=D+1).and.(i<=imax-D))then
            if((j>=jmax-D+1).and.(j<=jmax))then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         *((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                  &
                         + (DUDM(3,I,J)+DUDMB(3,I,J)))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J)))                          &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
            end if
         end if
    !***********************xlayer_bottom*********************************
         if ((i>=1).and.(i<=D))then
            if ((j>=D+1).and.(j<=jmax-D)) then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))                             &
                         * (dudm(2,i,j)+dudmB(2,i,j))*meshx(i)

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         *((DUDN(2,I,J)+DUDNB(2,I,J))                           &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))
           end if
         end if
    !***********************xlayer_top*********************************
        if((i>=imax-D+1).and.(i<=imax))then
            if((j>=D+1).and.(j<=jmax-D))then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))                             &
                         * (dudm(2,i,j)+dudmB(2,i,j))*meshx(i)

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         *((DUDN(2,I,J)+DUDNB(2,I,J))                           &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))                          &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))
           end if
         end if
    !***********************corners*********************************
    !***********************1*********************************
         if ((i>=1).and.(i<=D))then
            if((j>=jmax-D+1).and.(j<=jmax))then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))*meshx(i)


            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         * ((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                 &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
           end if
         end if
     !***********************2*********************************
         if ((i>=1).and.(i<=D)) then
            if ((j>=1).and.(j<=D)) then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))*meshx(i)

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         * ((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                 &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
           end if
         end if
!***********************3*********************************
         if ((i>=imax-D+1).and.(i<=imax)) then
            if ((j>=1).and.(j<=D)) then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))*meshx(i)

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         * ((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                 &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
            end if
         end if
!**********************4*********************************
         if ((i>=imax-D+1).and.(i<=imax)) then
            if ((j>=jmax-D+1).and.(j<=jmax)) then
            TAU(1,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudm(2,i,j)+dudmB(2,i,j))*meshx(i)

            TAU(2,I,J) =   (MI(I,J)+MIb(I,J))                                   &
                         * ((DUDN(2,I,J)+DUDNB(2,I,J))*meshy(j)                 &
                         + (DUDM(3,I,J)+DUDMB(3,I,J))*meshx(i))

            TAU(3,I,J) = - 2.0D0/3.0D0*(MI(I,J)+MIb(I,J))                       &
                         * ((DUDN(3,I,J)+DUDNB(3,I,J))*meshy(j)                 &
                         + (DUDM(2,I,J)+DUDMB(2,I,J))*meshx(i))                 &
                         + 2.0d0*(mi(i,j)+mib(i,j))*(dudn(3,i,j)+dudnb(3,i,j))*meshy(j)
           end if
         end if
!********************************************************
  end do
end do
!To make the comunication of the
!gosth point between the process
!location:mpiown.f90

call exchange_bandtau_x()
call exchange_bandtau_y()
end subroutine

!#####################
!#####################
!#####################

SUBROUTINE EULERNONLINEAR(DUDT,U)

IMPLICIT NONE                     
INTEGER::K,L,I,J 
REAL(KIND=IP),DIMENSION(1:4,I_S-G_P:D_S+G_P,S_S-G_P:E_S+G_P)::U,DUDT

!LOOP TO CALCULATE THE RHS OF THE EULER EQUATIONS (D/DT=RHS) WITH AND WITHOUT PML 

    DO J=S_S,E_S
      DO I=I_S,D_S
   !***IF TO CALCULATE THE INTERIOR DOMAIN, EULER EQUATION WITHOUT PML 
   !*********S(I,J) IS THE ACUSTIC FONT TO DEVELOPMENT THE INSTABILILITY IN DETERMINED POSITION AND IS NEGATIVE IN 
   !**********THE LHS
        IF( (I>=D+1).AND.(I<=IMAX-D))THEN
            IF((J>=D+1).AND.(J<=JMAX-D))THEN
               DUDT(1,I,J) = - (  A(1,1,I,J)*(DUDM(1,I,J)+DUDMB(1,I,J))          &
                                + A(1,2,I,J)*(DUDM(2,I,J)+DUDMb(2,I,J))          &
                                + B(1,1,I,J)*(DUDN(1,I,J)+DUDNb(1,I,J))          &
                                + B(1,3,I,J)*(DUDN(3,I,J)+DUDNb(3,I,J))          )  

               DUDT(2,I,J) = - (  A(2,2,I,J)*(DUDM(2,I,J)+DUDMb(2,I,J))          &
                                + A(2,4,I,J)*(DUDM(4,I,J)+DUDMB(4,I,J))          &
                                + B(2,2,I,J)*(DUDN(2,I,J)+DUDNb(2,I,J))          &
                                - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                   &
                                * (DTAUDM(1,I,J) + DTAUDN(2,I,J))                )

               DUDT(3,I,J) = - (  A(3,3,I,J)*(DUDM(3,I,J)+DUDMb(3,I,J))          &
                                + B(3,3,I,J)*(DUDN(3,I,J)+DUDNb(3,I,J))          &
                                + B(3,4,I,J)*(DUDN(4,I,J)+DUDNB(4,I,J))          &
                                - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                   &
                                * (DTAUDM(2,I,J) + DTAUDN(3,I,J))                )            

               DUDT(4,I,J) = - (  A(4,2,I,J)*(DUDM(2,I,J)+DUDMb(2,I,J))          &
                                + A(4,4,I,J)*(DUDM(4,I,J)+DUDMB(4,I,J))          &
                                + B(4,3,I,J)*(DUDN(3,I,J)+DUDNb(3,I,J))          &
                                + B(4,4,I,J)*(DUDN(4,I,J)+DUDNB(4,I,J))          &
                                - S_ACOU(I,J)                                    &
                                - (gamma-1.0D0)/Re/(U(1,I,J)+UB(1,I,J))*  (      &
				                        + (DUDM(2,I,J)+DUDMb(2,I,J))*TAU(1,I,J)          &
                                + (DUDM(3,I,J)+DUDMb(3,I,J))*TAU(2,I,J)          &
                                + (DUDN(2,I,J)+DUDNb(2,I,J))*TAU(2,I,J)          &
                                + (DUDM(3,I,J)+DUDMb(3,I,J))*TAU(3,I,J)          &
				                        - 1.0D0/RE/PR * (                                &
                                  D2TEMPDM(i,j) +                                &
                                  D2TEMPDN(i,j)             )                   ))    
            END IF
        END IF

!***If to calculate the pml domain, euler equation with PML 
    
!***********************ylayer_left***********************************
    
         if((i>=D+1).and.(i<=imax-D))then
            if((j>=1).and.(j<=D))then 
             call strechingy(U,dUdt,i,j)
            end if 
         end if 
    
!***********************ylayer_right***********************************
    
         if((i>=D+1).and.(i<=imax-D))then
            if((j>=jmax-D+1).and.(j<=jmax))then
             call strechingy(U,dUdt,i,j)
            end if 
         end if 
    
!***********************xlayer_bottom********************************* 
    
         if ((i>=1).and.(i<=D))then
            if ((j>=D+1).and.(j<=jmax-D)) then
             call strechingx(U,dUdt,i,j)
           end if 
         end if 
    
!***********************xlayer_top********************************* 
    
         if((i>=imax-D+1).and.(i<=imax))then
            if((j>=D+1).and.(j<=jmax-D))then
             call strechingx(U,dUdt,i,j)
           end if 
         end if 
    
!***********************corners********************************* 
    
!***********************1********************************* 
    
         if ((i>=1).and.(i<=D))then
            if((j>=jmax-D+1).and.(j<=jmax))then
             call strechingxy(U,dUdt,i,j)
           end if 
         end if 
    
!***********************2********************************* 
    
         if ((i>=1).and.(i<=D)) then
            if ((j>=1).and.(j<=D)) then
             call strechingxy(U,dUdt,i,j)
           end if 
         end if 
        
    
!***********************3********************************* 
    
    
         if ((i>=imax-D+1).and.(i<=imax)) then
            if ((j>=1).and.(j<=D)) then
             call strechingxy(U,dUdt,i,j)
           end if 
         end if 
    
!**********************4********************************* 
    
         if ((i>=imax-D+1).and.(i<=imax)) then
            if ((j>=jmax-D+1).and.(j<=jmax)) then
             call strechingxy(U,dUdt,i,j)
           end if 
         end if 
    
!******************************************************** 
       end do
    end do
   
end subroutine 

subroutine strechingx(u,dE,i,j)

implicit none                     
integer::k,l,i,j 
real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U,ub,dE

dE(1,i,j) = -(  A(1,1,i,j)*(dUdm(1,i,j)+dUdmB(1,i,j))*meshx(i)      &
              + A(1,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)      &
              + B(1,1,i,j)*(dUdn(1,i,j)+dUdnb(1,i,j))               & 
              + B(1,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))               )
!AQUI
dE(2,i,j) = -( A(2,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)       & 
             + A(2,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))*meshx(i)       &
             + B(2,2,i,j)*(dUdn(2,i,j)+dUdnb(2,i,j))                &
             - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                         &
             * (DTAUDM(1,I,J)*meshx(i) + DTAUDN(2,I,J))             )

dE(3,i,j) = -( A(3,3,i,j)*(dUdm(3,i,j)+dUdmb(3,i,j))*meshx(i)       &
             + B(3,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))                &
             + B(3,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))                &
             - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                         &
             * (DTAUDM(2,I,J)*meshx(i) + DTAUDN(3,i,j))             )

dE(4,i,j) = -( A(4,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)       &
             + A(4,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))*meshx(i)       &
             + B(4,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))                &
             + B(4,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))                &
             - s_acou(i,j)                                          &
             - (gamma-1.0D0)/Re/(U(1,I,J)+UB(1,I,J)) *         (    &
             + (DUDM(2,I,J)+DUDMb(2,I,J))*TAU(1,I,J)*meshx(i) 	    &
             + (DUDM(3,I,J)+DUDMb(3,I,J))*TAU(2,I,J)*meshx(i)       &
             + (DUDN(2,I,J)+DUDNb(2,I,J))*TAU(2,I,J)                &
             + (DUDN(3,I,J)+DUDNb(3,I,J))*TAU(3,I,J)                & 
				     - 1.0D0/RE/PR * (                                      &
               D2TEMPDM(i,j)*meshx(i)+D2TEMPDN(i,j) )              ))  

end subroutine 

subroutine strechingy(u,dE,i,j)

implicit none                     
integer::k,l,i,j 
real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U,ub,dE



dE(1,i,j) = -( A(1,1,i,j)*(dUdm(1,i,j)+dUdmB(1,i,j))                 &
             + A(1,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))                 &
             + B(1,1,i,j)*(dUdn(1,i,j)+dUdnb(1,i,j))*meshy(j)        &
             + B(1,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)        )
       
dE(2,i,j) = -( A(2,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))                 &
             + A(2,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))                 &
             + B(2,2,i,j)*(dUdn(2,i,j)+dUdnb(2,i,j))*meshy(j)        &
             - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                          &
             * (DTAUDM(1,i,j)+DTAUDN(2,i,j)*meshy(j))                )

dE(3,i,j) = -( A(3,3,i,j)*(dUdm(3,i,j)+dUdmb(3,i,j))                 &
             + B(3,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)        &
             + B(3,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))*meshy(j)        & 
             - 1.d0/Re/(U(1,I,J)+UB(1,I,J))                          &
             * (DTAUDM(2,i,j)+DTAUDN(3,i,j)*meshy(j))                )

dE(4,i,j) = -( A(4,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))                 &
             + A(4,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))                 &
             + B(4,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)        &
             + B(4,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))*meshy(j)        &
             - s_acou(i,j)                                           &
             - (gamma-1.0D0)/Re/(U(1,I,J)+UB(1,I,J)) * (             &
             + (DUDM(2,I,J)+DUDMb(2,I,J))*TAU(1,I,J)                 & 
             + (DUDM(3,I,J)+DUDMb(3,I,J))*TAU(2,I,J)                 &
             + (DUDN(2,I,J)+DUDNb(2,I,J))*TAU(2,I,J)*meshy(j)        &
             + (DUDN(3,I,J)+DUDNb(3,I,J))*TAU(3,I,J)*meshy(j)        &
				     - 1.0D0/RE/PR * (                                       &
               D2TEMPDM(i,j)+D2TEMPDN(i,j)*meshy(j) )               ))     

end subroutine 

subroutine strechingxy(u,dE,i,j)

implicit none                     
integer::k,l,i,j 
real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U,ub,dE

dE(1,i,j) = - (  A(1,1,i,j)*(dUdm(1,i,j)+dUdmB(1,i,j))*meshx(i)       &
               + A(1,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)       &
               + B(1,1,i,j)*(dUdn(1,i,j)+dUdnb(1,i,j))*meshy(j)       &
               + B(1,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)       )

dE(2,i,j) = - (  A(2,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)       &
               + A(2,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))*meshx(i)       &
               + B(2,2,i,j)*(dUdn(2,i,j)+dUdnb(2,i,j))*meshy(j)       &
               - 1.d0/Re/(U(1,I,J)+UB(1,I,J))*(DTAUDM(1,i,j)          &
               * meshx(i) + DTAUDN(2,i,j)*meshy(j))                   )
  
dE(3,i,j) = - (  A(3,3,i,j)*(dUdm(3,i,j)+dUdmb(3,i,j))*meshx(i)       &
               + B(3,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)       &
               + B(3,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))*meshy(j)       &
               - 1.d0/Re/(U(1,I,J)+UB(1,I,J))*(DTAUDM(2,i,j)          &
               * meshx(i) + DTAUDN(3,i,j)*meshy(j))                   )

dE(4,i,j) = - (  A(4,2,i,j)*(dUdm(2,i,j)+dUdmb(2,i,j))*meshx(i)       &
               + A(4,4,i,j)*(dUdm(4,i,j)+dUdmB(4,i,j))*meshx(i)       &
               + B(4,3,i,j)*(dUdn(3,i,j)+dUdnb(3,i,j))*meshy(j)       &
               + B(4,4,i,j)*(dUdn(4,i,j)+dUdnB(4,i,j))*meshy(j)       &
               - s_acou(i,j)                           &
               - (gamma-1.0D0)/Re/(U(1,I,J)+UB(1,I,J)) * (            &
               + (DUDM(2,I,J)+DUDMb(2,I,J))*TAU(1,I,J)*meshx(i)       &
               + (DUDM(3,I,J)+DUDMb(3,I,J))*TAU(2,I,J)*meshx(i)       &
               + (DUDN(2,I,J)+DUDNb(2,I,J))*TAU(2,I,J)*meshy(j)       &
               + (DUDN(3,I,J)+DUDNb(3,I,J))*TAU(3,I,J)*meshy(j)       &
				       - 1.0D0/RE/PR * (                                      &
                 D2TEMPDM(i,j)*meshx(i)+D2TEMPDN(i,j)*meshy(j) )     ))  


end subroutine 

end module
