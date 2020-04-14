module rhs
use global
use derivs
use eqeuler
!use mpi 
!USE MPIOWN
USE EMPIRICAL_COEFFICIENTS
contains   

subroutine rhs_euler(dUdt,U,nk)

   implicit none
     
   integer::k,i,j,nk
   real(kind=ip),dimension(1:nk,1:imax, 1:jmax)::dUdt,U
   real(kind=ip),dimension(1,1:imax, 1:jmax)::T_inst

!...Escomento instantaneo,varia com o tempo 
    A(1,1,1:imax, 1:jmax)  = U(2,1:imax, 1:jmax) !+ Ub(2,1:imax, 1:jmax)
    A(1,2,1:imax, 1:jmax)  = U(1,1:imax, 1:jmax) !+ Ub(1,1:imax, 1:jmax)
    !A(1,3,1:imax, 1:jmax)  = 0.0d0
    !A(1,4,1:imax, 1:jmax)  = 0.0d0
    
    !A(2,1,1:imax, 1:jmax)  = 0.d0 
    A(2,2,1:imax, 1:jmax)  = U(2,1:imax, 1:jmax) !+ Ub(2,1:imax, 1:jmax)
    !A(2,3,1:imax, 1:jmax)  = 0.d0
    A(2,4,1:imax, 1:jmax)  = 1.d0/(U(1,1:imax, 1:jmax))!+Ub(1,1:imax, 1:jmax))
    
    !A(3,1,1:imax, 1:jmax)  = 0.d0 
    !A(3,2,1:imax, 1:jmax)  = 0.d0
    A(3,3,1:imax, 1:jmax)  = U(2,1:imax, 1:jmax) !+ Ub(2,1:imax, 1:jmax)
    !A(3,4,1:imax, 1:jmax)  = 0.d0
    
    !A(4,1,1:imax, 1:jmax)  = 0.d0 
    A(4,2,1:imax, 1:jmax)  = gamma*U(4,1:imax, 1:jmax) !+ gamma*Ub(4,1:imax, 1:jmax)
    !A(4,3,1:imax, 1:jmax)  = 0.d0
    A(4,4,1:imax, 1:jmax)  = U(2,1:imax, 1:jmax) !+ Ub(2,1:imax, 1:jmax)
    
    
    B(1,1,1:imax, 1:jmax)  = U(3,1:imax, 1:jmax) !+ Ub(3,1:imax, 1:jmax)
    !B(1,2,1:imax, 1:jmax)  = 0.d0
    B(1,3,1:imax, 1:jmax)  = U(1,1:imax, 1:jmax) !+ Ub(1,1:imax, 1:jmax)
    !B(1,4,1:imax, 1:jmax)  = 0.d0
    
    !B(2,1,1:imax, 1:jmax)  = 0.d0 
    B(2,2,1:imax, 1:jmax)  = U(3,1:imax, 1:jmax) !+ Ub(3,1:imax, 1:jmax)
    !B(2,3,1:imax, 1:jmax)  = 0.d0
    !B(2,4,1:imax, 1:jmax)  = 0.d0
    
    !B(3,1,1:imax, 1:jmax)  = 0.d0 
    !B(3,2,1:imax, 1:jmax)  = 0.d0
    B(3,3,1:imax, 1:jmax)  = U(3,1:imax, 1:jmax) !+ Ub(3,1:imax, 1:jmax)
    B(3,4,1:imax, 1:jmax)  = 1.d0/(U(1,1:imax, 1:jmax))! +Ub(1,1:imax, 1:jmax))
    
    !B(4,1,1:imax, 1:jmax)  = 0.d0 
    !B(4,2,1:imax, 1:jmax)  = 0.d0
    B(4,3,1:imax, 1:jmax)  = gamma*U(4,1:imax, 1:jmax) !+ gamma*Ub(4,1:imax, 1:jmax)
    B(4,4,1:imax, 1:jmax)  = U(3,1:imax, 1:jmax) !+ Ub(3,1:imax, 1:jmax)

    !T = gamma(P+Pbase)/(Rho+Rhobase)-Tbase  
    !gammaB = U(1,1:imax, 1:jmax)/U(4,1:imax, 1:jmax)

    !TEMP(1:imax, 1:jmax)    = (gamma)*(U(4,1:imax, 1:jmax)+Ub(4,1:imax, 1:jmax))/ &
    !(U(1,1:imax, 1:jmax)+Ub(1,1:imax, 1:jmax))  -TEMPB(1:imax, 1:jmax) 

    TEMP(1,1:imax, 1:jmax) = (gamma)*(U(4,1:imax, 1:jmax))/(U(1,1:imax, 1:jmax)) 
    
    T_inst(:,:,:)= Temp(:,:,:)

    CALL MOLECULAR_MI(T_inst)
    CALL K_COEFFICIENT(T_inst)


!.........Derivadas das variaveis usadas no rhs da equa��o Euler,
!.........O Vetor de derivadas foi passado usando o modulo global
!.........location: derivs.f90 .......

    CALL deriv(U,nk) 
!...Create the tau varibles and make the exchanged with other process
!...Location:eqeuler.f90
!...Use the variable tau that it send by the global.

    CALL TAU_IJ()

!... RHS(right hand side) Equações de Euler with and without PML.
!... Using the variables U, necessary to calculate the rhs of euler equations
!... dUdt have the minus because this is the rhs of the euler equation dUdt=-deuler.
!... location: eqeuler.f90 .......

    call eulernonlinear(dUdt,U,size(U,1))

end subroutine



end module

