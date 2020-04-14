module diff
use global
use omp_lib
use openacc

contains
!.................................................................


subroutine deropm(dU,U,ii,fi,ij,fj,nk)

    !.......................................................................
    !     subrotina para calculo da derivada primeira de uma funcao u
    !     utilizando diferencas centradas de quarta ordem nos pontos
    !     internos do dominio, diferencas centradas de segunda ordem
    !     nos pontos vizinhos 'a fronteira e diferencas unilateral de
    !     segunda ordem nos pontos de fronteira.
    !.......................................................................
        implicit none
        
        integer  ::   k, i, j, ii, fi, ij, fj, nk
        real(kind=ip),dimension(nk,1:imax,1:jmax)   ::   U,dU


! !$acc kernels
        do k=1,nk
              
            ! !$omp parallel do	private (i)

            do j=1,jmax
                i=1
                du(k,i,j) = dx6(u,k,i,j,1,7,1,nk)
                i=2
                du(k,i,j) = dx5(u,k,i,j,1,7,1,nk)
                i=3
                du(k,i,j) = dx4(u,k,i,j,1,7,1,nk)

            ! !$acc parallel loop
            do i=4,imax-3
                    du(k,i,j)  = dudx(u,k,i,j,nk)
		end do

                i=imax-2
                du(k,i,j) = dx4(u,k,i,j,7,1,-1,nk)
                i=imax-1
                du(k,i,j) = dx5(u,k,i,j,7,1,-1,nk)
                i=imax
                du(k,i,j) = dx6(u,k,i,j,7,1,-1,nk)
    !....   ...................................................................
            end do 
            !$omp end parallel do
    !....   ...................................................................
     end do 
     
     call metricsX(du,size(du,1))
end subroutine
! !$acc end kernels

!$acc routine seq
real(kind=ip) function dudx(u,k,i,j,nk)

    integer  ::  k,i,j,z,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
    real(kind=ip),dimension(3) :: a

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudx=0.d0

    do z=1,3
         dudx= dudx + (a(z)*(u(k,i+z,j)-u(k,i-z,j)))/dm
    end do

end function dudx

!$acc routine seq
real(kind=ip) function dx4(u,k,i,j,iz,fz,p,nk)

    integer  ::  k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dx4=0.d0
        
       do z=iz,fz,p
          dx4= dx4+a(z)*u(k,i+z-3,j)/dm
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dx4=0.d0
   
       do z=iz,fz,p
          dx4= dx4+a(z)*u(k,i+z-5,j)/dm
       end do 

    end if

end function dx4


!$acc routine seq
real(kind=ip) function dx5(u,k,i,j,iz,fz,p,nk)

  integer :: k,i,j,z,iz,fz,p,nk
  real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
  real(kind=ip),dimension(7) :: a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dx5=0.d0
        
       do z=iz,fz,p
          dx5= dx5+a(z)*u(k,i+z-2,j)/dm
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dx5=0.d0
   
       do z=iz,fz,p
          dx5= dx5+a(z)*u(k,i+z-6,j)/dm
       end do 

    end if

end function dx5


!$acc routine seq
real(kind=ip) function dx6(u,k,i,j,iz,fz,p,nk)

    integer  ::  k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
    real(kind=ip),dimension(7) :: a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dx6=0.d0
        
       do z=iz,fz,p
          dx6= dx6+a(z)*u(k,i+z-1,j)/dm
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dx6=0.d0
   
       do z=iz,fz,p
          dx6= dx6+a(z)*u(k,i+z-7,j)/dm
       end do  

    end if

end function dx6

!.......................................................................
!.......................................................................
!.................................................................

subroutine deropn(dU,U,ii,fi,ij,fj, nk)
      implicit none
      integer  ::   k, i, j, ii, fi, ij, fj, nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: dU, U

      do k=1,nk
    !....   ...................................................................
            !$omp parallel do private (j)
            do i=ii,fi
    !....   ...................................................................
                j=1
                du(k,i,j) = dy6(u,k,i,j,1,7,1,nk)
                j=2
                du(k,i,j) = dy5(u,k,i,j,1,7,1,nk)
                j=3
                du(k,i,j) = dy4(u,k,i,j,1,7,1,nk)

                do j=4,jmax-3 
                    du(k,i,j)  = dudy(u,k,i,j,nk)
                end do

                j=jmax-2
                du(k,i,j) = dy4(u,k,i,j,7,1,-1,nk)
                j=jmax-1
                du(k,i,j) = dy5(u,k,i,j,7,1,-1,nk)
                j=jmax
                du(k,i,j) = dy6(u,k,i,j,7,1,-1,nk)
    !....   ...................................................................
            end do 
            !$omp end parallel do
    !....   ...................................................................
     end do 


	call metricsY(du,size(du,1))

end subroutine


!$acc routine seq
real(kind=ip) function dudy(u,k,i,j,nk)

    integer :: k,i,j,z,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
    real(kind=ip),dimension(3) :: a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudy=0.d0

    do z=1,3
         dudy= dudy + (a(z)*(u(k,i,j+z)-u(k,i,j-z)))/dn
    end do

end function dudy


!$acc routine seq
real(kind=ip) function dy4(u,k,i,j,iz,fz,p,nk)

    integer  ::  k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dy4=0.d0
        
       do z=iz,fz,p
          dy4= dy4+a(z)*u(k,i,j+z-3)/dn
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dy4=0.d0
   
       do z=iz,fz,p
          dy4= dy4+a(z)*u(k,i,j+z-5)/dn
       end do 

    end if

end function dy4

!$acc routine seq
real(kind=ip) function dy5(u,k,i,j,iz,fz,p,nk)

    integer  ::  k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dy5=0.d0
        
       do z=iz,fz,p
          dy5= dy5+a(z)*u(k,i,j+z-2)/dn
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dy5=0.d0
   
       do z=iz,fz,p
          dy5= dy5+a(z)*u(k,i,j+z-6)/dn
       end do 

    end if

end function dy5

!$acc routine seq
real(kind=ip) function dy6(u,k,i,j,iz,fz,p,nk)

    integer  ::  k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dy6=0.d0
        
       do z=iz,fz,p
          dy6= dy6+a(z)*u(k,i,j+z-1)/dn
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dy6=0.d0
   
       do z=iz,fz,p
          dy6= dy6+a(z)*u(k,i,j+z-7)/dn
       end do 

    end if

end function dy6

!.......................................................................

!derivada em y de uma variavel 
subroutine deropnb(du,u,ij,fj)
!.......................................................................
      implicit none
      integer  ::   k, i, j, ii, fi, ij, fj
      real(kind=ip),dimension(1:jmax) :: u,du
!....   ...................................................................
               j=1
               du(j) = dy6b(u,j,1,7,1)
               j=2
               du(j) = dy5b(u,j,1,7,1)
               j=3
               du(j) = dy4b(u,j,1,7,1)
               do j=4,jmax-3      
                   du(j)  = dudyb(u,j)
               end do
               j=jmax-2
               du(j) = dy4b(u,j,7,1,-1)
               j=jmax-1
               du(j) = dy5b(u,j,7,1,-1)
               j=jmax
               du(j) = dy6b(u,j,7,1,-1)

end subroutine

real(kind=ip) function dudyb(u,j)

    integer  ::  k,i,j,z
    real(kind=ip),dimension(1:jmax) :: u
    real(kind=ip),dimension(3)  ::  a   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudyb=0.d0

    do z=1,3
         dudyb= dudyb + (a(z)*(u(j+z)-u(j-z)))/dn
    end do

end function dudyb

real(kind=ip) function dy4b(u,j,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(1:jmax) :: u
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dy4b=0.d0
        
       do z=iz,fz,p
          dy4b= dy4b+a(z)*u(j+z-3)/dn
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dy4b=0.d0
   
       do z=iz,fz,p
          dy4b= dy4b+a(z)*u(j+z-5)/dn
       end do 

    end if

end function dy4b

real(kind=ip) function dy5b(u,j,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(1:jmax) :: u
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dy5b=0.d0
        
       do z=iz,fz,p
          dy5b= dy5b+a(z)*u(j+z-2)/dn
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dy5b=0.d0
   
       do z=iz,fz,p
          dy5b= dy5b+a(z)*u(j+z-6)/dn
       end do 

    end if

end function dy5b

real(kind=ip) function dy6b(u,j,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(1:jmax)  ::   u
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dy6b=0.d0
        
       do z=iz,fz,p
          dy6b= dy6b+a(z)*u(j+z-1)/dn
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dy6b=0.d0
   
       do z=iz,fz,p
          dy6b= dy6b+a(z)*u(j+z-7)/dn
       end do 

    end if

end function dy6b
!.......................................................................
 
!DERIVADA EM X DE UMA VARIAVEL
SUBROUTINE DEROPMB(DU,U,II,FI)
!.......................................................................
      IMPLICIT NONE
      INTEGER  ::   K, I, J, II, FI, IJ, FJ
      REAL(KIND=IP),DIMENSION(I_S-G_P:D_S+G_P)  ::   U,DU
!.......................................................................
                I=1
                DU(I) = DX6B(U,I,1,7,1)
                I=2
                DU(I) = DX5B(U,I,1,7,1)
                I=3
                DU(I) = DX4B(U,I,1,7,1)
                DO I=4,imax-3      
                    DU(I)  = DUDXB(U,I)
                END DO
                I=imax-2
                DU(I) = DX4B(U,I,7,1,-1)
                I=imax-1
                DU(I) = DX5B(U,I,7,1,-1)
                I=imax
                DU(I) = DX6B(U,I,7,1,-1)

END SUBROUTINE


real(kind=ip) function dudxb(u,i)

    integer  ::  k,i,j,z
    real(kind=ip),dimension(i_s-g_p:d_s+g_p)  ::   U,dU
    real(kind=ip),dimension(3)  ::  a 

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudxb=0.d0

    do z=1,3
         dudxb= dudxb + (a(z)*(u(i+z)-u(i-z)))/dm
    end do

end function dudxb

real(kind=ip) function dx4b(u,i,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(i_s-g_p:d_s+g_p)  ::   U,dU
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dx4b=0.d0
        
       do z=iz,fz,p
          dx4b= dx4b+a(z)*u(i+z-3)/dm
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dx4b=0.d0
   
       do z=iz,fz,p
          dx4b= dx4b+a(z)*u(i+z-5)/dm
       end do 

    end if

end function dx4b

real(kind=ip) function dx5b(u,i,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(i_s-g_p:d_s+g_p)  ::   U,dU
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dx5b=0.d0
        
       do z=iz,fz,p
          dx5b= dx5b+a(z)*u(i+z-2)/dm
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dx5b=0.d0
   
       do z=iz,fz,p
          dx5b= dx5b+a(z)*u(i+z-6)/dm
       end do 

    end if

end function dx5b

real(kind=ip) function dx6b(u,i,iz,fz,p)

    integer  ::  k,i,j,z,iz,fz,p
    real(kind=ip),dimension(i_s-g_p:d_s+g_p)  ::   U,dU
    real(kind=ip),dimension(7)  ::  a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dx6b=0.d0
        
       do z=iz,fz,p
          dx6b= dx6b+a(z)*u(i+z-1)/dm
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dx6b=0.d0
   
       do z=iz,fz,p
          dx6b= dx6b+a(z)*u(i+z-7)/dm
       end do 

    end if

end function dx6b


!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................
!.......................................................................

subroutine dermb(dudm,u)

      implicit none
      
      integer  ::   i, j
      
      real(kind=ip)  ::   dv4th, dv2nd, onesp, onesm
      real(kind=ip)  ::   ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(1:imax)  ::   u, dudm
     
    
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del


!.................. Dentro do Dominio ..................................
      do i=3,imax-2
         dudm(i) = dv4th(u(i+2),u(i+1),u(i-1),u(i-2),dm)
      enddo

!..........................Linha Especiais .............................
      i=2
         dudm(i) = dv2nd(u(i+1),u(i-1),dm)

!.......................................................................
      i=imax-1
         dudm(i) = dv2nd(u(i+1),u(i-1),dm)

!.......................................................................
      i=imax
         dudm(i) = onesm(u(i),u(i-1),u(i-2),dm)

!.......................................................................
      i=1
         dudm(i) = onesp(u(i),u(i+1),u(i+2),dm)
      
end subroutine

!Inserts transformations metrics
subroutine metricsY(dU,nk)
  real(kind=ip),dimension(nk,1:imax,1:jmax),intent(inout)   ::   dU
  integer :: i , j, k, nk
do k=1,nk
  do j = 1,jmax
    do i = 1,imax
      !***********************ylayer_left***********************************
         if((i>=D+1).and.(i<=imax-D))then
            if((j>=1).and.(j<=D))then
             du(k,i,j) = du(k,i,j)*meshy(j)
            end if
         end if
      !***********************ylayer_right***********************************
           if((i>=D+1).and.(i<=imax-D))then
              if((j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshy(j)
             end if
           end if
      !***********************corners********************************* 
    !***********************1********************************* 
           if ((i>=1).and.(i<=D))then
              if((j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshy(j)
             end if
           end if
  !***********************2********************************* 
           if ((i>=1).and.(i<=D)) then
              if ((j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
             end if
           end if
  !***********************3********************************* 
           if ((i>=imax-D+1).and.(i<=imax)) then
              if ((j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
             end if
           end if
  !**********************4********************************* 
           if ((i>=imax-D+1).and.(i<=imax)) then
              if ((j>=jmax-D+1).and.(j<=jmax)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
             end if
           end if
  !******************************************************** 
    end do
  end do
end do
end subroutine


!Inserts transformations metrics to temperature derivatives 
subroutine metricsX(dU,nk)
  real(kind=ip),dimension(nk,1:imax,1:jmax),intent(inout)   ::   dU
  integer :: i , j, k, nk
do k=1,nk
  do j = 1,jmax
    do i = 1,imax
 !***********************xlayer_bottom********************************* 
           if ((i>=1).and.(i<=D))then
              if ((j>=D+1).and.(j<=jmax-D)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
 !***********************xlayer_top********************************* 
           if((i>=imax-D+1).and.(i<=imax))then
              if((j>=D+1).and.(j<=jmax-D))then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
 !***********************corners********************************* 
 !***********************1********************************* 
           if ((i>=1).and.(i<=D))then
              if((j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
  !***********************2********************************* 
           if ((i>=1).and.(i<=D)) then
              if ((j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
  !***********************3********************************* 
           if ((i>=imax-D+1).and.(i<=imax)) then
              if ((j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
  !**********************4********************************* 
           if ((i>=imax-D+1).and.(i<=imax)) then
              if ((j>=jmax-D+1).and.(j<=jmax)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
             end if
           end if
  !******************************************************** 
    end do
  end do
end do
end subroutine





end module
