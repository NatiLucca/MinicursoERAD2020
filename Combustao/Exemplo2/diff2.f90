module diff
use global
contains
!.................................................................
subroutine deropm(dU,U,ii,fi,ij,fj)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj
      
      real(kind=ip),dimension(4,im,jm):: U,dU

      do k=1,4
!....   ...................................................................
         do j=ij,fj
!....   ...................................................................
           do i=ii+3,fi-3      
               du(k,i,j)  = dudx(u,k,i,j)
           end do
           i=ii+2
           du(k,i,j) = dx4(u,k,i,j,1,7,1)
           i=fi-2
           du(k,i,j) = dx4(u,k,i,j,7,1,-1)
!....     ...................................................................
           i=ii+1
           du(k,i,j) = dx5(u,k,i,j,1,7,1)
           i=fi-1
           du(k,i,j) = dx5(u,k,i,j,7,1,-1)
!....     ...................................................................
           i=ii
           du(k,i,j) = dx6(u,k,i,j,1,7,1)
           i=fi
           du(k,i,j) = dx6(u,k,i,j,7,1,-1)
!....   ...................................................................
         end do 
!....   ...................................................................
     end do 
end subroutine

real function dudx(u,k,i,j)

    integer::k,i,j,z
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(3)::a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudx=0.d0

    do z=1,3
         dudx= dudx + (a(z)*(u(k,i+z,j)-u(k,i-z,j)))/dm
    end do

end function dudx

real function dx4(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

real function dx5(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

real function dx6(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

subroutine deropn(dU,U,ii,fi,ij,fj)

!.......................................................................
!
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj
      
      real(kind=ip),dimension(4,im,jm):: U,dU

      do k=1,4
!....   ...................................................................
         do i=ii,fi
!....   ...................................................................
           do j=ij+3,fj-3      
               du(k,i,j)  = dudy(u,k,i,j)
           end do
           j=ij+2
           du(k,i,j) = dy4(u,k,i,j,1,7,1)
           j=fj-2
           du(k,i,j) = dy4(u,k,i,j,7,1,-1)
!....     ...................................................................
           j=ij+1
           du(k,i,j) = dy5(u,k,i,j,1,7,1)
           j=fj-1
           du(k,i,j) = dy5(u,k,i,j,7,1,-1)
!....     ...................................................................
           j=ij
           du(k,i,j) = dy6(u,k,i,j,1,7,1)
           j=fj
           du(k,i,j) = dy6(u,k,i,j,7,1,-1)
!....   ...................................................................
         end do 
!....   ...................................................................
     end do 

end subroutine

real function dudy(u,k,i,j)

    integer::k,i,j,z
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(3)::a

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

real function dy4(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

real function dy5(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

real function dy6(u,k,i,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(4,im,jm):: U
    real(kind=ip),dimension(7)::a

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

subroutine deropnb(du,u,ij,fj)

!.......................................................................
!
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj
      
      real(kind=ip),dimension(jm):: u,du

!....   ...................................................................
           do j=ij+3,fj-3      
               du(j)  = dudyb(u,j)
           end do
           j=ij+2
           du(j) = dy4b(u,j,1,7,1)
           j=fj-2
           du(j) = dy4b(u,j,7,1,-1)
!....     ...............................................................
           j=ij+1
           du(j) = dy5b(u,j,1,7,1)
           j=fj-1
           du(j) = dy5b(u,j,7,1,-1)
!....     ...............................................................
           j=ij
           du(j) = dy6b(u,j,1,7,1)
           j=fj
           du(j) = dy6b(u,j,7,1,-1)
!....   ...................................................................

end subroutine

real function dudyb(u,j)

    integer::k,i,j,z
    real(kind=ip),dimension(jm):: u
    real(kind=ip),dimension(3)::a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudyb=0.d0

    do z=1,3
         dudyb= dudyb + (a(z)*(u(j+z)-u(j-z)))/dn
    end do

end function dudyb

real function dy4b(u,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(jm):: u
    real(kind=ip),dimension(7)::a

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

real function dy5b(u,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(jm):: u
    real(kind=ip),dimension(7)::a

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

real function dy6b(u,j,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(jm):: u
    real(kind=ip),dimension(7)::a

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

subroutine deropmb(dU,U,ii,fi)

!.......................................................................
!
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj
      
      real(kind=ip),dimension(im):: U,dU

!....   ...................................................................
           do i=ii+3,fi-3      
               du(i)  = dudxb(u,i)
           end do
           i=ii+2
           du(i) = dx4b(u,i,1,7,1)
           i=fi-2
           du(i) = dx4b(u,i,7,1,-1)
!....     ...............................................................
           i=ii+1
           du(i) = dx5b(u,i,1,7,1)
           i=fi-1
           du(i) = dx5b(u,i,7,1,-1)
!....     ...............................................................
           i=ii
           du(i) = dx6b(u,i,1,7,1)
           i=fi
           du(i) = dx6b(u,i,7,1,-1)
!....   ...................................................................

end subroutine

real function dudxb(u,i)

    integer::k,i,j,z
    real(kind=ip),dimension(im):: U
    real(kind=ip),dimension(3)::a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudxb=0.d0

    do z=1,3
         dudxb= dudxb + (a(z)*(u(i+z)-u(i-z)))/dm
    end do

end function dudxb

real function dx4b(u,i,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(im):: U
    real(kind=ip),dimension(7)::a

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

real function dx5b(u,i,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(im):: U
    real(kind=ip),dimension(7)::a

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

real function dx6b(u,i,iz,fz,p)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(im):: U
    real(kind=ip),dimension(7)::a

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
subroutine dernb(dudn,u)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: i, j
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(jm):: u, dudn
     
    
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del


!.................. Dentro do Dominio ..................................
      do j=3,jmax-2
         dudn(j) = dv4th(u(j+2),u(j+1),u(j-1),u(j-2),dn)
      enddo

!..........................Linha Especiais .............................
      j=2
         dudn(j) = dv2nd(u(j+1),u(j-1),dn)

!.......................................................................
      j=jmax-1
         dudn(j) = dv2nd(u(j+1),u(j-1),dn)

!.......................................................................
      j=jmax
         dudn(j) = onesm(u(j),u(j-1),u(j-2),dn)

!.......................................................................
      j=1
         dudn(j) = onesp(u(j),u(j+1),u(j+2),dn)
      
end subroutine

!.................................................................
subroutine dermb(dudm,u)

      implicit none
      
      integer:: i, j
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(jm):: u, dudm
     
    
      
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



subroutine derm(dU,U,ii,fi,ij,fj)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii,fi,ij,fj
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(4,im,jm):: U,dU
    
     
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del

do k=1,4
!.................. Dentro do Dominio ..................................
      do j=ij,fj
      do i=ii+2,fi-2      
         du(k,i,j)  = dv4th(u(k,i+2,j),u(k,i+1,j),u(k,i-1,j),u(k,i-2,j),dm)
      enddo
      enddo
!.......................................................................
      i=ii
      do j=ij,fj
         du(k,i,j) = onesp(u(k,i,j),u(k,i+1,j),u(k,i+2,j),dm)
      enddo
!......................................................................      
      
      i=ii+1
      do j=ij,fj
         du(k,i,j) = dv2nd(u(k,i+1,j),u(k,i-1,j),dm)
      enddo
!.......................................................................     
      
      i=fi
      do j=ij,fj
         du(k,i,j) = onesm(u(k,i,j),u(k,i-1,j),u(k,i-2,j),dm)
      enddo

!.......................................................................
      i=fi-1
      do j=ij,fj
         du(k,i,j) = dv2nd(u(k,i+1,j),u(k,i-1,j),dm)
      enddo

!.......................................................................
end do      
end subroutine

subroutine dern(dU,U,ii,fi,ij,fj)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: ii,fi,ij,fj,k, i,j
      
      real(kind=ip):: dv4th, dv2nd, onesp, onesm
      real(kind=ip):: ap2, ap1, am1, am2, del, a 
      real(kind=ip),dimension(4,im,jm):: U,dU
     
     
      
!......................................................function statement
      dv4th(ap2,ap1,am1,am2,del) =                                      &
     &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
      dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
      onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
      onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del

do  k=1,4
!.................. Dentro do Dominio ..................................
      do j=ij+2,fj-2
      do i=ii,fi
         du(k,i,j) = dv4th(u(k,i,j+2),u(k,i,j+1),u(k,i,j-1),u(k,i,j-2),dn)
      enddo
      enddo

!..........................Linha Especiais .............................
      j=ij+1
      do i=ii,fi
         du(k,i,j) = dv2nd(u(k,i,j+1),u(k,i,j-1),dn)
      enddo

!.......................................................................
      j=fj-1
      do i=ii,fi
         du(k,i,j) = dv2nd(u(k,i,j+1),u(k,i,j-1),dn)
      enddo

!.......................................................................
      j=fj
      do i=ii,fi
         du(k,i,j) = onesm(u(k,i,j),u(k,i,j-1),u(k,i,j-2),dn)
      enddo

!.......................................................................
      j=ij
      do i=ii,fi
         du(k,i,j) = onesp(u(k,i,j),u(k,i,j+1),u(k,i,j+2),dn)
      enddo
!.......................................................................
end do      
end subroutine

end module
