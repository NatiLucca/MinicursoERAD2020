module mpiown 

use global
!use mpi 

contains   

subroutine exchange_result_b(U_procs,U)

!**********************************************************************************
! Subroutine to join the different matrix of each process with 
! the master process.
! In this subroutine is created a complete matrix with the different 
! matrix  for each process. 

! Input 
! matrix of each process
! U_procs,
! Other information is passed by the global.

! Output 
! U, a completed matrix.

! Created: Jhonatan Aguirre 
! Date   : 08-04-2017 
! Update : 08-04-2017 
!**********************************************************************************
!*****************************************************************************80
  use mpi

  implicit none
  real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U_procs
  real(kind=ip),dimension(1:4,1:imax,1:jmax)::U
  !real(kind=ip),dimension(1:4,i_s:d_s+3,s_s:e_s+3)::U_aux
  real(kind=ip),allocatable,dimension(:,:,:)::U_aux
  integer:: tag,k,column,sender,i,j,tama,source
  integer:: ss,es,is,ds 
  integer:: cont
  integer:: status(MPI_STATUS_SIZE)
  
      if(procs_id/=master)then 

         tama=4*(d_s-i_s+1)*(e_s-s_s+1)

         call  MPI_Send (U_procs(1:4,i_s:d_s,s_s:e_s),tama,MPI_DOUBLE_PRECISION,master,&
                & procs_id,MPI_COMM_WORLD,error)
      else 
 
      do k=1,num_procs-1 
       
          sender  = k 

          is      = com_vector(sender*5+2) 
          ds      = com_vector(sender*5+3) 
          ss      = com_vector(sender*5+4) 
          es      = com_vector(sender*5+5) 

         tama=4*(ds-is+1)*(es-ss+1)

          allocate(U_aux(1:4,is:ds,ss:es))

          call MPI_Recv (U_aux(1:4,is:ds,ss:es),tama,MPI_DOUBLE_PRECISION,&
                      !&    MPI_ANY_SOURCE,MPI_ANY_TAG, & 
                      &    k,k, & 
                      &    MPI_COMM_WORLD, status,error )
          !write(*,*)status(mpi_source)

          U(1:4,is:ds,ss:es)=U_aux(1:4,:,:)
          
          deallocate(U_aux)

      end do 


      U(1:4,i_s:d_s,s_s:e_s)=U_procs(1:4,i_s:d_s,s_s:e_s)

     end if

end subroutine

subroutine size_vector()
!**********************************************************************************
! Subroutine to interchange the size of the differents process with 
! the master process.
! In it is created a vector with the different size of the x and y size for each 
! process and its send to the mastes, where this dates are save in the a vector. 

! Input 
! s_s,e_s,d_s,i_s
! All information is passed by the global.

! Output 
! All information is passed by the global.

! Created: Jhonatan Aguirre 
! Date   : 08-04-2017 
! Update : 08-04-2017 
!**********************************************************************************
  
  implicit none
  !vector to safe the  size of the different procs 
  integer,dimension(1:5):: com_procs 
  integer::tama1, tama2 

  !vector to comunicate the size of the different procs 

  com_procs(1) =procs_id
  com_procs(2) =i_s
  com_procs(3) =d_s
  com_procs(4) =s_s
  com_procs(5) =e_s



 !the size of the information thar is send
 !each processs
  tama1 = 5

  call MPI_GATHER(com_procs(1:5), tama1, MPI_INTEGER, com_vector,&
      &           tama1, MPI_INTEGER, master, mpi_comm_world, error)
  

   if(procs_id==master)then
  !   write(*,*)'com', com_vector
   end if
 
end subroutine


subroutine decomp_band (m, n, r, s, e)

!*****************************************************************************80
!
!! DECOMP_BAND determines the part of the work to be done by this process.
!
!  Discussion:
!
!    This routine is a replacement for the routine MPE_Decomp1d, and
!    is an informed guess as to how it works.
!
!    Given M things in a row, divide them up in consecutive chunks among 
!    N processes so that the number assigned to each process is nearly
!    equal.
!
!    If M is exactly divisible by N, then each process gets M/N consecutive
!    tasks.
!
!    Otherwise, if there is a remainder of P, then the first P processes get
!    M/N+1 tasks, and the last (N-P) processes get M/N tasks.
!
!  Example:
!
!    M = 27, N = 10
!
!    R    S   E
!   --   --  --
!    0:   1   3
!    1:   4   6
!    2:   7   9
!    3:  10  12
!    4:  13  15
!    5:  16  18
!    6:  19  21
!    7:  22  23
!    8:  24  25
!    9:  26  27
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of tasks to be divided among processes.
!
!    Input, integer N, the number of processes.
!
!    Input, integer R, the rank of a particular process.  Because of
!    peculiarities of MPI, R will be between 0 and N-1.
!
!    Output, integer S, E, the first and last tasks to be given to process R.
!    S will be 1 or greater, E will be M or less.
  implicit none

  integer e
  integer m
  integer n
  integer part
  integer r
  integer s
  integer whole

  whole = m / n
  part = m - whole * n

  if ( r + 1 <= part ) then
    s =    r             * ( whole + 1 ) + 1
    e =  ( r + 1 )       * ( whole + 1 )
  else
    s =   r             * whole + 1 + part
    e = ( r + 1 )       * whole     + part
  end if

  return
end subroutine

!call decomp_band (dim_procs, cartesian_id, i_s, d_s, s_s, e_s)

subroutine decomp_matrix ()

!*****************************************************************************80
!
!! DECOMP_BAND determines the part of the work to be done by this process.
!
!  Discussion:
!
!    This routine is a replacement for the routine MPE_Decomp1d, and
!    is an informed guess as to how it works.
!
!    Given M things in a row, divide them up in consecutive chunks among 
!    N processes so that the number assigned to each process is nearly
!    equal.
!
!    If M is exactly divisible by N, then each process gets M/N consecutive
!    tasks.
!
!    Otherwise, if there is a remainder of P, then the first P processes get
!    M/N+1 tasks, and the last (N-P) processes get M/N tasks.
!
!  Example:
!
!    M = 27, N = 10
!
!    R    S   E
!   --   --  --
!    0:   1   3
!    1:   4   6
!    2:   7   9
!    3:  10  12
!    4:  13  15
!    5:  16  18
!    6:  19  21
!    7:  22  23
!    8:  24  25
!    9:  26  27
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 September 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of tasks to be divided among processes.
!
!    Input, integer N, the number of processes.
!
!    Input, integer R, the rank of a particular process.  Because of
!    peculiarities of MPI, R will be between 0 and N-1.
!
!    Output, integer S, E, the first and last tasks to be given to process R.
!    S will be 1 or greater, E will be M or less.
  implicit none

  integer:: part_x, part_y
  integer:: whole_x, whole_y
  integer:: r

  whole_x = imax / dims(1)
  whole_y = jmax / dims(2)

  part_x = imax - whole_x * dims(1)
  part_y = jmax - whole_y * dims(2)


     if ( coords(1) + 1 <= part_x ) then
       i_s =  coords(1)       * ( whole_x + 1 ) + 1
       d_s = (coords(1) + 1 ) * ( whole_x + 1 )
     else
       i_s =  coords(1)      *  whole_x  + 1 + part_x 
       d_s = (coords(1) + 1 )*  whole_x      + part_x 
     end if


     if ( coords(2) + 1 <= part_y ) then

       s_s =  coords(2)      * ( whole_y + 1 ) + 1
       e_s = (coords(2) + 1 )* ( whole_y + 1 )
     else
       s_s =  coords(2)      * whole_y  + 1 + part_y 
       e_s = (coords(2) + 1 )* whole_y      + part_y 
     end if

end subroutine

subroutine exchange_bandtau_x()
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

      req  = MPI_REQUEST_NULL
      tama = 3*3*(e_s-(s_s)+1)
      tag  = 1

      call MPI_SendRecv ( tau(1:3,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          tau(1:3,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          cartesian_comm, req(1), ierr)
      tag = 2
      call MPI_SendRecv ( tau(1:3,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          tau(1:3,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          cartesian_comm, req(1), ierr)
  return
end subroutine

subroutine exchange_bandtau_y()
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

    req = MPI_REQUEST_NULL
    tama=3*3*(d_s-(i_s)+1)

      tag = 1
      call MPI_SendRecv ( tau(1:3,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          tau(1:3,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
      tag = 2
      call MPI_SendRecv ( tau(1:3,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          tau(1:3,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)
      return
end subroutine

subroutine exchange_band_dudm(U)
   implicit none
   real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

      req = MPI_REQUEST_NULL
   
      tama=4*3*(e_s-(s_s)+1)

! left to right 
!**********************************************************
      tag = 1
      call MPI_SendRecv ( U(1:4,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          U(1:4,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                       cartesian_comm, req(1), ierr)
!**********************************************************
!right to left 
      tag = 2
      call MPI_SendRecv ( U(1:4,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          U(1:4,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                       cartesian_comm, req(1), ierr)
  return
end subroutine

subroutine exchange_band_dudmb(Ub)
   implicit none
   real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::Ub
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

      req = MPI_REQUEST_NULL
   
      tama=4*3*(e_s-(s_s)+1)

! left to right 
!**********************************************************
      tag = 1
      call MPI_SendRecv ( Ub(1:4,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          Ub(1:4,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                       cartesian_comm, req(1), ierr)
!**********************************************************
!right to left 
      tag = 2
      call MPI_SendRecv ( Ub(1:4,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          Ub(1:4,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                       cartesian_comm, req(1), ierr)
  return
end subroutine

subroutine exchange_band_dudn(U)
   implicit none
   real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::dUdt,U
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   !integer status_array(MPI_STATUS_SIZE,48), req(48), ierr
   integer:: tag,tag2,k,tama
    req = MPI_REQUEST_NULL
    tama=4*3*(d_s-(i_s)+1)
!**********************************************************
!  ! Bottom to Top
      tag = 1
      call MPI_SendRecv ( U(1:4,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          U(1:4,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
!**********************************************************
!  !  Top to Bottom 
      tag = 2!+10*(k-1)
      call MPI_SendRecv ( U(1:4,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          U(1:4,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)
end subroutine

subroutine exchange_band_dudnb(Ub)
   implicit none
   real(kind=ip),dimension(1:4,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::dUdt,Ub
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   !integer status_array(MPI_STATUS_SIZE,48), req(48), ierr
   integer:: tag,tag2,k,tama
    req = MPI_REQUEST_NULL
    tama=4*3*(d_s-(i_s)+1)
!**********************************************************
!  ! Bottom to Top
      tag = 1
      call MPI_SendRecv ( Ub(1:4,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          Ub(1:4,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
!**********************************************************
!  !  Top to Bottom 
      tag = 2!+10*(k-1)
      call MPI_SendRecv ( Ub(1:4,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          Ub(1:4,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)
end subroutine


!**********************************************************************************
subroutine exchange_bandtemp_y()
!SUBROUTINE FOR EXCHANGE VALUES OF THE TEMPERATURE BETWEEN NEIGHBORING PROCESS - J DIRECTION
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   !integer status_array(MPI_STATUS_SIZE,48), req(48), ierr
   integer:: tag,tag2,k,tama
   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
    req = MPI_REQUEST_NULL
    tama=1*3*(d_s-(i_s)+1)
!**********************************************************
!  ! Bottom to Top
      tag = 1
      call MPI_SendRecv ( temp(1,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          temp(1,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
!**********************************************************
!  !  Top to Bottom 
      tag = 2
      call MPI_SendRecv ( temp(1,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          temp(1,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)
end subroutine
!**********************************************************************************



!**********************************************************************************
subroutine exchange_bandtemp_x()
!SUBROUTINE FOR EXCHANGE VALUES OF THE TEMPERATURE BETWEEN NEIGHBORING PROCESS - I DIRECTION
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
      req = MPI_REQUEST_NULL
   
   !This is the size of the matrix that will be 
   !exchaneg
   !3(three varibles)*3(goths points)*(size of the band in y direction) 
      tama=1*3*(e_s-(s_s)+1)

! left to right 
!**********************************************************

      tag = 1
      

      call MPI_SendRecv ( temp(1,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          temp(1,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                       cartesian_comm, req(1), ierr)

!**********************************************************
!right to left 
      tag = 2
      call MPI_SendRecv ( temp(1,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          temp(1,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                       cartesian_comm, req(1), ierr)

  return

end subroutine




!**********************************************************************************
subroutine exchange_dtempdn()
!SUBROUTINE FOR EXCHANGE VALUES OF THE TEMPERATURE DERIVATIVE BETWEEN NEIGHBORING PROCESS - J DIRECTION
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   !integer status_array(MPI_STATUS_SIZE,48), req(48), ierr
   integer:: tag,tag2,k,tama

   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
    req = MPI_REQUEST_NULL
 
    tama=1*3*(d_s-(i_s)+1)

!**********************************************************
!  ! Bottom to Top
      tag = 1


      call MPI_SendRecv ( DTEMPDN(1,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          DTEMPDN(1,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
!**********************************************************
!  !  Top to Bottom 

      tag = 2
      call MPI_SendRecv ( DTEMPDN(1,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          DTEMPDN(1,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)

end subroutine
!**********************************************************************************





!**********************************************************************************
subroutine exchange_dtempdm()
!SUBROUTINE FOR EXCHANGE VALUES OF THE TEMPERATURE DERIVATIVE BETWEEN NEIGHBORING PROCESS - I DIRECTION
   implicit none
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   integer:: tag,tag2,k,tama

   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
      req = MPI_REQUEST_NULL
   
   !This is the size of the matrix that will be 
   !exchaneg
   !3(three varibles)*3(goths points)*(size of the band in y direction) 
      tama=1*3*(e_s-(s_s)+1)

! left to right 
!**********************************************************

      tag = 1
      

      call MPI_SendRecv ( dtempdm(1,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          dtempdm(1,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                       cartesian_comm, req(1), ierr)

!**********************************************************
!right to left 
      tag = 2
      call MPI_SendRecv ( dtempdm(1,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          dtempdm(1,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                       cartesian_comm, req(1), ierr)

  return

end subroutine



subroutine exchange_band_x(U,nk)
!*****************************************************************************
!
!! EXCHANGE_BAND swaps interface information between neighboring processes.
!
!  Discussion:
!
!    The standard MPI send and receive routines are used.  These are
!    blocking routines, which means that control will not return from
!    a call to MPI_Send until the data has been buffered or has reached
!    the receiving process.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 Abril 2017
!
!  Author:
!
!   Jhonatan Andres Aguirre manco 
!
!  Parameters:
!
!    Input, integer NX, the number of (interior) grid lines in the X direction.
!
!    Input, integer S_size, E_size, indicates the first and last indices of Y data
!    to be controlled by this process.
!
!    Input, integer CARTESIAN_COMM, the identifier for the Cartesian 
!    communicator.
!
!    Input, integer BOT, TOP, the identifiers of the processes "below" and
!    "above" this process, with which it must exchange information.
!
!    Input/output, real A(0:NX+1,S_size-1:E_size+1), an array of data.  Column
!    S must be sent to the "bottom" process and "E" to the top process.
!    Values for S-1 and E+1 will be received from the bottom and top
!    processes.
!
   implicit none
   integer:: tag,tag2,k,tama, nk
   real(kind=ip),dimension(1:nk,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::U
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 

   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
      req = MPI_REQUEST_NULL
   
      tama=nk*3*(e_s-(s_s)+1)

! left to right 
!**********************************************************

      tag = 1
      

      call MPI_SendRecv ( U(1:nk,d_s-2:d_s  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                          U(1:nk,i_s-3:i_s-1,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                       cartesian_comm, req(1), ierr)

!**********************************************************
!right to left 
      tag = 2
      call MPI_SendRecv ( U(1:nk,i_s:i_s+2  ,s_s:e_s), tama, MPI_DOUBLE_PRECISION, left , tag,&
                          U(1:nk,d_s+1:d_s+3,s_s:e_s), tama, MPI_DOUBLE_PRECISION, right, tag,&
                       cartesian_comm, req(1), ierr)

  return

end subroutine

subroutine exchange_band_y(U,nk)

!*****************************************************************************80
!
!! EXCHANGE_BAND swaps interface information between neighboring processes.
!
!  Discussion:
!
!    The standard MPI send and receive routines are used.  These are
!    non blocking routines !

!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!   02 Abril 2017 
!  Author:
!
!    Jhonatan Andres Aguirre Manco
!
!  Parameters:
!
!    Input/output, real U(1:4,i_size:d_size,x), an array of data. 
!     Column S must be sent to the "bottom" process and "E" to the top process.
!    Values for S-1 and E+1 will be received from the bottom and top
!    processes.

   implicit none
   integer:: tag,tag2,k,tama,nk
   real(kind=ip),dimension(1:nk,i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::dUdt,U
   integer:: status_array(MPI_STATUS_SIZE,12), ierr
   integer,dimension(12)::req 
   !integer status_array(MPI_STATUS_SIZE,48), req(48), ierr

   !its needed because no all processor exchange dates and is no clear whas it 
   !the reques values
    req = MPI_REQUEST_NULL
 
    tama=nk*3*(d_s-(i_s)+1)

!**********************************************************
!  ! Bottom to Top
      tag = 1


      call MPI_SendRecv ( U(1:nk,i_s:d_s,e_s-2:e_s   ), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          U(1:nk,i_s:d_s,s_s-3:s_s-1 ), tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          cartesian_comm, req(1), ierr)
!**********************************************************
!  !  Top to Bottom 

      tag = 2!+10*(k-1)
      call MPI_SendRecv ( U(1:nk,i_s:d_s,s_s:s_s+2)  , tama, MPI_DOUBLE_PRECISION, bot, tag,&
                          U(1:nk,i_s:d_s,e_s+1:e_s+3), tama, MPI_DOUBLE_PRECISION, top, tag,&
                          cartesian_comm, req(1), ierr)

end subroutine


SUBROUTINE EXCHANGEB_BAND_X(U)
!*****************************************************************************80
!
!! EXCHANGE_BAND SWAPS INTERFACE INFORMATION BETWEEN NEIGHBORING PROCESSES.
!
!  DISCUSSION:
!
!    THE STANDARD MPI SEND AND RECEIVE ROUTINES ARE USED.  THESE ARE
!    BLOCKING ROUTINES, WHICH MEANS THAT CONTROL WILL NOT RETURN FROM
!    A CALL TO MPI_SEND UNTIL THE DATA HAS BEEN BUFFERED OR HAS REACHED
!    THE RECEIVING PROCESS.
!
!  LICENSING:
!
!    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
!
!  MODIFIED:
!
!    27 SEPTEMBER 2017
!
!     BY: MATHEUS SILVA
!
!  AUTHOR:
!
!   JHONATAN ANDRES AGUIRRE MANCO 
!
!  PARAMETERS:
!
!    INPUT, INTEGER NX, THE NUMBER OF (INTERIOR) GRID LINES IN THE X DIRECTION.
!
!    INPUT, INTEGER S_SIZE, E_SIZE, INDICATES THE FIRST AND LAST INDICES OF Y DATA
!    TO BE CONTROLLED BY THIS PROCESS.
!
!    INPUT, INTEGER CARTESIAN_COMM, THE IDENTIFIER FOR THE CARTESIAN 
!    COMMUNICATOR.
!
!    INPUT, INTEGER BOT, TOP, THE IDENTIFIERS OF THE PROCESSES "BELOW" AND
!    "ABOVE" THIS PROCESS, WITH WHICH IT MUST EXCHANGE INFORMATION.
!
!    INPUT/OUTPUT, REAL A(0:NX+1,S_SIZE-1:E_SIZE+1), AN ARRAY OF DATA.  COLUMN
!    S MUST BE SENT TO THE "BOTTOM" PROCESS AND "E" TO THE TOP PROCESS.
!    VALUES FOR S-1 AND E+1 WILL BE RECEIVED FROM THE BOTTOM AND TOP
!    PROCESSES.
!
   IMPLICIT NONE
   REAL(KIND=IP),DIMENSION(I_S-G_P:D_S+G_P,S_S-G_P:E_S+G_P)::U
   INTEGER:: STATUS_ARRAY(MPI_STATUS_SIZE,12), IERR
   INTEGER,DIMENSION(12)::REQ 
   INTEGER:: TAG,TAG2,K,TAMA

   !ITS NEEDED BECAUSE NO ALL PROCESSOR EXCHANGE DATES AND IS NO CLEAR WHAS IT 
   !THE REQUES VALUES
      REQ = MPI_REQUEST_NULL
   
      TAMA=4*3*(E_S-(S_S)+1) 

! LEFT TO RIGHT 
!**********************************************************

      TAG = 1
      

      CALL MPI_SENDRECV ( U(D_S-2:D_S  ,S_S:E_S), TAMA, MPI_DOUBLE_PRECISION, RIGHT, TAG,&
                          U(I_S-3:I_S-1,S_S:E_S), TAMA, MPI_DOUBLE_PRECISION, LEFT , TAG,&
                       CARTESIAN_COMM, REQ(1), IERR)

!**********************************************************
!RIGHT TO LEFT 
      TAG = 2
      CALL MPI_SENDRECV ( U(I_S:I_S+2  ,S_S:E_S), TAMA, MPI_DOUBLE_PRECISION, LEFT , TAG,&
                          U(D_S+1:D_S+3,S_S:E_S), TAMA, MPI_DOUBLE_PRECISION, RIGHT, TAG,&
                       CARTESIAN_COMM, REQ(1), IERR)

  RETURN

END SUBROUTINE





SUBROUTINE EXCHANGEB_BAND_Y(U)

!*****************************************************************************80
!
!! EXCHANGE_BAND SWAPS INTERFACE INFORMATION BETWEEN NEIGHBORING PROCESSES.
!
!  DISCUSSION:
!
!    THE STANDARD MPI SEND AND RECEIVE ROUTINES ARE USED.  THESE ARE
!    NON BLOCKING ROUTINES !

!  LICENSING:
!
!    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE.
!
!  MODIFIED:
!
!    27 SEPTEMBER 2017
!
!     BY: MATHEUS SILVA
! 
!  AUTHOR:
!
!    JHONATAN ANDRES AGUIRRE MANCO
!
!  PARAMETERS:
!
!    INPUT/OUTPUT, REAL U(1:4,I_SIZE:D_SIZE,X), AN ARRAY OF DATA. 
!     COLUMN S MUST BE SENT TO THE "BOTTOM" PROCESS AND "E" TO THE TOP PROCESS.
!    VALUES FOR S-1 AND E+1 WILL BE RECEIVED FROM THE BOTTOM AND TOP
!    PROCESSES.

   IMPLICIT NONE
   REAL(KIND=IP),DIMENSION(I_S-G_P:D_S+G_P,S_S-G_P:E_S+G_P)::DUDT,U
   INTEGER:: STATUS_ARRAY(MPI_STATUS_SIZE,12), IERR
   INTEGER,DIMENSION(12)::REQ 
   !INTEGER STATUS_ARRAY(MPI_STATUS_SIZE,48), REQ(48), IERR
   INTEGER:: TAG,TAG2,K,TAMA

   !ITS NEEDED BECAUSE NO ALL PROCESSOR EXCHANGE DATES AND IS NO CLEAR WHAS IT 
   !THE REQUES VALUES
    REQ = MPI_REQUEST_NULL
 
    TAMA=4*3*(D_S-(I_S)+1)

!**********************************************************
!  ! BOTTOM TO TOP
      TAG = 1


      CALL MPI_SENDRECV ( U(I_S:D_S,E_S-2:E_S   ), TAMA, MPI_DOUBLE_PRECISION, TOP, TAG,&
                          U(I_S:D_S,S_S-3:S_S-1 ), TAMA, MPI_DOUBLE_PRECISION, BOT, TAG,&
                          CARTESIAN_COMM, REQ(1), IERR)
!**********************************************************
!  !  TOP TO BOTTOM 

      TAG = 2!+10*(K-1)
      CALL MPI_SENDRECV ( U(I_S:D_S,S_S:S_S+2)  , TAMA, MPI_DOUBLE_PRECISION, BOT, TAG,&
                          U(I_S:D_S,E_S+1:E_S+3), TAMA, MPI_DOUBLE_PRECISION, TOP, TAG,&
                          CARTESIAN_COMM, REQ(1), IERR)

END SUBROUTINE


end module
