!Author: Matheus Silva
!Last Update: 09/11/2018
!
!Use this program to generate:
!- Mesh
!- Transversal Profiles
!- Longitudinal Profiles
!From a data file
program grid
IMPLICIT NONE 
INTEGER						:: i,j,ii,jj, ni, nj	 
REAL(8) 					:: GAMMA, pi
parameter (ni=521, nj=481, pi=3.1415926535897932)
REAL(8), DIMENSION	(1:ni,1:nj)	:: x,y,aa,bb,cc,dd,ee,ff,gg,meshx,meshy
REAL(8), DIMENSION	(:,:), allocatable	:: aaa,bbb,ccc,ddd,eee,fff,ggg
real(8) :: RhoMin, RhoMax, Umin, Umax, Vmin, Vmax, Pmin, Pmax, VorMin, VorMax, SOMin, SOMax
integer :: X1, X2, X3, X4, X5, X6, Y1, Y2, Y3, Y4, Y5, XLimitInf, XLimitSup, YLimitSup, YLimitInf
character*100 DataFile

DataFile = trim("contorno_15.dat")

!Interest domain limits to get data
XlimitInf = 0
XlimitSup = 14
YLimitInf = -1
YLimitSup = 1

!Used in the temperature
GAMMA = 1.4d0




RhoMin = 100.0d0
RhoMax = -100.0d0
Umin = 100.0d0
Umax = -100.0d0
Vmin = 100.0d0
Vmax = -100.0d0
Pmin = 100.0d0
Pmax = -100.0d0
VorMin = 100.0d0
VorMax = -100.0d0
SOMin = 100.0d0
SOMax = -100.0d0

!Read Data File
open(2,file="contorno_15.dat",status='old')
    do i=1,ni
        do j=1,nj
            read(2,*) meshx(i,j), meshy(i,j) ,aa(i,j) ,bb(i,j) ,cc(i,j) ,dd(i,j) ,ee(i,j), gg(i,j)
        end do
    end do
close(2)

!Inserting values to x(i,j) and y(i,j)
do j = 1,nj
    do i =1,ni
    	x(i,j) = meshx(i,j)
    end do
end do
do j = 1,nj
    do i =1,ni
	y(i,j) = meshy(1,j)
    end do
end do

!Plotting mesh DataFile used to plot in GnuPlot
open (1,file="mesh.dat")
do j=1,nj,2 !plotar na direção de i
    do i=1,ni,2
    write(1,*) x(i,j), y(i,j)
    end do
    if (j<nj) then
    jj=j+1
    do i=ni,1,-2
    write(1,*) x(i,jj), y(i,jj)
    end do
    end if
end do
do i=1,ni,2 !plotar na direção de j
    do j=nj,1, -2
    write(1,*) x(i,j), y(i,j)
    end do
    if (i<ni) then
    ii=i+1
    do j=1, nj,2
    write(1,*) x(ii,j), y(ii,j)
    end do
    end if
end do
close(1)

!X1, X2...Xn; Y1,...Yn are pre-defined values 
!Where the profiles will be generated
do i = 1, ni
	if (x(i,1) >= XLimitInf) then
		X1 = i
		exit
	end if
end do
do i = 1, ni
	if (x(i,1) >=3.0d0) then
		X2 = i
		exit
	end if
end do
do i = 1, ni
	if (x(i,1) >=5.0d0) then
		X3 = i
		exit
	end if
end do
do i = 1, ni
	if (x(i,1) >=8.0d0) then
		X4 = i
		exit
	end if
end do
do i = 1, ni
	if (x(i,1) >=12.0d0) then
		X5 = i
		exit
	end if
end do
do i = 1, ni
	if (x(i,1) >= XLimitSup) then
		X6 = i
		exit
	end if
end do
do j = 1, nj
	if (y(1,j) >= YLimitInf) then
		Y1 = j
		exit
	end if
end do
do j = 1, nj
	if (y(1,j) >=-0.50d0) then
		Y2 = j
		exit
	end if
end do
do j = 1, nj
	if (y(1,j) >=0.0d0) then
		Y3 = j
		exit
	end if
end do
do j = 1, nj
	if (y(1,j) >=0.50d0) then
		Y4 = j
		exit
	end if
end do
do j = 1, nj
	if (y(1,j) >= YLimitSup) then
		Y5 = j
		exit
	end if
end do
open (1,file="TProfile1.dat")
	do j=Y1, Y5
		write(1,*) x(X1,j), y(X1,j), aa(X1,j), bb(X1,j), cc(X1,j), dd(X1,j), ee(X1,j), GAMMA*dd(X1,j)/aa(X1,j)
	end do
close(1)
open (2,file="TProfile2.dat")
	do j=Y1, Y5
		write(2,*) x(X2,j), y(X2,j), aa(X2,j), bb(X2,j), cc(X2,j), dd(X2,j), ee(X2,j), GAMMA*dd(X2,j)/aa(X2,j)
	end do
close(2)
open (3,file="TProfile3.dat")
	do j=Y1, Y5
		write(3,*) x(X3,j), y(X3,j), aa(X3,j), bb(X3,j), cc(X3,j), dd(X3,j), ee(X3,j), GAMMA*dd(X3,j)/aa(X3,j)
	end do
close(3)
open (4,file="TProfile4.dat")
	do j=Y1, Y5
		write(4,*) x(X4,j), y(X4,j), aa(X4,j), bb(X4,j), cc(X4,j), dd(X4,j), ee(X4,j), GAMMA*dd(X4,j)/aa(X4,j)
	end do
close(4)
open (5,file="TProfile5.dat")
	do j=Y1, Y5
		write(5,*) x(X5,j), y(X5,j), aa(X5,j), bb(X5,j), cc(X5,j), dd(X5,j), ee(X5,j), GAMMA*dd(X5,j)/aa(X5,j)
	end do
close(5)
open (6,file="TProfile6.dat")
	do j=Y1, Y5
		write(6,*) x(X6,j), y(X6,j), aa(X6,j), bb(X6,j), cc(X6,j), dd(X6,j), ee(X6,j), GAMMA*dd(X6,j)/aa(X6,j)
	end do
close(6)
open (7,file="LProfile1.dat")
	do i=X1, X6
		write(7,*) x(i,Y1), y(i,Y1), aa(i,Y1), bb(i,Y1), cc(i,Y1), dd(i,Y1), ee(i,Y1), GAMMA*dd(i,Y1)/aa(i,Y1)
	end do
close(7)
open (7,file="LProfile5.dat")
	do i=X1, X6
		write(7,*) x(i,Y2), y(i,Y2), aa(i,Y2), bb(i,Y2), cc(i,Y2), dd(i,Y2), ee(i,Y2), GAMMA*dd(i,Y2)/aa(i,Y2)
	end do
close(7)
open (7,file="LProfile6.dat")
	do i=X1, X6
		write(7,*) x(i,Y3), y(i,Y3), aa(i,Y3), bb(i,Y3), cc(i,Y3), dd(i,Y3), ee(i,Y3), dd(i,Y3)/aa(i,Y3)
	end do
close(7)

allocate(aaa(X1:X6, Y1:Y5))
allocate(bbb(X1:X6, Y1:Y5))
allocate(ccc(X1:X6, Y1:Y5))
allocate(ddd(X1:X6, Y1:Y5))
allocate(eee(X1:X6, Y1:Y5))
allocate(fff(X1:X6, Y1:Y5))
allocate(ggg(X1:X6, Y1:Y5))

aaa=0.0d0
bbb=0.0d0
ccc=0.0d0
ddd=0.0d0
eee=0.0d0
fff=0.0d0
ggg=0.0d0

do j = Y1,Y5
    do i = X1,X6
	aaa(i,j) = aa(i,j)
	bbb(i,j) = bb(i,j)
	ccc(i,j) = cc(i,j)
	ddd(i,j) = dd(i,j)
	eee(i,j) = ee(i,j)
	fff(i,j) = GAMMA*dd(i,j)/aa(i,j)
	ggg(i,j) = gg(i,j)
    end do
end do

open (1,file="minmax.dat")
	write(1,*) minval(aaa), minval(bbb), minval(ccc), minval(ddd), minval(eee) , minval(fff) , minval(ggg)
	write(1,*) maxval(aaa), maxval(bbb), maxval(ccc), maxval(ddd), maxval(eee) , maxval(fff) , maxval(ggg)
close(1)


end program
