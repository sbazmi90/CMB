program genus
implicit none

integer i,j,k,num,numc,numr,kmax,kmin,io,ncr,np,kn
Real*8  mean,var,dx,dummy,dmin,dmax,dk,w1,w2,w3,w4
Real*8 ,allocatable :: d(:,:),gv(:),gh(:),g(:)
Real*8 ,allocatable :: upnv(:,:),dwxv(:,:),dwnv(:,:),upxv(:,:),nupnv(:,:),ndwxv(:,:),ndwnv(:,:),nupxv(:,:)
Real*8 ,allocatable :: upnh(:,:),dwxh(:,:),dwnh(:,:),upxh(:,:),nupnh(:,:),ndwxh(:,:),ndwnh(:,:),nupxh(:,:)
 Character*100 nam,snam
 Character*1000000 cdummy

!Print*, 'What is the file name?'
Read*, nam
! Call getarg(1,nam)
open(1,file=nam)

!Constants
dmax=-10000;dmin=10000
np=100

!Counting rows
numr=0
do
    read (1,*,end=101) dummy
    numr=numr+1
end do
101 rewind(1)

numc=0
do
read (1,'(A)',advance='NO',iostat=io) cdummy
if(io<0) exit
if(io>0) stop 'problem in reading'
end do
rewind(1)

do i=2,1000000
if (cdummy(i:i).eq.' '.and.cdummy(i-1:i-1).ne.' ')  numc=numc+1
end do
!print*,numc,numr

Allocate(d(numc,numr))

!Reading
do j=1,numr
    read (1,*) (d(i,j),i=1,numc)
end do
d=d
!Making data standard
if (.true.) then
mean=0; var=0
do j=1,numr
  do i=1,numc
    mean=mean+d(i,j)
  end do
end do
d=d-mean/(numr*numc)
do j=1,numr
  do i=1,numc
    var=var+d(i,j)**2
  end do
end do
d=d/sqrt(var/(numr*numc))
end if

!Finding minimum and maximum of data
do j=1,numr
  do i=1,numc
    if(dmin>d(i,j)) dmin=d(i,j) 
    if(dmax<d(i,j)) dmax=d(i,j)
  end do
end do

 close(1)

dx=(dmax-dmin)/real(np)
kmin=floor(dmin/dx)
kmax=floor(dmax/dx)
num=kmax-kmin+1
!print*,dx

allocate(upnv(kmin:kmax,numr),dwxv(kmin:kmax,numr),dwnv(kmin:kmax,numr),upxv(kmin:kmax,numr))
!allocate(nupnv(num,numr),ndwxv(num,numr),ndwnv(num,numr),nupxv(num,numr))
allocate(upnh(kmin:kmax,numc),dwxh(kmin:kmax,numc),dwnh(kmin:kmax,numc),upxh(kmin:kmax,numc))
!allocate(nupnh(num,numc),ndwxh(num,numc),ndwnh(num,numc),nupxh(num,numc))
allocate(gv(kmin:kmax),gh(kmin:kmax),g(kmin:kmax))

upnv=0;dwxv=0;dwnv=0;upxv=0 !;nupnv=0;ndwxv=0;ndwnv=0;nupxv=0
upnh=0;dwxh=0;dwnh=0;upxh=0 !;nupnh=0;ndwxh=0;ndwnh=0;nupxh=0
gv=0;gh=0;g=0

! Vertical 
dummy=0
do i=1,numc-1
  do j=2,numr-1
    dk=floor(d(i+1,j)/dx)-floor(d(i,j)/dx)
    dk=abs(dk)/dk
    do k=floor(d(i,j)/dx),floor(d(i+1,j)/dx),int(dk)

! Minima
      if (floor(d(i,j+1)/dx)>floor(d(i,j)/dx).and.floor(d(i,j-1)/dx)>floor(d(i,j)/dx)) then
! Up
        if (dk>0) upnv(k,j)=upnv(k,j)+1
! Down
        if (dk<0) dwnv(k,j)=dwnv(k,j)+1
      end if
! Maxima
      if (floor(d(i,j+1)/dx)<floor(d(i,j)/dx).and.floor(d(i,j-1)/dx)<floor(d(i,j)/dx)) then
! Up
        if (dk>0) upxv(k,j)=upxv(k,j)+1
! Down
        if (dk<0) dwxv(k,j)=dwxv(k,j)+1
      end if

    end do
  end do
end do

! Horizontal 
do j=1,numr-1
  do i=2,numc-1
    dk=floor(d(i,j+1)/dx)-floor(d(i,j)/dx)
    dk=abs(dk)/dk
    do k=floor(d(i,j)/dx),floor(d(i,j+1)/dx),int(dk)

! Minima
      if (floor(d(i+1,j)/dx)>floor(d(i,j)/dx).and.floor(d(i-1,j)/dx)>floor(d(i,j)/dx)) then
! Up
        if (dk>0) upnh(k,i)=upnh(k,i)+1
! Down
        if (dk<0) dwnh(k,i)=dwnh(k,i)+1
      end if
! Maxima
      if (floor(d(i+1,j)/dx)<floor(d(i,j)/dx).and.floor(d(i-1,j)/dx)<floor(d(i,j)/dx)) then
! Up
        if (dk>0) upxh(k,i)=upxh(k,i)+1
! Down
        if (dk<0) dwxh(k,i)=dwxh(k,i)+1
      end if

    end do
  end do
end do

! Averaging Vertical Genus
do k=kmin,kmax
  do j=2,numr-1
    gv(k)=gv(k)+upnv(k,j)-dwxv(k,j)+dwnv(k,j)-upxv(k,j)
  end do
end do
gv=gv/real(numr-2)

! Averaging Horizontal Genus
do k=kmin,kmax
  do i=2,numc-1
    gh(k)=gh(k)+upnh(k,j)-dwxh(k,j)+dwnh(k,j)-upxh(k,j)
  end do
end do
gh=gh/real(numc-2)

! Directional Averaging
do k=kmin,kmax
    g(k)=(gh(k)+gv(k))/4.0
end do

if (.false.) then
do k=kmin,kmax
  w1=0;w2=0;w3=0;w4=0
  do j=2,numr-1
    w1=w1+upnv(k,j)
    w2=w2+dwxv(k,j)
    w3=w3+dwnv(k,j)
    w4=w4+upxv(k,j)
  end do
  do i=2,numc-1
    w1=w1+upnh(k,j)
    w2=w2+dwxh(k,j)
    w3=w3+dwnh(k,j)
    w4=w4+upxh(k,j)
  end do
    print*, k,w1,w2,w3,w4
end do
end if

snam=trim(adjustl(nam))//'_genus'
open(2,file=snam)
do k=kmin,kmax
  write(2,*) (k+0.5)*dx,-g(k)
end do
 close(2)


end
