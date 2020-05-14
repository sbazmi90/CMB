Program PDF
implicit none

integer i,j,k,hn,num,numc,numr,io
Real*8  mean,var,skew,kur,dummy
Real*8 ,allocatable :: d(:,:)
 Character*100 nam
 Character*1000000 cdummy

!Print*, 'What is the file name?'
Read*, nam
open(1,file=nam)

!Constants
hn=10000000
num=10000

!Counting rows
numr=0
do i=1,hn
    read (1,*,end=101) dummy
    numr=numr+1
end do
101 rewind(1)

numc=0
do
read (1,'(A)',advance='NO',iostat=io) cdummy
if(io<0) exit
if(io>0) stop 'problem reading'
end do
rewind(1)

do i=2,1000000
if (cdummy(i:i).eq.' '.and.cdummy(i-1:i-1).ne.' ')  numc=numc+1
end do
!print*, 'Size:',numc,numr
Allocate(d(numc,numr))

!Reading
do j=1,numr
    read (1,*) (d(i,j),i=1,numc)
end do

 close(1)

!making standard data
mean=0;var=0;skew=0;kur=0
do i=1,numc
  do j=1,numr
    mean=mean+d(i,j)
  end do
end do
mean=mean/(numc*numr)
d=d-mean

!variance
do i=1,numc
  do j=1,numr
    var=var+d(i,j)**2
  end do
end do
var=sqrt(var/(numr*numc))
d=d/var

!Skewness
do i=1,numc
  do j=1,numr
    skew=skew+d(i,j)**3
  end do
end do
skew=skew/(numr*numc)

!Kurtosis
do i=1,numc
  do j=1,numr
    kur=kur+d(i,j)**4
  end do
end do
kur=kur/(numr*numc)

print*, mean,var,skew,kur

nam=trim(adjustl(nam))//'_pdf'
!open(2,file=nam)

!do k=kmin,kmax
!  write(2,*) k*dx,p(k)/(numc*numr*dx)
!end do
! close(2)


end
