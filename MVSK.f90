program vddm
implicit none

integer io,i,j,hn,l,w,npix
Real*8 pi,dummy,mean,s,s_1,s_2,s_12,e_1,e_2,k_12,skew,kur
Real*8 s1,ss1,ss2,k,k1,k2,k3
Real*8, allocatable :: map(:,:),d1(:,:),d2(:,:),d11(:,:),d22(:,:),d12(:,:)
 Character*100 nam
 Character*1000000 cdummy

!Print*, 'What is the file name?'
Read*, nam
open(1,file=nam)

!Constants
pi=4.0*atan(1.0)
hn=10000000

!Counting rows
w=0
do i=1,hn
    read (1,*,end=101) dummy
    w=w+1
end do
101 rewind(1)

l=0
do
read (1,'(A)',advance='NO',iostat=io) cdummy
if(io<0) exit
if(io>0) stop 'problem reading'
end do
rewind(1)

do i=2,1000000
if (cdummy(i:i).eq.' '.and.cdummy(i-1:i-1).ne.' ')  l=l+1
end do

allocate(map(l,w),d1(l,w),d2(l,w),d11(l,w),d22(l,w),d12(l,w))

!Reading
do j=1,w
    read (1,*) (map(i,j),i=1,l)
end do
 close(1)

mean=0;s=0
do i=1,l
  do j=1,w
    mean=map(i,j)+mean
  end do
end do
mean=mean/(l*w)
do i=1,l
  do j=1,w
    s=s+(map(i,j)-mean)**2
  end do
end do 
s=sqrt(s/(l*w))

do i=1,l
  do j=1,w
    skew=skew+((map(i,j)-mean)/s)**3
  end do
end do 
skew=skew/(l*w)

do i=1,l
  do j=1,w
    kur=kur+((map(i,j)-mean)/s)**4
  end do
end do 
kur=kur/(l*w)

 call vdd(map,l,w,d1,d2,s_1,s_2)
 call vdd(d1,l,w,d11,d12,e_1,s_12)
 call vdd(d2,l,w,d12,d22,s_12,e_2)

s1=0
do i=1,l
  do j=1,w
    s1=s1+d11(i,j)**2
  end do
end do 
s1=sqrt(s1/(l*w))

ss1=0
do i=1,l
  do j=1,w
    ss1=ss1+(-3.)*map(i,j)**2*d11(i,j)*s/(4.*s1**2)
  end do
end do 
ss1=sqrt(ss1/(l*w))

ss2=0
do i=1,l
  do j=1,w
    ss2=ss2+(-3.)*d1(i,j)**2*d11(i,j)*s**3/(s1**4)
  end do
end do 
ss2=sqrt(ss2/(l*w))

k1=0
do i=1,l
  do j=1,w
    k1=k1+(map(i,j)**3*d11(i,j))/(s**4*s1**2)
  end do
end do 
k1=sqrt(k1/(l*w))

k2=0
do i=1,l
  do j=1,w
    k2=k2+(2.*map(i,j)*d1(i,j)**2*d11(i,j)+d1(i,j)**4)/(s**2*s1**4)
  end do
end do 
k2=sqrt(k2/(l*w))

k3=0
do i=1,l
  do j=1,w
    k3=k3+(d1(i,j)**4)/(2.*s**2*s1**4)
  end do
end do 
k3=sqrt(k3/(l*w))

print*, skew/s,ss1,ss2,kur/s**6,k1,k2,k3

end 

!####### subroutines #######

subroutine vdd(map,l,w,d1,d2,s_1,s_2)

integer i,j,l,w,npix,nind
Real*8 s_1,s_2,m1,m2,map(l,w),d1(l,w),d2(l,w)
!Real*8, allocatable :: map(:,:),d1(:,:),d2(:,:)

!allocate(map(l,w),d1(l-1,w),d2(l,w-1))
do i=1,l
  nind=mod(i,l)+1
  d1(i,:)=map(nind,:)-map(i,:)
end do
do j=1,w
  nind=mod(j,w)+1
  d2(:,j)=map(:,nind)-map(:,j)
end do

m1=0
s_1=0
m2=0
s_2=0

do i=1,l
  do j=1,w
    m1=d1(i,j)+m1
  end do
end do
m1=m1/((l)*w)

do i=1,l
  do j=1,w
    s_1=(d1(i,j)-m1)**2+s_1
  end do
end do
s_1=sqrt(s_1/((l)*w))

do i=1,l
  do j=1,w
    m2=d2(i,j)+m2
  end do
end do
m2=m2/((w)*l)

do i=1,l
  do j=1,w
    s_2=(d2(i,j)-m2)**2+s_2
  end do
end do
s_2=sqrt(s_2/((w)*l))

end subroutine

