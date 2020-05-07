program cross
implicit none

integer method,i,j,k,hn,num,numc,numr,kmax,kmin,io,ncr,np,k1,k33,nn1,l,kn,ty,tz
Real*8  mean,var,dx,dummy,dmin,dmax,next,past,p(-1000:1000),saddel
integer, allocatable:: xc(:),yc(:),nc(:),ncv(:),xcv(:),ycv(:)
Real*8 ,allocatable :: d(:,:),n(:,:),nt(:),v_nt(:),nv(:,:),hh(:,:),hv(:,:)
 Character*100 nam,snam,thrs
 Character*1000000 cdummy
!open(3,file='horizontal')
!open(4,file='vertical')

!Print*, 'What is the file name?'
Read*, nam
open(1,file=nam)
!Print*, 'Which method do you want?'
Read*, method

!Constants
hn=10000000
dmax=-10000;dmin=10000
kmax=-10000;kmin=10000
np=200

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
if(io>0) stop 'problem in reading'
end do
rewind(1)

do i=2,1000000
if (cdummy(i:i).eq.' '.and.cdummy(i-1:i-1).ne.' ')  numc=numc+1
end do
!print*,numc,numr

Allocate(d(numc,numr),xc(numr*numc),yc(numr*numc),xcv(numr*numc),ycv(numr*numc))

!Reading
!numc=600
!numr=600
!do i=1,numc
do j=1,numr
    read (1,*) (d(i,j),i=1,numc)
! read(1,*)ty,tz,d(i,j)
!enddo
end do

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
num=floor((dmax-dmin)/dx)
!print*,dx
!*******************PDF making
kmax=-10000;kmin=10000
do j=1,numr
  do i=1,numc 
   k1=floor(d(i,j)/dx) 
   p(k1)=p(k1)+1
   if(kmin>k1) kmin=k1
   if(kmax<k1) kmax=k1
  enddo
enddo
p=p/(numr*numc)

do k=kmin,kmax
write(50,*)k*dx,p(k)/dx
enddo

!*******************************





allocate(n(numr,-num:num),nv(numc,-num:num),nt(-num:num),hh(numc,numr),hv(numc,numr))
allocate(v_nt(-num:num),nc(-num:num),ncv(-num:num))
n=0;nt=0;v_nt=0;nc=0;ncv=0;hh=0;hv=0



!**************************
kmax=-10000;kmin=10000



do j=1,numr
  do i=1,numc-1
   
   k1=floor(d(i,j)/dx)  
   k33=floor(d(i+1,j)/dx)
   if(k1.lt.k33)then
    nn1=k33-k1
    do l=1,nn1
      kn=k1+(l-1)
       
    !if(d(i+1,j)>=d(i,j)+dx) then

    !  do k=floor(d(i,j)/dx),floor(d(i+1,j)/dx)

        if (method==1) then 
          n(j,kn)=n(j,kn)+1
          nc(kn)=nc(kn)+1
          xc(nc(kn))=i
          yc(nc(kn))=j
         if(kn==0) hh(i,j)=1
           !hh(i,j)=1
          if(kmin>kn) kmin=kn
          if(kmax<kn) kmax=kn
        end if

        if (method==2) then
      !   next= floor(d(i,j+1)/dx)-floor(d(i,j)/dx)  
         next= floor(d(i,j+1)/dx)-floor(d(i,j-1)/dx)
         saddel=  floor(d(i,j+1)/dx)+floor(d(i,j-1)/dx)-2*floor(d(i,j)/dx)        
! next=d(i,j+1)-d(i,j)
          past=floor(d(i,j-1)/dx)-floor(d(i,j)/dx) 
         ! past=d(i,j-1)-d(i,j)
         ! if (next*past>0.or.next**2+past**2==0) then
!          if (next==0.0.and. past==0.0) then
          if (next==0.0.and. saddel.ne.0.0) then
          ! if (next==0.0) then
! print*,next
                n(j,kn)=n(j,kn)+1
            nc(kn)=nc(kn)+1
            xc(nc(kn))=i
            yc(nc(kn))=j
           if(kn==-1) hh(i,j)=1
            !hh(i,j)=1
            if(kmin>kn) kmin=kn
            if(kmax<kn) kmax=kn
          end if 
        end if

      end do
    end if
  end do
end do
!n=n/(numc-1)
!Position of each cross
do k=kmin,kmax
write(thrs,'(I4)') k
snam=trim(adjustl(nam))//'_t='//trim(adjustl(thrs))
!open(2,file=snam)
  do i=1,nc(k)
!    write(2,*) xc(i),yc(i),d(xc(i),yc(i))*var+mean
  end do
  close(2)
end do

!Averaging
do k=kmin,kmax
  ncr=0
  do j=1,numr
    !if (n(j,k)/=0) ncr=ncr+1
    nt(k)=nt(k)+n(j,k)
  end do
  !if (ncr/=0) nt(k)=nt(k)/ncr
nt(k)=nt(k)/(numr*(numc-1))
end do

!print*,333
!***************vertical
kmin=1000000
kmax=-1000000
do i=1,numc
  do j=1,numr-1
   
   k1=floor(d(i,j)/dx)  
   k33=floor(d(i,j+1)/dx)
   if(k1.lt.k33)then
   !print*,22   
     nn1=k33-k1
    do l=1,nn1
      kn=k1+(l-1)
       
    !if(d(i+1,j)>=d(i,j)+dx) then

    !  do k=floor(d(i,j)/dx),floor(d(i+1,j)/dx)

        if (method==1) then 
        !print*,i,kn,num         
          nv(i,kn)=nv(i,kn)+1
          ncv(kn)=ncv(kn)+1
          xcv(ncv(kn))=i
          ycv(ncv(kn))=j
         if(kn==0) hv(i,j)=1
          
          !print*,kn
          if(kmin>kn) kmin=kn
          if(kmax<kn) kmax=kn
        end if

        if (method==2) then
        ! next= floor(d(i+1,j)/dx)-floor(d(i,j)/dx)  
        next= floor(d(i+1,j)/dx)-floor(d(i-1,j)/dx)
         saddel=  floor(d(i+1,j)/dx)+floor(d(i-1,j)/dx)-2*floor(d(i,j)/dx)         
 ! next=d(i,j+1)-d(i,j)
          past=floor(d(i-1,j)/dx)-floor(d(i,j)/dx) 
         ! past=d(i,j-1)-d(i,j)
         ! if (next*past>0.or.next**2+past**2==0) then
         ! if (next==0.0.and. past==0.0) then
         !if (next==0.0) then
          if (next==0.0.and. saddel.ne.0.0) then
           ! print*,next
                nv(i,kn)=nv(i,kn)+1
            ncv(kn)=ncv(kn)+1
            xcv(ncv(kn))=i
            ycv(ncv(kn))=j
            if(kn==-1) hv(i,j)=1
            if(kmin>kn) kmin=kn
            if(kmax<kn) kmax=kn
          end if 
        end if

      end do
    end if
  end do
end do



do k=kmin,kmax
  ncr=0
  do i=1,numc
   ! if (n(i,k)/=0) ncr=ncr+1
v_nt(k)=v_nt(k)+nv(i,k)    
!v_nt(k)=v_nt(k)+(n(i,k)-nt(k))**2
  end do
 !v_nt(k)=sqrt(v_nt(k)/ncr)
v_nt(k)=v_nt(k)/((numr-1)*numc)
end do

snam=trim(adjustl(nam))//'_cross'
!open(2,file=snam)

do k=kmin,kmax
  if((nt(k)+v_nt(k))/2.0.ge.1e-5) write(20,*)  k*dx,(nt(k)+v_nt(k))/2.0
end do

!do i=1,numc
!write(3,'(10000f30.15)')(hh(i,j),j=1,numr)
!write(4,'(10000f30.15)')(hv(i,j),j=1,numr)
!end do




 close(2)


end
