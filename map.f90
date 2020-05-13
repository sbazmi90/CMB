program fgn
implicit none

integer ( kind = 4 )  i,j,l,m,lensav,ier,ldim,lenwrk,dir,model
real*8 pi,s_0,k_c,eta_x,eta_y,h,rand,mean,var,gam,skk
real*8 ,allocatable ::WSAVE(:),WORK(:),map(:,:)
integer, dimension(:), allocatable :: seed
integer :: values(1:8), kk
 complex ( kind = 8 ),allocatable :: c(:,:,:)
 complex (8)   ic
 character*100 nam

 !CALL RANDOM_SEED()
pi=4.0*atan(1.0)
ic=(0,1)

 
 eta_x=1
 eta_y=1
 s_0=1
 k_c=10
 h=0.8


l=1024;
m=1024;

ldim=l
lensav=2*(L+M) + INT(LOG(REAL(L))/log(2.0))+ INT(LOG(REAL(M))/log(2.0)) + 8
lenwrk=2*L*M+10

allocate(c(4,l,m),WSAVE(lensav),WORK(lenwrk),map(ldim,m))
 call random_seed(size=kk)
allocate(seed(1:kk))
seed(:) = values(8)
 call random_seed(put=seed)

 c=0;
do dir=1,4

!i is x and j is y component
do i=1,l
  do j=1,m
    call random_number(rand)
    skk=((s_0*k_c**h)**2)*(4*pi*eta_x*eta_y/(k_c**2+(eta_x*real(i))**2+(eta_y*real(j))**2)**(h+1))
    c(dir,i,j)=exp(ic*rand*2*pi)*sqrt(skk)
    !c(dir,i,j)=exp(ic*rand*2*pi)*2.0*(s_0*k_c**h)* &
!& sqrt(pi*eta_x*eta_y/(k_c**2+(eta_x*real(i))**2+(eta_y*real(j))**2)**(h+1))
  end do 
end do

 call cfft2i (l,m,wsave,lensav,ier)
 call cfft2b (ldim,l,m,c(dir,:,:),wsave,lensav,work,lenwrk,ier )
 
end do

mean=0;
do i=1,l
  do j=1,m
    map(i,j)=0.5*(real(c(1,i,j)+c(2,l-i+1,j)+c(3,i,m-j+1)+c(4,l-i+1,m-j+1)))
    mean=mean+map(i,j)
  end do 
end do
mean=mean/(l*m)
 
map=map-mean


var=0;
do i=1,l
  do j=1,m
    var=var+map(i,j)**2
  end do 
end do


map=map/sqrt(var/(m*l))

mean=0;
do i=1,l
  do j=1,m
    mean=mean+map(i,j)
  end do 
end do
mean=mean/(l*m)
print*,mean

  nam='map'
  open (1,file=nam)
    do i=1,l
      write (1,*) (map(j,i),j=1,m)
    end do 
    close(1)

deallocate(c,WSAVE,WORK,map,seed)
end 
