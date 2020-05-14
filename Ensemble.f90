program avg
implicit none

integer c,lc,uc,kx,k,i,method,kmax,kmin,num
real(8) ds,dummy,kp,xp
real(8),allocatable:: x(:),e(:),nu(:)
!real , allocatable ::
 character(256) :: nam,bcn,ecn,str,cs

kmax=-10000
kmin=10000
num=10000

allocate(x(-num:num),e(-num:num),nu(-num:num))
x=0; e=0; nu=0

!Print*, 'Enter the beginning common name:'
!read*, bcn
bcn='pdf'
!Print*, 'Enter the ending common name:'
!read*, ecn
ecn=''
!Print*, 'Enter the number of files:'
read*, uc
!Print*, 'Enter ds:'
!read*, ds
ds=0.05

lc=1

do c=lc,uc

  write(cs,'(i4)') c
  str=trim(adjustl(bcn))//trim(adjustl(cs))//trim(adjustl(ecn))
  str=trim(adjustl(str))
  open(1,file=str)


    do 
      read(1,*,end=10) kp,xp

      if (xp/xp==1.0) then 
        k=floor(kp/ds)
        x(k)=x(k)+xp
        nu(k)=nu(k)+1

        if (k>kmax) then
          kmax=k
        end if
        if (kmin>k) then
          kmin=k
        end if

      end if

    enddo 

10  rewind(1)
    do 
      read(1,*,end=11) kp,xp

      if (xp/xp==1.0) then 
        k=floor(kp/ds)
        e(k)=e(k)+(xp-x(k)/nu(k))**2
      end if

    enddo 

11  close(1)

end do

do k=kmin,kmax
  if (nu(k)*e(k)/=0.0) then
    x(k)=x(k)/nu(k)
    e(k)=sqrt(e(k)/nu(k))
  end if
enddo

nam=trim(adjustl(bcn))//trim(adjustl(ecn))//'_avg'
open(2,file=nam)
do k=kmin,kmax
  if (x(k)/=0.0) write(2,*)  (k+0.5)*ds,x(k),e(k)
!  write(2,*)  (k+0.5)*ds,x(k),e(k)
end do

end program



