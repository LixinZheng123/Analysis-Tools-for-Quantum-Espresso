program average

integer           :: i,j
real              :: r, gofr(8)
character(len=40) :: filename
character(len=1)  :: i_ch

namelist /input/ filename
read(*,input)

do i=0,7
  write(i_ch,'(i1)') i
  j=i+1
  open(unit=j, file=(trim(filename)//'_'//trim(i_ch)//'.gOO'), status='old')
enddo
!
open(unit=9, file=(trim(filename)//'.gOO'), status='unknown')
!
do j=1,199
  do i=1,8
    read(i,*) r, gofr(i)
  enddo
  write(9,*) r, (gofr(8)+gofr(1)+gofr(2)+gofr(3)+gofr(4)+gofr(5)+gofr(6)+gofr(7))/8
enddo

do i=1,8
  close(i)
enddo

end program average
