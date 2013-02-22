program matrix_tst
implicit none
real                                :: start, finish
integer                             :: xsize, ysize, x, y, cont, veses
Integer, Dimension(:,:),Allocatable :: mat2d1, mat2d2, mat2dres

ysize = 2000
xsize = 2000

allocate(mat2d1(0:ysize,0:xsize))
cont=0
do x=0, xsize, 1
    do y=0, xsize, 1
        cont=cont+1
        mat2d1(y,x)=cont
    end do
end do

allocate(mat2d2(0:ysize,0:xsize))
mat2d2=transpose(mat2d1)

allocate(mat2dres(0:ysize,0:xsize))
cont=0
do x=0, xsize, 1
    do y=0, xsize, 1
        cont=cont+1
        mat2dres(y,x)=cont
    end do
end do

call cpu_time(start)
write(*,*) "start =", start

do x=0, xsize, 1
    do y=0, xsize, 1
        mat2dres(y, x) = mat2d1(y, x) + mat2d2(y, x)
    end do
end do
call cpu_time(finish)
write(*,*) "finish =", finish
write(*,*) "Time =", finish-start

end program matrix_tst
