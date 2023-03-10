
module  SortSearch_module
use point_module
use polygon_module
!use Qtree_module

implicit none


contains


recursive function binarySearch (a, value) result (index)
!use point_module
implicit none
type(point), intent(in) :: a(:)
type(point), intent(in) :: value
integer          :: index, mid, n

n = size(a)
mid = size(a)/2 + 1
if (size(a) == 0) then
index = 0
else if (a(mid) > value) then
index = binarySearch(a(:mid-1), value)
else if (.not.(a(mid) > value)) then
index = binarySearch(a(mid+1:), value)
if (index .ne. 0) then
index = mid + index
end if
else
index = mid
end if
end function binarySearch

subroutine binarySearch_2 (a, value, index)
!use point_module
type(point), intent(in) :: a(:)
type(point), intent(in) :: value
integer, intent (out) :: index

integer :: startIndex, stopIndex, middle

startIndex = 1
stopIndex = size(a)
middle = Floor(0.5*(stopIndex + startIndex))

do while( .not.(a(middle) == value) .and. (startIndex < stopIndex))

if( a(middle) > value ) then

stopIndex = middle - 1

else if ( .not. (a(middle) > value) ) then

startIndex = middle + 1

end if

middle = Floor(0.5*(stopIndex + startIndex))

end do

if(a(middle) == value) then

index = middle
return

else
index = 0
return
end if

end subroutine

!subroutine nearestPolygonsSR(Qtr,polygons,list,status)
!
!    use Qtree_module
!    use point_module
!    use polygon_module
!
!    implicit none
!    type (Qtree), intent(in), pointer :: Qtr
!    type (polygon), intent(in) :: polygons(:)
!    integer, intent(inout) :: list(5), status
!
!    integer :: startIndex, stopIndex, middle, counter
!    real*8 :: tol, distQP, rQ, rP
!
!    tol = 1d-10
!    counter = 1
!    startIndex = 1
!    stopIndex = size(Polygons)
!    middle = Floor(0.5*(stopIndex + startIndex))
!
!
!    !rQ = pt2pt_dist(Qtr%center, Qtr%boundary(1))
!    distQP = pt2pt_dist(Qtr%center, Polygons(middle)%center)
!
!
!        if (Qtr%center > Polygons(middle)%center) then
!
!            if (distQP > (rQ+rP)) then
!
!                startIndex = middle + 1
!
!            else
!
!                list(counter) = middle
!                counter = counter + 1
!                startIndex = middle
!
!            end if
!
!        else if (Qtr%center < Polygons(middle)%center) then
!
!
!            if (distQP > (rQ+rP)) then
!
!                stopIndex = middle - 1
!
!            else
!
!                list(counter) = middle
!                counter = counter + 1
!                stopIndex = middle
!
!            end if
!
!
!        else if (Qtr%center == Polygons(middle)%center) then
!
!
!
!        else
!
!
!        end if
!
!        middle = Floor(0.5*(stopIndex + startIndex))
!
!
!end subroutine


subroutine Merge(A,NA,B,NB,C,NC)

!use point_module
implicit none
integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
type(point), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
type(point), intent(in)     :: B(NB)
type(point), intent(in out) :: C(NC)

integer :: I,J,K

I = 1; J = 1; K = 1;
do while(I <= NA .and. J <= NB)
if (A(I) <= B(J)) then
C(K) = A(I)
I = I+1
else
C(K) = B(J)
J = J+1
endif
K = K + 1
enddo
do while (I <= NA)
C(K) = A(I)
I = I + 1
K = K + 1
enddo
return

end subroutine merge

recursive subroutine MergeSort(A,N,T)

!use point_module
implicit none
integer, intent(in) :: N
type(point), dimension(N), intent(in out) :: A
type(point), dimension((N+1)/2), intent (out) :: T
type(point) :: V
integer :: NA,NB

if (N < 2) return
if (N == 2) then
if (A(1) > A(2)) then
V = A(1)
A(1) = A(2)
A(2) = V
endif
return
endif
NA=(N+1)/2
NB=N-NA

call MergeSort(A,NA,T)
call MergeSort(A(NA+1),NB,T)

if (A(NA) > A(NA+1)) then
T(1:NA)=A(1:NA)
call Merge(T,NA,A(NA+1),NB,A,N)
endif
return

end subroutine MergeSort

subroutine MergeSortSR(n,a)
!use point_module
implicit none
integer, intent(in) :: n
type(point),intent(in out) :: a(n)
type(point) :: temp((n+1)/2)

call MergeSort(a,n,temp)

!write(*,*) 'Sorted!'
!do i=1,n
! write(56,'(20f12.4)') A(i)%x, A(i)%y
!end do

end subroutine MergeSortSR



subroutine Merge_2 (A,NA,B,NB,C,NC)

!use point_module
!use polygon_module

implicit none
integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
type(polygon), intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
type(polygon), intent(in)     :: B(NB)
type(polygon), intent(in out) :: C(NC)

integer :: I,J,K

I = 1; J = 1; K = 1;
do while(I <= NA .and. J <= NB)
if (A(I)%center <= B(J)%center) then
C(K) = A(I)
I = I+1
else
C(K) = B(J)
J = J+1
endif
K = K + 1
enddo
do while (I <= NA)
C(K) = A(I)
I = I + 1
K = K + 1
enddo
return

end subroutine merge_2

recursive subroutine MergeSort_2(A,N,T)

!use point_module
!use polygon_module

implicit none
integer, intent(in) :: N
type(polygon), dimension(N), intent(in out) :: A
type(polygon), dimension((N+1)/2), intent (out) :: T
type(polygon) :: V
integer :: NA,NB

if (N < 2) return
if (N == 2) then
if (A(1)%center > A(2)%center) then
V = A(1)
A(1) = A(2)
A(2) = V
endif
return
endif
NA=(N+1)/2
NB=N-NA

call MergeSort_2(A,NA,T)
call MergeSort_2(A(NA+1),NB,T)

if (A(NA)%center > A(NA+1)%center) then
T(1:NA)=A(1:NA)
call Merge_2(T,NA,A(NA+1),NB,A,N)
endif
return

end subroutine MergeSort_2

subroutine sortPolygonsSR(n,a)
!use point_module
!use polygon_module

implicit none
integer, intent(in) :: n
type(polygon),intent(in out) :: a(n)
type(polygon) :: temp((n+1)/2)

call MergeSort_2(a,n,temp)

!write(*,*) 'Sorted!'
!do i=1,n
! write(56,'(20f12.4)') A(i)%x, A(i)%y
!end do

end subroutine sortPolygonsSR



end module
