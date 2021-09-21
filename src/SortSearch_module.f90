
module  SortSearch_module
use point_module

implicit none

contains

subroutine binarySearch (a, value, index)
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
end subroutine binarySearch

subroutine Merge(A,NA,B,NB,C,NC)
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
end subroutine MergeSort

subroutine MergeSortSR(n,a)
    implicit none
    integer, intent(in) :: n
    type(point),intent(in out) :: a(n)
    type(point) :: temp((n+1)/2)

    call MergeSort(a,n,temp)

end subroutine MergeSortSR

end module
