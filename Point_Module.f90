module point_module
implicit none

Private
Public :: point, assignment(=), operator(==), operator(>), operator(<), operator(<=), operator(>=)
     
Public ::pt2pt_dist, orientation 

Type :: point 
    real*8 :: x
    real*8 :: y
end Type
  
Interface assignment(=)
    module procedure point_to_array
    module procedure array_to_point
    !module procedure array_to_pointArray
    !module procedure pointArray_to_array
end Interface


Interface operator (>)
    module procedure greater_than
end Interface

Interface operator (>=)
    module procedure greater_equal_than
end Interface

Interface operator (<)
    module procedure less_than
end Interface

Interface operator (<=)
    module procedure less_equal_than
end Interface

Interface operator (==)
    module procedure equal_to
end Interface

contains
    
    subroutine  point_to_array(array_result, pt)
        real*8, dimension(2), intent(out) :: array_result
        type(point), intent(in) :: pt
        array_result(1) = pt%x
        array_result(2) = pt%y  
    end subroutine point_to_array

    subroutine  array_to_point(pt_result, array)
        type(point), intent(out) :: pt_result
        real*8, dimension(2), intent(in) :: array
        pt_result%x = array(1)
        pt_result%y = array(2)
    end subroutine array_to_point
    
    !subroutine  array_to_pointArray(pts_result, array)
    !    type (point), intent(out) :: pts_result(:)
    !    real*8, intent(in) :: array(:)
    !    integer ::i, n
    !    
    !    n = size(array)/2
    !    
    !    do i = 1,n
    !    pts_result(i) = array(2*i-1:2*i)
    !    end do
    !    
    !end subroutine array_to_pointArray
    
    !subroutine  pointArray_to_array(array_result, pts)
    !    type (point), intent(in) :: pts(:)
    !    real*8, intent(out) :: array_result(:)
    !    integer ::i, n
    !    
    !    n = size(pts)
    !    
    !    do i = 1,n
    !     array_result(2*i-1:2*i) = pts(i)
    !    end do
    !    
    !end subroutine pointArray_to_array
    

    function equal_to (pt_1, pt_2)
        logical :: equal_to
        type (point), intent(in) :: pt_1, pt_2
        real*8 :: tol
        
        tol = 1.d-10
        
        if((abs(pt_1%y - pt_2%y) .le. tol).and.(abs(pt_1%x - pt_2%x) .le. tol)) then
            equal_to = .true.
            return
        else 
            equal_to = .false.
        end if 
    end function
    
    function greater_than(pt_1, pt_2)
        logical :: greater_than
        type (point), intent(in) :: pt_1, pt_2
        real*8 :: tol
    
        tol = 1.d-10
    
        if (abs(pt_1%y - pt_2%y) .le. tol) then
            if (pt_1%x .gt. pt_2%x) then
                greater_than = .true.
                return
            else 
                greater_than = .false.
                return
            end if 
        else if (pt_1%y .gt. pt_2%y) then
            greater_than = .true.
            return
        else 
           greater_than = .false.
            return 
        end if 
    end function
    
    function greater_equal_than(pt_1, pt_2)
        logical :: greater_equal_than
        type (point), intent(in) :: pt_1, pt_2
        
        if (pt_1 > pt_2 .or. pt_1 == pt_2) then
            greater_equal_than = .true.
        else
            greater_equal_than = .false.
        end if 
        
    end function
    
    function less_than(pt_1, pt_2)
        logical :: less_than
        type (point), intent(in) :: pt_1, pt_2
        real*8 :: tol
    
        tol = 1.d-10
    
        if (abs(pt_1%y - pt_2%y) .le. tol) then
            if (pt_1%x .lt. pt_2%x) then
                less_than = .true.
                return
            else 
                less_than = .false.
                return
            end if 
        else if (pt_1%y .lt. pt_2%y) then
            less_than = .true.
            return
        else 
           less_than = .false.
            return 
        end if 
    end function
    
    function less_equal_than(pt_1, pt_2)
        logical :: less_equal_than
        type (point), intent(in) :: pt_1, pt_2
        
        if (pt_1 < pt_2 .or. pt_1 == pt_2) then
            less_equal_than = .true.
        else
            less_equal_than = .false.
        end if 
        
    end function
    
    function pt2pt_dist(pt_1, pt_2)
        real* 8 :: pt2pt_dist
        type(point), intent(in) :: pt_1
        type(point), intent(in) :: pt_2
       
        pt2pt_dist = sqrt((pt_2%x-pt_1%x)**2+(pt_2%y-pt_1%y)**2)
        
    end function pt2pt_dist
    
    ! Trinagle orientation 
    ! orientation = 1 -> CCW
    ! orientation = -1 -> CW
    ! orintation = 0 -> A,B,C lie in a line
    function orientation(A,B,C)
        integer :: orientation 
        type(point), intent(in) :: A
        type(point), intent(in) :: B
        type(point), intent(in) :: C
        real *8 :: d, tol
        
        tol = 1.d-10
        
         d = (B%x-A%x)*(C%y-A%y)-(B%y-A%y)*(C%x-A%x)
        
        if  ( (d .lt. 0d0)  .and. (abs(d) .gt. tol) ) then 
            orientation = -1
        else if ( (d .gt. 0d0) .and. (abs(d) .gt. tol) ) then
            orientation = 1
        else 
           orientation = 0
        end if 
        
    
    end function orientation 
    
end module    

module seed_point_module
use point_module
implicit none 

type :: seed_point
    type(point) :: pos
    integer :: wpoly
end type 
    

end module

module intrsc_point_module
use point_module
implicit none

type :: intrsc_point
    type(point) :: pos
    integer :: wpoly
    integer :: wside
end type


end module





