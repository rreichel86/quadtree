module Qtree_input
    use polygon_module
    use point_module
    logical :: data_seeds
    integer :: num_poly, num_seeds, num_mat_sets
    type(polygon), allocatable :: polygons(:)
    type(point), allocatable :: seeds(:)
    integer :: level_min, max_seed_Q
    
    contains
        subroutine QtrInputReset()
            implicit none
            integer :: i, istat, np
            np = num_poly
            num_poly = 0
            num_seeds = 0
            level_min = 0
            max_seed_Q = 0
            num_mat_sets = 0
            
            if ( allocated(seeds) ) deallocate(seeds, stat=istat)
            
            do i = 1, np
                call polygon_delete(polygons(i))
            end do
            
            deallocate(polygons, stat=istat)
            
        end subroutine
    
end module

