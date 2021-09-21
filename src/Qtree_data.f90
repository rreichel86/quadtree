module Qtree_data
    use point_module
    use polygon_module
    implicit none
    integer :: num_node, num_elem, counter
    integer :: num_intrsc_pts
    integer :: debug_zhl
    integer :: mate_zhl
    integer, allocatable :: elm_typ_ma(:,:)
    !integer :: num_mat_sets
    type(point), allocatable :: nodes(:), temp_nodes(:), sc_nodes(:)
    integer, allocatable :: elements(:,:)
    logical, allocatable :: nodes_mask(:)
    
    
    contains
        subroutine QtrDataReset()
            implicit none
            integer :: istat
            num_node = 0
            num_elem = 0
            num_intrsc_pts = 0
            counter = 0
            mate_zhl = 0
            !!elm_typ = 0
            
            deallocate(elm_typ_ma, nodes, elements, temp_nodes, &
                       nodes_mask,  stat=istat)
             
        end subroutine
    
end module
