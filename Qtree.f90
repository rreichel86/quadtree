subroutine QtreeSR(numPolygons, polygons, numSeeds, seeds)

    use point_module
    use seed_point_module
    use Qtree_module
    use polygon_module
    use segment_module
    use SortSearch_module
    use Qtree_data   
    use Qtree_input, only: level_min, max_seed_Q, num_mat_sets
    
    implicit none 
    
   integer, intent(in) :: numPolygons
    type(polygon), intent(inout) :: polygons(numPolygons)
     integer, intent(in) :: numSeeds
    type(point), intent(in) :: seeds(numSeeds)
    type (Qtree), pointer :: root => null()
    type (QtreeList), pointer :: QtrList => null()
    integer :: numVertices, numTotalPoints
    type (seed_point), allocatable :: totalPointsArr(:)

    real*8 :: pi, start, finish, xinp, yinp
    integer :: istat, i,j, n, idx, zhl
    
    pi = 4.*atan(1.)
    
    call cpu_time(start)

    numVertices = 0 
    do i = 1, numPolygons
        numVertices = numVertices + polygons(i)%num_vertices
    end do 
    
    numTotalPoints = numSeeds + numVertices
    Allocate(totalPointsArr(numTotalPoints), stat=istat)
    
    zhl = 0
    do i = 1, numPolygons
        n = polygons(i)%num_vertices
        do j = 1, n
            zhl = zhl + 1
            totalPointsArr(zhl)%pos = polygons(i)%vertices(j)
            totalPointsArr(zhl)%wpoly = i
        end do
    end do

    if (numSeeds .gt. 0) then
        do i = 1, numSeeds
            totalPointsArr(numVertices+i)%pos = seeds(i)
            totalPointsArr(numVertices+i)%wpoly = 0
        end do 
    end if 
    
    ! Quadtree decomposition
    ! Initilize
    call QtrInit(root,polygons(1)%num_vertices,polygons(1)%vertices)
    
    ! Subdivide
    call QtrSubdivide(root,level_min,numTotalPoints,totalPointsArr,max_seed_Q)

    allocate(QtrList)
    call Qtr2List(root,QtrList)
        
    ! Balance
    call QtrBalance(QtrList)

    ! Compute intersections
    call QIntrsPts(root,root, numPolygons, polygons)
    call Qtr2List(root,QtrList)

    call QtrList%countPoints_(num_node)
    num_node = num_node + numVertices
    Deallocate(totalPointsArr, Stat=istat)
    Allocate (Temp_nodes(num_node),nodes_mask(num_node), Stat=istat)
    nodes_mask = .false.

    zhl = 0
    do i = 1, numPolygons
        n = polygons(i)%num_vertices
        do j = 1, n
            zhl = zhl + 1
            Temp_nodes(zhl) = polygons(i)%vertices(j)
        end do
    end do

    n = num_node - numVertices + 1
    call QtrList%savePoints_(n,Temp_nodes(zhl+1:num_node))

    ! sort Temp_nodes
    call MergeSortSR(num_node, Temp_nodes)

    ! filter  Temp_nodes
    do i = 1, numPolygons
        n = polygons(i)%num_vertices
        do j = 1, n
            call binarySearch_2 (Temp_nodes, polygons(i)%vertices(j), idx)
            if(idx .ne. 0)  nodes_mask(idx) = .true.
        end do
    end do
    
    call coordinates(root,root)
    
    num_node = count(nodes_mask)
    Allocate (nodes(num_node+num_elem), elements(num_elem,16), elm_typ_ma(12,num_mat_sets), Stat=istat)
    nodes(1:num_node) = pack(Temp_nodes,nodes_mask)
    
    call MergeSortSR(num_node, nodes(1:num_node))
    deallocate( Temp_nodes, nodes_mask, Stat = istat)
    ! end filter node coords
    
    open(unit=55, file='selm.txt', status='unknown')
        call connectivity(root,root)
    close(55)
    
    mate_zhl = 0
    do i = 1, 12
        do j = 1, num_mat_sets
            
           if ( elm_typ_ma(i,j) .ne. 0 ) then
                mate_zhl = mate_zhl + 1
                elm_typ_ma(i,j) = mate_zhl
           end if      
           
        end do 
    end do 
    
    open(unit=56, file='scor.txt', status='unknown')
        !write(56, '(i6)') num_node
        do i=1, num_node + num_elem    
            write(56,'(i6,2f32.16)')  i, nodes(i)         
        end do
    close(56)
    
    call cpu_time(finish)
    write(*,1999)
    write(*,2000) finish-start
    
    return
    
1000    format(i6','f32.16','f32.16)
1999    format(/4x,'Q T R E E  M E S H')
2000    format(5x,'Total Time            =',f10.3,' seconds')   
    
end subroutine    


subroutine HowManyPoints(bx,n,pts,azhl)
    use point_module
    implicit none
    type(point), intent(in) :: bx(4)
    integer, intent(in) :: n
    type(point), intent(in) :: pts(n)
    integer, intent(out) :: azhl
    integer :: i
    logical :: containsPoint, temp(n)
        
        azhl = 0
        do i=1,n
            temp(i) = containsPoint(bx,pts(i))
        end do
        azhl = count(temp,1)
        
end subroutine

subroutine HowManySeeds(bx,n,seeds,azhl)
    use point_module
    use seed_point_module
    implicit none
    type(point), intent(in) :: bx(4)
    integer, intent(in) :: n
    type(seed_point), intent(in) :: seeds(n)
    integer, intent(out) :: azhl
    integer :: i
    logical :: containsPoint, temp(n)
        
        azhl = 0
        do i=1,n
            temp(i) = containsPoint(bx,seeds(i)%pos)
        end do
        azhl = count(temp,1)
        
end subroutine

subroutine whichPoints(bx,n,pts,azhl,lpts)
    use point_module
    implicit none 
    type(point), intent(in) :: bx(4)
    integer, intent(in) :: n
    type(point), intent(in) :: pts(n)
    integer, intent(in) :: azhl
    type(point), intent(out) :: lpts(azhl)
    integer :: i
    logical :: containsPoint, temp(n)
    
        do i=1,n
            temp(i) = containsPoint(bx,pts(i))
        end do
        
        lpts = pack(pts,temp)

end subroutine

subroutine whichSeeds(bx,n,seeds,azhl,lseeds)
    use point_module
    use seed_point_module
    implicit none 
    type(point), intent(in) :: bx(4)
    integer, intent(in) :: n
    type(seed_point), intent(in) :: seeds(n)
    integer, intent(in) :: azhl
    type(seed_point), intent(out) :: lseeds(azhl)
    integer :: i
    logical :: containsPoint, temp(n)
    
        do i=1,n
            temp(i) = containsPoint(bx,seeds(i)%pos)
        end do
        
        lseeds = pack(seeds,temp)

end subroutine


logical function containsPoint(bx,pt)
        use point_module
        use polygon_module
        
        implicit none
        type(point), intent(in) :: bx(4)
        type(point), intent(in) :: pt
        integer :: sign
        logical :: point_in_polygon
        
        type(polygon) :: box
        
        call polygon_init(box, bx, 0) 
        call point_in_polygonSR (box, pt, sign, point_in_polygon, .true.)
      
        containsPoint = point_in_polygon
        
    end function
    
    subroutine write_seeds(arr)
    use point_module
    implicit none 
    
    real*8, intent(in) :: arr(2)
    
        open(unit=57, file='D:\seed.dat', status='unknown', position='append')
        
            write(57, '(2f32.16)') arr(1), arr(2) 
        
        close(57)
    end subroutine 
 
    subroutine read_seeds(n, pts_arr)
    use point_module 
    implicit none 
    
    integer, intent(in) :: n
    type (point) :: pts_arr(n)
    real*8 :: xinp, yinp
    integer :: i, istat
    
        open(unit=57, file='D:\seed.dat', status='old', action='Read', &
            iostat=istat)
       
            do i=1, n     
                read(57, '(2f32.16)', iostat=istat)   xinp, yinp  
                pts_arr(i) = point(xinp,yinp)
            end do
        
        close(57)
        
    end subroutine
    
    subroutine ellipse (a,b,pt_c,alpha,n,pts)
        use point_module
        implicit none
        type(point), intent(in) :: pt_c 
        real*8, intent(in) :: a, b, alpha
        integer, intent(in) :: n
        type(point), intent(inout) :: pts(n)
        
        real*8 :: theta, pi
        integer :: i 
        
        pi = 4.0d0*datan(1.0d0)
        
        theta = 0d0
        do i=1,n
           
            pts(i)%x = pt_c%x + a*cos(theta)*cosd(alpha) - b*sin(theta)*sind(alpha)
            pts(i)%y = pt_c%y + a*cos(theta)*sind(alpha) + b*sin(theta)*cosd(alpha)
            theta = theta + 2*pi/n
            
        end do
        
    end subroutine
    
    subroutine ucase (string)
        implicit none
        character(len=*), intent(inout) :: string
        
        integer :: i, length
        
        length = len(string)
        
        do i = 1, length
            if (lge(string(i:i),'a').and. lle(string(i:i),'z')) then
                string(i:i) = achar( ichar( string(i:i)) - 32 )
            end if
        
        end do
      
      end subroutine ucase
