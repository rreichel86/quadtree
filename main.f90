program quadtree_main


     use point_module
    use polygon_module
    use Qtree_module
    use Qtree_input
    
    implicit none
    
    integer i, j, k, istat
    real(kind=8) :: td(10)
    character(len=5) :: token
    
    integer :: num_vertices, num_helperPoints, mat_nro
    integer :: num_node, counter, idx, aidx, eidx
    integer, allocatable :: divisions(:)
    logical, allocatable :: aHelperPoints(:)
    type(point), allocatable :: vertices(:), helperPoints(:), tmpPoints(:)
    type(point) :: centerPoint
    real(kind=8) :: alpha, radius(2)
    real(kind=8) :: dxi, xi, xcoord, ycoord
    
    character(len=20) filename
    
    
    
    td = 0d0
    
    
    
    
    
    
    
    write(*,*) 'Please enter input file name'
    !read(*,*) filename 
    filename = 'iYeti.txt'
    write(*,1000) filename
    1000 format ('  ','The input file mame is: ', A) 
    
    open(unit=56, file=filename, status='old', action='Read', iostat=istat)
    
    if (istat == 0) then
        
        do 
            
            read(56, *) token   
            
            ! poly
            ! num_poly,num_mapa_sets
            ! // Definition einer Ellipse
            ! nro,typ,mapa_nro,num_vertices,num_helpers,a,b,xc,yc,alpha
            ! // Definition einer Ellipse
            ! ....
            ! // Definition eines beliebigen Polygonzuges
            ! nro,typ,mapa_nro,num_vertices,num_helpers
            ! //
            ! Vertex_nro,division,x-coord,y-coord
            ! .... // Vertices gegen den Uhrzeigersinn angeben!
            ! num_vertices
            ! //
            ! Vertex_nro,x-coord,y-coord
            ! .... //
            ! num_helpers
            ! // Definition eines beliebigen Polygonzuges
            ! ....
            ! num_poly
            
            if (token .eq. 'poly ') then
                
                read(56,*) td(1:2)
                num_poly = int(td(1))
                
                allocate(polygons(num_poly), stat=istat)
                
                ! Loop over all polygons
                do i = 1, num_poly 
                    
                    read(56,*) td(1:10)
                    mat_nro = int(td(3))
                    num_vertices = int(td(4))
                    num_helperPoints = int(td(5))
                    radius = td(6:7)
                    centerPoint = point(td(8),td(9))
                    alpha = td(10)
                    
                    ! Definition einer Ellipse
                    if (radius(1) .gt. 0d0) then 
                        
                        allocate(vertices(num_vertices), stat=istat)
                        vertices = point(0d0,0d0)
                        
                        call ellipse(radius(1),radius(2),centerPoint,alpha,num_vertices,vertices)
                        
                        call polygon_init(polygons(i), vertices, mat_nro)
                        
                        deallocate(vertices, stat=istat)
                        
                    ! Definition eines belibigen Polygonzuges
                    else 
                        
                        allocate(tmpPoints(num_vertices),divisions(num_vertices), &
                        & helperPoints(num_vertices), aHelperPoints(num_vertices), &
                        & stat = istat)
                        
                        
                        tmpPoints = point(0d0, 0d0)
                        divisions = 0
                        helperPoints = point(0d0, 0d0)
                        aHelperPoints = .false.
                        counter = 0
                        
                        ! Loop over vertices 
                        do j = 1, num_vertices
                            read(56,*) td(1:4)
                            divisions(j) = int(td(2))
                            tmpPoints(j) = point(td(3),td(4))
                            if (divisions(j) .gt. 0) counter = counter + divisions(j) - 1
                        end do
                        
                        ! Loop over helper points 
                        do k = 1, num_helperPoints
                            read(56,*) td(1:3)
                            idx = int(td(1))
                            helperPoints(idx) = point(td(2),td(3)) 
                            aHelperPoints(idx) = .true.
                        end do 
                        
                        allocate(vertices(num_vertices + counter), stat=istat)
                        
                        vertices = point(0d0,0d0)
                        
                        
                        counter = 0
                        do j = 1, num_vertices
                            counter = counter + 1
                            vertices(counter) = tmpPoints(j)
                            if (divisions(j) .eq. 0) cycle
                            dxi = 2.d0 / divisions(j) 
                            xi = -1.d0 + dxi
                            
                            aidx = i
                            eidx = i+1
                            if (i .eq. num_vertices) eidx = 1
                            
                            do k = 1, divisions(j) - 1
                                counter = counter + 1
                                
                                if (aHelperPoints(aidx)) then 
                                    
                                    xcoord =  (1/2.d0)*xi*(xi-1) * tmpPoints(aidx)%x + &
                                    &  (1/2.d0)*xi*(xi+1) * tmpPoints(eidx)%x + &
                                    &   1 - xi**2 * helperPoints(aidx)%x
                                    
                                    ycoord =  (1/2.d0)*xi*(xi-1) * tmpPoints(aidx)%y + &
                                    &  (1/2.d0)*xi*(xi+1) * tmpPoints(eidx)%y + &
                                    &   1 - xi**2 * helperPoints(aidx)%y
                                    
                                
                                else
                                    
                                    xcoord =  (1/2.d0)*(xi-1) * tmpPoints(aidx)%x + &
                                    &  (1/2.d0)*(xi+1) * tmpPoints(eidx)%x 
                                    
                                    ycoord =  (1/2.d0)*(xi-1) * tmpPoints(aidx)%y + &
                                    &  (1/2.d0)*(xi+1) * tmpPoints(eidx)%y 
                                    
                                    
                                end if 
                                vertices(counter) = point(xcoord,ycoord)
                                xi = xi + dxi
                                
                            end do
                            
                        end do
                        
                        deallocate(tmpPoints, divisions, helperPoints, aHelperPoints, &
                        & stat = istat)
                        
                        call polygon_init(polygons(i), vertices, mat_nro)
                        
                        deallocate(vertices, stat=istat)
                        
                        
                    end if 
                    
                end do  
                
                ! seed (Optional)
                ! num_seeds
                ! nro,division,x-coord,y-coord
                ! ....
                ! num_seeds
            
            else if (token .eq. 'seed ') then
                
                
                ! qtree
                ! level_min, max_seeds_cells
            
            else if (token .eq. 'qtree') then
                read(56,*) td(1:2)
                level_min = int(td(1))
                max_seed_q = 1
            
            else if (token .eq. 'end  ') then 
                
                exit
                
            end if 
            
        end do 
    
    else 
        write(*,1010) istat
        1010 format(' ','Error opening file: iostat =', i6) 
        
    end if 
    
    close(unit=56)
    
!    if (polygons(1)%init) then
!        if (allocated(seeds)) then
!            num_seeds = size(seeds)
!        call QtreeSR0(num_poly,polygons, num_seeds, seeds) 
!        else 
!            call QtreeSR1(num_poly,polygons)
!        end if  
!    end if 
!    
    
end program quadtree_main



