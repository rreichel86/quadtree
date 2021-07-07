program quadtree_main

    use point_module
    use polygon_module
    use Qtree_module
    use Qtree_input
    ! use Qtree_data
    
    implicit none
    
    integer i, j, k, istat, zhl
    
    character(len=20) filenameIn
    character(len=5) token
    character cc*4, yyy*80, xxx*80
    logical pcomp
    integer num_nodes_elem, region
    
    
    write(*,1000)
    1000 format (5x,'Please enter input file name') 
    ! read(*,'(A)') filenameIn
    filenameIn = 'Test.txt'
    write(*,1010) filenameIn
    1010 format (5x,'The input file mame is: ', A) 
    
    ! Read input file 
    call readInputFile(filenameIn)

    ! Generate Quadtree mesh
    if (polygons(1)%init) then
        if (allocated(seeds)) then
            num_seeds = size(seeds)
            call QtreeSR(num_poly,polygons, num_seeds, seeds) 
        else
            call QtreeSR(num_poly,polygons, 0, point(0d0,0d0))
        end if
    end if

    ! TODO: MENU 
    write(*,1080) 
    write(*,1081)
    write(*,1082)
    !read(*,1083) token
    token = 'exit'
    1080 format(5x,'Type COMMAND to continue:'/)
    1081 format(7x,'MPLOT'/ &
        &       7x,'PARV'/ &
        &       7x,'MFEM'/ &
        &       7x,'FEAP'/ &
        &       7x,'RFINE'/ &
        &       7x,'EXIT'/)
    1082 format(3x,'> ',$)
    1083 format(a5)
    if (pcomp(token,'mplot',5)) then
        write(*,*) '    mplot'
    else if (pcomp(token,'parv',4)) then
        write(*,*) '    parv'
    else if (pcomp(token,'mfem',4)) then
        write(*,*) '    nfem'
    else if (pcomp(token,'feap',4)) then
        write(*,*) '    feap'
    else if (pcomp(token,'rfine',5)) then
        write(*,*) '    rfine'
    else if (pcomp(token,'exit',4)) then
        write(*,*) '    exit'
    else 
        write(*,1100) 
        1100 format(5x,'Undefined command')
    end if 

    ! TODO: subroutine writeFeapInputFile()

end program quadtree_main

subroutine readInputFile(filenameIn)
    use codat
    use iofile
    use iosave
    use point_module
    use polygon_module
    use Qtree_input

    implicit none

    character(len=20), intent(in) ::  filenameIn

    real(kind=8) :: td(10)
    character cc*4, yyy*80, xxx*80
    logical pcomp
    
    integer :: num_vertices, num_helperPoints, mat_nro
    integer :: ipos,cumDiv, idx, aidx, eidx
    integer, allocatable :: divisions(:)
    logical, allocatable :: aHelperPoints(:)
    type(point), allocatable :: vertices(:), helperPoints(:), tmpPoints(:)
    type(point) :: centerPoint
    real(kind=8) :: alpha, radius(2)
    real(kind=8) :: dxi, xi, xcoord, ycoord
    
    integer i, j, k, istat
    
    ior = 59
    iow = 0
    lfile = 0
    lread = .false.
    lsave = .false.
    coflg = .false.

    td = 0d0

    open(unit=ior, file=filenameIn, status='old', iostat=istat)
    if (istat == 0) then

        do
1           call pintio(ior,yyy,8)
            read(yyy,1020,err=1) cc
            1020 format(a4)
            write(*,*) ' '

            ! [poly]
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
            
            if (pcomp(cc,'poly',4)) then
                write(*,*) '     <POLY>'
                call dinput(td,2)
                num_poly = int(td(1))
                num_mat_sets = int(td(2))

                allocate(polygons(num_poly), stat=istat)

                ! Loop over all polygons
                do i = 1, num_poly 

                call dinput(td,10)
                mat_nro = int(td(3))
                num_vertices = int(td(4))
                num_helperPoints = int(td(5))
                radius = td(6:7)
                centerPoint = point(td(8),td(9))
                alpha = td(10)

                if (num_mat_sets .lt. mat_nro) mat_nro = 0

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
                    cumDiv = 0

                    ! Loop over vertices
                    do j = 1, num_vertices
                        call dinput(td,4)
                        divisions(j) = int(td(2))
                        tmpPoints(j) = point(td(3),td(4))
                        if (divisions(j) .gt. 0) cumDiv = cumDiv + divisions(j) - 1
                    end do
                    ! end loop over vertices

                    ! Loop over helper points
                    do k = 1, num_helperPoints
                        call dinput(td, 3)
                        idx = int(td(1))
                        helperPoints(idx) = point(td(2),td(3)) 
                        aHelperPoints(idx) = .true.
                    end do 
                    ! end loop over helper points

                    allocate(vertices(num_vertices + cumDiv), stat=istat)

                    vertices = point(0d0,0d0)

                    ipos = 0
                    do j = 1, num_vertices
                        ipos = ipos + 1
                        vertices(ipos) = tmpPoints(j)
                        if (divisions(j) .eq. 0) cycle
                        dxi = 2.d0 / divisions(j)
                        xi = -1.d0 + dxi

                        aidx = i
                        eidx = i+1
                        if (i .eq. num_vertices) eidx = 1

                        do k = 1, divisions(j) - 1
                            ipos = ipos + 1

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
                            vertices(ipos) = point(xcoord,ycoord)
                            xi = xi + dxi

                        end do
                    end do

                    deallocate(tmpPoints, divisions, helperPoints, aHelperPoints, &
                        & stat = istat)

                    call polygon_init(polygons(i), vertices, mat_nro)

                    deallocate(vertices, stat=istat)

                end if
                end do
                ! end loop over all polygons

                write(*,*) '     </POLY>'

            ! [seed] (Optional)
            ! num_seeds
            ! nro,division,x-coord,y-coord
            ! ....
            ! num_seeds
            else if (pcomp(cc,'seed',4)) then
                
                write(*,*) '     <SEED>'
                write(*,*) '     </SEED>'
                
            ! [qtree]
            ! level_min, max_seeds_cells
            else if (pcomp(cc,'qtre',4)) then
                write(*,*) '     <QTRE>'
                call dinput(td,2)
                level_min = int(td(1))
                max_seed_q = 1
                write(*,*) '     </QTRE>'

            ! [end]
            else if (pcomp(cc,'end',3)) then
                write(*,*) '     </END>'
                
                exit
                
            end if
        end do
    else
        write(*,1070) istat
        1070 format(5x,'Error opening file: iostat =', i6) 
    end if
    
    close(unit=ior)

end subroutine


2000 format('feap')
3000 format(/,'coor')
3010 format(i6,',,',f22.16,',',f22.16)
4000 format(/,'elem')
4010 format(i6,',',i6,',')
4020 format(i6,',')
4030 format(i6)
5000 format(/,'mate')
6000 format(/,'end')
7000 format('inte')
8000 format('stop')

end subroutine
