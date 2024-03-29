program quadtree_main

    use point_module
    use polygon_module
    use Qtree_module
    use Qtree_input
    
    implicit none
    
    integer i, j, k, istat, zhl
    
    character(len=229) filenameIn
    character(len=5) token
    character cc*4, yyy*80, xxx*80
    logical pcomp
    integer num_nodes_elem, region
    
    integer :: cstat, estat
    character(len=100) :: cmsg
    type (Qtree), pointer :: Quadtree

    integer cnt, len, status
    character(len=128) commandArg
    logical exists, interactiveMode

    filenameIn = ' '
    interactiveMode = .false.

    ! get command line arguments
    cnt = command_argument_count ()

    ! ./Quadtree -h
    ! ./Quadtree --help
    ! ./Quadtree input.txt
    ! ./Quadtree -i input.txt
    ! ./Quadtree -interactive input.txt
    do i = 1, cnt
        call get_command_argument (i, commandArg, len, status)
        if (status .ne. 0) then 
            write(*,*) 'get_command_argument failed: status = ', status
            stop
        end if

        write(*,*)  

        if ( commandArg(1:len)  == '-h' .or. commandArg(1:len) == '--help' ) then
            write(*,*) 'A Quadtree based mesh generator'
            write(*,*) 'Usage: ./Quadtree [arguments] input.txt'
            write(*,*)
            write(*,*) 'Arguments:'
            write(*,*) ' -h or --help         Prints this help message and exit'
            write(*,*) ' -i, --iteractive     Iteractive mode'
            stop

        else if ( commandArg(1:len)  == '-i' .or. commandArg(1:len) == '--iteractive' ) then
            interactiveMode = .true.

        else if ( commandArg(len-3:len) == '.txt' ) then

            inquire (File = commandArg(1:len), EXIST = exists)
            if (.not. exists) then
                write(*,*) 'Specified input file: ', commandArg(1:len) 
                write(*,*) 'does not exits.'
                stop
            end if
            filenameIn = commandArg(1:len)
            write(*,1010) filenameIn
            1010 format (5x,'The input file mame is: ', A)
            exit

        else
            write(*,*) 'Unknown option argument = ', commandArg(1:len)
            write(*,*) 'More info with: ./Quadtree --help'
            stop
        end if
    end do

    if ( filenameIn == ' ' ) then
        write(*,1000)
        1000 format (5x,'Please enter input file name: ',$) 
        read(*,'(A)') filenameIn
        write(*,1010) filenameIn
    end if 

    ! Read input file 
    call readInputFile(filenameIn)

    ! Generate Quadtree mesh
    if (polygons(1)%init) then
        allocate(Quadtree)
        if (allocated(seeds)) then
            num_seeds = size(seeds)
            call QtreeMeshSR(Quadtree, num_poly,polygons, num_seeds, seeds) 
        else
            call QtreeMeshSR(Quadtree, num_poly,polygons, 0, point(0d0,0d0))
        end if
    end if

    if( .not. interactiveMode ) stop

    do

    write(*,1080) 
    write(*,1082)
    read(*,1083) token
    1080 format(5x,'Type available COMMAND to continue:'/ &
        &       5x,'For help type "HELP"'/)
    1082 format(3x,'> ',$)
    1083 format(a5)
    if (pcomp(token,'mplot',5)) then

        ! write(*,1084)
        ! 1084 format(5x,'mplot'/)
        call execute_command_line("matlab -nodesktop -nosplash -r 'QtreePlotMesh; exit'", exitstat=estat, &
            & cmdstat=cstat, cmdmsg=cmsg)

        if (cstat .gt. 0) then
            print *, "    Command execution failed with error ", trim(cmsg)
        else if (cstat .lt. 0) then 
            print *, "    Command execution not supported"
        else
            print *, "    Command completed with status ", estat
        end if
        

    else if (pcomp(token,'parv',4)) then

        ! write(*,1085)
        ! 1085 format(5x,'parv'/)

    else if (pcomp(token,'mfem',4)) then

        ! write(*,1086)
        ! 1086 format(5x,'mfem'/)
        call writeInputFileForMFEM()

    else if (pcomp(token,'feap',4)) then

        ! write(*,1087)
        ! 1087 format(5x,'feap'/)
        ! call writeFeapInputFile()

    else if (pcomp(token,'rfine',5)) then

        ! write(*,1088)
        ! 1088 format(5x,'rfine'/)
        call read_seeds()
        call QtreeMeshRefineSR(Quadtree, num_poly,polygons, num_seeds, seeds) 
        

    else if (pcomp(token,'help',4)) then
        
        write(*,1081)
        1081 format(/5x,'MPLOT',' - Plot Quadtree Mesh in MATLAB '/ &
            ! &        5x,'PARV ',' - Save Quadtree Mesh Data for visualization with the program PARAVIEW' / &
            &        5x,'MFEM ',' - Generate MFEM input files' / &
            ! &        5x,'FEAP ',' - Genarate FEAP input file' / &
            &        5x,'RFINE',' - Refine Quadtree Mesh' / &
            &        5x,'HELP ',' - List available COMMANDS' / &
            &        5x,'EXIT '/)

    else if (pcomp(token,'exit',4)) then
        
        ! write(*,1089)
        ! 1089 format(5x,'exit'/)
        exit

    else 
        write(*,1100) 
        1100 format(5x,'Undefined command')

    end if
    end do

end program quadtree_main

subroutine readInputFile(filenameIn)
    use codat
    use iofile
    use iosave
    use point_module
    use polygon_module
    use Qtree_input

    implicit none

    character(len=40), intent(in) ::  filenameIn

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
            ! write(*,*) ' '

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
                ! write(*,*) '     <POLY>'
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

                ! write(*,*) '     </POLY>'

            ! [seed] (Optional)
            ! num_seeds
            ! nro,division,x-coord,y-coord
            ! ....
            ! num_seeds
            else if (pcomp(cc,'seed',4)) then
                
                ! write(*,*) '     <SEED>'
                ! write(*,*) '     </SEED>'
                
            ! [qtree]
            ! level_min, max_seeds_cells
            else if (pcomp(cc,'qtre',4)) then
                ! write(*,*) '     <QTRE>'
                call dinput(td,2)
                level_min = int(td(1))
                max_seed_q = 1
                ! write(*,*) '     </QTRE>'

            ! [end]
            else if (pcomp(cc,'end',3)) then
                ! write(*,*) '     </END>'
                
                exit
                
            end if
        end do
    else
        write(*,1070) istat
        1070 format(5x,'Error opening file: iostat =', i6) 
        stop
    end if
    
    close(unit=ior)

end subroutine
 
subroutine read_seeds()

    use point_module
    use Qtree_input

    implicit none

    integer :: ior
    real*8 :: x, y
    integer :: i, istat

    ior = 60
    if ( allocated(seeds) ) deallocate(seeds)
    
    open(unit=60, file='./Output/mfem/seeds.txt', status='old', action='Read', &
        iostat=istat)

        read(60, '(i6)') num_seeds
        allocate(seeds(num_seeds), stat=istat)
        do i=1, num_seeds
            read(60, '(2f32.16)') x,y
            seeds(i) = point(x, y)
        end do

    close(unit=ior)
    
end subroutine

subroutine writeInputFileForMFEM()

    use Qtree_data

    implicit none
    
    integer :: iow
    integer i,j,k
    integer numsec, num_nodes_elem, region

    ! coor.txt
    ! [ number, x-coor, y-coor ]
    iow = 60
    open(unit=iow, file='./Output/mfem/coor.txt', status='unknown')
        do i=1, num_node + num_elem
            write(iow,2000) i, nodes(i)
        end do
    close(unit=iow)
    ! sections.txt
    ! [ isec, region, node_1,...,node_nsec ]
    iow = 60
    open(unit=iow, file='./Output/mfem/sections.txt', status='unknown')
        numsec  = 0
        do i = 1, num_elem
            num_nodes_elem = elements(i,1)
            region = elements(i,2)

            do j = 1, num_nodes_elem
                numsec = numsec + 1
                write(iow, 2010, advance='no') numsec, region
                if (j .ne. num_nodes_elem) write(iow,2020) elements(i,2+j), elements(i,3+j), num_node + i
                if (j .eq. num_nodes_elem) write(iow,2020) elements(i,2+j), elements(i,3), num_node + i
            end do
        end do
    close(unit=iow)
    ! ord.txt
    ! [ isec, pgrad, qgrad ]
    iow = 60
    open(unit=iow, file='./Output/mfem/ord.txt', status='unknown')
        do i=1, numsec
            write(iow,2020) i,1,1
        end do
    close(unit=iow)
2000 format(i6,2f32.16)
2010 format(2i6)
2020 format(3i6)

end subroutine

subroutine writeFeapInputFile()

    use Qtree_data

    implicit none
    
    integer :: iow
    integer i, j, k
    integer num_nodes_elem, region

    iow = 60

    open(unit=iow, file='./Output/feap/iTest.feap', status='unknown')

        !FEAP
        write(iow,2000)
        !COOR
        write(iow,3000)
        ! The indices (1:num_nodes) correspond to polygonal elements' vertices
        ! and the indices (num_nodes + 1 : num_nodes + num_elem) correspond to 
        ! the scaling center of the respective polygonal element.
        do i=1, num_node + num_elem
            write(iow, 3010) i, nodes(i)
        end do 
        !ELEM
        write(iow,4000)
        do i=1, num_elem
            num_nodes_elem = elements(i,1)
            region = elements(i,2)
            
            k = elm_typ_ma(num_nodes_elem, region)
            write(iow, 4010, advance='no') i, k
            do j= 1, num_nodes_elem
                write(iow, 4020, advance='no') elements(i,2+j)
            end do
            write(iow, 4030) num_node + i
        end do
        !MATE
        write(iow,5000)
        !END 
        write(iow,6000)
        !INTE
        write(iow,7000)
        !STOP
        write(iow,8000)

    close(unit=iow)

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
