module  Qtree_module
use point_module
use seed_point_module
use intrsc_point_module


interface binaryTransformation
    module procedure binaryTransformation_2
    module procedure binaryTransformation_4
end interface


type Qtree
    integer :: level = 0
    integer :: ref(36)
    real(8), allocatable :: vertices(:,:)
    integer, allocatable :: idxVertices(:)
    integer, allocatable :: signos(:)
    real(8), allocatable :: seeds(:,:)
    type (Qtree), pointer :: NW => null() !1
    type (Qtree), pointer :: SW => null() !2
    type (Qtree), pointer :: NE => null() !3
    type (Qtree), pointer :: SE => null() !4
    type (Qtree), pointer :: father => null()
    contains
        Procedure, Pass :: containsPoint_
end type

type QtreePtr
    type (Qtree), Pointer :: Q
end type

type QtreeNode
    type(Qtree), pointer :: Q
    type(QtreeNode), pointer :: next => null()
end type

type QtreeList
    type(QtreeNode), pointer :: HEAD => null()
    type(QtreeNode), pointer :: TAIL => null()
    contains
        Procedure, Pass :: append_
        Procedure, Pass :: pop_
        Procedure, Pass :: printPoints_
        Procedure, Pass :: countPoints_
        Procedure, Pass :: savePoints_
        Procedure, PAss :: isEmpty_
end type

contains
    
    ! QtreeList method
    ! Append a QtreeNode to the QtreeList
    subroutine append_(this, Q)
        implicit none
        class(QtreeList) :: this
        type(Qtree), pointer :: Q
        
        if ( .not. associated(this%HEAD) ) then
            Allocate(this%HEAD)
            this%TAIL => this%HEAD
            this%TAIL%next => null()
            this%TAIL%Q => Q
        else
            Allocate(this%TAIL%next)
            this%TAIL => this%TAIL%next
            this%TAIL%next => null()
            this%TAIL%Q => Q

        end if
    end subroutine
    
    ! QtreeList method
    ! Check if a QtreeList is empty
    logical function isEmpty_(this)
        implicit none
        class(QtreeList) :: this

        isEmpty_ = .false.
        if (  .not. associated(this%HEAD) .and. .not. associated(this%TAIL) ) then
            isEmpty_ = .true.
        end if
    end function
    
    ! QtreeList method
    ! delete/return first QtreeNode from QtreeList
    subroutine pop_(this, Q)
        implicit none
        class(QtreeList) :: this
        type(Qtree), pointer  :: Q
        type(QtreeNode), pointer :: Qnode
        
        if ( this%isEmpty_() ) then
            Q => null()
            return
        end if
        allocate(Qnode)
        Qnode => this%HEAD
        this%HEAD => This%HEAD%next
        if ( .not. associated(this%HEAD) ) this%TAIL => null()
        Q => Qnode%Q
        deallocate(Qnode)
    end subroutine

    ! QtreeList method
    ! Store in an array the points stored in each QtreeNode
    subroutine savePoints_(this,numPoints,pointsArr)
        implicit none
        class(QtreeList) :: this
        integer :: numPoints
        type(point) :: pointsArr(numPoints)
        type(QtreeNode) , pointer :: Qnode

        integer i, zhl, numSeeds, numIntrscPoints

        zhl = 0
        Qnode => this%HEAD
        do
            if ( .not. associated(Qnode) ) exit
            ! print Quad's boundary nodes
            do i = 1, 4
                zhl = zhl + 1
                pointsArr(zhl) = Qnode%Q%Boundary(i)
            end do
            ! print intrsc_points
            numIntrscPoints = Qnode%Q%num_intrsc_points
            if (numIntrscPoints .ne. 0) then
                do i = 1, numIntrscPoints
                    zhl = zhl + 1
                    pointsArr(zhl) = Qnode%Q%intrsc_points(i)%pos
                end do
            end if
            Qnode => Qnode%next
        end do
    end subroutine

    ! QtreeList method
    ! count total number of points stored in each QtreeNode
    subroutine countPoints_(this,numPoints)
        implicit none
        class(QtreeList) :: this
        integer :: numPoints
        type(QtreeNode) , pointer :: Qnode

        numPoints = 0
        Qnode => this%HEAD
        do
            if ( .not. associated(Qnode) ) exit
            ! Quad's boundary nodes
            numPoints = numPoints + 4
            ! intrsc_points
            numPoints = numPoints + Qnode%Q%num_intrsc_points
            Qnode => Qnode%next
        end do

    end subroutine

    recursive subroutine Qtr2List(Qtr,QtrList)
        implicit none
        type (Qtree), pointer :: Qtr
        type (QtreeList), pointer :: QtrList

        if (.not. associated(Qtr%NW)) then
            call QtrList%append_(Qtr)
        else
            call Qtr2List(Qtr%NW,QtrList)
            call Qtr2List(Qtr%SW,QtrList)
            call Qtr2List(Qtr%NE,QtrList)
            call Qtr2List(Qtr%SE,QtrList)
        end if

    end subroutine

    recursive subroutine QtrPrintLeaves(iow, Qtr)
        implicit none
        integer :: iow
        type (Qtree), pointer :: Qtr
        type (QtreeList), pointer :: QtrList

        integer :: i

        if (.not. associated(Qtr%NW)) then
            do i = 1, 4
                if (iow .gt. 0) then
                    write(iow,2000) Qtr%Boundary(i)
                else
                    write(*,2000) Qtr%Boundary(i)
                end if
            end do
        else
            call QtrPrintLeaves(iow, Qtr%NW)
            call QtrPrintLeaves(iow, Qtr%SW)
            call QtrPrintLeaves(iow, Qtr%NE)
            call QtrPrintLeaves(iow, Qtr%SE)
        end if
2000 format(2f22.16)

    end subroutine

    recursive subroutine QtrSubdivideNS(Qtr,level_min)
        use point_module
        use seed_point_module
        use Qtree_data

        implicit none
        Type (Qtree), pointer :: Qtr
        integer, intent(in) :: level_min

        Type(point) :: v(9), childrenBoundary(4,4)
        integer :: i, p, istat, level, azhl, numseeds
        Type(QtreePtr), pointer :: children(:)

        if (Qtr%level .ge. 18) return
        if (Qtr%level .eq. 0) call subdivideQ (Qtr)
        if (.not. associated(Qtr%NW)) then
           call subdivideQ (Qtr)
        end if
        if (associated(Qtr%NW)) then

            allocate(children(4), stat=istat)
            children(1)%Q => Qtr%NW
            children(2)%Q => Qtr%SW
            children(3)%Q => Qtr%NE
            children(4)%Q => Qtr%SE

            do i = 1, 4
                if( children(i)%Q%level.lt.level_min ) call QtrSubdivide (children(i)%Q,level_min)
            end do

        end if
        if ( associated(children) ) deallocate (children, stat=istat)
    end subroutine

    subroutine subdivideQ (Q)
        use point_module

        implicit none
        Type(Qtree), pointer :: Q

        Type(point) :: v(9)
        Type(point) :: childrenBoundary(4,4)
        integer :: i, p, level, istat
        Type(QtreePtr), pointer :: children(:)

        if ( associated(Q%NW) ) return

        do i = 1,4
            v(i) = Q%Boundary(i)
        end do

        v(5)%x = ( v(1)%x + v(2)%x )/2d0
        v(5)%y = ( v(1)%y + v(2)%y )/2d0
        v(6)%x = ( v(2)%x + v(3)%x )/2d0
        v(6)%y = ( v(2)%y + v(3)%y )/2d0
        v(7)%x = ( v(3)%x + v(4)%x )/2d0
        v(7)%y = ( v(3)%y + v(4)%y )/2d0
        v(8)%x = ( v(4)%x + v(1)%x )/2d0
        v(8)%y = ( v(4)%y + v(1)%y )/2d0
        v(9)%x = ( v(1)%x + v(2)%x + v(3)%x + v(4)%x )/4d0
        v(9)%y = ( v(1)%y + v(2)%y + v(3)%y + v(4)%y )/4d0

        childrenBoundary(1,:) = [v(8),v(9),v(7),v(4)]
        childrenBoundary(2,:) = [v(1),v(5),v(9),v(8)]
        childrenBoundary(3,:) = [v(9),v(6),v(3),v(7)]
        childrenBoundary(4,:) = [v(5),v(2),v(6),v(9)]

        allocate (Q%NW, Q%SW, Q%NE, Q%SE, Stat=istat)
        allocate(children(4), stat=istat)

        do i = 1, 4
            children(i)%Q => Null()
        end do

        children(1)%Q => Q%NW
        children(2)%Q => Q%SW
        children(3)%Q => Q%NE
        children(4)%Q => Q%SE

        level = Q%level + 1
        do i = 1, 4
            p = i - 1
            children(i)%Q%level = level
            children(i)%Q%ref = Q%ref
            children(i)%Q%ref(2*level-1) = int(p/2) + 1
            children(i)%Q%ref(2*level) = mod(p,2) + 1
            children(i)%Q%father => Q
            children(i)%Q%Boundary = childrenBoundary(i,:)
        end do

        if ( allocated (Q%seeds) ) then
            do i = 1, 4
                if ( children(i)%Q%containsPoint_(Q%seeds(1)%pos) ) then
                    allocate(children(i)%Q%seeds(1))
                    children(i)%Q%seeds(1) = Q%seeds(1)
                end if
            end do
            deallocate (Q%seeds)
        end if
        if ( associated(children) ) deallocate (children, stat=istat)
    end subroutine

    logical function containsPoint_(this,pt)
        use point_module

        implicit none
        Class(Qtree) :: this
        type(point), intent(in) :: pt

        real*8 :: xmin, ymin, xmax, ymax

        xmin = this%Boundary(1)%x
        ymin = this%Boundary(1)%y
        xmax = this%Boundary(3)%x
        ymax = this%Boundary(3)%y

        containsPoint_ = .false.

        if (xmin .gt. pt%x) return
        if (xmax .lt. pt%x) return
        if (ymin .gt. pt%y) return
        if (ymax .lt. pt%y) return

        containsPoint_ = .true.
    end function

    logical function isQ_in(Qtr)
        implicit none
        type(qtree), pointer :: Qtr
        logical :: in_cond(4), cq0
        integer :: ms, mwp, mip

         cq0 = .false.
         ms = sum(Qtr%signo)
         mwp = sum(Qtr%wpoly)
         mip = Qtr%num_intrsc_points

         if (ms .eq. 0) cq0 =  (Qtr%wpoly(1) .eq. Qtr%wpoly(2) ) &
                               & .and. ( Qtr%wpoly(2) .eq. Qtr%wpoly(3) ) &
                               & .and. ( Qtr%wpoly(3) .eq. Qtr%wpoly(4) )

        in_cond(1) = (ms .gt. -1) &
                     .and. (.not.(ms .eq. 0 .and.  mwp .gt. 4  .and. cq0 .and. mip .eq. 0))
        in_cond(2) = (ms .eq. -1 .and. mip .ne. 0)
        in_cond(3) = (ms .eq. -2 .and. mip .ne. 0)
        in_cond(4) = (ms .eq. -4 .and. mwp .gt. 4 .and. mip .ne. 0)

        if( count(in_cond) .eq.0 ) then
           isQ_in = .false.
        else
           isQ_in = .true.
        end if


     end function

     recursive subroutine QsDelete(Qtr)
        implicit none
        Type (Qtree), pointer :: Qtr
        logical :: Qchildren(4)
        integer :: istat
        
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
        
        
        if (count(QChildren).eq.0) then
                
               if(allocated(Qtr%seeds)) deallocate(Qtr%seeds, stat=istat)
               deallocate(Qtr, stat=istat)
               
                
        else if (count(QChildren).gt.0) then
        
            if (QChildren(1)) call QsDelete(Qtr%NW)
            if (QChildren(2)) call QsDelete(Qtr%SW)
            if (QChildren(3)) call QsDelete(Qtr%NE)
            if (QChildren(4)) call QsDelete(Qtr%SE)
            
        end if
    
    end subroutine

     subroutine QtrInit(Qtr,n,Boundary)
        type(Qtree), pointer, intent(inout) :: Qtr 
        integer, intent(in) :: n
        type(point), intent(in) :: Boundary(n)
        real*8 :: xmin, xmax, ymin, ymax, mx

        integer :: i, istat

        if ( .not. associated(Qtr) ) Allocate(Qtr, Stat=istat)

        xmin = minval(Boundary(1:n)%x)
        xmax = maxval(Boundary(1:n)%x)
        ymin = minval(Boundary(1:n)%y)
        ymax = maxval(Boundary(1:n)%y)

        mx = max(xmax-xmin,ymax-ymin)

        Qtr%Boundary(1) = point(xmin,ymin)
        Qtr%Boundary(2) = point(xmax,ymin)
        Qtr%Boundary(3) = point(xmax,ymax)
        Qtr%Boundary(4) = point(xmin,ymax)

    end subroutine 
    
    subroutine QtrDelete(Qtr)
        type(Qtree), pointer, intent(inout) :: Qtr 
       
    end subroutine
    
    subroutine QtrBalance(QtrList)
        implicit none
        type(QtreeList), pointer :: QtrList
        type(Qtree), pointer :: Q, AQ, NQ
        type(QtreePtr), pointer :: neighbours(:)

        integer :: i, counter, dir, level, levelNQ, refNQ(36)
        logical :: existNQRef, isSplit

        allocate(neighbours(4))
        do i = 1, 4
            neighbours(i)%Q => null()
        end do

        do
            if ( QtrList%isEmpty_() ) exit
            call QtrList%pop_(Q)

            isSplit = .false.
            counter = 0
            ! Search for the current Quad's Neighbours
            ! Loop over directions:
            ! 1 - West
            ! 2 - South
            ! 3 - East
            ! 4 - North

            do dir = 1, 4

                level = Q%level
                refNQ = Q%ref
                existNQRef = .false.

                ! Possible Quad's Neighbour at current direction
                call QtrRefNeighbourQ(level,refNQ(1:2*level),dir,0,existNQRef)
                if (existNQRef) then

                    call nearestCommonAncestor(Q,level,refNQ(1:2*level),AQ)
                    call QtrSearchNeighbourQ(AQ,level,refNQ(1:2*level),NQ)

                    counter = counter + 1
                    neighbours(counter)%Q => NQ

                    if (isSplit) continue
                    ! Check if current Quad has to be split
                    if ( splitQ(Q, dir, NQ) ) then 
                        ! split current Quad 
                        call subdivideQ (Q)

                        isSplit = .true.
                        call QtrList%append_(Q%NW)
                        call QtrList%append_(Q%SW)
                        call QtrList%append_(Q%NE)
                        call QtrList%append_(Q%SE)
                    end if
                end if
            end do

            ! if curret Quad was split
            ! Check if current Quad's Neighbours need to be split 
            if (isSplit) then
                do i = 1, counter
                    NQ => neighbours(i)%Q
                    levelNQ = NQ%level
                    ! Check if Quad's Neighbour is a leaf
                    if ( .not. associated(NQ%NW) ) then
                        ! Check if Quad's Neighbour is larger than
                        ! the current Quad
                        if ( (level+1) - levelNQ > 1 ) then
                            call QtrList%append_(NQ)
                        end if
                    end if
                end do
            end if
        end do
    end subroutine

    logical function splitQ(Q,dir,NQ)
        implicit none 
        type(Qtree), intent(in), pointer :: Q
        integer, intent(in) :: dir 
        type(Qtree), intent(in), pointer :: NQ
        type(Qtree), pointer :: NQchild1, NQchild2
        
        splitQ = .false.
        ! check if current Quad is a leaf
        if ( associated(Q%NW) ) return

        ! check if Quad's neighbour is a leaf
        if ( .not. associated(NQ%NW) ) return

        ! neighbour Quad's children to be checked
        if (dir .eq. 1) then ! West NQ
            NQchild1 => NQ%NE
            NQchild2 => NQ%SE
        else if (dir .eq. 3) then ! East NQ
            NQchild1 => NQ%NW
            NQchild2 => NQ%SW
        else if (dir .eq. 2) then ! South NQ
            NQchild1 => NQ%NW
            NQchild2 => NQ%NE
        else if (dir .eq. 4) then ! North NQ
            NQchild1 => NQ%SW
            NQchild2 => NQ%SE
        end if

        ! check if neighbour Quad's children are leaves
        if ( associated(NQchild1%NW) ) splitQ = .true.
        if ( associated(NQchild2%NW) ) splitQ = .true.

    end function

    subroutine nearestCommonAncestor(Q, levelNQ, refNQ, AQ)
        implicit none
        Type(Qtree), intent(in), pointer :: Q
        integer, intent(in) :: levelNQ
        integer, intent(in) :: refNQ(2*levelNQ)
        type(Qtree), pointer :: AQ

        ! local variables
        integer l, NCA

        NCA = 0
        do l = levelNQ, 1, -1
            if ( (Q%ref(2*l-1) .eq. refNQ(2*l-1)) .and. (Q%ref(2*l) .eq. refNQ(2*l)) ) then 
                NCA = l
                exit 
            end if 
        end do
        
        AQ => Q%father

        do l = levelNQ, NCA, -1
            if (AQ%level .gt. NCA) then
                AQ => AQ%father
            else
                exit
            end if
        end do

    end subroutine 

    subroutine QtrSearchNeighbourQ(AQ,levelNQ,refNQ,NQ)
        implicit none
        type(Qtree), intent(in), pointer :: AQ
        integer, intent(in) :: levelNQ
        integer, intent(in) :: refNQ(2*levelNQ)
        Type(Qtree), pointer :: NQ

        ! local variables
        integer l, NCA, pos(2)

        NQ => AQ

        NCA = AQ%level
        do l = (NCA+1), levelNQ, 1
            ! has children ?
            if (associated(NQ%NW)) then ! yes
                pos = [refNQ(2*l-1), refNQ(2*l)]
                if ( (pos(1) .eq. 1) .and. (pos(2) .eq. 1) ) then
                    NQ => NQ%NW
                else if ( (pos(1) .eq. 1) .and. (pos(2) .eq. 2) ) then
                    NQ => NQ%SW
                else if ( (pos(1) .eq. 2) .and. (pos(2) .eq. 1) ) then
                    NQ => NQ%NE
                else if ( (pos(1) .eq. 2) .and. (pos(2) .eq. 2) ) then
                    NQ => NQ%SE
                end if
            else ! no
                exit
            end if
        end do

    end subroutine

    subroutine QtrRefNeighbourQ(level,ref,dir,Ntyp,status)

        implicit none 
        integer, intent(in) :: level
        integer, intent(inout) :: ref(2*level)
        integer, intent(in) :: dir, Ntyp
        logical, intent(out) :: status

        integer :: N(level,2), lim(2)
        integer :: l

        status = .true.
        !lim(1:2) = 1
        
        ! covert into integer form N_x and N_y
        N(1:level,1) = ref(1:2*level:2)
        N(1:level,2) = ref(2:2*level:2)

        do l = 1, level 
           N(l,1) = N(l,1) - 1
           N(l,2) = N(l,2) - 1
        end do

        ! compute lim_x and lim_y 
        lim(1) = calcLim(level, N(1:level,1))
        lim(2) = calcLim(level, N(1:level,2))

        if (Ntyp == 0) then
            call refEdgeNeighbourQ(dir, level, lim, N, status)
        else if (Ntyp == 1) then
            call refCornerNeighbourQ(dir, level, lim, N, status)
        end if
        
        ! determine reference of searched neighbour by interweaving new N_x and N_y
        if (status) then 
            do l = 1, level 
               N(l,1) = N(l,1) + 1
               N(l,2) = N(l,2) + 1
            end do
            
            ref(1:2*level:2) = N(1:level,1)
            ref(2:2*level:2) = N(1:level,2)
        end if 

    end subroutine

    subroutine refEdgeNeighbourQ(dir, level, lim, N, status)
    ! refEdgeNeighbourQ: perform binary transformation to obtain the
    ! 4 possible edge neighbours:
    ! 1 - West
    ! 2 - South
    ! 3 - East
    ! 4 - North

        implicit none 
        integer, intent(in) :: dir
        integer, intent(in) :: level
        integer, intent(in) :: lim(2)
        integer, intent(inout) :: N(level,2)
        logical, intent(out) :: status
        
        ! local variable 
        logical :: has_same_father 

        has_same_father = .true.
        status = .true.

        ! Check if the neighbour (dir) we are looking for has the saem father
        ! as the current Quad
        if ( (N(level,1) .eq. 0) .and. (N(level,2) .eq. 0) ) then ! NW Quad

            if ( (dir .eq. 1) .or. (dir .eq. 4) )  then 
                has_same_father = .false.
            end if 

        else if ( (N(level,1) .eq. 0) .and. (N(level,2) .eq. 1) ) then ! SW Quad

            if ( (dir .eq. 1) .or. (dir .eq. 2) )  then 
                has_same_father = .false.
            end if 

        else if ( (N(level,1) .eq. 1) .and. (N(level,2) .eq. 0) ) then ! NE Quad

            if ( (dir .eq. 3) .or. (dir .eq. 4) )  then 
                has_same_father = .false.
            end if 

        else if ( (N(level,1) .eq. 1) .and. (N(level,2) .eq. 1) ) then ! SE Quad

            if ( (dir .eq. 2) .or. (dir .eq. 3) )  then
                has_same_father = .false.
            end if 

        end if

        !has the same father ?
        if (has_same_father) then ! yes
        ! perform binary operation on N(level)
            if ( (dir .eq. 1) .or. (dir .eq. 3) )  then
                call binaryTransformation(level, N(1:level,1))
            else if ( (dir .eq. 2) .or. (dir .eq. 4) )  then
                call binaryTransformation(level, N(1:level,2))
            end if 
        else ! no
        ! perform binary operation on N(level) to N(lim-1)
            if ( (dir .eq. 1) .or. (dir .eq. 3) )  then 
                call binaryTransformation(lim(1), level, N(1:level,1), status)
            else if ( (dir .eq. 2) .or. (dir .eq. 4) )  then 
                call binaryTransformation(lim(2), level, N(1:level,2), status)
            end if 
        end if 

    end subroutine

    subroutine refCornerNeighbourQ(dir, level, lim, N, existNQ)
    ! refEdgeNeighbourQ: perform binary transformation to obtain the
    ! 4 possible edge neighbours:
    ! 1 - South West
    ! 2 - South East
    ! 3 - North East
    ! 4 - North West

        implicit none 
        integer, intent(in) :: dir
        integer, intent(in) :: level
        integer, intent(in) :: lim(2)
        integer, intent(inout) :: N(level,2)
        logical, intent(out) :: existNQ
        
        ! local variable 
        logical ::  existNQ_1, existNQ_2

        existNQ_1 = .true.
        existNQ_2 = .true.

        existNQ = .true.


        if ( N(level,2) .eq. 0 ) then  ! North Quads
            if ( (dir .eq. 1) .or. (dir .eq. 2) )  call binaryTransformation(level, N(1:level,2))
            if ( (dir .eq. 3) .or. (dir .eq. 4) )  call binaryTransformation(lim(2), level, N(1:level,2), existNQ_1)
        else if ( N(level,2) .eq. 1 ) then  ! South Quads
            if ( (dir .eq. 1) .or. (dir .eq. 2) )  call binaryTransformation(lim(2), level, N(1:level,2), existNQ_1)
            if ( (dir .eq. 3) .or. (dir .eq. 4) )  call binaryTransformation(level, N(1:level,2))
        end if 

        if ( N(level,1) .eq. 0 ) then  ! West Quads
            if ( (dir .eq. 2) .or. (dir .eq. 3) )  call binaryTransformation(level, N(1:level,1))
            if ( (dir .eq. 1) .or. (dir .eq. 4) )  call binaryTransformation(lim(1), level, N(1:level,1), existNQ_2)
        else if ( N(level,1) .eq. 1 ) then  ! East Quads
            if ( (dir .eq. 2) .or. (dir .eq. 3) )  call binaryTransformation(lim(1), level, N(1:level,1), existNQ_2)
            if ( (dir .eq. 1) .or. (dir .eq. 4) )  call binaryTransformation(level, N(1:level,1))
        end if 

        existNQ =  existNQ_1 .and. existNQ_2

    end subroutine

    integer function calcLim(level, N)
    ! calcLim: calculate lim working backwards through N
    ! lim is teh level at which N(i) first becomes not equal to N(i-1)

        implicit none 
        integer, intent(in) :: level
        integer, intent(in) :: N(level)

        !local variable
        integer :: l

        do l = level,1,-1
            if (l .eq. 1) then
                calcLim = l
                exit
            else if ( N(l) .ne. N(l-1) ) then
                calcLim = l
                exit
            end if
        end do

    end function
    
    subroutine binaryTransformation_2(level, N)
    ! binaryTransformation: perform binary operation on N

        implicit none 
        integer, intent(in) :: level
        integer, intent(inout) :: N(level) 

        N(level) = 1 - N(level)

    end subroutine

    subroutine binaryTransformation_4(lim, level, N, status)
    ! binaryTransformation: perform binary operation on N

        implicit none 
        integer, intent(in) :: lim
        integer, intent(in) :: level
        integer, intent(inout) :: N(level) 
        logical, intent(out) :: status

        !local variable 
        integer :: l

        status = .true.

        if (lim .eq. 1) then 
            status = .false.
            return 
        end if 

        ! N(i) = 1 -> N(i) = 0
        ! N(i) = 0 -> N(i) = 1

        do l = level, lim-1, -1
            N(l) = 1 - N(l)
        end do 

    end subroutine

end module



