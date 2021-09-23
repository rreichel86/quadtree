module  Qtree_module

use point_module

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
        Procedure, Pass :: isEmpty_
        Procedure, Pass :: countPoints_
        Procedure, Pass :: countNodes_
        Procedure, Pass :: savePoints_
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

        integer i, zhl

        zhl = 0
        Qnode => this%HEAD
        do
            if ( .not. associated(Qnode) ) exit
            do i = 1, 4
                zhl = zhl + 1
                pointsArr(zhl) = Qnode%Q%vertices(i,:)
            end do
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
            numPoints = numPoints + 4
            Qnode => Qnode%next
        end do

    end subroutine
    
    subroutine countNodes_(this,numNodes)
        implicit none
        class(QtreeList) :: this
        integer :: numNodes
        type(QtreeNode) , pointer :: Qnode

        numNodes = 0
        Qnode => this%HEAD
        do
            if ( .not. associated(Qnode) ) exit
            numNodes = numNodes + 1
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

        integer :: i

        if (.not. associated(Qtr%NW) .and. allocated(Qtr%vertices) ) then
            do i = 1, 4
                if (iow .gt. 0) then
                    write(iow,2000) Qtr%vertices(i,:)
                else
                    write(*,2000) Qtr%vertices(i,:)
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

        implicit none
        Type (Qtree), pointer :: Qtr
        integer, intent(in) :: level_min

        integer :: i, istat, level
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
                if( children(i)%Q%level.lt.level_min ) call QtrSubdivideNS (children(i)%Q,level_min)
            end do

        end if
        if ( associated(children) ) deallocate (children, stat=istat)
    end subroutine

    subroutine subdivideQ (Q)
        use point_module

        implicit none
        Type(Qtree), pointer :: Q

        real(8) :: v(9,2)
        integer :: i, level, istat
        integer :: idx(4,4)
        Type(QtreePtr), pointer :: children(:)

        if ( associated(Q%NW) ) return

        do i = 1,4
            v(i,:) = Q%vertices(i,:)
        end do

        v(5,1) = ( v(1,1) + v(2,1) )/2d0
        v(5,2) = ( v(1,2) + v(2,2) )/2d0
        v(6,1) = ( v(2,1) + v(3,1) )/2d0
        v(6,2) = ( v(2,2) + v(3,2) )/2d0
        v(7,1) = ( v(3,1) + v(4,1) )/2d0
        v(7,2) = ( v(3,2) + v(4,2) )/2d0
        v(8,1) = ( v(4,1) + v(1,1) )/2d0
        v(8,2) = ( v(4,2) + v(1,2) )/2d0
        v(9,1) = ( v(1,1) + v(2,1) + v(3,1) + v(4,1) )/4d0
        v(9,2) = ( v(1,2) + v(2,2) + v(3,2) + v(4,2) )/4d0

        idx(1,:) = [8,9,7,4]
        idx(2,:) = [1,5,9,8]
        idx(3,:) = [9,6,3,7]
        idx(4,:) = [5,2,6,9]

        allocate (Q%NW, Q%SW, Q%NE, Q%SE, Stat=istat)
        allocate(children(4), stat=istat)

        children(1)%Q => Q%NW
        children(2)%Q => Q%SW
        children(3)%Q => Q%NE
        children(4)%Q => Q%SE

        level = Q%level + 1
        do i = 1, 4
            children(i)%Q%level = level
            children(i)%Q%ref = Q%ref
            children(i)%Q%ref(level) = i
            children(i)%Q%father => Q
            children(i)%Q%vertices(1,:) = v(idx(i,1),:)
            children(i)%Q%vertices(2,:) = v(idx(i,2),:)
            children(i)%Q%vertices(3,:) = v(idx(i,3),:)
            children(i)%Q%vertices(4,:) = v(idx(i,4),:)
        end do

        if ( allocated (Q%seeds) ) then
            do i = 1, 4
                if ( children(i)%Q%containsPoint_(Q%seeds(1,:)) ) then
                    allocate(children(i)%Q%seeds(1,2))
                    children(i)%Q%seeds(1,:) = Q%seeds(1,:)
                end if
            end do
            deallocate (Q%seeds)
        end if
        if ( associated(children) ) deallocate (children, stat=istat)
    end subroutine

    logical function containsPoint_(this,point)
        implicit none

        Class(Qtree) :: this
        real(8), intent(in) :: point(2)

        real*8 :: xmin, ymin, xmax, ymax

        xmin = this%vertices(1,1)
        ymin = this%vertices(1,2)
        
        xmax = this%vertices(3,1)
        ymax = this%vertices(3,2)
        
        containsPoint_ = .false.
        
        if ( point(1) .lt. xmin ) then
            return
        else if ( point(1) .gt. xmax ) then
            return
        else if ( point(2) .lt. ymin ) then
            return
        else if ( point(2) .gt. ymax ) then
            return
        end if

        containsPoint_ = .true.

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

    subroutine QtrInit(Qtr,n,points)
        implicit none
        type(Qtree), pointer, intent(inout) :: Qtr 
        integer, intent(in) :: n
        real(8), intent(in) :: points(n,2)

        ! local variables
        integer :: istat
        real*8 :: xmin, xmax, ymin, ymax

        if ( .not. associated(Qtr) ) Allocate(Qtr, Stat=istat)

        xmin = minval(points(:,1))
        ymin = minval(points(:,2))

        xmax = maxval(points(:,1))
        ymax = maxval(points(:,2))

        allocate( Qtr%vertices(4,2) )
        Qtr%vertices(1,1:2) = point(xmin,ymin)
        Qtr%vertices(2,1:2) = point(xmax,ymin)
        Qtr%vertices(3,1:2) = point(xmax,ymax)
        Qtr%vertices(4,1:2) = point(xmin,ymax)

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



