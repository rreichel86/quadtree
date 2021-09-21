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
    ! write out the points stored in each QtreeNode
    subroutine printPoints_(this,iow)
        use Qtree_data, only: num_node
        implicit none
        class(QtreeList) :: this
        integer :: iow
        type(QtreeNode) , pointer :: Qnode

        integer i, numSeeds, numIntrscPoints

        num_node = 0
        Qnode => this%HEAD
        do
            if ( .not. associated(Qnode) ) exit
            ! print Quad's boundary nodes
            do i = 1, 4
                write(iow,'(2f32.16)') Qnode%Q%Boundary(i)
            end do
            num_node = num_node + 4
            ! print intrsc_points
            numIntrscPoints = Qnode%Q%num_intrsc_points
            if (numIntrscPoints .ne. 0) then
                do i = 1, numIntrscPoints
                    write(iow,'(2f32.16)') Qnode%Q%intrsc_points(i)%pos
                end do
                num_node = num_node + numIntrscPoints
            end if
            ! print seeding_points
            if ( allocated(Qnode%Q%seeds) ) then
                numSeeds = size(Qnode%Q%seeds)
                do i = 1,  numSeeds
                    write(iow,'(2f32.16)') Qnode%Q%seeds(i)%pos
                end do
                num_node = num_node + numSeeds
            end if
            Qnode => Qnode%next
        end do
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

    subroutine searchPoint(Qtr,pt,CQ)
        use point_module

        implicit none
        Type(Qtree), pointer :: Qtr, CQ
        Type(point), intent(in) :: pt

        CQ => Qtr
        do
            if (CQ%level .eq. 18) return
            if ( .not. associated(CQ%NW) ) exit
            if ( CQ%NW%containsPoint_(pt) ) then
                CQ => CQ%NW
            else if ( CQ%SW%containsPoint_(pt) ) then
                CQ => CQ%SW
            else if ( CQ%NE%containsPoint_(pt) ) then
                CQ => CQ%NE
            else if ( CQ%SE%containsPoint_(pt) ) then
                CQ => CQ%SE
            end if
        end do
    end subroutine

    subroutine QtrRefine(Qtr,seed)
        use point_module

        implicit none
        Type(Qtree), pointer :: Qtr, Q
        Type(point), intent(in) :: seed

        Q => Qtr
        do
            if (Q%level .eq. 18) return
            if (Q%level .eq. 0) call subdivideQ (Q)
            if ( .not. associated(Q%NW) ) then
                call subdivideQ (Q)
                exit
            end if
            if ( Q%NW%containsPoint_(seed) ) then
                Q => Q%NW
            else if ( Q%SW%containsPoint_(seed) ) then
                Q => Q%SW
            else if ( Q%NE%containsPoint_(seed) ) then
                Q => Q%NE
            else if ( Q%SE%containsPoint_(seed) ) then
                Q => Q%SE
            end if
        end do
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
    
    subroutine connect_h_neighbor(ctrl,Qtr,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
        use Qtree_data
        use point_module
        
        implicit none
        
        integer, intent(in) :: ctrl
        Type (Qtree), pointer :: Qtr, NQ, AQ
        integer, intent(in) ::  itr_nro,  level, ref(2*level)
        integer, intent(inout) :: zhl, wpoly(14), signo(14)
        type (point), intent(inout) :: temp_coor(14)
       
        logical :: existNQRef, QOchildren(4)
        logical :: cond0(2), cond1(2), cond2(2)
        logical :: cond3
        integer :: i,signo_temp, wpoly_temp, istat, refNQ(2*level)
        type(point) :: pt_temp
        type(QtreePtr), pointer :: children(:)
        
        character*2 :: nb(4)
        integer :: dir(4)
        
        existNQRef = .false.
        nb = ['SS','EE','NN','WW']
        dir = [2,3,4,1]
        pt_temp = point(0d0,0d0)
        refNQ = ref

        call QtrRefNeighbourQ(level,refNQ,dir(itr_nro),0,existNQRef)
                                        
        if (existNQRef) then ! exist a ref?
                        
            ! ref -> Q
            call nearestCommonAncestor(Qtr,level,refNQ,AQ)
            call QtrSearchNeighbourQ(AQ,level,refNQ,NQ)

            QOchildren = [associated(NQ%NW),associated(NQ%SW),associated(NQ%NE),associated(NQ%SE)]
            
            if (count(QOchildren).gt.0) then
            
                Allocate(children(4), Stat=istat)
                do i=1,4
                    children(i)%Q => Null()
                end do
                                
                children(1)%Q => NQ%NW
                children(2)%Q => NQ%SW
                children(3)%Q => NQ%NE 
                children(4)%Q => NQ%SE
                
            end if
            
            cond0(1) =      isQ_in(NQ)
            cond0(2) = .not.isQ_in(NQ) .and. NQ%num_mat_sets .eq.1
            
            !if ( (NQ%level .eq. level) .and. cond0(1) ) then
            if ( (NQ%level .eq. level) .and. count(QOchildren).gt.0) then  ! has same level as Neighbor and has Children
                                    
                !if (count(QOchildren).eq.0) then  ! children?
                                     
                !else if (count(QOchildren).gt.0) then ! children?
                                        
                    if (nb(itr_nro) .eq. 'SS') then
                                        
                        cond1(1) = QOchildren(1) .and. children(1)%Q%signo(3).ge.0 .and. children(1)%Q%wpoly(3).ge.1
                        cond2(1) = QOchildren(3) .and. children(3)%Q%signo(4).ge.0 .and. children(3)%Q%wpoly(4).ge.1
                        
                        !condition point in \omega_#>1
                        cond1(2) = QOchildren(1) .and. children(1)%Q%signo(3).lt.0 .and. children(1)%Q%wpoly(3).gt.1
                        cond2(2) = QOchildren(3) .and. children(3)%Q%signo(4).lt.0 .and. children(3)%Q%wpoly(4).gt.1
                                            
                        if ( cond1(ctrl) ) then
                            pt_temp = children(1)%Q%Boundary(3)   
                            signo_temp = children(1)%Q%signo(3)
                            wpoly_temp = children(1)%Q%wpoly(3)                         
                        else if ( cond2(ctrl) ) then
                            pt_temp = children(3)%Q%Boundary(4)
                            signo_temp = children(3)%Q%signo(4)
                            wpoly_temp = children(3)%Q%wpoly(4)
                        end if 
                                            
                        cond3 = temp_coor(zhl)%x .gt. pt_temp%x
                                            
                    else if (nb(itr_nro) .eq. 'EE') then
                                            
                        cond1(1) = QOchildren(1) .and. children(1)%Q%signo(1).ge.0 .and. children(1)%Q%wpoly(1).ge.1
                        cond2(1) = QOchildren(2) .and. children(2)%Q%signo(4).ge.0 .and. children(2)%Q%wpoly(4).ge.1
                        
                         !condition point in \omega_#>1
                        cond1(2) = QOchildren(1) .and. children(1)%Q%signo(1).lt.0 .and. children(1)%Q%wpoly(1).gt.1
                        cond2(2) = QOchildren(2) .and. children(2)%Q%signo(4).lt.0 .and. children(2)%Q%wpoly(4).gt.1
                        
                                        
                        if ( cond1(ctrl) ) then
                            pt_temp = children(1)%Q%Boundary(1)
                            signo_temp = children(1)%Q%signo(1)
                            wpoly_temp = children(1)%Q%wpoly(1)
                                                 
                        else if ( cond2(ctrl) ) then
                            pt_temp = children(2)%Q%Boundary(4)
                            signo_temp = children(2)%Q%signo(4)
                            wpoly_temp = children(2)%Q%wpoly(4)
                        end if   
                                            
                        cond3 = temp_coor(zhl)%y .gt. pt_temp%y
                                                   
                    else if (nb(itr_nro) .eq. 'NN') then
                                            
                        cond1(1) = QOchildren(2) .and. children(2)%Q%signo(2).ge.0 .and. children(2)%Q%wpoly(2).ge.1
                        cond2(1) = QOchildren(4) .and. children(4)%Q%signo(1).ge.0 .and. children(4)%Q%wpoly(1).ge.1
                        
                         !condition point in \omega_#>1
                        cond1(2) = QOchildren(2) .and. children(2)%Q%signo(2).lt.0 .and. children(2)%Q%wpoly(2).gt.1
                        cond2(2) = QOchildren(4) .and. children(4)%Q%signo(1).lt.0 .and. children(4)%Q%wpoly(1).gt.1
                                        
                        if ( cond1(ctrl) ) then
                            pt_temp = children(2)%Q%Boundary(2)
                            signo_temp = children(2)%Q%signo(2)
                            wpoly_temp = children(2)%Q%wpoly(2)
                                                
                        else if ( cond2(ctrl) ) then
                            pt_temp = children(4)%Q%Boundary(1)
                            signo_temp = children(4)%Q%signo(1)
                            wpoly_temp = children(4)%Q%wpoly(1)
                        end if     
                                            
                        cond3 = temp_coor(zhl)%x .lt. pt_temp%x                               
                                            
                    else if (nb(itr_nro) .eq. 'WW') then
                                            
                        cond1(1) = QOchildren(3) .and. children(3)%Q%signo(2).ge.0 .and. children(3)%Q%wpoly(2).ge.1
                        cond2(1) = QOchildren(4) .and. children(4)%Q%signo(3).ge.0 .and. children(4)%Q%wpoly(3).ge.1
                        
                        !condition point in \omega_#>1
                        cond1(2) = QOchildren(3) .and. children(3)%Q%signo(2).lt.0 .and. children(3)%Q%wpoly(2).gt.1
                        cond2(2) = QOchildren(4) .and. children(4)%Q%signo(3).lt.0 .and. children(4)%Q%wpoly(3).gt.1
                                        
                        if ( cond1(ctrl) ) then
                            pt_temp = children(3)%Q%Boundary(2)
                            signo_temp = children(3)%Q%signo(2)
                            wpoly_temp = children(3)%Q%wpoly(2)
                                                 
                        else if ( cond2(ctrl) ) then 
                            pt_temp = children(4)%Q%Boundary(3)
                            signo_temp = children(4)%Q%signo(3)
                            wpoly_temp = children(4)%Q%wpoly(3)
                        end if 
                                            
                        cond3 = temp_coor(zhl)%y .lt. pt_temp%y
                                            
                    end if  
                                        
                    if ( cond1(ctrl) .or. cond2(ctrl) ) then
                                                
                        if (temp_coor(zhl) .eq. pt_temp) then
                        else if ( cond3 ) then
                            zhl = zhl + 1
                            temp_coor(zhl) = temp_coor(zhl-1)
                            temp_coor(zhl-1) = pt_temp
                            
                            signo(zhl) = signo(zhl-1)
                            signo(zhl-1) = signo_temp
                            
                            wpoly(zhl) = wpoly(zhl-1)
                            wpoly(zhl-1) = wpoly_temp
                            
                        else
                            zhl = zhl + 1
                            temp_coor(zhl) = pt_temp 
                            signo(zhl) = signo_temp
                            wpoly(zhl) = wpoly_temp
                        end if 
                                                                         
                    end if
                                        
                !end if ! children?
                                    
            else    ! same level as neighbor?
            end if  ! same level as neighbor?
                                       
        end if  ! exist ref?
        
        if ( associated(children) ) deallocate (children, stat=istat)
        
    end subroutine 
    
    
    subroutine connect_h_seeds(ctrl,Qtr,zhl,temp_coor,wpoly,signo)
        use Qtree_data
        use point_module
        use segment_module
        
        
        implicit none
        
        integer, intent(in) :: ctrl
        Type (Qtree), pointer :: Qtr
        integer, intent(inout) :: zhl, wpoly(14), signo(14)
        type (point), intent(inout) :: temp_coor(14)
        logical :: cond(3,2)
        integer :: i,j,a,b,k, nseeds, eng, c_b, c_f
        type(point) :: temp(14) !,pt_a, pt_b
        !type(segment) :: S_temp
        
        temp = point(0d0,0d0)
        
        nseeds = size(Qtr%seeds)   
        
        do j=1, nseeds
        do i=1, zhl
        
            if(i .lt. zhl) then
                a = i
                b = i+1
            else 
                a = i
                b = 1
            end if 
            
            if ( b .eq. zhl ) then
                c_f = 1
            else
                c_f = b+1
            end if
            
            if ( a .eq. 1 ) then
                c_b = zhl
            else
                c_b = a-1
            end if
            
            
                       
            cond(:,1)  = [signo(a) .eq. 0 .and. signo(b) .eq. 0 &
                     &, signo(a) .eq. 1 .and. signo(b) .eq. 0 &
                     &, signo(a) .eq. 0 .and. signo(b) .eq. 1]
            
            cond(:,2)  = [signo(a) .eq. 0 .and. signo(b) .eq. 0 &
                     &, (signo(a) .eq. -1 .and. wpoly(a) .gt. 1) .and. signo(b) .eq. 0 &
                     &, signo(a) .eq. 0 .and. (signo(b) .eq. -1 .and. wpoly(b) .gt. 1)]      
            
            
                    ! first: condition I
                    !( signo(a) .eq. 0 .and. signo(b) .eq. 0 ) .and.  ( wpoly(a) .eq. Qtr%seeds(j)%wpoly .and. wpoly(b) .eq. Qtr%seeds(j)%wpoly )
                    
            
                       
                  !pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(a), temp_coor(c_b)) )  --> ! delete a
                  !pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(b), temp_coor(c_f)) )  --> ! delete b
            
                     
                  ! first: condition II
                  ! pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(a), temp_coor(b)) )
                   
                   !signo(a) .eq. 0 
                  
                  !( signo(a) .eq. 0 .and. signo(c_b) .eq. 0 ) .and.  ( wpoly(a) .eq. Qtr%seeds(j)%wpoly .and. wpoly(c_b) .eq. Qtr%seeds(j)%wpoly )
                  ! delete a
                  !signo(b) .eq. 0 
                  
                  !( signo(b) .eq. 0 .and. signo(c_f) .eq. 0 ) .and.  ( wpoly(b) .eq. Qtr%seeds(j)%wpoly .and. wpoly(c_f) .eq. Qtr%seeds(j)%wpoly )
                  ! delete b
                  
                  
                  
            if ( count( cond(:,ctrl) ) .gt. 0 ) then         
                            
                !pt_a = temp_coor(a)
                !pt_b = temp_coor(b)
                !S_temp = segment(pt_a, pt_b)
                
                
                !do j=1, nseeds
                    eng = 0        
                    k = 0
                    if(i .lt. zhl)  temp(b:zhl) = temp_coor(b:zhl)
                                
                    if (Qtr%seeds(j)%wpoly == 0) cycle
                    if ( temp_coor(a) == Qtr%seeds(j)%pos .or. temp_coor(b) == Qtr%seeds(j)%pos) cycle
                           
                       ! -- condition I --
                        if ( (signo(a) .eq. 0 .and. signo(b) .eq. 0) &
                             & .and. (wpoly(a) .eq. Qtr%seeds(j)%wpoly .and. wpoly(b) .eq. Qtr%seeds(j)%wpoly) ) then
                                    
                            
                            if (pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(a), temp_coor(c_b)) )) then     !delete a
                                    
                                temp_coor(a) = Qtr%seeds(j)%pos
                                exit 
                                
                            else if (pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(b), temp_coor(c_f)) )) then  !delete b
                                
                                temp_coor(b) = Qtr%seeds(j)%pos
                                exit 
                                
                            else 
                                 k = k + 1
                                 temp_coor(a+k) = Qtr%seeds(j)%pos 
                                 eng = eng + 1
                            end if    
                                
                        ! -- condition II --    
                        else if (pt_in_segment( Qtr%seeds(j)%pos, segment(temp_coor(a), temp_coor(b)) )) then
                                       
                                
                             if ( (signo(a) .eq. 0 .and. signo(c_b) .eq. 0) &
                                  & .and. (wpoly(a) .eq. Qtr%seeds(j)%wpoly .and. wpoly(c_b) .eq. Qtr%seeds(j)%wpoly) ) then  !delete a
                                 
                                  temp_coor(a) = Qtr%seeds(j)%pos
                                  exit 
                             
                             else if ( (signo(b) .eq. 0 .and. signo(c_f) .eq. 0) &
                                       & .and. (wpoly(b) .eq. Qtr%seeds(j)%wpoly .and. wpoly(c_f) .eq. Qtr%seeds(j)%wpoly) ) then !delte b
                                      
                                   temp_coor(b) = Qtr%seeds(j)%pos    
                                   exit 
                             else 
                                  k = k + 1
                                  temp_coor(a+k) = Qtr%seeds(j)%pos    
                                  eng = eng + 1
                             end if     
                             
                                        
                        else
                        
                            cycle
                                        
                        end if 
                                    
                        if(i .lt. zhl) temp_coor(b+k:zhl+k) = temp(b:zhl)     
                        zhl = zhl + k
                !end do 
                            
            end if   
            if (eng .eq. 1) exit
            
        end do  
        end do
    
    
    
    end subroutine
    
    recursive subroutine connectivity(Qtr)

        use Qtree_data
        use Qtree_input 
        use SortSearch_module
        use polygon_module
        use segment_module
        use point_module
        
        implicit none
        
        Type (Qtree), pointer :: Qtr
        logical :: Qchildren(4), cc, cc_2
        integer :: itr_nro
        
        !var for Q*
        integer ::  zhl, signo(14), wpoly(14), temp_elem(14), signo_temp, wpoly_temp, nmat
        type(point) :: temp_coor(14), x0
        
        !vars for Q'
        integer ::  zhl_2, signo_2(14), wpoly_2(14), temp_elem_2(14), signo_temp_2, wpoly_temp_2, nmat_2
        type(point) :: temp_coor_2(14)
        
        integer ::  sch, level, ref(36)
        integer ::  Q(5), SAS, SES, SANR, SENR
        
        integer :: i, index, istat
        
        real*8 :: tol 
    
        tol = 1.d-10
 
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
         
        if (count(Qchildren).eq.0) then
                
         ! In condition 
         if ( isQ_in(Qtr) .and. Qtr%num_mat_sets .eq. 1 ) then
         
               temp_coor = point(0d0,0d0)
               temp_elem = 0
               
               signo = 0
                    
                sch = 0
                zhl = 1
                Q = [1,2,3,4,1]
                do itr_nro = 1,4
                
                        ref = Qtr%ref
                        level = Qtr%level
                    
                        SAS = Qtr%signo( Q(itr_nro) )
                        SES = Qtr%signo( Q(itr_nro+1) )
                        SANR = Qtr%wpoly( Q(itr_nro) )
                        SENR = Qtr%wpoly( Q(itr_nro+1) )
                        
                        
                        if ( (SAS .ge. 0 .and. SANR .ge. 1) ) then
                            temp_coor(zhl) = Qtr%Boundary( Q(itr_nro) )  
                            signo(zhl) = SAS 
                            wpoly(zhl) = SANR
                        end if 
                        
                        if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) cycle
                        if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .eq. SENR)) ) cycle
                        
                        if ( ( SAS .eq. 1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) then    ! Case 1
                            
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                         
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. 1 .and. SENR .eq. 1) ) then  ! Case 1
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                        
                        else if ( ( SAS .ge. 0 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) ) then  ! Case 2
                            
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                            
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .ge. 0 .and. SENR .eq. 1) ) then  ! Case 2
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly    
                              
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) then  ! Case 3
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                        
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) ) then ! Case 3
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                        else if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR)) ) then ! Case 4
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                        
                        else if ( ( SAS .eq. 0 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR) ) then  ! Case 5
                            zhl = zhl + 1
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                            
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) .and. (SANR .ne. SENR) ) then  ! Case 5
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly      
                            
                        !else if ( ( SAS .eq. 1 .and. SANR .eq. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) ) then  ! Case 6
                        !    
                        !    if ( (Qtr%num_intrsc_points .ne. 0) .and. (sch .lt. Qtr%num_intrsc_points) ) then
                        !        zhl = zhl + 1
                        !        sch = sch + 1
                        !        temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                        !        signo(zhl) = 0
                        !        wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                        !    end if 
                        !else if ( ( SAS .eq. 0 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) ) then  ! Case 6.5   
                        !    if ( (Qtr%num_intrsc_points .ne. 0) .and. sch .ge. 1 ) then
                        !        zhl = zhl -1
                        !    end if     
                        !else if ( ( SAS .eq. 0 .and. SANR .gt. 1 )  .and. ( SES .eq. 1 .and. SENR .eq. 1) ) then  ! Case 6
                        !
                        !    if ( (Qtr%num_intrsc_points .ne. 0) .and. (sch .lt. Qtr%num_intrsc_points) ) then
                        !        sch = sch + 1
                        !            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                        !            signo(zhl) = 0
                        !            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly  
                        !    end if 
                            
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. 0 .and. SENR .eq. 1) ) then
                            zhl = zhl - 1 
                            
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) ) then
                            zhl = zhl - 1
                        
                        end if    
                           
                        call connect_h_neighbor(1,Qtr,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                        zhl = zhl + 1
                        
                end do
                
                zhl = zhl - 1
                counter = counter + 1
                
                if ( allocated(Qtr%seeds) ) call connect_h_seeds(1,Qtr,zhl,temp_coor,wpoly,signo)
                
                !------------
                do i = 1, zhl
                    call binarySearch (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
                end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
                
                nmat = Qtr%mat_nros(1)
                
                if (elm_typ_ma(zhl,nmat).eq.0)  then
                    !mate_zhl = mate_zhl + 1
                    elm_typ_ma(zhl,nmat) = 1
                end if 
                
                elements(counter,1) = zhl 
                elements(counter,2) = nmat
                elements(counter,3:16) = temp_elem
                
                write(55,'(17i6)') counter, zhl, nmat, temp_elem
                 !------------
                
         else if ( isQ_in(Qtr) .and. Qtr%num_mat_sets .gt.1 ) then  
                
               !Q*
               temp_coor = point(0d0,0d0)
               temp_elem = 0
               wpoly = 0
               signo = 0
               zhl = 1
               cc = .true.
                
               !Q'
               temp_coor_2 = point(0d0,0d0)
               temp_elem_2 = 0
               wpoly_2 = 0
               signo_2 = 0
               zhl_2 = 1
               cc_2 = .true.
               
               sch = 0
                
                Q = [1,2,3,4,1]
                do itr_nro = 1,4
                
                        ref = Qtr%ref
                        level = Qtr%level
                    
                        SAS = Qtr%signo( Q(itr_nro) )
                        SES = Qtr%signo( Q(itr_nro+1) )
                        SANR = Qtr%wpoly( Q(itr_nro) )
                        SENR = Qtr%wpoly( Q(itr_nro+1) )
                        
                        
                        if ( (SAS .ge. 0 .and. SANR .ge. 1) ) then
                            temp_coor(zhl) = Qtr%Boundary( Q(itr_nro) )  
                            signo(zhl) = SAS 
                            wpoly(zhl) = SANR
                        end if 
                        
                        ! out*
                        ! on*
                        if ( (SAS .le. 0 .and. SANR .gt. 1) ) then
                            temp_coor_2(zhl_2) = Qtr%Boundary( Q(itr_nro) )  
                            signo_2(zhl_2) = SAS 
                            wpoly_2(zhl_2) = SANR
                        end if
                        
                        if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) cycle
                        !if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .eq. SENR)) ) cycle
                        
                        ! in - out  
                        if ( ( SAS .eq. 1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) then    ! Case 1
                            
                            sch = sch + 1
                            
                            zhl = zhl + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl_2 = zhl_2 - 1
                            cc_2 = .false.
                            
                         
                        ! out - in  
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. 1 .and. SENR .eq. 1) ) then  ! Case 1
                        
                            sch = sch + 1
                            
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl_2 = zhl_2 - 1
                            cc_2 = .false.
                            
                        ! in - out*  
                        else if ( ( SAS .ge. 0 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) ) then  ! Case 2
                            
                            
                            sch = sch + 1
                            
                            zhl = zhl + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                            
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
                            
                        ! out* - in     
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .ge. 0 .and. SENR .eq. 1) ) then  ! Case 2
                        
                            sch = sch + 1
                            
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly    
                            
                            zhl_2 = zhl_2 + 1
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
                         
                        ! out* - out   
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1) ) then  ! Case 3
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            zhl_2 = zhl_2 + 1
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
                            
                            
                            sch = sch + 1
                            zhl = zhl + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                        ! out - out*
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) ) then ! Case 3
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            sch = sch + 1
                            zhl = zhl + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
                            
                        ! out* - out**    ???
                        else if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR)) ) then ! Case 4
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                            
                            sch = sch + 1
                            zhl = zhl + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly
                        
                        ! on** - out*  
                        else if ( ( SAS .eq. 0 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR) ) then  ! Case 5
                            
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                            
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
 
                            
                        ! out* - on**    
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) .and. (SANR .ne. SENR) ) then  ! Case 5
                        
                            sch = sch + 1
                            temp_coor(zhl) = Qtr%intrsc_points(sch)%pos
                            signo(zhl) = 0
                            wpoly(zhl) = Qtr%intrsc_points(sch)%wpoly   
                            
                            zhl_2 = zhl_2 + 1
                            temp_coor_2(zhl_2) = Qtr%intrsc_points(sch)%pos
                            signo_2(zhl_2) = 0
                            wpoly_2(zhl_2) = Qtr%intrsc_points(sch)%wpoly
                            
                        ! in - in  
                        ! on - in   
                        ! in - on
                        ! on - on    
                        else if ( ( SAS .ge. 0 .and. SANR .eq. 1 )  .and. ( SES .ge. 0 .and. SENR .eq. 1) ) then
                        
                            zhl_2 = zhl_2 - 1
                            cc_2 = .false.
                          
                        ! in - on*    
                        else if ( ( SAS .ge. 0 .and. SANR .eq. 1 )  .and. ( SES .ge. 0 .and. SENR .gt. 1) ) then
                        
                            zhl_2 = zhl_2 - 1
                            cc_2 = .false.
                          
                         ! on* - in   
                        else if ( ( SAS .ge. 0 .and. SANR .gt. 1 )  .and. ( SES .ge. 0 .and. SENR .eq. 1) ) then
                        
                            !zhl_2 = zhl_2 - 1
                            cc_2 = .false.
                                     
                         ! out - on
                        else if ( ( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. 0 .and. SENR .eq. 1) ) then
                            zhl = zhl - 1 
                            zhl_2 = zhl_2 - 1
                            
                        ! out* - on*    
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) ) then
                            zhl = zhl - 1
                            cc = .false.
                            !zhl_2 = zhl_2 - 1
                            
                        ! out* - out*   
                        else if ( ( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .eq. SENR) ) then
                             cc = .false.
                             zhl = zhl - 1
                            !zhl_2 = zhl_2 - 1
                        end if    
                        
                           
                        if (cc) call connect_h_neighbor(1,Qtr,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                         zhl = zhl + 1
                        cc = .true.
                        
                        if (cc_2) call connect_h_neighbor(2,Qtr,itr_nro,ref,level,zhl_2,temp_coor_2,wpoly_2,signo_2) ! ??
                        zhl_2 = zhl_2 + 1 
                        cc_2 = .true.
                        
                end do
                
                zhl = zhl - 1
                counter = counter + 1
                
                zhl_2 = zhl_2 - 1
                
                if ( allocated(Qtr%seeds) ) call connect_h_seeds(1,Qtr,zhl,temp_coor,wpoly,signo)
                
                !------------
                do i = 1, zhl
                    call binarySearch (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
                end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
                
                nmat = Qtr%mat_nros(1)
                
                if (elm_typ_ma(zhl,nmat).eq.0)  then
                    !mate_zhl = mate_zhl + 1
                    elm_typ_ma(zhl,nmat) = 1
                end if 
                
                elements(counter,1) = zhl 
                elements(counter,2) = nmat
                elements(counter,3:16) = temp_elem

                write(55,'(17i6)') counter, zhl, nmat, temp_elem
                !------------
                
                counter = counter + 1
                
                if ( allocated(Qtr%seeds) ) call connect_h_seeds(2,Qtr,zhl_2,temp_coor_2,wpoly_2,signo_2)
                
                !------------
                do i = 1, zhl_2
                    call binarySearch (nodes(1:num_node), temp_coor_2(i), index)
                    temp_elem_2(i) = index
                    if (index == 0) write(*,*) counter, temp_coor_2(i), zhl_2
                end do 
                
                call polykernel(zhl_2,temp_coor_2(1:zhl_2),x0)
                nodes(num_node + counter) = x0
                
                nmat = Qtr%mat_nros(2)
                
                if (elm_typ_ma(zhl_2,nmat).eq.0)  then
                    !mate_zhl = mate_zhl + 1
                    elm_typ_ma(zhl_2,nmat) = 1
                end if 
                
                elements(counter,1) = zhl_2 
                elements(counter,2) = nmat
                elements(counter,3:16) = temp_elem_2
                
                write(55,'(17i6)') counter, zhl_2, nmat, temp_elem_2
                !------------
                
         else if ( (.not.isQ_in(Qtr) .and. Qtr%num_mat_sets .eq.1) )then
                
               temp_coor = point(0d0,0d0)
               temp_elem = 0
               
               signo = 0
                    
               !sch = 0
               zhl = 1
        
               do itr_nro = 1, 4
                 
                ref = Qtr%ref
                level = Qtr%level
                
                temp_coor(zhl) = Qtr%Boundary(itr_nro)
                signo(zhl) = Qtr%signo(itr_nro) 
                wpoly(zhl) = Qtr%wpoly(itr_nro)
                            
                call connect_h_neighbor(2,Qtr,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                zhl = zhl + 1
               
               end do
               
               zhl = zhl - 1
               counter = counter + 1
               
               !------------
               do i = 1, zhl
                    call binarySearch (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
               end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
               
                nmat = Qtr%mat_nros(1)
                
               if (elm_typ_ma(zhl,nmat).eq.0)  then
                   !mate_zhl = mate_zhl + 1
                    elm_typ_ma(zhl,nmat) = 1
                end if 
                
                elements(counter,1) = zhl 
                elements(counter,2) = nmat
                elements(counter,3:16) = temp_elem
                
                write(55,'(17i6)') counter, zhl, nmat, temp_elem
                !--------------
                
         end if 
        
            
        else if (count(Qchildren).gt.0) then
        
            if (QChildren(1)) call connectivity(Qtr%NW)
            if (QChildren(2)) call connectivity(Qtr%SW)
            if (QChildren(3)) call connectivity(Qtr%NE)
            if (QChildren(4)) call connectivity(Qtr%SE)
            
        end if 
   return 
1000    format(i6','i6','i6','i6,','i6','i6','i6','i6,','i6','i6','i6','i6,','i6','i6','i6','i6)
        
    end subroutine
    
    
   recursive subroutine QIntrsPts(Qtr,num_poly,polygons)
        use Qtree_data
        !use  Qtree_input
        use SortSearch_module
        use polygon_module
        use segment_module
        
        implicit none
        
        Type (Qtree), pointer :: Qtr
        integer, intent(in) :: num_poly
        Type (polygon), intent(in) :: polygons(num_poly)
        Type (segment) :: SQ
        type (point) :: pt1, pt2
        logical :: Qchildren(4), intrsc1, intrsc2, status(4), collision
        integer :: i,j, zhl, signo, Q(5), ms, mwp, nclist, clist(4), wpoly, nr1, nr2, SAS, SES, SANR, SENR
        integer :: mat_nro, test, list(4), nlist 
        
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
        
        if ( count(QChildren).eq.0 ) then !children
        
            debug_zhl = debug_zhl + 1
            
            ! polygons w +1
            do i = 1, 4
                call point_in_polygonSR (polygons(1), Qtr%Boundary(i), Qtr%signo(i), status(i), .true.)
            end do 
            
            ! polygons w -1
            nclist = 0
            clist = 0
            
           
            do j = 2,num_poly
            
                collision = rectRectCollision(Qtr,polygons(j))
                if (.not. collision) cycle
                nclist = nclist + 1
                clist(nclist) = j
                
           end do
           
            status = .true.
            nlist = 0
            list = 0
            
           do j = 1, nclist
                wpoly = clist(j)
                do i = 1, 4
                    if ( (Qtr%signo(i) .eq. -1 .or. Qtr%signo(i) .eq. 0) .and. (Qtr%wpoly(i) .ge.1) ) cycle 
                    if ( ( Qtr%signo(i) .eq. 1 ) .and. ( Qtr%wpoly(i) .eq.1 ) ) then
                    
                        call point_in_polygonSR (polygons(wpoly), Qtr%Boundary(i), Qtr%signo(i), status(i), .false.)
                    
                        if (Qtr%signo(i) .eq. -1) then 
                            Qtr%wpoly(i) = wpoly
                        else if (Qtr%signo(i) .eq. 0) then
                            Qtr%wpoly(i) = wpoly
                        end if
                        
                    end if     
                end do  
                test = count(status)
                if (count(status) .le. 4) then
                    nlist = nlist + 1
                    list(nlist) = clist(j)
                end if 
           end do
                
                ms = sum(Qtr%signo)
                mwp = sum(Qtr%wpoly)
                
                
                if ( (ms .eq.  4) .or.  (ms .eq. 3) ) then
                
                    Qtr%mat_nros(1) = polygons(1)%mat_nro
                    Qtr%num_mat_sets = 1
            
                ! condition to compute intrsc_points
                else if ( (ms .gt. -3) .and. (ms .lt. 3) &
                        & .or. ( ms .eq. -4 .and. mwp .gt.4 ) ) then

                    zhl = 0
                    Q = [1,2,3,4,1]
                    do j = 1, 4 
                    
                        SAS = Qtr%signo( Q(j) )
                        SES = Qtr%signo( Q(j+1) )
                        SANR = Qtr%wpoly( Q(j) )
                        SENR = Qtr%wpoly( Q(j+1) )
                        
                        ! CASE 1: in-out || out-in
                        if ( (( SAS .eq. 1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1))  & 
                            & .or. (( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. 1 .and. SENR .eq. 1))  ) then  
                            
                            SQ = segment( Qtr%Boundary( Q(j) ), Qtr%Boundary( Q(j+1) ) )
                            
                            if ( SAS .eq. -1 .and. SANR .eq. 1 )  then 
                                nr1 = SANR
                            else if ( SES .eq. -1 .and. SENR .eq. 1 ) then 
                                nr1 = SENR
                            end if     
                            
                            call intrsc_segment_polygon (polygons(nr1), SQ, pt1, intrsc1)
                           
                            if (intrsc1) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt1
                                Qtr%intrsc_points(zhl)%wpoly = nr1
                            end if   
                        
                        ! Case 2: in-out* || out*-in
                        !         on-out* || out*-on    
                        else if ( (( SAS .ge. 0 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1))  & 
                            & .or. (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .ge. 0 .and. SENR .eq. 1))  ) then
                            
                            SQ = segment( Qtr%Boundary( Q(j) ), Qtr%Boundary( Q(j+1) ) )
                            
                            if ( SAS .eq. -1 .and. SANR .gt. 1 )  then 
                                nr1 = SANR
                            else if ( SES .eq. -1 .and. SENR .gt. 1 ) then 
                                nr1 = SENR
                            end if     
                            
                            call intrsc_segment_polygon (polygons(nr1), SQ, pt1, intrsc1)
                           
                            if (intrsc1) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt1
                                Qtr%intrsc_points(zhl)%wpoly = nr1
                            end if 
                        
                        !Case 3: out-out* || out*-out
                                 
                        else if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .eq. 1))  & 
                            & .or. (( SAS .eq. -1 .and. SANR .eq. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1))  ) then
                            
                            SQ = segment( Qtr%Boundary( Q(j) ), Qtr%Boundary( Q(j+1) ) )   
                            
                            if ( SAS .eq. -1 .and. SANR .gt. 1 ) then
                                nr1 = SANR
                                nr2 = SENR
                            else if ( SAS .eq. -1 .and. SANR .eq. 1 ) then
                                nr1 = SANR
                                nr2 = SENR
                            end if 
                            
                            call intrsc_segment_polygon (polygons(nr1), SQ, pt1, intrsc1)
                            call intrsc_segment_polygon (polygons(nr2), SQ, pt2, intrsc2)


                            ! 8.01.2017   pt1 == pt2 ?
                            if (intrsc1) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt1
                                Qtr%intrsc_points(zhl)%wpoly = nr1
                            end if
                            
                            if (intrsc2) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt2
                                Qtr%intrsc_points(zhl)%wpoly = nr2
                            end if
                            
                        !Case 4: out**-out* || out*-out**    
                        else if ( (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR)) ) then
                            
                            SQ = segment( Qtr%Boundary( Q(j) ), Qtr%Boundary( Q(j+1) ) )   
                            
                            !if ( SAS .eq. -1 .and. SANR .gt. 1 ) then
                                nr1 = SANR
                                nr2 = SENR
                            !else if ( SAS .eq. -1 .and. SANR .eq. 1 ) then
                                !nr1 = SANR
                                !nr2 = SENR
                            !end if 
                            
                            call intrsc_segment_polygon (polygons(nr1), SQ, pt1, intrsc1)
                            call intrsc_segment_polygon (polygons(nr2), SQ, pt2, intrsc2)

                            ! 8.01.2017   pt1 == pt2 ?
                            if ( intrsc1 ) then
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt1
                                Qtr%intrsc_points(zhl)%wpoly = nr1
                            end if
                            
                            if ( intrsc2) then
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt2
                                Qtr%intrsc_points(zhl)%wpoly = nr2
                            end if
                            
                        ! Case 5: on*-out** || out**-on*
                        !         on**-out* || out*-on**    
                        else if ( (( SAS .eq. 0 .and. SANR .gt. 1 )  .and. ( SES .eq. -1 .and. SENR .gt. 1) .and. (SANR .ne. SENR))  & 
                            & .or. (( SAS .eq. -1 .and. SANR .gt. 1 )  .and. ( SES .eq. 0 .and. SENR .gt. 1) .and. (SANR .ne. SENR))  ) then
                            
                            SQ = segment( Qtr%Boundary( Q(j) ), Qtr%Boundary( Q(j+1) ) )
                            
                            if ( SAS .eq. -1 .and. SANR .gt. 1 )  then 
                                nr1 = SANR
                            else if ( SES .eq. -1 .and. SENR .gt. 1 ) then 
                                nr1 = SENR
                            end if     
                            
                            call intrsc_segment_polygon_2 (polygons(nr1), SQ, pt1, intrsc1)
                           
                            if (intrsc1) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt1
                                Qtr%intrsc_points(zhl)%wpoly = nr1
                            end if 
                         
                        !! Case 6: in - on* || on* - in
                        !else if ( (( sas .eq. 1 .and. sanr .eq. 1 )  .and. ( ses .eq. 0 .and. senr .gt. 1)) &
                        !    & .or. (( sas .eq. 0 .and. sanr .gt. 1 )  .and. ( ses .eq. 1 .and. senr .eq. 1))) then
                        !    
                        !    sq = segment( qtr%boundary( q(j) ), qtr%boundary( q(j+1) ) )
                        !    
                        !    if ( sas .eq. 0 .and. sanr .gt. 1 )  then 
                        !        nr1 = sanr
                        !    else if ( ses .eq. 0 .and. senr .gt. 1 ) then 
                        !        nr1 = senr
                        !    end if     
                        !    
                        !    call intrsc_segment_polygon_2 (polygons(nr1), sq, pt1, intrsc1)
                        !   
                        !    if (intrsc1 ) then  
                        !        zhl = zhl + 1
                        !        qtr%intrsc_points(zhl)%pos = pt1
                        !        qtr%intrsc_points(zhl)%wpoly = nr1
                        !        write(56,'(2f32.16)') pt1%x, pt1%y   
                        !    end if 
                            
                        end if 
                             
                    end do 
                
                    if (zhl .gt. 0) then
                        num_intrsc_pts = num_intrsc_pts + zhl
                        Qtr%num_intrsc_points = zhl
                    end  if      
                    
                    ! Condition more than 1 mat_set
                    
                    if ( (Qtr%num_intrsc_points .gt. 0) .and.  (sum(Qtr%wpoly) .eq. 4) ) then
                    
                        Qtr%mat_nros(1) = polygons(1)%mat_nro
                        Qtr%num_mat_sets = 1
                     
                    else if ( (Qtr%num_intrsc_points .gt. 0) .and.  (sum(Qtr%wpoly) .gt. 4) ) then
                        
                        Qtr%num_mat_sets = 1
                        Qtr%mat_nros(1) = polygons(1)%mat_nro
                        mat_nro = polygons(list(1))%mat_nro
                        
                        if(mat_nro .gt. 0) then
                            Qtr%num_mat_sets = 2
                            Qtr%mat_nros(2) = mat_nro
                        end if 
                    
                    else if ( .not. isQ_in(Qtr) ) then
                        
                        mat_nro = polygons(list(1))%mat_nro
                        
                        if(mat_nro .gt. 0) then
                            Qtr%num_mat_sets = 1
                            Qtr%mat_nros(1) = mat_nro
                        end if
                        
                    else if ( isQ_in(Qtr) ) then    
                    
                        Qtr%mat_nros(1) = polygons(1)%mat_nro
                        Qtr%num_mat_sets = 1   
                        
                    end if 
                    
                end if
            
        else if (count(QChildren).gt.0) then  ! no children 
            if (QChildren(1)) call QIntrsPts(Qtr%NW,num_poly,polygons)
            if (QChildren(2)) call QIntrsPts(Qtr%SW,num_poly,polygons)
            if (QChildren(3)) call QIntrsPts(Qtr%NE,num_poly,polygons)
            if (QChildren(4)) call QIntrsPts(Qtr%SE,num_poly,polygons)
        end if 
    end subroutine 
     
    recursive subroutine coordinates(Qtr)    ! später umbenennen
        use Qtree_data
        use SortSearch_module
        use polygon_module
        
        implicit none
        
        Type (Qtree), pointer :: Qtr
        logical :: Qchildren(4)
        integer :: i, index  !, ms, mwp, mip
        
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
        if (count(QChildren).eq.0) then
            
            !ms = sum(Qtr%signo)
            !mwp = sum(Qtr%wpoly)
            !mip = Qtr%num_intrsc_points
            
            ! in Condition
            !if ( ms .gt. -2 .or. (ms .eq. -2 .and. mip .ne. 0) &
            !   & .or. ( ms .eq. -4 .and. mwp .gt.4 .and. mip .ne. 0 ) ) then
            
            if( isQ_in(Qtr) )then   
                
                num_elem = num_elem + 1
                
                    do i = 1, 4
                         if(Qtr%signo(i).ge.0) then
                            call binarySearch (Temp_nodes, Qtr%Boundary(i), index)
                            nodes_mask(index) = .true. 
                        end if     
                    end do 
                    
                    if(Qtr%num_intrsc_points .gt. 0) then
                    
                        do i = 1, Qtr%num_intrsc_points
                            call binarySearch (Temp_nodes, Qtr%intrsc_points(i)%pos, index)
                            nodes_mask(index) = .true.     
                        end do
                            
                    end if
            end if 
            
            if ( (isQ_in(Qtr) .and. Qtr%num_mat_sets .gt.1)  .or. (.not.isQ_in(Qtr) .and. Qtr%num_mat_sets .eq.1)  ) then
            
                num_elem = num_elem + 1
                
                     !do j = 1, nwp
                        do i = 1, 4
                    
                            if ( Qtr%signo (i) .gt. 0) cycle
                            if ( (Qtr%signo(i) .eq. -1 .and. Qtr%wpoly(i) .gt. 1) ) then
                                
                                call binarySearch (Temp_nodes, Qtr%Boundary(i), index)
                                nodes_mask(index) = .true. 
                                
                            else if ( (Qtr%signo(i) .eq. 0 .and. Qtr%wpoly(i) .gt. 1) ) then
                                
                                call binarySearch (Temp_nodes, Qtr%Boundary(i), index)
                                nodes_mask(index) = .true. 
                                
                            end if   
                           
                        end do    
                     !end do   
            end if
            
        else if (count(QChildren).gt.0) then
            if (QChildren(1)) call coordinates(Qtr%NW)
            if (QChildren(2)) call coordinates(Qtr%SW)
            if (QChildren(3)) call coordinates(Qtr%NE)
            if (QChildren(4)) call coordinates(Qtr%SE)
        end if 
    end subroutine

    recursive subroutine filter3(root,Qtr)    ! später umbenennen
        use Qtree_data
        use SortSearch_module
        use polygon_module
        
        implicit none
        
        Type (Qtree), pointer :: root, Qtr
        !integer, intent(in) :: nwp, wpolys(nwp)
        logical :: Qchildren(4)
        integer :: i,j, index
        
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
        if (count(QChildren).eq.0) then
            
            ! Condition
            if ( (isQ_in(Qtr) .and. Qtr%num_mat_sets .gt.1)  .or. (.not.isQ_in(Qtr) .and. Qtr%num_mat_sets .eq.1)  ) then
            
            !if ( (.not.isQ_in(Qtr) .and. Qtr%num_mat_sets .eq.1)  ) 
            num_elem = num_elem + 1
                
                     !do j = 1, nwp
                        do i = 1, 4
                    
                            if ( Qtr%signo (i) .ge. 0) cycle
                            if ( Qtr%signo(i) .eq. -1 .and. Qtr%wpoly(i) .gt. 1) then
                                call binarySearch (Temp_nodes, Qtr%Boundary(i), index)
                                nodes_mask(index) = .true. 
                            end if  
                        end do    
                     !end do   
            end if
        else if (count(QChildren).gt.0) then
            if (QChildren(1)) call filter3(root,Qtr%NW)
            if (QChildren(2)) call filter3(root,Qtr%SW)
            if (QChildren(3)) call filter3(root,Qtr%NE)
            if (QChildren(4)) call filter3(root,Qtr%SE)
        end if 
    end subroutine
    
    subroutine QtrReset(QtrList)
        implicit none
        type(QtreeList), pointer :: QtrList
        type(Qtree), pointer :: Q


        do
            if ( QtrList%isEmpty_() ) exit
            call QtrList%pop_(Q)
            
            Q%num_intrsc_points = 0
            Q%signo = 1
            Q%wpoly = 1
            Q%num_mat_sets = 0
            Q%mat_nros = [0,0,0]


        end do
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

    recursive subroutine QtrSubdivide(Qtr,level_min,n,seeds,max_seed_Q)
        use point_module
        use seed_point_module
        use Qtree_data

        implicit none
        Type (Qtree), pointer :: Qtr
        integer, intent(in) :: level_min
        integer, intent(in), optional  :: n
        type(seed_point), intent(in), optional :: seeds(*)
        integer, intent(in), optional :: max_seed_Q

        real*8 :: xmin, xmax, ymin, ymax, xmid, ymid
        real*8 :: boundary(8)
        integer :: i, istat, level, azhl, numseeds
        logical :: containsPoint

        if (Qtr%level .ge. 18) return
        if (.not. associated(Qtr%NW)) then

            Allocate (Qtr%NW, Qtr%SW, Qtr%NE, Qtr%SE, Stat=istat)

            xmin = minval(Qtr%Boundary(1:4)%x)
            xmax = maxval(Qtr%Boundary(1:4)%x)
            ymin = minval(Qtr%Boundary(1:4)%y)
            ymax = maxval(Qtr%Boundary(1:4)%y)
            xmid = sum(Qtr%Boundary(1:4)%x)/4.d0
            ymid = sum(Qtr%Boundary(1:4)%y)/4.d0

            level = Qtr%level + 1

            ! NW Child
            Qtr%NW%level = level
            Qtr%NW%ref = Qtr%ref
            Qtr%NW%ref(2*level-1)= 1
            Qtr%NW%ref(2*level)= 1
            Qtr%NW%father => Qtr
            boundary= [xmin,ymid,xmid,ymid,xmid,ymax,xmin,ymax]

            do i = 1,4
                Qtr%NW%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do

            ! SW Child
            Qtr%SW%level = level
            Qtr%SW%ref = Qtr%ref
            Qtr%SW%ref(2*level-1)= 1
            Qtr%SW%ref(2*level)= 2
            Qtr%SW%father => Qtr
            boundary= [xmin,ymin,xmid,ymin,xmid,ymid,xmin,ymid]

            do i = 1,4
                Qtr%SW%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do

            ! NE Child
            Qtr%NE%level = level
            Qtr%NE%ref = Qtr%ref
            Qtr%NE%ref(2*level-1)= 2
            Qtr%NE%ref(2*level)= 1
            Qtr%NE%father => Qtr
            boundary= [xmid,ymid,xmax,ymid,xmax,ymax,xmid,ymax]

            do i = 1,4
                Qtr%NE%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do

            ! SE Child
            Qtr%SE%level = level
            Qtr%SE%ref = Qtr%ref
            Qtr%SE%ref(2*level-1)= 2
            Qtr%SE%ref(2*level)= 2
            Qtr%SE%father => Qtr
            boundary= [xmid,ymin,xmax,ymin,xmax,ymid,xmid,ymid]

            do i = 1,4
                Qtr%SE%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do

            if (present(n) .and. present(seeds).and. present(max_seed_Q)) then

                call HowManySeeds(Qtr%NW%Boundary,n,seeds,azhl)
                if (azhl.gt.0) then
                    allocate(Qtr%NW%seeds(azhl))
                    call whichSeeds(Qtr%NW%Boundary,n,seeds,azhl,Qtr%NW%seeds)

                    if(azhl.eq.2) then
                        if (Qtr%NW%seeds(1)%pos == Qtr%NW%seeds(2)%pos) then 
                            azhl = azhl - 1
                        end if 
                    end if

                    if((azhl.gt.max_seed_Q).or.(Qtr%NW%level.lt.level_min)) then
                        call QtrSubdivide(Qtr%NW,level_min,azhl,Qtr%NW%seeds,max_seed_Q)
                        if (allocated (Qtr%NW%seeds).and.(azhl.gt.max_seed_Q)) deallocate (Qtr%NW%seeds)
                    end if

                else if ((azhl.eq.0).and.(Qtr%NW%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%NW,level_min)
                end if


                call HowManySeeds(Qtr%SW%Boundary,n,seeds,azhl)
                if (azhl.gt.0) then
                    allocate(Qtr%SW%seeds(azhl))
                    call whichSeeds(Qtr%SW%Boundary,n,seeds,azhl,Qtr%SW%seeds)

                    if(azhl.eq.2) then
                        if (Qtr%SW%seeds(1)%pos == Qtr%SW%seeds(2)%pos) then
                            azhl = azhl - 1
                        end if
                    end if

                    if((azhl.gt.max_seed_Q).or.(Qtr%SW%level.lt.level_min)) then
                        call QtrSubdivide(Qtr%SW,level_min,azhl,Qtr%SW%seeds,max_seed_Q)
                        if (allocated (Qtr%SW%seeds).and.(azhl.gt.max_seed_Q)) deallocate (Qtr%SW%seeds)
                    end if

                else if ((azhl.eq.0).and.(Qtr%SW%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%SW,level_min)
                end if

                call HowManySeeds(Qtr%NE%Boundary,n,seeds,azhl)
                if (azhl.gt.0) then
                    allocate(Qtr%NE%seeds(azhl))
                    call whichSeeds(Qtr%NE%Boundary,n,seeds,azhl,Qtr%NE%seeds)

                    if(azhl.eq.2) then
                        if (Qtr%NE%seeds(1)%pos == Qtr%NE%seeds(2)%pos) then 
                            azhl = azhl - 1
                        end if
                    end if

                    if((azhl.gt.max_seed_Q).or.(Qtr%NE%level.lt.level_min)) then
                        call QtrSubdivide(Qtr%NE,level_min,azhl,Qtr%NE%seeds,max_seed_Q)
                        if (allocated (Qtr%NE%seeds).and.(azhl.gt.max_seed_Q)) deallocate (Qtr%NE%seeds)
                    end if

                else if ((azhl.eq.0).and.(Qtr%NE%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%NE,level_min)
                end if 

                call HowManySeeds(Qtr%SE%Boundary,n,seeds,azhl)
                if (azhl.gt.0) then
                    allocate(Qtr%SE%seeds(azhl))
                    call whichSeeds(Qtr%SE%Boundary,n,seeds,azhl,Qtr%SE%seeds)

                    if(azhl.eq.2) then
                        if (Qtr%SE%seeds(1)%pos == Qtr%SE%seeds(2)%pos) then
                            azhl = azhl - 1
                        end if 
                    end if 

                    if((azhl.gt.max_seed_Q).or.(Qtr%SE%level.lt.level_min)) then
                        call QtrSubdivide(Qtr%SE,level_min,azhl,Qtr%SE%seeds,max_seed_Q)
                        if (allocated (Qtr%SE%seeds).and.(azhl.gt.max_seed_Q)) deallocate (Qtr%SE%seeds) 
                    end if

                else if ((azhl.eq.0).and.(Qtr%SE%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%SE,level_min)
                end if 

            else 

                if (allocated (Qtr%seeds)) then

                    numseeds = size(Qtr%seeds)

                    call HowManySeeds(Qtr%NW%Boundary,numseeds,Qtr%seeds,azhl)
                    if (azhl.gt.0) then
                        allocate(Qtr%NW%seeds(azhl))
                        call whichSeeds(Qtr%NW%Boundary,numseeds,Qtr%seeds,azhl,Qtr%NW%seeds)
                    end if

                    call HowManySeeds(Qtr%SW%Boundary,numseeds,Qtr%seeds,azhl)
                    if (azhl .gt.0) then
                        allocate(Qtr%SW%seeds(azhl))
                        call whichSeeds(Qtr%SW%Boundary,numseeds,Qtr%seeds,azhl,Qtr%SW%seeds)
                    end if

                    call HowManySeeds(Qtr%NE%Boundary,numseeds,Qtr%seeds,azhl)
                    if (azhl .gt.0) then
                        allocate(Qtr%NE%seeds(azhl))
                        call whichSeeds(Qtr%NE%Boundary,numseeds,Qtr%seeds,azhl,Qtr%NE%seeds)
                    end if

                    call HowManySeeds(Qtr%SE%Boundary,numseeds,Qtr%seeds,azhl)
                    if (azhl .gt.0) then
                        allocate(Qtr%SE%seeds(azhl))
                        call whichSeeds(Qtr%SE%Boundary,numseeds,Qtr%seeds,azhl,Qtr%SE%seeds)
                    end if


                    deallocate (Qtr%seeds) 


                end if



                if(Qtr%NW%level.lt.level_min) then
                    call QtrSubdivide(Qtr%NW,level_min)
                end if

                if(Qtr%SW%level.lt.level_min) then
                    call QtrSubdivide(Qtr%SW,level_min)
                end if

                if(Qtr%NE%level.lt.level_min) then
                    call QtrSubdivide(Qtr%NE,level_min)
                end if

                if(Qtr%SE%level.lt.level_min) then
                    call QtrSubdivide(Qtr%SE,level_min)
                end if

            end if

        else if (associated(Qtr%NW)) then

            if (present(n) .and. present(seeds).and. present(max_seed_Q)) then
                call HowManySeeds(Qtr%NW%Boundary,n,seeds,azhl)
                if((azhl.gt.max_seed_Q).or.(Qtr%NW%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%NW,level_min,n,seeds,max_seed_Q)
                end if

                call HowManySeeds(Qtr%SW%Boundary,n,seeds,azhl)
                if((azhl.gt.max_seed_Q).or.(Qtr%SW%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%SW,level_min,n,seeds,max_seed_Q)
                end if

                call HowManySeeds(Qtr%NE%Boundary,n,seeds,azhl)
                if((azhl.gt.max_seed_Q).or.(Qtr%NE%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%NE,level_min,n,seeds,max_seed_Q)
                end if

                call HowManySeeds(Qtr%SE%Boundary,n,seeds,azhl)
                if((azhl.gt.max_seed_Q).or.(Qtr%SE%level.lt.level_min)) then
                    call QtrSubdivide(Qtr%SE,level_min,n,seeds,max_seed_Q)
                end if
            else
                if(Qtr%NW%level.lt.level_min) then
                    call QtrSubdivide(Qtr%NW,level_min)
                end if

                if(Qtr%SW%level.lt.level_min) then
                    call QtrSubdivide(Qtr%SW,level_min)
                end if

                if(Qtr%NE%level.lt.level_min) then
                    call QtrSubdivide(Qtr%NE,level_min)
                end if

                if(Qtr%SE%level.lt.level_min) then
                    call QtrSubdivide(Qtr%SE,level_min)
                end if
            end if

        end if

    end subroutine


    logical function rectRectCollision(Qtr,polygon_1)
    use polygon_module
    !use Qtree_module

    implicit none
    type(Qtree), pointer, intent(in) :: Qtr
    type(polygon), intent(in) :: polygon_1

    real*8 :: Q_xmin, Q_xmax, Q_ymin, Q_ymax
    real*8 :: P_xmin, P_xmax, P_ymin, P_ymax


    Q_xmin =  Qtr%Boundary(1)%x
    Q_ymin =  Qtr%Boundary(1)%y
    Q_xmax =  Qtr%Boundary(3)%x
    Q_ymax =  Qtr%Boundary(3)%y

    P_xmin =  polygon_1%xmin
    P_ymin =  polygon_1%ymin
    P_xmax =  polygon_1%xmax
    P_ymax =  polygon_1%ymax

    if ((Q_ymin .gt. P_ymax .or. Q_ymax .lt. P_ymin) .or. (Q_xmin .gt. P_xmax .or. Q_xmax .lt. P_xmin)) then
        rectRectCollision = .false.
        return
    else
        rectRectCollision = .true.
        return
    end if
    end function

end module



