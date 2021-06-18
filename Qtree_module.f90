
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
            
            deallocate(seeds, stat=istat)
            
            do i = 1, np
                call polygon_delete(polygons(i))          
            end do 
            
            deallocate(polygons, stat=istat)
            
        end subroutine
    
    
    
    
end module


module  Qtree_module
use point_module
use seed_point_module
use intrsc_point_module

type Qtree 

    type(point) :: Boundary(4) 
    integer, allocatable :: elem(:)
    type(seed_point), allocatable :: seeds(:)
    type(point) :: center
    type(intrsc_point) :: intrsc_points(4)
    integer :: num_intrsc_points = 0
    integer :: level = 0 
    integer :: ref(36)
    integer :: signo(4) = 1
    integer :: wpoly(4) = 1
    logical :: status(4) = .true.   ! wird abgeschaft RR
    integer :: num_mat_sets = 0
    integer :: mat_nros(3)  = [0,0,0]
    type (Qtree), pointer :: NW => null()
    type (Qtree), pointer :: SW => null()
    type (Qtree), pointer :: NE => null()
    type (Qtree), pointer :: SE => null()
    !type (Qtrees), pointer :: children(:)
    
end type

type Qtrees
    type (Qtree), Pointer :: Q
end type

contains
    
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
                
               if(allocated(Qtr%elem)) deallocate(Qtr%elem, stat=istat)
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
        
        nullify(Qtr)
        Allocate(Qtr, Stat=istat)
        !Allocate(Qtr%children(4), Stat=istat)
        !do i=1,4
        !    Qtr%children(i)%Q => Null()
        !end do
        
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
    
    subroutine connect_h_neighbor(ctrl,root,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
        use Qtree_data
        use point_module
        
        implicit none
        
        integer, intent(in) :: ctrl
        Type (Qtree), pointer :: root, Qtr, Qtr_out
        integer, intent(in) ::  itr_nro,  level, ref(2*level)
        integer, intent(inout) :: zhl, wpoly(14), signo(14)
        type (point), intent(inout) :: temp_coor(14)
       
        logical :: status, QOchildren(4)
        logical :: cond0(2), cond1(2), cond2(2)
        logical :: cond3
        integer :: i,signo_temp, wpoly_temp, istat, ref_tmp(2*level)
        type(point) :: pt_temp
        type(Qtrees), pointer :: children(:)
        
        character*2 :: nb(4)
        
        nb = ['SS','EE','NN','WW']
        pt_temp = point(0d0,0d0)
        ref_tmp = ref
        
        call QtrFindNeighborQ(level,ref_tmp,nb(itr_nro),status) ! Neighbors in a given direction -> ref
                                        
        if (status) then ! exist a ref?
                        
            ! ref -> Q
            call QtrSearchQ(root,level,ref_tmp,1,Qtr_out)
            
            QOchildren = [associated(Qtr_out%NW),associated(Qtr_out%SW),associated(Qtr_out%NE),associated(Qtr_out%SE)]
            
            if (count(QOchildren).gt.0) then
            
                Allocate(children(4), Stat=istat)
                do i=1,4
                    children(i)%Q => Null()
                end do
                                
                children(1)%Q => Qtr_out%NW
                children(2)%Q => Qtr_out%SW
                children(3)%Q => Qtr_out%NE 
                children(4)%Q => Qtr_out%SE
                
            end if                                     
            
            cond0(1) =      isQ_in(Qtr_out)
            cond0(2) = .not.isQ_in(Qtr_out) .and. Qtr_out%num_mat_sets .eq.1
            
            !if ( (Qtr_out%level .eq. level) .and. cond0(1) ) then
            if ( (Qtr_out%level .eq. level) .and. count(QOchildren).gt.0) then  ! has same level as Neighbor and has Children
                                    
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
    
    recursive subroutine connectivity(root,Qtr)

        use Qtree_data
        use Qtree_input 
        use SortSearch_module
        use polygon_module
        use segment_module
        use point_module
        
        implicit none
        
        Type (Qtree), pointer :: root, Qtr
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
                           
                        call connect_h_neighbor(1,root,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                        zhl = zhl + 1
                        
                end do
                
                zhl = zhl - 1
                counter = counter + 1
                
                if ( allocated(Qtr%seeds) ) call connect_h_seeds(1,Qtr,zhl,temp_coor,wpoly,signo)
                
                !------------
                do i = 1, zhl
                    call binarySearch_2 (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
                end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
                
                allocate (Qtr%elem(zhl), Stat = istat)
                Qtr%elem = temp_elem(1:zhl)
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
                        
                           
                        if (cc) call connect_h_neighbor(1,root,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                         zhl = zhl + 1
                        cc = .true.
                        
                        if (cc_2) call connect_h_neighbor(2,root,itr_nro,ref,level,zhl_2,temp_coor_2,wpoly_2,signo_2) ! ??
                        zhl_2 = zhl_2 + 1 
                        cc_2 = .true.
                        
                end do
                
                zhl = zhl - 1
                counter = counter + 1
                
                zhl_2 = zhl_2 - 1
                
                if ( allocated(Qtr%seeds) ) call connect_h_seeds(1,Qtr,zhl,temp_coor,wpoly,signo)
                
                !------------
                do i = 1, zhl
                    call binarySearch_2 (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
                end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
                
                allocate (Qtr%elem(zhl), Stat = istat)
                Qtr%elem = temp_elem(1:zhl)
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
                    call binarySearch_2 (nodes(1:num_node), temp_coor_2(i), index)
                    temp_elem_2(i) = index
                    if (index == 0) write(*,*) counter, temp_coor_2(i), zhl_2
                end do 
                
                call polykernel(zhl_2,temp_coor_2(1:zhl_2),x0)
                nodes(num_node + counter) = x0
                
                !allocate (Qtr%elem(zhl), Stat = istat)
                !Qtr%elem = temp_elem(1:zhl)
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
                            
                call connect_h_neighbor(2,root,itr_nro,ref,level,zhl,temp_coor,wpoly,signo)
                zhl = zhl + 1
               
               end do
               
               zhl = zhl - 1
               counter = counter + 1
               
               !------------
               do i = 1, zhl
                    call binarySearch_2 (nodes(1:num_node), temp_coor(i), index)
                    temp_elem(i) = index
                    if (index == 0) write(*,*) counter, temp_coor(i), index
               end do 
                
                call polykernel(zhl,temp_coor(1:zhl),x0)
                nodes(num_node + counter) = x0
               
                allocate (Qtr%elem(zhl), Stat = istat)
                Qtr%elem = temp_elem(1:zhl)
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
        
            if (QChildren(1)) call connectivity(root,Qtr%NW)
            if (QChildren(2)) call connectivity(root,Qtr%SW)
            if (QChildren(3)) call connectivity(root,Qtr%NE)
            if (QChildren(4)) call connectivity(root,Qtr%SE)
            
        end if 
   return 
1000    format(i6','i6','i6','i6,','i6','i6','i6','i6,','i6','i6','i6','i6,','i6','i6','i6','i6)
        
    end subroutine
    
    
   recursive subroutine QIntrsPts(root,Qtr,num_poly,polygons)
        use Qtree_data
        !use  Qtree_input
        use SortSearch_module
        use polygon_module
        use segment_module
        
        implicit none
        
        Type (Qtree), pointer :: root, Qtr
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
                                write(56,'(2f32.16)') pt1%x, pt1%y   
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
                                write(56,'(2f32.16)') pt1%x, pt1%y   
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
                                write(56,'(2f32.16)') pt1%x, pt1%y   
                            end if
                            
                            if (intrsc2) then  
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt2
                                Qtr%intrsc_points(zhl)%wpoly = nr2
                                write(56,'(2f32.16)') pt2%x, pt2%y   
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
                                write(56,'(2f32.16)') pt1%x, pt1%y   
                            end if
                            
                            if ( intrsc2) then
                                zhl = zhl + 1
                                Qtr%intrsc_points(zhl)%pos = pt2
                                Qtr%intrsc_points(zhl)%wpoly = nr2
                                write(56,'(2f32.16)') pt2%x, pt2%y   
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
                                write(56,'(2f32.16)') pt1%x, pt1%y   
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
            if (QChildren(1)) call QIntrsPts(root,Qtr%NW,num_poly,polygons)
            if (QChildren(2)) call QIntrsPts(root,Qtr%SW,num_poly,polygons)
            if (QChildren(3)) call QIntrsPts(root,Qtr%NE,num_poly,polygons)
            if (QChildren(4)) call QIntrsPts(root,Qtr%SE,num_poly,polygons)
        end if 
    end subroutine 
     
    recursive subroutine coordinates(root,Qtr)    ! später umbenennen
        use Qtree_data
        use SortSearch_module
        use polygon_module
        
        implicit none
        
        Type (Qtree), pointer :: root, Qtr
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
                            call binarySearch_2 (Temp_nodes, Qtr%Boundary(i), index)
                            nodes_mask(index) = .true. 
                        end if     
                    end do 
                    
                    if(Qtr%num_intrsc_points .gt. 0) then
                    
                        do i = 1, Qtr%num_intrsc_points
                            call binarySearch_2 (Temp_nodes, Qtr%intrsc_points(i)%pos, index)
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
                                
                                call binarySearch_2 (Temp_nodes, Qtr%Boundary(i), index)
                                nodes_mask(index) = .true. 
                                
                            else if ( (Qtr%signo(i) .eq. 0 .and. Qtr%wpoly(i) .gt. 1) ) then
                                
                                call binarySearch_2 (Temp_nodes, Qtr%Boundary(i), index)
                                nodes_mask(index) = .true. 
                                
                            end if   
                           
                        end do    
                     !end do   
            end if
            
        else if (count(QChildren).gt.0) then
            if (QChildren(1)) call coordinates(root,Qtr%NW)
            if (QChildren(2)) call coordinates(root,Qtr%SW)
            if (QChildren(3)) call coordinates(root,Qtr%NE)
            if (QChildren(4)) call coordinates(root,Qtr%SE)
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
                                call binarySearch_2 (Temp_nodes, Qtr%Boundary(i), index)
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
    
    
    recursive subroutine QtrBalance(root,Qtr,level_min)
        
        implicit none
        Type (Qtree), pointer :: root
        Type (Qtree), pointer :: Qtr
        Type (Qtree), pointer :: Qtr_out
        integer, intent(in) :: level_min
        logical :: status
        integer :: i, ref(36), level
        character*2 :: nb(8)
        
        nb = ['NN','NW','WW','SW','SS','SE','EE','NE']
        
        if (.not. associated(Qtr%NW)) then
        
             if (Qtr%level.gt.level_min) then
                       
                do i = 1,8
                    ref = Qtr%ref
                    level = Qtr%level
                    call QtrFindNeighborQ(level,ref,nb(i),status)
                        if (status) then
                            call QtrSearchQ(root,level,ref,1,Qtr_out)
                        
                            if (Qtr_out%level .lt. level-1) then
                                call QtrSubdivide(Qtr_out,Qtr_out%level+1)
                            end if
                            
                        end if 
                end do
                      
                
             end if
         
        else if (associated(Qtr%NW)) then
        
            call QtrBalance(root,Qtr%NW,level_min)
            call QtrBalance(root,Qtr%SW,level_min)
            call QtrBalance(root,Qtr%NE,level_min)
            call QtrBalance(root,Qtr%SE,level_min)
            
        end if 
        
    end subroutine 
    
    recursive subroutine QtrSearchQ(Qtr,level,ref,i,Qtr_out)
        implicit none 
        Type (Qtree), intent(in), pointer :: Qtr
        integer, intent(in) :: level, ref(36), i
        Type (Qtree), intent(out), pointer :: Qtr_out
        logical :: Qchildren(4)
        integer :: nx, ny
        
        Qchildren = [associated(Qtr%NW),associated(Qtr%SW),associated(Qtr%NE),associated(Qtr%SE)]
        
        if (count(QChildren).gt.0) then
            nx = ref(2*i-1)
            ny = ref(2*i) 
     
            if (i.gt.level) then
                Qtr_out => Qtr
                return 
            end if 
     
            if ((nx.eq.1) .and. (ny.eq.1) .and. QChildren(1)) then      
            
                call QtrSearchQ(Qtr%NW,level,ref,i+1,Qtr_out)
                
            else if ((nx.eq.1) .and. (ny.eq.2) .and. QChildren(2)) then
            
                call QtrSearchQ(Qtr%SW,level,ref,i+1,Qtr_out)
                
            else if ((nx.eq.2) .and. (ny.eq.1) .and. QChildren(3)) then
            
                call QtrSearchQ(Qtr%NE,level,ref,i+1,Qtr_out)
                
            else if ((nx.eq.2) .and. (ny.eq.2) .and. QChildren(4)) then
            
                call QtrSearchQ(Qtr%SE,level,ref,i+1,Qtr_out)
                
            else     
               Qtr_out => Qtr
               return  
            end if 
            
        else if (count(QChildren).eq.0) then
            Qtr_out => Qtr
            return
        end if 
    end subroutine 
    
    subroutine QtrFindNeighborQ(level,ref,dir,status)
        implicit none 
        integer, intent(in) :: level
        integer, intent(inout) :: ref(2*level)
        logical, intent(out) :: status
        character*2 :: dir, pos 
        integer :: Nx(level), Ny(level)
        integer :: i, lim_x, lim_y
    
        status = .true.
        lim_x = 1
        lim_y = 1

        
        Nx = ref(1:2*level:2)
        Ny = ref(2:2*level:2)
    
        do  i = 1, level
            Nx(i) = Nx(i) - 1
            Ny(i) = Ny(i) - 1
        end do
   
        do i = level, 1, -1
            if (i .gt. 1) then
                if (Nx(i) .ne. Nx(i-1)) then 
                    lim_x = i-1
                    exit 
                end if 
            end if 
        end do 
    
        do i = level, 1, -1
            if (i.gt. 1) then
                if (Ny(i) .ne. Ny(i-1)) then 
                    lim_y = i-1
                    exit 
                end if
            end if     
        end do
    
        if ((Nx(level).eq.0) .and. (Ny(level).eq.0)) then      
            pos = 'NW'
        else if ((Nx(level).eq.0) .and. (Ny(level).eq.1)) then
            pos = 'SW'
        else if ((Nx(level).eq.1) .and. (Ny(level).eq.0)) then
            pos = 'NE'
        else if ((Nx(level).eq.1) .and. (Ny(level).eq.1)) then
            pos = 'SE'
        end if 
    
        call ucase(dir)

!NN        
        if (dir == 'NN') then
            if (pos == 'NW') then
            
                do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
                end do 
               
            else if (pos == 'SW') then   
                Ny(level) = abs(Ny(level) - 1)
            else if (pos == 'NE') then
            
                do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
                end do
               
            else if (pos == 'SE') then
                Ny(level) = abs(Ny(level) - 1)
            end if 
            
!NW            
        else if (dir == 'NW') then
            if (pos == 'NW') then
            
                do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
                end do      
               
            else if (pos == 'SW') then
            
               do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
               do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
               Ny(level) = abs(Ny(level)-1) 
               
            else if (pos == 'NE') then
            
               do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do 
               Nx(level) = abs(Nx(level)-1)    
               
            else if (pos == 'SE') then
            
                Nx(level) = abs(Nx(level) - 1)
                Ny(level) = abs(Ny(level) - 1)
                
            end if
            
!WW            
        else if (dir == 'WW') then
        
            if (pos == 'NW') then
            
                do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
            
            else if (pos == 'SW') then
            
                do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
               
            else if (pos == 'NE') then
            
                Nx(level) = abs(Nx(level) - 1)
                
            else if (pos == 'SE') then
            
                Nx(level) = abs(Nx(level) - 1)
                
            end if

!SW            
        else if (dir == 'SW') then
            if (pos == 'NW') then
            
              do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                Ny(level) = abs(Ny(level)-1)
               
            else if (pos == 'SW') then
            
                do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
                end do
                
            else if (pos == 'NE') then
            
                Nx(level) = abs(Nx(level) - 1)
                Ny(level) = abs(Ny(level) - 1) 
                
            else if (pos == 'SE') then
            
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Nx(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
               do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do 
               Nx(level) = abs(Nx(level)-1)
            end if

!SS            
        else if (dir == 'SS') then
            if (pos == 'NW') then
            
                Ny(level) = abs(Ny(level)-1)
                
            else if (pos == 'SW') then
            
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
            
               do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do     
               
            else if (pos == 'NE') then
            
                Ny(level) = abs(Ny(level)-1)
                
            else if (pos == 'SE') then
            
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
                
               do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do
            end if
            
!SE            
        else if (dir == 'SE') then
            if (pos == 'NW') then
            
                Nx(level) = abs(Nx(level) - 1)
                Ny(level) = abs(Ny(level) - 1)
                
            else if (pos == 'SW') then
            
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                 do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do 
               Nx(level) = abs(Nx(level)-1)
               
            else if (pos == 'NE') then
            
                do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
               
               do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
               Ny(level) = abs(Ny(level) - 1)
               
            else if (pos == 'SE') then
            
                do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Ny(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
                end do  
                      
            end if
            
!EE            
        else if (dir == 'EE') then
            if (pos == 'NW') then
            
                Nx(level) = abs(Nx(level) - 1)
                
            else if (pos == 'SW') then
                Nx(level) = abs(Nx(level) - 1)
                
            else if (pos == 'NE') then
            
               do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
               
            else if (pos == 'SE') then
                
                do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
                   
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                    
            end if
            
!NE            
        else if (dir == 'NE') then
            if (pos == 'NW') then
            
                do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do 
               Nx(level) = abs(Nx(level)-1)
               
            else if (pos == 'SW') then
            
               Nx(level) = abs(Nx(level) - 1)
               Ny(level) = abs(Ny(level) - 1)  
               
            else if (pos == 'NE') then
            
                do i = 1, level
                    if (Ny(i).eq.1) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                
                do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do
                
                if (status .eqv. .false.) return
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
               end do
               do i = lim_y, level
                    Ny(i) = abs(Ny(i) - 1)
               end do
               
            else if (pos == 'SE') then
            
                do i = 1, level
                    if (Nx(i).eq.0) then
                        status = .true.    
                        exit    
                    end if 
                    status = .false.
                end do 
                
                if (status .eqv. .false.) return
                
                do i = lim_x, level
                    Nx(i) = abs(Nx(i) - 1)
                end do
                Ny(level) = abs(Ny(level)-1)
                     
            end if
            
        end if
        
        do  i = 1, level
           Nx(i) = Nx(i) + 1
           Ny(i) = Ny(i) + 1
        end do
        
        ref(1:2*level:2) = Nx
        ref(2:2*level:2) = Ny
        
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
        
        real*8 :: xmin, xmax, ymin, ymax, xmid, ymid, x0, y0
        real*8 :: boundary(8)
        integer :: i, istat, level, azhl, numseeds
        logical :: containsPoint
      
        if (Qtr%level .ge. 18) return
        if (.not. associated(Qtr%NW)) then
            
            Allocate (Qtr%NW, Qtr%SW, Qtr%NE, Qtr%SE, Stat=istat)
            
            !nullify(Qtr%NW%NW, Qtr%NW%SW, Qtr%NW%NE, Qtr%NW%SE)
            !nullify(Qtr%SW%NW, Qtr%SW%SW, Qtr%SW%NE, Qtr%SW%SE)
            !nullify(Qtr%NE%NW, Qtr%NE%SW, Qtr%NE%NE, Qtr%NE%SE)
            !nullify(Qtr%SE%NW, Qtr%SE%SW, Qtr%SE%NE, Qtr%SE%SE)
            
            xmin = minval(Qtr%Boundary(1:4)%x)
            xmax = maxval(Qtr%Boundary(1:4)%x)
            ymin = minval(Qtr%Boundary(1:4)%y)
            ymax = maxval(Qtr%Boundary(1:4)%y)
            xmid = sum(Qtr%Boundary(1:4)%x)/4.d0
            ymid = sum(Qtr%Boundary(1:4)%y)/4.d0
            
            Qtr%center = point(xmid,ymid)
            
            if (Qtr%level.eq.0) num_node = num_node + 9
            if (Qtr%level.gt.0) num_node = num_node + 5
            if (Qtr%level.eq.0) write(56,'(2f32.16)') xmin, ymin
            if (Qtr%level.eq.0) write(56,'(2f32.16)') xmax, ymin
            if (Qtr%level.eq.0) write(56,'(2f32.16)') xmax, ymax
            if (Qtr%level.eq.0) write(56,'(2f32.16)') xmin, ymax
            write(56,'(2f32.16)') xmid, ymid
            write(56,'(2f32.16)') xmid, ymin
            write(56,'(2f32.16)') xmax, ymid
            write(56,'(2f32.16)') xmid, ymax
            write(56,'(2f32.16)') xmin, ymid
            
            
            
            
            level = Qtr%level + 1
            Qtr%NW%level = level
            Qtr%NW%ref = Qtr%ref
            Qtr%NW%ref(2*level-1)= 1
            Qtr%NW%ref(2*level)= 1
            boundary= [xmin,ymid,xmid,ymid,xmid,ymax,xmin,ymax]
            
            do i = 1,4
                Qtr%NW%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do
            
            x0 = sum(Qtr%NW%Boundary(1:4)%x)/4.d0
            y0 = sum(Qtr%NW%Boundary(1:4)%y)/4.d0
            Qtr%NW%center = point(x0,y0)
            
            Qtr%SW%level = level
            Qtr%SW%ref = Qtr%ref
            Qtr%SW%ref(2*level-1)= 1
            Qtr%SW%ref(2*level)= 2
            boundary= [xmin,ymin,xmid,ymin,xmid,ymid,xmin,ymid]
            
            do i = 1,4
                Qtr%SW%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do
            
            x0 = sum(Qtr%SW%Boundary(1:4)%x)/4.d0
            y0 = sum(Qtr%SW%Boundary(1:4)%y)/4.d0
            Qtr%SW%center = point(x0,y0)
            
            Qtr%NE%level = level
            Qtr%NE%ref = Qtr%ref
            Qtr%NE%ref(2*level-1)= 2
            Qtr%NE%ref(2*level)= 1
            boundary= [xmid,ymid,xmax,ymid,xmax,ymax,xmid,ymax]
            
            
            do i = 1,4
                Qtr%NE%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do
            
            x0 = sum(Qtr%NE%Boundary(1:4)%x)/4.d0
            y0 = sum(Qtr%NE%Boundary(1:4)%y)/4.d0
            Qtr%NE%center = point(x0,y0)
            
            Qtr%SE%level = level
            Qtr%SE%ref = Qtr%ref
            Qtr%SE%ref(2*level-1)= 2
            Qtr%SE%ref(2*level)= 2
            boundary= [xmid,ymin,xmax,ymin,xmax,ymid,xmid,ymid]
            
            do i = 1,4
                Qtr%SE%Boundary(i) = point(boundary(2*i-1), boundary(2*i))
            end do
            
            x0 = sum(Qtr%SE%Boundary(1:4)%x)/4.d0
            y0 = sum(Qtr%SE%Boundary(1:4)%y)/4.d0
            Qtr%SE%center = point(x0,y0)
     
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
                !else 
                
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
                !else 
                
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
                !else
                
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
                !else
                
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



