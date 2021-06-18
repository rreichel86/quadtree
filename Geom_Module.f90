module segment_module  ! line segment module
    use point_module
    use vector_module
     
    Type :: segment 
        type(point) :: pt_1
        type(point) :: pt_2
    end Type 
    
contains
    
    function pt_in_segment (pt, segm)
        implicit none
        logical pt_in_segment
        type(point), intent(in) :: pt
        type(segment), intent(in) :: segm
        
        real*8 :: t0_1, t0_2, tol
        
        tol = 1.d-10
        pt_in_segment = .false.
        
        if(segm%pt_1 == pt) then
            pt_in_segment = .true.
            return
        end if 
        
        if(segm%pt_2 == pt) then
            pt_in_segment = .true.
            return
        end if
        
        if (abs(segm%pt_2%x-segm%pt_1%x) .le. tol) then
        
            if (abs(segm%pt_1%x-pt%x) .le. tol) then
               
                t0_1 = (pt%y - segm%pt_1%y) / (segm%pt_2%y-segm%pt_1%y)
                if ((t0_1 < 0.d0 .or. t0_1 > 1.d0)) return
            else
               return
            end if 
            
        else if (abs(segm%pt_2%y-segm%pt_1%y) .le. tol) then   
        
            if (abs(segm%pt_1%y-pt%y) .le. tol) then
            
                t0_1 = (pt%x - segm%pt_1%x) / (segm%pt_2%x-segm%pt_1%x)
                if ((t0_1 < 0.d0 .or. t0_1 > 1.d0)) return
            else
                return
            end if 
        else     
        
           t0_1 = (pt%x - segm%pt_1%x) / (segm%pt_2%x-segm%pt_1%x)
           t0_2 = (pt%y - segm%pt_1%y) / (segm%pt_2%y-segm%pt_1%y)
           !if (.not.(abs(t0_1 - t0_2) .le. tol)) return 
           
           if ((t0_1 < 0 .or. t0_1 > 1) .or. (t0_2 < 0 .or. t0_2 > 1)) return
           
        end if 
        
            pt_in_segment = .true.
            return 
    
    end function 
    
    function pt2segm_dist (segm, pt)
        use vector_module
        use point_module
        implicit none
        real*8 :: pt2segm_dist
        type(point), intent(in) :: pt
        type(segment), intent(in) :: segm
        
        type(vector) :: v
        type(vector) :: w0
        type(point) :: pt_s
        real*8 :: c1, c2, b
        
        v = segm%pt_2 - segm%pt_1
        w0 = pt - segm%pt_1
        
        c1 = w0 .dot. v
        c2 = v .dot. v
        
        
        if (c1 <= 0 ) then 
        
            pt2segm_dist = pt2pt_dist(pt, segm%pt_1)
            return
            
        else if (c2 <= c2 ) then 
            
            pt2segm_dist = pt2pt_dist(pt, segm%pt_2)
            return
            
        else 
            b = c1/c2
            pt_s = segm%pt_1 + (b*v)
            pt2segm_dist = pt2pt_dist(pt, pt_s)
            return

        end if 
        
    end function pt2segm_dist
    
    subroutine intrsc_pt(segm_1, segm_2, Pt0, intrsc )
        implicit none
        type(segment), intent(in) :: segm_1, segm_2
        type(point), intent(out)  :: Pt0
        logical, intent(out) :: intrsc
        
        type(point) :: P1,P2,Q1,Q2
        type(vector) :: v1, v2
        real*8 :: tol, A, As, At, s0, t0
        logical :: cond1, cond2
        
        intrsc = .false.
        Pt0 = point(0d0,0d0)
        
        tol = 1.d-10
    
        P1 = segm_1%pt_1
        P2 = segm_1%pt_2
        v1 = P2 - P1 
        
        Q1 = segm_2%pt_1
        Q2 = segm_2%pt_2
        v2 = Q2 - Q1
        
        A =  -(P2%x - P1%x)*(Q2%y - Q1%y) + (P2%y - P1%y)*(Q2%x - Q1%x)
        
        if (abs(A) .le. tol) return
        
        As = -(Q1%x - P1%x)*(Q2%y - Q1%y) + (Q1%y - P1%y)*(Q2%x - Q1%x)
        At = (P2%x - P1%x)*(Q1%y - P1%y) - (P2%y - P1%y)*(Q1%x - P1%x)
        
        s0 = As/A  !   s0 >= 0  .and. s0 <= 1
        t0 = At/A  !   t0 >= 0  .and. t0 <= 1
        
        if ( (abs(s0) .gt. tol) .and. (s0 .lt. 0) ) return
        if ( (abs(t0) .gt. tol) .and. (t0 .lt. 0) ) return
        
        
        !if ((s0 .lt. 0) .or. (t0 .lt. 0)) return
        !if ((s0 .gt. 1) .or. (t0 .gt. 1)) return
        
        
        cond1 = ( (s0 .le. tol) .or. (s0 .gt. 0)) .and. ((abs(s0-1) .le. tol) .or. (s0 .lt. 1))
        cond2 = ( (t0 .le. tol) .or. (t0 .gt. 0)) .and. ((abs(t0-1) .le. tol) .or. (t0 .lt. 1))
        
 
        if ( cond1 .and. cond2) then  
            intrsc = .true.
        else 
            intrsc = .false.
        end if 
        !
        !
        if(intrsc) then
             Pt0 = P1 + s0 * v1 
        !    !Pt2 = Q1 + t0 * v2     
        end if     
    
    end subroutine 
    
    subroutine intrsc_pt2(ray_1, segm_2, Pt0, intrsc )
        implicit none
        type(segment), intent(in) :: ray_1, segm_2
        type(point), intent(out)  :: Pt0
        logical, intent(out) :: intrsc
        
        type(point) :: P1,P2,Q1,Q2
        type(vector) :: v1, v2
        real*8 :: tol, A, As, At, s0, t0
        logical :: cond1
        
        intrsc = .false.
        Pt0 = point(0d0,0d0)
        
        tol = 1.d-10
    
        P1 = ray_1%pt_1
        P2 = ray_1%pt_2
        v1 = P2 - P1 
        
        Q1 = segm_2%pt_1
        Q2 = segm_2%pt_2
        v2 = Q2 - Q1
        
        A =  -(P2%x - P1%x)*(Q2%y - Q1%y) + (P2%y - P1%y)*(Q2%x - Q1%x)
        
        if (abs(A) .le. tol) return
        
        As = -(Q1%x - P1%x)*(Q2%y - Q1%y) + (Q1%y - P1%y)*(Q2%x - Q1%x)
        At = (P2%x - P1%x)*(Q1%y - P1%y) - (P2%y - P1%y)*(Q1%x - P1%x)
        
        s0 = As/A  
        t0 = At/A  !   t0 >= 0  .and. t0 <= 1
        
        if ( (abs(t0) .gt. tol) .and. (t0 .lt. 0) ) return
        
        
        cond1 = ( (t0 .le. tol) .or. (t0 .gt. 0)) .and. ((abs(t0-1) .le. tol) .or. (t0 .lt. 1))
        
 
        if (cond1) then  
            intrsc = .true.
        else 
            intrsc = .false.
        end if 
        !
        !
        if(intrsc) then
             !Pt2 = P1 + s0 * v1 
             Pt0 = Q1 + t0 * v2     
        end if     
    
    end subroutine 
    
    
end module


module polygon_module
    use point_module
    use vector_module
    use segment_module
    
    private 
    public :: polygon, point_in_polygonSR, polygon_init, polygon_dist, &
        intrsc_segment_polygon, polygon_CCW_test, polygon_delete, polykernel, &
         intrsc_segment_polygon_2
    
    
    Type :: polygon
        integer :: num_vertices 
        type(point), allocatable :: vertices(:)
        type(point), allocatable :: kernel(:)
        type(segment), allocatable :: sides(:)
        type(point) :: center
        real*8 :: xmin, xmax, ymin, ymax
        integer :: mat_nro 
        logical :: init = .false.
    end Type    
    
contains 
    
    subroutine polygon_delete(polygon_1)
        implicit none
        type(polygon) :: polygon_1
        integer :: istat
        
        if(polygon_1%init) then
            polygon_1%num_vertices = 0
            polygon_1%mat_nro = 0
            polygon_1%center = point(0.d0,0.d0)
            polygon_1%xmin = 0.d0
            polygon_1%xmax = 0.d0
            polygon_1%ymin = 0.d0
            polygon_1%ymax = 0.d0
    
            deallocate(polygon_1%vertices, polygon_1%sides, stat=istat)
            if (allocated(polygon_1%kernel))  deallocate(polygon_1%kernel, stat=istat)
            
            polygon_1%init = .false.
            
        end if 
    
    end subroutine
    
    function polygon_CCW_test(polygon_1) 
        implicit none
        logical :: polygon_CCW_test
        Type(polygon), intent(in) :: polygon_1
        integer :: i, n
        real*8 :: area 
        
        n = polygon_1%num_vertices
        
        area = polygon_1%vertices(n)%x * polygon_1%vertices(1)%y - &
               polygon_1%vertices(1)%x * polygon_1%vertices(n)%y
        do i = 1, n-1
    
            area = area + polygon_1%vertices(i)%x * polygon_1%vertices(i+1)%y
            area = area - polygon_1%vertices(i+1)%x * polygon_1%vertices(i)%y
      
        end do   
        area = 0.5*area 
        
        if( area > 0) then
            polygon_CCW_test = .true.
            return 
        else
            polygon_CCW_test = .false.
            return 
        end if     
        
        
    end function 

    subroutine intrsc_segment_polygon(polygon_1, segm_2, pt0, intrsc)
        implicit none
        type(polygon), intent(in) :: polygon_1
        type(segment), intent(in) :: segm_2
        type(point), intent(out) :: pt0
        logical, intent(out) :: intrsc
        
        type(point):: pt
        integer :: i
        
        intrsc = .false.
        do i= 1, polygon_1%num_vertices
            
            call intrsc_pt(polygon_1%sides(i), segm_2, pt, intrsc )
            if (intrsc) then
                pt0 = pt
                exit
            end if 
        end do 
        
    end subroutine
    
    
     subroutine intrsc_segment_polygon_2(polygon_1, segm_2, pt0, intrsc)
        implicit none
        type(polygon), intent(in) :: polygon_1
        type(segment), intent(in) :: segm_2
        type(point), intent(out) :: pt0
        logical, intent(out) :: intrsc
        
        type(point):: pt
        integer :: i
        
        intrsc = .false.
        do i= 1, polygon_1%num_vertices
            
            call intrsc_pt(polygon_1%sides(i), segm_2, pt, intrsc )
            
            if ( intrsc .and. (.not.(segm_2%pt_1 .eq. pt)) .and. (.not.(segm_2%pt_2 .eq. pt)) ) then
                pt0 = pt
                return
            end if 
        end do 
        
        intrsc = .false.
        pt0 = point(0d0, 0d0)
        
    end subroutine
    
    
    subroutine polygon_dist(polygon_1, pt, min_dist, side_nro)
        implicit none
        type(polygon), intent(in) :: polygon_1
        type(point), intent(in) :: pt
        real*8, intent(out) :: min_dist
        integer, intent(out) :: side_nro
        
        real*8 , dimension(polygon_1%num_vertices) :: distances
        integer :: i
        
        
        do i= 1, polygon_1%num_vertices
            
            distances(i) = pt2segm_dist (polygon_1%sides(i), pt)
        
        end do 
        
        min_dist = minval(distances)
        side_nro = minloc(distances,1)
        
    
    end subroutine 
    
    subroutine polygon_init(polygon_1, pts, mat_nro)
        implicit none
        type(polygon), intent(out) :: polygon_1
        type(point), intent(in) :: pts(:)
        integer, intent(in) :: mat_nro
        
        real*8 :: xmid, ymid
        integer :: n, i, istat
        
        n = size(pts)
        polygon_1%mat_nro = mat_nro
        polygon_1%num_vertices = n
        allocate(polygon_1%vertices(n), polygon_1%sides(n), stat=istat)
        
        polygon_1%vertices = pts
        
        polygon_1%xmin = minval(polygon_1%vertices(:)%x)
        polygon_1%xmax = maxval(polygon_1%vertices(:)%x)
        polygon_1%ymin = minval(polygon_1%vertices(:)%y)
        polygon_1%ymax = maxval(polygon_1%vertices(:)%y)
        
        !xmid = sum(polygon_1%vertices(:)%x)/n
        !ymid = sum(polygon_1%vertices(:)%y)/n
        
        xmid = (polygon_1%xmin + polygon_1%xmax) / 2.d0
        ymid = (polygon_1%ymin + polygon_1%ymax) / 2.d0
        
        polygon_1%center= point(xmid,ymid)
        
        
        do i = 1, n-1
            polygon_1%sides(i) = segment(pts(i),pts(i+1))
        end do 
        polygon_1%sides(n) = segment(pts(n),pts(1))
        
        polygon_1%init = .true.
        
    end subroutine polygon_init 

    subroutine vertex_concavity (polygon_1, n, list)
        implicit none 
        
        type(polygon), intent(in) :: polygon_1
        integer, intent(in) :: n
        logical, intent(out) :: list(n)
    
        integer :: num_vertices, i
        type(point) :: A, B, C 
        num_vertices = polygon_1%num_vertices
        
        do i = 1, num_vertices
            
            if (i .eq. 1) then
                A = polygon_1%vertices(num_vertices)
            else 
                A = polygon_1%vertices(i-1)
            end if 
            
                B = polygon_1%vertices(i) 
                
            if (i .eq. num_vertices ) then
                C = polygon_1%vertices(1)  
            else
                C = polygon_1%vertices(i+1)
            end if 
            
            
            if (orientation(A,B,C) .eq. -1) then
                list(i) = .true.
            end if     
            
        end do
        
    
    end subroutine 
    
    
    subroutine polygon_kernel(polygon_1)
      implicit none
      type(polygon), intent(inout) :: polygon_1
      
      logical, allocatable :: mask_concave_vertices(:)
      integer, allocatable :: idx_concave_vertices(:)
      logical :: intrsc
      type(point) :: A, B, D, C1, C2, U
      integer :: i,ii,idx,j,k, zhl, status, num_vertices, num_concave_vertices, O1, O2
      type(point), allocatable :: kernel(:)
      logical, allocatable :: mask_kernel(:)
      !integer :: nkernel 
      
      num_vertices = polygon_1%num_vertices 
      
      allocate(mask_concave_vertices(num_vertices), polygon_1%kernel(num_vertices), stat=status)
      
      !Determine the concavity of each vertex of polygon_1
      call vertex_concavity(polygon_1,num_vertices,mask_concave_vertices)
      polygon_1%kernel = polygon_1%vertices 
      
      num_concave_vertices = count(mask_concave_vertices)
      
      if (num_concave_vertices .eq.0) then
        deallocate(mask_concave_vertices, stat=status)
        return
      else
          
          allocate (idx_concave_vertices(num_concave_vertices), stat=status)
          
          zhl = 0
          do i=num_vertices,1,-1
              if (mask_concave_vertices(i)) then
                  zhl = zhl + 1
                  idx_concave_vertices(zhl) = i
              end if     
          end do
          
          deallocate(mask_concave_vertices, stat=status)
          
      end if 
      
      
     do ii = 1, num_concave_vertices 
         
         ! A = v_(i-1)    vertex adjacent the concave vertex i 
         ! B = v_(i)      concave vertex i
         ! D = v_(i+1)    vertex adjacent the concave vertex i 
         
         if (idx_concave_vertices(ii) .eq. 1) then
             A = polygon_1%vertices(num_vertices)
         else 
             A = polygon_1%vertices(idx_concave_vertices(ii)-1)
         end if     
         
            B = polygon_1%vertices(idx_concave_vertices(ii))
         
         if (idx_concave_vertices(ii) .eq. num_vertices) then
             D = polygon_1%vertices(1)
         else 
             D = polygon_1%vertices(idx_concave_vertices(ii)+1)
         end if
         
         do i = 1, size(polygon_1%kernel) 
             !A-B-C1 triangle 
             !A-B_C2 triangle
             C1 = polygon_1%kernel(i)
             if ( i .eq. size(polygon_1%kernel) ) then
                 idx = 1
             else
                 idx = i+1
             end if     
             C2 = polygon_1%kernel(idx)
             
             ! compute orientation of the triangles
             O1 = orientation(A,B,C1)
             O2 = orientation(A,B,C2)
             
             if (O1 .eq. 0 ) cycle
             if (O2 .eq. 0 ) cycle
             
             ! if orientation of the triangles are opposite.
             ! compute intersection U of ray A-B and polygon edge C1-C2
             if (O1 .eq. -O2) then  
                 call intrsc_pt2(segment(A,B), segment(C1,C2), U, intrsc)
                 
                 if (intrsc) allocate(kernel( size(polygon_1%kernel)+1 ), stat=status)
                 
                 ! insert U in Kernel between C1-C2
                 if (idx .eq. 1) then 
                     kernel(1:size(polygon_1%kernel)) = polygon_1%kernel(1:size(polygon_1%kernel))
                     kernel(size(polygon_1%kernel)+1) = U
                     call move_alloc(kernel,polygon_1%kernel)
                 else 
                    kernel(1:i) = polygon_1%kernel(1:i)
                    kernel(i+1) = U
                    kernel(idx+1:size(polygon_1%kernel)+1) = polygon_1%kernel(idx:size(polygon_1%kernel))
                    call move_alloc(kernel,polygon_1%kernel)    
                 end if
                     
             end if 
             
         end do 
         
         do i = 1, size(polygon_1%kernel)
             !B-D-C1 triangle     
             !B-D_C2 triangle 
             C1 = polygon_1%kernel(i)
             if ( i .eq. size(polygon_1%kernel) ) then
                 idx = 1
             else
                 idx = i+1
             end if     
             C2 = polygon_1%kernel(idx)
             
             ! compute orientation of the triangles
             O1 = orientation(B,D,C1)
             O2 = orientation(B,D,C2)
             
             if (O1 .eq. 0 ) cycle
             if (O2 .eq. 0 ) cycle
             
             ! if orientation of the triangles are opposite.
             ! compute intersection U of ray A-B and polygon edge C1-C2
             if (O1 .eq. -O2) then
                 
                 
                 call intrsc_pt2(segment(B,D), segment(C1,C2), U, intrsc)
                 
                 
                 if (intrsc) allocate(kernel(size(polygon_1%kernel)+1), stat=status)
                 
                 ! insert U in Kernel between C1-C2
                 if (idx .eq. 1) then 
                     kernel(1:size(polygon_1%kernel)) = polygon_1%kernel(1:size(polygon_1%kernel))
                     kernel(size(polygon_1%kernel)+1) = U
                     call move_alloc(kernel,polygon_1%kernel)
                 else 
                    kernel(1:i) = polygon_1%kernel(1:i)
                    kernel(i+1) = U
                    kernel(idx+1:size(polygon_1%kernel)+1) = polygon_1%kernel(idx:size(polygon_1%kernel))
                    call move_alloc(kernel,polygon_1%kernel)   
                 end if
                     
             end if 
               
         end do 
         
         
         allocate(mask_kernel(size(polygon_1%kernel)), stat = status)
         mask_kernel = .true.
         do k = 1, size(polygon_1%kernel)
             C1 = polygon_1%kernel(k)
             
             O1 = orientation(A,B,C1)
             O2 = orientation(B,D,C1)
             
             !if orientation of triangle A-B-C1 or triangle  B-D-C1 is CW ...
             ! Store C1s 
             if ( (O1 .eq. -1) .or. (O2 .eq. -1) ) then
                mask_kernel(k) = .false.
             end if 
             
         end do
         
         
         if (count(mask_kernel) .eq. 0) then
             deallocate(polygon_1%kernel, stat = status)
             exit 
         end if 
         
         ! Delete C1s from the Kernel
         allocate(kernel(count(mask_kernel)), stat = status)
         kernel = pack(polygon_1%kernel,mask_kernel)
         call move_alloc(kernel, polygon_1%kernel)
         deallocate(mask_kernel, stat = status) 
     end do 
      
     
    end subroutine 
    
    
    subroutine polykernel(num_points, list_points, x0)
      implicit none
      integer, intent(in) :: num_points 
      type(point), intent(in) :: list_points(num_points)
      type(point), intent(out) :: x0
      
      type(polygon) :: polygon_1
      integer :: n
      
      call polygon_init(polygon_1, list_points, 1)
      call polygon_kernel(polygon_1)
      
        n = size(polygon_1%kernel)
      
        x0%x = sum(polygon_1%kernel(:)%x)/n
        x0%y = sum(polygon_1%kernel(:)%y)/n
      
      
      call polygon_delete(polygon_1)
      
    end subroutine 
    
    
    
    
    
    subroutine point_in_polygonSR (polygon_1, pt_q, sign, point_in_polygon, normal)
        implicit none
         
        type(polygon), intent(in) :: polygon_1
        logical, intent(out) :: point_in_polygon
        type(point), intent(in) :: pt_q
        integer,intent(out) :: sign
        logical,intent(in) :: normal
        
        integer :: num_sides, i
        num_sides = polygon_1%num_vertices
        
        sign = -1
        
        do i=1, num_sides
           sign = sign * right_cross(polygon_1%sides(i), pt_q)
           
           if (sign == 0) exit 
        end do 
        
          if (normal) then
            if(sign .eq. 1)  point_in_polygon = .true.
            if(sign .eq. 0)  point_in_polygon = .true.
            if(sign .eq. -1) point_in_polygon = .false.
          end if 
          
          if (.not. normal) then
          
            if(sign .eq. 1) then
                sign = -1
            else if(sign .eq. -1) then
                sign = 1
            end if     
            
            if(sign .eq. 1)  point_in_polygon = .true.
            if(sign .eq. 0)  point_in_polygon = .true.
            if(sign .eq. -1) point_in_polygon = .false.
          
          end if 
          
          
    end subroutine 
    
    function right_cross(side, pt)
        implicit none
        integer :: right_cross
        type(segment) :: side
        type(point) :: pt 
        
        type(point) :: A, B, C 
        real*8 :: tol, delta
        
        tol = 1.d-10
        A = pt
        B = side%Pt_1
        C = side%Pt_2
        
        if((abs(A%y - B%y) .le. tol).and.(abs(A%y - C%y) .le. tol)) then
        
            if (((A%x < B%x).and.(A%x < C%x )) .or. ((A%x > B%x).and.(A%x > C%x ))) then
            
                right_cross = 1
                return
                
            else
            
                right_cross = 0
                return
                
            end if 
        end if 
        
        if ((A == B) .or. (A == C)) then
        
            right_cross = 0
            return
            
        end if 
        if (B%y > C%y) then
            B = side%Pt_2
            C = side%Pt_1
        end if 
        
        if ((A%y <= B%y).or.(A%y > C%y)) then
        
            right_cross = 1
            return
            
        end if 
        delta = (B%x-A%x)*(C%y-A%y)-(B%y-A%y)*(C%x-A%x)
        
        if ((delta > 0) .and. (delta > tol)) then
        
            right_cross = -1
            return
            
        else if ((delta < 0) .and. (abs(delta) > tol)) then
        
            right_cross = 1
            return
            
        else 
            if (pt_in_segment (pt, side)) then
                right_cross = 0
                return
            else 
                right_cross = 1
            end if 
            
        end if 
    
    end function 
    
end module