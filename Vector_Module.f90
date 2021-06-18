module vector_module
    use point_module
    implicit none 
    Private
    Public :: vector, assignment(=), operator(+), operator(-), &
              operator(*), operator(/), operator(.DOT.)
    
    
    Type :: vector 
        real*8 :: x
        real*8 :: y
    end Type 
    
    Interface assignment (=)
        module procedure array_to_vector
        module procedure vector_to_array
        module procedure point_to_vector
        module procedure vector_to_point
    end Interface 
    
    Interface operator (+)
         module procedure vector_add_vv
         module procedure vector_add_vp
         module procedure vector_add_pv
    end Interface 
    
    Interface operator (-)
        module procedure vector_subtract_vv
        module procedure vector_subtract_pp
    end Interface
    
    Interface operator (*)
        module procedure vector_times_real
        module procedure real_times_vector
        module procedure vector_times_int
        module procedure int_times_vector
    end Interface
    
    Interface operator (/)
        module procedure vector_div_real
        module procedure vector_div_int
    end Interface
    
    Interface operator (.DOT.)
        module procedure dot_product
    end Interface
    
    contains
        subroutine  array_to_vector(vec_result, array)
            Type (vector), intent(out) :: vec_result
            real*8, dimension(2), intent(in) :: array
            
            vec_result%x = array(1)
            vec_result%y = array(2)
            
        end subroutine array_to_vector
        
        subroutine  vector_to_array(array_result, vec)
            real*8, dimension(2), intent(out) :: array_result
            Type (vector), intent(in) :: vec
            
            array_result(1) = vec%x
            array_result(2) = vec%y
            
        end subroutine vector_to_array
        
        subroutine  point_to_vector(vec_result, pt)
            Type (vector), intent(out) :: vec_result
            Type (point), intent(in) :: pt
            
            vec_result%x = pt%x
            vec_result%y = pt%y
            
        end subroutine point_to_vector
        
        subroutine  vector_to_point(pt_result, vec)
            Type (point), intent(out) :: pt_result
            Type (vector), intent(in) :: vec
            
            
             pt_result%x = vec%x
             pt_result%y = vec%y
            
        end subroutine vector_to_point
            
        function vector_add_vv(vec_1, vec_2)
            type (vector) :: vector_add_vv
            type (vector), intent(in) :: vec_1, vec_2
            
            vector_add_vv%x = vec_1%x + vec_2%x
            vector_add_vv%y = vec_1%y + vec_2%y
            
        end function vector_add_vv
        
        function vector_add_vp(vec_1, pt_2)
            type (vector) :: vector_add_vp
            type (vector), intent(in) :: vec_1
            type (point), intent(in) :: pt_2
            
            vector_add_vp%x = vec_1%x + pt_2%x
            vector_add_vp%y = vec_1%y + pt_2%y
            
        end function vector_add_vp
        
        function vector_add_pv(pt_1, vec_2)
            type (vector) :: vector_add_pv
            type (point), intent(in) :: pt_1
            type (vector), intent(in) :: vec_2
            
            vector_add_pv%x = pt_1%x + vec_2%x
            vector_add_pv%y = pt_1%y + vec_2%y
            
        end function vector_add_pv
        
        
        
        function vector_subtract_vv(vec_1, vec_2)
            type (vector) :: vector_subtract_vv
            type (vector), intent(in) :: vec_1, vec_2
            
            vector_subtract_vv%x = vec_1%x - vec_2%x
            vector_subtract_vv%y = vec_1%y - vec_2%y
            
        end function vector_subtract_vv
        
        
        function vector_subtract_pp(pt_1, pt_2)
            type (vector) :: vector_subtract_pp
            type (point), intent(in) :: pt_1, pt_2
            
            vector_subtract_pp%x = pt_1%x - pt_2%x
            vector_subtract_pp%y = pt_1%y - pt_2%y
            
        end function vector_subtract_pp
        
        
        
        
        
        
        function vector_times_real(vec_1, real_2)
            type (vector) :: vector_times_real
            type (vector), intent(in) :: vec_1
            real*8, intent(in) :: real_2
            
            vector_times_real%x = vec_1%x * real_2
            vector_times_real%y = vec_1%y * real_2
        
        end function vector_times_real
        
        function real_times_vector(real_1, vec_2)
            type (vector) :: real_times_vector
            type (vector), intent(in) :: vec_2
            real*8, intent(in) :: real_1
            
            real_times_vector%x = real_1 * vec_2%x 
            real_times_vector%y = real_1 * vec_2%y 
        
        end function real_times_vector
        
        function vector_times_int(vec_1, int_2)
            type (vector) :: vector_times_int
            type (vector), intent(in) :: vec_1
            integer, intent(in) :: int_2
            
            vector_times_int%x = vec_1%x * int_2
            vector_times_int%y = vec_1%y * int_2
        
        end function vector_times_int
        
        function int_times_vector(int_1, vec_2)
            type (vector) :: int_times_vector
            type (vector), intent(in) :: vec_2
            integer, intent(in) :: int_1
            
            int_times_vector%x = int_1 * vec_2%x 
            int_times_vector%y = int_1 * vec_2%y 
        
        end function int_times_vector
        
       function vector_div_real(vec_1, real_2)
            type (vector) :: vector_div_real
            type (vector), intent(in) :: vec_1
            real*8, intent(in) :: real_2
            
            vector_div_real%x = vec_1%x / real_2
            vector_div_real%y = vec_1%y / real_2
            
       end function vector_div_real
       
       
       function vector_div_int(vec_1, int_2)
            type (vector) :: vector_div_int
            type (vector), intent(in) :: vec_1
            integer, intent(in) :: int_2
            
            vector_div_int%x = vec_1%x / Real(int_2)
            vector_div_int%y = vec_1%y / Real(int_2)
       
       end function vector_div_int
    
       function dot_product(vec_1, vec_2)
            real*8 :: dot_product 
            type (vector), intent(in) :: vec_1, vec_2
            
            dot_product = vec_1%x*vec_2%x + vec_1%y*vec_2%y
       end function dot_product
        
       
       
end module     