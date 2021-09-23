program quadtree_main

    use point_module
    use Qtree_module
    
    implicit none
    type(Qtree), pointer :: Q3 => null()
    type(QtreeList), pointer :: QtrList => null()
    real(8) :: pointsArray(4,2)
    real(8) :: seedsArray(2,2)

    ! A square of side length 10
    pointsArray(1,:) = [ 0, 0]
    pointsArray(2,:) = [10, 0]
    pointsArray(3,:) = [10,10]
    pointsArray(4,:) = [ 0,10]

    seedsArray(1,:) = [ 2, 2]
    seedsArray(2,:) = [ 3, 3]

    ! <QUADTREE>
    ! Quadtree decomposition
    ! Initilize
    allocate(Q3)
    call QtrInit(Q3, 4, pointsArray)

    open(unit=55, file='./mtlb/vertices0.txt', status='unknown')
        call QtrPrintLeaves(55, Q3)
    close(unit=55)
    
end program quadtree_main
