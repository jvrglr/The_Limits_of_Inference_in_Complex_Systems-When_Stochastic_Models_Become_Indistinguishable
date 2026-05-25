module declarations_module
    !Module for public variables to be used in the rest of the modules
    implicit none
    public
    
    integer*4 :: ios
    double precision :: mu, K, D, x, t
  contains
  
    subroutine assingments()
      implicit none
        mu = 1d0
        K = 2.5d0
        D = 0.2d0
    end subroutine assingments
  
end module declarations_module
