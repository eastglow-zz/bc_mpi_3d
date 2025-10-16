      module input
      implicit none
      integer(4) :: imori,jmori,kmori,np ! variables for global indices
      integer(4) :: iml, jml, kml
      integer(4) :: dimen
      integer(4) :: dimx, dimy, dimz

      real(8) :: DXL=0.1
      real(8) :: DYL=0.1
      real(8) :: DZL=0.1

      !Variables for MPI parallelization
      integer(4) :: NPX, NPY, NPZ
      integer(4) :: nprocs, myrank
      real(8) :: ttime, dttime

      !File output path
      character(len=300) :: opath = './outputs'

      logical :: error_flag = .false. 
      character(len=1000) :: error_log = 'No errors.'

      contains 

      ! ----------------------------------------------------------------

      subroutine get_dimension(im, jm, km)
      implicit none 

      integer(4),intent(in) :: im, jm, km

      dimx = 1 
      dimy = 1
      dimz = 1
      if(im .le. 1) dimx = 0
      if(jm .le. 1) dimy = 0
      if(km .le. 1) dimz = 0 
      dimen = dimx + dimy + dimz
      ! write(*,*)'Dimension:', dimen
      ! write(*,*)'dimx, dimy, dimz:', dimx, dimy, dimz
      ! write(*,*)'imori, jmori, kmori:', imori, jmori, kmori
      if(dimen .eq. 0) then
        write(*,*)'Error, get_dimension() in grobal_vars.for'
        write(*,*)'  No dimension specified'
        write(*,*)'  Aborting...'
        call abort()
      end if

      return 
      end subroutine get_dimension

      ! ----------------------------------------------------------------

      end module input