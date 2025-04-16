      module input
      implicit none
      integer(4) :: imori,jmori,kmori,np ! variables for global indices
      integer(4) :: iml, jml, kml

      real(8) :: DXL=0.1
      real(8) :: DYL=0.1
      real(8) :: DZL=0.1

      !Variables for MPI parallelization
      integer(4) :: NPX, NPY, NPZ
      integer(4) :: nprocs, myrank
      real(8) :: ttime, dttime

      !File output path
      character(len=300) :: opath = './outputs'

      end module input