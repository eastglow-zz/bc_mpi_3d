module input
implicit none
integer(4) :: im,jm,km,np ! variables for global indices
integer(4) :: iml,jml,kml ! variables for local indices' bound

real(8) :: DXL=1.0
real(8) :: DYL=1.0
real(8) :: DZL=1.0

!Variables for MPI parallelization
integer(4) :: NPX, NPY, NPZ
integer(4) :: nprocs, myrank

!File output path
character(len=300) :: opath = './outputs'

end module input