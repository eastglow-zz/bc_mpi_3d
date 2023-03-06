      program bctest_par
      use input
      use boundary_conditions
      use draw
      use outputs
      implicit none
      include 'mpif.h'

      integer(4) :: istatus(MPI_STATUS_SIZE)
      integer(4),allocatable :: a(:,:,:,:)
      real(8),allocatable :: b(:,:,:,:)
      integer(4) :: numghost
      integer(4) :: il,iu, jl,ju, kl,ku
      integer(4) :: ilb, iub, jlb, jub, klb, kub

      integer(4) :: i,j,k,nn, ll
      integer(4) :: ig, jg, kg
      integer(4) :: ierr

      NPX = 2
      NPY = 2
      NPZ = 1

      im = 100
      jm = 100
      km = 2
      np = 4

      numghost = 2

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

      call get_ijkrange(il,iu,jl,ju,kl,ku,myrank)
      iml=iu-il+1
      jml=ju-jl+1
      kml=ku-kl+1

      ilb = 1-numghost
      iub = iml+numghost
      jlb = 1-numghost
      jub = jml+numghost
      klb = 1-numghost
      kub = kml+numghost


      allocate(a(ilb:iub,jlb:jub,klb:kub,np))
      allocate(b(ilb:iub,jlb:jub,klb:kub,np))

      !initialization
      a(:,:,:,:) = 0
      b(:,:,:,:) = 0.0

      do k=1,kml
      do j=1,jml
      do i=1,iml
            a(i,j,k,1)=myrank
      enddo
      enddo
      enddo

      call put_sphere_i4(a(:,:,:,1),1,1,1, 20,10, numghost,'adiabatic'
     &                           ,myrank)

      call put_sphere_i4(a(:,:,:,1),40,30,1,25,10,numghost,'adiabatic'
     &                           ,myrank)


      call output_i4_raw(a(:,:,:,1),numghost,'a1',0,myrank,'test data')

      call output_i4_pvtr(a(:,:,:,1),numghost,'a1',0,myrank 
     &                     ,"test pvtr")


      deallocate(a)
      deallocate(b)
      call MPI_FINALIZE(ierr)

      end program bctest_par


      subroutine print_aijk(a,ng)
      use input
      implicit none

      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      integer(4) :: i,j,k

      write(*,*) 'Subroutine print_aijk() called'

      !print
      do k=1,kml
      do j=1,jml
      do i=1,iml
            write(*,*)a(i,j,k)
      enddo
      enddo
      enddo

      end subroutine print_aijk
