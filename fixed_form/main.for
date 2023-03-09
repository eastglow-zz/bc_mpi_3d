      program bctest_par
      use input
      use boundary_conditions
      use outputs
      implicit none
      include 'mpif.h'

      integer(4) :: istatus(MPI_STATUS_SIZE)
      integer(4),allocatable :: a(:,:,:,:)
      real(8),allocatable :: b(:,:,:,:)
      real(8),allocatable :: c(:,:,:,:)
      integer(4) :: numghost
      integer(4) :: il,iu, jl,ju, kl,ku
      integer(4) :: ilb, iub, jlb, jub, klb, kub

      integer(4) :: i,j,k,nn, ll
      integer(4) :: ig, jg, kg
      integer(4) :: ierr
      integer(4) :: pidx, pidy, pidz

      NPX = 4
      NPY = 4
      NPZ = 1

      im = 64
      jm = 64
      km = 4
      np = 5

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
      allocate(c(ilb:iub,jlb:jub,klb:kub,0:1))

      !initialization
      a(:,:,:,:) = 0
      b(:,:,:,:) = 0.0
      c(:,:,:,:) = 0.0
      call distance_from_center(c(:,:,:,told(1)),numghost,1)

      call init_bc_call_count(0)

      do i=1,1000
      call init_bc_call_count(0)

      call boundary_condition_r8_par('periodic',c(:,:,:,told(i))
     &                                                    ,numghost)

      call simple_diffusion(c(:,:,:,tnew(i)),c(:,:,:,told(i)), 
     &                                                  numghost, i)

      if(mod(i,100).eq.0)then
      call boundary_condition_r8_par('periodic',c(:,:,:,tnew(i))
     &                                                    ,numghost)
      call output_r8_pvtr(c(:,:,:,tnew(i)),numghost,'c',i,myrank 
     &                     ,"test pvtr")

            if(myrank.eq.0)then
                  call verbose_bc_call_count()
                  call verbose_bc_max_mpi_tag()
            endif
      endif
 
      enddo !do i=1,1000


      deallocate(a)
      deallocate(b)
      deallocate(c)
      call MPI_FINALIZE(ierr)


      contains


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


      subroutine distance_from_center(a,ng,timestep)
      use input
      use boundary_conditions
      implicit none
      include 'mpif.h'

      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng, timestep

      integer(4) :: i,j,k, ig,jg,kg
      real(8) :: rsq, r

      do k=1,kml
      do j=1,jml
      do i=1,iml
            call get_global_ijk(ig,jg,kg,i,j,k,myrank)
            rsq=(ig-IM/2)**2+(jg-JM/2)**2+(kg-KM/2)**2
            r=sqrt(rsq) + timestep
            a(i,j,k)=r
      enddo
      enddo
      enddo

      end subroutine distance_from_center



      subroutine simple_diffusion(a,ao,ng,timestep)
      use input
      use boundary_conditions
      implicit none
      include 'mpif.h'

      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      real(8),intent(inout) :: ao(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng, timestep

      real(8) :: D, dx, dt

      integer(4) :: i,j,k

      D=1.0
      dx=1.0
      dt=0.8*0.15*DXL*DXL/D

      do k=1,kml
      do j=1,jml
      do i=1,iml
            a(i,j,k)=ao(i,j,k) 
     &     + D*dt*( (ao(i-1,j,k)-2.0*ao(i,j,k)+ao(i+1,j,k))/DXL/DXL 
     &             +(ao(i,j-1,k)-2.0*ao(i,j,k)+ao(i,j+1,k))/DYL/DYL  
     &             +(ao(i,j,k-1)-2.0*ao(i,j,k)+ao(i,j,k+1))/DZL/DZL ) 
      enddo
      enddo
      enddo

      end subroutine simple_diffusion


      integer(4) function told(n)
      implicit none
      integer(4),intent(in) :: n

      told=int(mod(n,2))
      return
      end function told


      integer(4) function tnew(n)
      implicit none
      integer(4),intent(in) :: n

      tnew=int(mod(n+1,2))
      return
      end function tnew

      end program bctest_par