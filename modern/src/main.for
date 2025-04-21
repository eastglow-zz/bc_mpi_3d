      program bctest_par

      use mpi_f08

      use input
      use boundary_conditions
      use pvoutputs
      implicit none

      type(MPI_Status) :: istatus
      type(MPI_Comm) :: comm 
      integer(4),allocatable :: a(:,:,:,:)
      real(8),allocatable :: b(:,:,:,:)
      real(8),allocatable :: c(:,:,:,:)
      integer(4) :: numghost
      integer(4) :: ngx, ngy, ngz
      integer(4) :: il,iu, jl,ju, kl,ku
      integer(4) :: ilb, iub, jlb, jub, klb, kub

      integer(4) :: i,j,k,nn, ll
      integer(4) :: ig, jg, kg
      integer(4) :: ierr
      integer(4) :: pidx, pidy, pidz
      real(8) :: wtime_start, wtime_end, wtime 

      character(len=100) :: bctype 

      bctype = 'ADIABATIC'

      NPX = 2
      NPY = 2
      NPZ = 1

      imori = 32
      jmori = 32
      kmori = 1
      np = 5

      call get_dimension(imori, jmori, kmori)

      numghost = 2
      ngx = numghost * dimx
      ngy = numghost * dimy
      ngz = numghost * dimz

      comm = MPI_COMM_WORLD
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(comm,nprocs,ierr)
      call MPI_COMM_RANK(comm,myrank,ierr)

      call init_bc(NPX,NPY,NPZ,IMORI,JMORI,KMORI,myrank, comm)
      iml = bc_iml 
      jml = bc_jml 
      kml = bc_kml 
      call init_pvoutputs(IMORI,JMORI,KMORI,bc_iml,bc_jml,bc_kml,
     &                         bc_il, bc_iu, bc_jl, bc_ju, bc_kl, bc_ku, 
     &                0.d0, 0.d0, 0.d0, dble(DXL), dble(DYL), dble(DZL),
     &                  nprocs, bc_pidx, bc_pidy,bc_pidz, bc_pidg, comm)

      ! call get_ijkrange(il,iu,jl,ju,kl,ku,myrank)
      ! iml=iu-il+1
      ! jml=ju-jl+1
      ! kml=ku-kl+1

      ilb = 1-ngx
      iub = iml+ngx
      jlb = 1-ngy
      jub = jml+ngy
      klb = 1-ngz
      kub = kml+ngz


      allocate(a(ilb:iub,jlb:jub,klb:kub,np))
      allocate(b(ilb:iub,jlb:jub,klb:kub,np))
      allocate(c(ilb:iub,jlb:jub,klb:kub,0:1))

      !initialization
      a(:,:,:,:) = 0
      b(:,:,:,:) = 0.0
      c(:,:,:,:) = 0.0
      call distance_from_center(c(:,:,:,told(0)),ngx,ngy,ngz,1)

      wtime_start = MPI_WTIME()

      call init_bc_call_count(0)

      ttime = 0.0d0
      do i=0,1000
        call init_bc_call_count(0)
        !call boundary_condition_r8('PERIODIC',c(:,:,:,told(i))          ! wtime_mpif08.csv
        call boundary_condition_r8(bctype,c(:,:,:,told(i)), numghost)    ! wtime_modern.csv

        call simple_diffusion(c(:,:,:,tnew(i)),c(:,:,:,told(i)), 
     &                                                  ngx,ngy,ngz, i)

        if(mod(i,10).eq.0)then
          !call boundary_condition_r8('PERIODIC',c(:,:,:,tnew(i))        ! wtime_mpif08.csv
          call boundary_condition_r8(bctype,c(:,:,:,tnew(i)),numghost)  ! wtime_modern.csv
    !       call output_r8_pvtr_2d(c(:,:,:,tnew(i)),numghost,'c',i,ttime,  
    !  &                                              myrank ,"test pvtr")
    !       call output_r8_pvtr(c(:,:,:,tnew(i)),numghost,'c',i,ttime,  
    !  &                                              myrank ,"test pvtr")
          call output_r8_pvtr(c(:,:,:,tnew(i)),numghost,'c',  
     &                                     i,ttime, myrank ,"test pvtr")

          if(myrank.eq.0)then
            call verbose_bc_call_count()
            call verbose_bc_max_mpi_tag()
          endif
        endif
 
        ttime = ttime + dttime
      enddo !do i=0,10000

      wtime_end = MPI_WTIME()

      wtime = wtime_end - wtime_start

      if (myrank.eq.0) then 
        write(*,*)'Wall time spent (s): ', wtime
      endif 

      deallocate(a)
      deallocate(b)
      deallocate(c)

      call finalize_bc()
      call finalize_pvoutputs()
      call MPI_FINALIZE(ierr)


      contains

      ! ----------------------------------------------------------------

      subroutine print_aijk(a,ngx,ngy,ngz)
      use input
      implicit none

      integer(4),intent(in) :: ngx,ngy,ngz
      integer(4),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)

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

      ! ----------------------------------------------------------------

      subroutine distance_from_center(a,ngx,ngy,ngz,timestep)
      use input
      use boundary_conditions
      implicit none

      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: timestep

      integer(4) :: i,j,k, ig,jg,kg
      real(8) :: rsq, r

      do k=1,kml
      do j=1,jml
      do i=1,iml
            call get_global_ijk(ig,jg,kg,i,j,k,myrank)
            rsq=(ig-IMORI/2)**2+(jg-JMORI/2)**2+(kg-KMORI/2)**2
            r=sqrt(rsq) + timestep
            a(i,j,k)=r
      enddo
      enddo
      enddo

      end subroutine distance_from_center

      ! ----------------------------------------------------------------

      subroutine simple_diffusion(a,ao,ngx,ngy,ngz,timestep)
      use input
      use boundary_conditions
      implicit none

      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      real(8),intent(inout) :: 
     &                     ao(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: timestep

      real(8) :: D, dx, dt
      real(8) :: axx, ayy, azz 

      integer(4) :: i,j,k

      D=1.0
      dx=1.0
      dt=0.8*0.15*DXL*DXL/D
      dttime = dt

      do k=1,kml
      do j=1,jml
      do i=1,iml
        axx = 0.0d0 
        ayy = 0.0d0
        azz = 0.0d0
        if(ngx.gt.0)
     &             axx = (ao(i-1,j,k)-2.0*ao(i,j,k)+ao(i+1,j,k))/DXL/DXL
        if(ngy.gt.0)
     &             ayy = (ao(i,j-1,k)-2.0*ao(i,j,k)+ao(i,j+1,k))/DYL/DYL
        if(ngz.gt.0) 
     &             azz = (ao(i,j,k-1)-2.0*ao(i,j,k)+ao(i,j,k+1))/DZL/DZL
        a(i,j,k)=ao(i,j,k) + D*dt*( axx + ayy + azz ) 
      enddo
      enddo
      enddo

      end subroutine simple_diffusion

      ! ----------------------------------------------------------------

      integer(4) function told(n)
      implicit none
      integer(4),intent(in) :: n

      told=int(mod(n,2))
      return
      end function told

      ! ----------------------------------------------------------------

      integer(4) function tnew(n)
      implicit none
      integer(4),intent(in) :: n

      tnew=int(mod(n+1,2))
      return
      end function tnew

      ! ----------------------------------------------------------------

      end program bctest_par