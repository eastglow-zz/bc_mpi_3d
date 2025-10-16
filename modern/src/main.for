      program bctest_par

      use mpi_f08

      use input
      use boundary_conditions
      use diffusion_equation
      use pvoutputs

      implicit none

      type(MPI_Status) :: istatus
      type(MPI_Comm) :: comm 
      integer(4),allocatable :: a(:,:,:,:)
      real(8),allocatable :: b(:,:,:,:)
      real(8),allocatable :: c(:,:,:,:)
      integer(4),allocatable :: pf(:,:,:) 
      real(8),allocatable :: mob(:,:,:)
      real(8),allocatable :: mob_op(:,:,:)
      real(8),allocatable :: flow_factor(:,:,:)
      integer(4) :: numghost
      integer(4) :: ngx, ngy, ngz
      integer(4) :: il,iu, jl,ju, kl,ku
      integer(4) :: ilb, iub, jlb, jub, klb, kub

      integer(4) :: i,j,k,nn, ll
      integer(4) :: ig, jg, kg
      integer(4) :: ierr
      integer(4) :: pidx, pidy, pidz
      real(8) :: wtime_start, wtime_end, wtime 

      real(8) :: total_c 

      character(len=100) :: bctype 

      bctype = 'ADIABATIC'

      NPX = 1
      NPY = 5
      NPZ = 1

      imori = 100
      jmori = 100
      kmori = 1
      np = 1

      DXL = 0.5
      DYL = 0.5 
      DZL = 0.5

      call get_dimension(imori, jmori, kmori)

      numghost = 1
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
      allocate(pf(ilb:iub,jlb:jub,klb:kub))
      allocate(mob(ilb:iub,jlb:jub,klb:kub))
      allocate(mob_op(ilb:iub,jlb:jub,klb:kub))
      allocate(flow_factor(ilb:iub,jlb:jub,klb:kub))

      !initialization
      a(:,:,:,:) = 0
      b(:,:,:,:) = 0.0
      c(:,:,:,:) = 0.0
      pf(:,:,:) = 1 ! Compressible background phase
      mob(:,:,:) = Mobility
      flow_factor(:,:,:) = 1.0d0  ! Uniform flow factor (full flow)

      call phaseid_ic_box(pf,ngx,ngy,ngz, 10, imori, 1, jmori, 1, 1, 2) ! Incompressible liquid phase 
      call phaseid_ic_box(pf,ngx,ngy,ngz, imori/2-5, imori/2+5, 
     &                                               1, jmori, 1, 1, 3) ! Incompressible liquid phase

      !call distance_from_center(c(:,:,:,told(0)),ngx,ngy,ngz,1)

    !   call diffu_ic_stepfunction_x(c(:,:,:,told(0)),ngx,ngy,ngz, 
    !  &                                                    1, 5, 100.0d0)
    !   call diffu_ic_stepfunction_x(c(:,:,:,tnew(0)),ngx,ngy,ngz, 
    !  &                                                    1, 5, 100.0d0)

      call diffu_ic_by_phase(c(:,:,:,told(0)),pf(:,:,:),ngx,ngy,ngz,
     &                                          0.5d0, 0.15d0, 0.25d0)
      call diffu_ic_by_phase(c(:,:,:,tnew(0)),pf(:,:,:),ngx,ngy,ngz,
     &                                          0.5d0, 0.15d0, 0.25d0)

      wtime_start = MPI_WTIME()

      call init_bc_call_count(0)

      ttime = 0.0d0
      do i=0,1000000
        call init_bc_call_count(0)
        !call boundary_condition_r8('PERIODIC',c(:,:,:,told(i))          ! wtime_mpif08.csv
        call boundary_condition_i4(bctype,pf(:,:,:), numghost)
        call boundary_condition_r8(bctype,c(:,:,:,told(i)), numghost)
        call boundary_condition_r8(bctype,mob(:,:,:), numghost)
        
        call calc_diffusion_dt(dttime, mob, DXL)
        call calc_diffusion_mobility(mob_op,flow_factor,mob,
     &                          c(:,:,:,told(i)),pf,ngx,ngy,ngz, dttime)
        call boundary_condition_r8(bctype,mob_op(:,:,:), numghost) 

    !     call simple_diffusion(c(:,:,:,tnew(i)),c(:,:,:,told(i)), 
    !  &                                                  ngx,ngy,ngz, i)
    !     call diffusion_sourced(c(:,:,:,tnew(i)),c(:,:,:,told(i)), 
    !  &                                    ngx,ngy,ngz, i,  0.0d0, 0.3d0)

    !     call diffusion_biphase_box(c(:,:,:,tnew(i)),c(:,:,:,told(i)),
    !  &                                                    ngx,ngy,ngz,i, 
    !  &                 imori/2-5, imori/2+5, jmori/2-5, jmori/2+5, 1, 1, 
    !  &                                                   0.3d0, 100.0d0)

        call diffusion_by_flux(c(:,:,:,tnew(i)),c(:,:,:,told(i)),
     &                                  mob_op,pf,ngx,ngy,ngz,i, dttime)

        if(mod(i,1000).eq.0 .or. error_flag)then
          !call boundary_condition_r8('PERIODIC',c(:,:,:,tnew(i))        ! wtime_mpif08.csv
          call boundary_condition_r8(bctype,c(:,:,:,tnew(i)),numghost)  ! wtime_modern.csv
    !       call output_r8_pvtr_2d(c(:,:,:,tnew(i)),numghost,'c',i,ttime,  
    !  &                                              myrank ,"test pvtr")
    !       call output_r8_pvtr(c(:,:,:,tnew(i)),numghost,'c',i,ttime,  
    !  &                                              myrank ,"test pvtr")
          call output_r8_pvtr(c(:,:,:,tnew(i)),numghost,'c',  
     &                                     i,ttime, myrank ,"test pvtr")
          call output_i4_pvtr(pf(:,:,:),numghost,'pf',i,ttime,myrank,
     &                                                        "phaseID")
          call output_r8_pvtr(mob_op(:,:,:),numghost,'mob',i,ttime,
     &                                                myrank,"mobility")

          total_c = total_amount_c(c(:,:,:,tnew(i)),ngx,ngy,ngz)

          if(myrank.eq.0)then
            !write(*,*)'Time step = ', i, 'Output made' 
            write(*,*)'Total concentration at time step ', i, ' = ', 
     &                                                           total_c
            call verbose_bc_call_count()
            call verbose_bc_max_mpi_tag()
          endif
        endif

        if(error_flag) then
          if(myrank.eq.0) then
            write(*,*)'Error detected at time step ', i
            write(*,*)'  Error log message:'
            write(*,*)'  !!! '//trim(error_log)
            write(*,*)'  Output made up to the error point.'
            write(*,*)'  Aborting...'
          endif
          exit
        endif
 
        ttime = ttime + dttime
      enddo !do i=0,100000

      wtime_end = MPI_WTIME()

      wtime = wtime_end - wtime_start

      if (myrank.eq.0) then 
        write(*,*)'Wall time spent (s): ', wtime
      endif 

      deallocate(a)
      deallocate(b)
      deallocate(c)
      deallocate(pf)
      deallocate(mob)
      deallocate(mob_op)

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