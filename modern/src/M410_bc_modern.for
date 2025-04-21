!      ! Author: Dong-Uk Kim (GitHub: eastglow-zz)

!      ! ========= Example usage of boundary_conditions module: =========
!      
!      program example_boundary_conditions
!      use mpi_f08  ! Necessary dependency
!      use boundary_conditions 
!      implicit none  
!
!      integer(4) :: npx, npy, npz 
!      integer(4) :: imori, jmori, kmori 
!      integer(4) :: nprocs
!      integer(4) :: myrank  
!      integer(4) :: numghost 
!      integer(4) :: ierr 
!      type(MPI_Comm) :: comm 
!      integer(4),allocatable :: integer4_array(:,:,:)
!      real(8),allocatable :: real8_array(:,:,:)
!
!      npx = 2  ! # of parallel partitioning along x-axis
!      npy = 2  ! # of parallel partitioning along y-axis
!      npz = 1  ! # of parallal partitioning along z-axis
!      imori = 32  ! # of global grid along x-axis 
!      jmori = 32  ! # of global grid along y-axis
!      kmori = 2   ! # of global grid along z-axis 
!      
!      numghost = 1 ! # of ghost layer, maximum: 2
!
!      ! Initialization of MPI routines 
!      comm = MPI_COMM_WORLD
!      call MPI_INIT(ierr)
!      call MPI_COMM_SIZE(comm,nprocs,ierr)
!      call MPI_COMM_RANK(comm,myrank,ierr)
!
!      ! Initialization of boundary_conditions module (necessary)
!      call init_bc(NPX,NPY,NPZ,IMORI,JMORI,KMORI,myrank, comm)
!      ! Now you can use bc_iml, bc_jml, bc_kml, and bc_dimx, bc_dimy, bc_dimz to allocate the local data arrays
!      ! bc_iml, bc_jml, bc_kml: the number of local grid along x,y,z axis respectively 
!      ! bc_dimx, bc_dimy, bc_dimz: dimension on/off switch integer.
!      !   in 3D case, bc_dimx=1, bc_dimy=1, bc_dimz=1.
!      !   if 2D in xy plane, bc_dimx=1, bc_dimy=1, bc_dimz=0, 
!      !     and etc...
!
!      ! Local data array allocations (necessary)
!      allocate(
!        integer4_array(1-numghost*bc_dimx:bc_iml+numghost*bc_dimx,
!     &                 1-numghost*bc_dimy:bc_jml+numghost*bc_dimy,
!     &                 1-numghost*bc_dimz:bc_kml+numghost*bc_dimz ) )
!
!      allocate(
!        real8_array(   1-numghost*bc_dimx:bc_iml+numghost*bc_dimx,
!     &                 1-numghost*bc_dimy:bc_jml+numghost*bc_dimy,
!     &                 1-numghost*bc_dimz:bc_kml+numghost*bc_dimz ) )
!
!      integer4_array(:,:,:) = myrank 
!      real8_array(:,:,:) = myrank
!
!      ! Use of boundary condition routines (necessary)
!      call boundary_condition_i4('ADIABATIC',
!     &                                   integer4_array(:,:,:),numghost)
!
!      call boundary_condition_r8('ADIABATIC',
!     &                                      real8_array(:,:,:),numghost)
!
!      ! reinitialization of MPI tag values (necessary)
!      call bc_call_count(0)
!
!      ! Finalizing boundary_conditions module (necessary)
!      call finalize_bc()
!
!      call MPI_FINALIZE(ierr)
!
!      end program example_boundary_conditions

!      ! ================================================================
!      ! ========= Available boundary condition routines:       =========

!      boundary_condition_r8(bctype, array3d_name, num_ghost_grid) ! For real(8) arrays
!      boundary_condition_r4(bctype, array3d_name, num_ghost_grid) ! For real(8) arrays 
!      boundary_condition_i4(bctype, array3d_name, num_ghost_grid) ! For integer(4) arrays

!      ! ================================================================
!      ! ========= Available boundary condition types (all UPPER CASE LETTERS): 

!      PERIODIC
!      ADIABATIC
!      ADIABATIC_X
!      ADIABATIC_Y
!      ADIABATIC_Z
!      ADIABATIC_XY
!      ADIABATIC_YZ
!      ADIABATIC_ZX

!      ! ================================================================
      
      module boundary_conditions
      use mpi_f08
      implicit none

      integer(4) :: bc_call_count
      integer(4) :: bc_max_mpi_tag = 0
      integer(4),parameter :: bc_mpitag_jumpstep = 13
      integer(4),parameter :: bc_ng_max = 2

      integer(4) :: bc_dimen 
      integer(4) :: bc_dimx, bc_dimy, bc_dimz

      integer(4) :: bc_npx=1, bc_npy=1, bc_npz=1 !Number of partitioning along x,y, and z axis
      integer(4) :: bc_iml, bc_jml, bc_kml 
      integer(4) :: bc_imori, bc_jmori, bc_kmori
      integer(4) :: bc_il, bc_iu, bc_jl, bc_ju, bc_kl, bc_ku
      integer(4) :: bc_pidg
      integer(4) :: bc_pidg_xm, bc_pidg_xp
      integer(4) :: bc_pidg_ym, bc_pidg_yp
      integer(4) :: bc_pidg_zm, bc_pidg_zp
      integer(4) :: bc_pidx, bc_pidy, bc_pidz

      type(MPI_Comm) :: bc_comm

      type(MPI_Datatype) :: bc_xslab_r8(1:bc_ng_max)
      type(MPI_Datatype) :: bc_yslab_r8(1:bc_ng_max)
      type(MPI_Datatype) :: bc_zslab_r8(1:bc_ng_max)
      type(MPI_Datatype) :: bc_xslab_r4(1:bc_ng_max)
      type(MPI_Datatype) :: bc_yslab_r4(1:bc_ng_max)
      type(MPI_Datatype) :: bc_zslab_r4(1:bc_ng_max)
      type(MPI_Datatype) :: bc_xslab_i4(1:bc_ng_max)
      type(MPI_Datatype) :: bc_yslab_i4(1:bc_ng_max)
      type(MPI_Datatype) :: bc_zslab_i4(1:bc_ng_max)
      integer(4) :: bc_idx_start(3) = (/0, 0, 0/)
      integer(4) :: bc_arrsize(3,1:bc_ng_max)
      integer(4) :: bc_subsize_x(3,1:bc_ng_max)
      integer(4) :: bc_subsize_y(3,1:bc_ng_max)
      integer(4) :: bc_subsize_z(3,1:bc_ng_max)


      logical :: bc_initialized = .false.
      

      contains

      ! ----------------------------------------------------------------

      subroutine init_bc(inpx,inpy,inpz,i_all,j_all,k_all,irank,comm)

      ! inpx: The number of parallel partitioning along x-axis
      ! inpy: The number of parallel partitioning along y-axis 
      ! inpz: The number of parallel partitioning along z-axis 
      ! i_all: The number of grid before the parallel partitioning, along x-axis 
      ! j_all: The number of grid before the parallel partitioning, along y-axis
      ! k_all: The number of grid before the parallel partitioning, along z-axis 
      ! irank: Rank of this process 
      ! comm: MPI communicator handle

      implicit none
      integer(4),intent(in) :: inpx,inpy,inpz !saved in bc_npx, bc_npy, bc_npz
      integer(4),intent(in) :: i_all, j_all, k_all !saved in bc_imori,bc_jmori,bc_kmori
      integer(4),intent(in) :: irank !saved in bc_pidg
      type(MPI_Comm),intent(in) :: comm

      integer(4) :: ilb,iub,jlb,jub,klb,kub
      integer(4) :: pidx, pidy, pidz 
      integer(4) :: ng
      integer(4) :: ierr

      bc_npx = inpx
      bc_npy = inpy
      bc_npz = inpz

      bc_imori = i_all
      bc_jmori = j_all
      bc_kmori = k_all

      call bc_get_dimension(i_all,j_all,k_all)
      
      ! Saving the neighbor process ids, assuming the periodicity
      bc_pidg = irank
      call get_pidxyz(pidx, pidy, pidz, bc_pidg)
      bc_pidx = pidx 
      bc_pidy = pidy 
      bc_pidz = pidz

      bc_initialized = .true. ! This statement should be here. Not coming after get_pid()

      bc_pidg_xm = get_pid(bc_pidx-1,bc_pidy,bc_pidz)
      bc_pidg_xp = get_pid(bc_pidx+1,bc_pidy,bc_pidz)
      bc_pidg_ym = get_pid(bc_pidx,bc_pidy-1,bc_pidz)
      bc_pidg_yp = get_pid(bc_pidx,bc_pidy+1,bc_pidz)
      bc_pidg_zm = get_pid(bc_pidx,bc_pidy,bc_pidz-1)
      bc_pidg_zp = get_pid(bc_pidx,bc_pidy,bc_pidz+1)

      bc_comm = comm

      call get_ijkrange(ilb,iub,jlb,jub,klb,kub,bc_pidg)
      bc_iml=iub-ilb+1
      bc_jml=jub-jlb+1
      bc_kml=kub-klb+1
      bc_il = ilb 
      bc_iu = iub 
      bc_jl = jlb 
      bc_ju = jub 
      bc_kl = klb 
      bc_ku = kub 

      do ng=1,bc_ng_max
        bc_arrsize(1,ng) = 2*ng + bc_iml
        bc_arrsize(2,ng) = 2*ng + bc_jml
        bc_arrsize(3,ng) = 2*ng + bc_kml

        ! For 2D (or 1D - not guaranteed)
        if(bc_dimx .eq. 0) bc_arrsize(1,ng) = 1
        if(bc_dimy .eq. 0) bc_arrsize(2,ng) = 1
        if(bc_dimz .eq. 0) bc_arrsize(3,ng) = 1

        bc_subsize_x(1,ng) = ng
        bc_subsize_x(2,ng) = bc_arrsize(2,ng)
        bc_subsize_x(3,ng) = bc_arrsize(3,ng)

        bc_subsize_y(1,ng) = bc_arrsize(1,ng)
        bc_subsize_y(2,ng) = ng
        bc_subsize_y(3,ng) = bc_arrsize(3,ng)

        bc_subsize_z(1,ng) = bc_arrsize(1,ng)
        bc_subsize_z(2,ng) = bc_arrsize(2,ng)
        bc_subsize_z(3,ng) = ng

        ! For real(8) type arrays
        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_x(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                      MPI_DOUBLE_PRECISION, bc_xslab_r8(ng), ierr)
        call mpi_type_commit(bc_xslab_r8(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_y(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                      MPI_DOUBLE_PRECISION, bc_yslab_r8(ng), ierr)
        call mpi_type_commit(bc_yslab_r8(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_z(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                      MPI_DOUBLE_PRECISION, bc_zslab_r8(ng), ierr)
        call mpi_type_commit(bc_zslab_r8(ng), ierr)

        ! For real(4) type arrays
        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_x(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                                  MPI_REAL, bc_xslab_r4(ng), ierr)
        call mpi_type_commit(bc_xslab_r4(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_y(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                                  MPI_REAL, bc_yslab_r4(ng), ierr)
        call mpi_type_commit(bc_yslab_r4(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_z(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                                  MPI_REAL, bc_zslab_r4(ng), ierr)
        call mpi_type_commit(bc_zslab_r4(ng), ierr)

        ! For integer(4) type arrays
        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_x(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                               MPI_INTEGER, bc_xslab_i4(ng), ierr)
        call mpi_type_commit(bc_xslab_i4(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_y(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                               MPI_INTEGER, bc_yslab_i4(ng), ierr)
        call mpi_type_commit(bc_yslab_i4(ng), ierr)

        call mpi_type_create_subarray(3, bc_arrsize(:,ng), 
     &             bc_subsize_z(:,ng), bc_idx_start, MPI_ORDER_FORTRAN,
     &                               MPI_INTEGER, bc_zslab_i4(ng), ierr)
        call mpi_type_commit(bc_zslab_i4(ng), ierr)
      enddo

      call init_bc_call_count(0)

      end subroutine init_bc

      ! ----------------------------------------------------------------

      subroutine finalize_bc()
      implicit none 

      integer(4) :: ng
      integer(4) :: ierr

      do ng=1,bc_ng_max
        call mpi_type_free(bc_xslab_r8(ng),ierr)
        call mpi_type_free(bc_yslab_r8(ng),ierr)
        call mpi_type_free(bc_zslab_r8(ng),ierr)
        call mpi_type_free(bc_xslab_r4(ng),ierr)
        call mpi_type_free(bc_yslab_r4(ng),ierr)
        call mpi_type_free(bc_zslab_r4(ng),ierr)
        call mpi_type_free(bc_xslab_i4(ng),ierr)
        call mpi_type_free(bc_yslab_i4(ng),ierr)
        call mpi_type_free(bc_zslab_i4(ng),ierr)
      enddo 

      return 
      end subroutine finalize_bc

      ! ----------------------------------------------------------------

      subroutine bc_get_dimension(im, jm, km)
      implicit none 

      integer(4),intent(in) :: im, jm, km

      bc_dimx = 1 
      bc_dimy = 1
      bc_dimz = 1
      if(im .le. 1) bc_dimx = 0
      if(jm .le. 1) bc_dimy = 0
      if(km .le. 1) bc_dimz = 0 
      bc_dimen = bc_dimx + bc_dimy + bc_dimz
      ! write(*,*)'Dimension:', dimen
      ! write(*,*)'dimx, dimy, dimz:', dimx, dimy, dimz
      ! write(*,*)'imori, jmori, kmori:', imori, jmori, kmori
      if(bc_dimen .eq. 0) then
        write(*,*)'Error, bc_get_dimension() in M410_bc_modern.for'
        write(*,*)'  No dimension specified'
        write(*,*)'  Aborting...'
        call abort()
      end if

      return 
      end subroutine bc_get_dimension

      ! ----------------------------------------------------------------

      subroutine init_bc_call_count(val)
      implicit none
      integer(4),intent(in) :: val

       bc_call_count = val

      end subroutine init_bc_call_count

      ! ----------------------------------------------------------------


      subroutine verbose_bc_call_count()
      implicit none
      
       write(*,*) 'bc_call_count:', bc_call_count

      end subroutine verbose_bc_call_count

      ! ----------------------------------------------------------------

      subroutine verbose_bc_max_mpi_tag()
      implicit none

      write(*,*) 'bc_max_mpi_tag:', bc_max_mpi_tag

      end subroutine verbose_bc_max_mpi_tag

      ! ----------------------------------------------------------------

      subroutine monitor_bc_max_mpi_tag()
      implicit none
      !include 'mpif.h'
      if(bc_max_mpi_tag.ge.32767)then
            write(*,*)'Reminder: bc_max_mpi_tag exceeded 32767.',
     &           bc_max_mpi_tag
            if(bc_max_mpi_tag.ge.MPI_TAG_UB)then
            write(*,*)'ERROR: bc_max_mpi_tag exceeded MPI_TAG_UB.',
     &        bc_max_mpi_tag
            write(*,*)'MPI_TAG_UB for this system:',MPI_TAG_UB
            write(*,*)'Aborting due to high bc_max_mpi_tag.'
            call abort
            endif
      endif
      end subroutine monitor_bc_max_mpi_tag

      ! ----------------------------------------------------------------

      subroutine bound_check(flagid, i,j,k, ng)
      implicit none
      character(len=*),intent(in) :: flagid
      integer(4),intent(in) :: i,j,k, ng

      logical :: go_nogo
      character(len=100) :: err_sign

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      err_sign=" "
      go_nogo = .true.

      if(i.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(iL)"
      endif
      if(i.gt.bc_iml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(iU)"
      endif
      if(j.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jL)"
      endif
      if(j.gt.bc_jml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jU)"
      endif
      if(k.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kL)"
      endif
      if(k.gt.bc_kml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kU)"
      endif

      if(go_nogo.eqv. .false.) then
            write(*,*)"Tried to access the out of array bound:"
            write(*,*)"flagid=",trim(flagid),", i,j,k = ",i,j,k
            write(*,*)"err_sign: ",trim(err_sign)," ,rank=",bc_pidg
            call abort
      endif

      end subroutine bound_check

      ! ----------------------------------------------------------------

      integer(4) function get_pid(pidx, pidy, pidz)
      implicit none
      integer(4),intent(in) :: pidx, pidy, pidz

      integer(4) :: pidxw, pidyw, pidzw

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      pidxw = pidx
      if(pidxw.gt.bc_npx-1)pidxw=0
      if(pidxw.lt.0)pidxw=bc_npx-1

      pidyw = pidy
      if(pidyw.gt.bc_npy-1)pidyw=0
      if(pidyw.lt.0)pidyw=bc_npy-1

      pidzw = pidz
      if(pidzw.gt.bc_npz-1)pidzw=0
      if(pidzw.lt.0)pidzw=bc_npz-1

      get_pid = pidxw + pidyw*bc_npx + pidzw*bc_npx*bc_npy
      return
      end function get_pid

      ! ----------------------------------------------------------------

      subroutine get_pidxyz(pidx, pidy, pidz, pid)
      implicit none
      integer(4),intent(inout) :: pidx, pidy, pidz
      integer(4),intent(in) :: pid

      pidz = int(pid/bc_npx/bc_npy)
      pidy = int((pid - pidz*bc_npx*bc_npy)/bc_npx)
      pidx = pid - pidz*bc_npx*bc_npy - pidy*bc_npx

      return
      end subroutine get_pidxyz

      ! ----------------------------------------------------------------

      subroutine get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      implicit none
      integer(4),intent(inout) :: il,iu,jl,ju,kl,ku
      integer(4),intent(in) :: pid

      integer(4) :: pidx, pidy, pidz


      call get_pidxyz(pidx, pidy, pidz, pid)

      ! x bounds
      call para_range(1,bc_imori,bc_npx,pidx,il,iu)
      ! y bounds
      call para_range(1,bc_jmori,bc_npy,pidy,jl,ju)
      ! z bounds
      call para_range(1,bc_kmori,bc_npz,pidz,kl,ku)

      return
      end subroutine get_ijkrange

      ! ----------------------------------------------------------------

      subroutine get_global_ijk(ig,jg,kg, ilocal,jlocal,klocal, pid)
      implicit none
      integer(4),intent(inout) :: ig, jg, kg !global ijk index
      integer(4),intent(in) :: ilocal, jlocal, klocal !local ijk index
      integer(4),intent(in) :: pid

      integer(4) :: il,iu, jl,ju, kl,ku

      call get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      ig = il + (ilocal-1)  !for spatial index starts from 1, not 0
      jg = jl + (jlocal-1)  !for spatial index starts from 1, not 0
      kg = kl + (klocal-1)  !for spatial index starts from 1, not 0

      return
      end subroutine get_global_ijk

      ! ----------------------------------------------------------------

      SUBROUTINE para_range(n1,n2,npart,irank,ista,iend)
      ! Deviding system size bc_kmori into KM = bc_kmori/(porcess number)
      ! npart: total number of porcesses (npart=48 in my computer)
      ! irank: current process name (irank=1~47)
      ! ista and iend: first and last grids (J-direction) on each process in the total grid system
      ! use: CALL para_range(1,bc_kmori,nprocs,myrank,ista,iend)
      implicit none
      integer(4),intent(in) :: n1,n2,npart
      integer(4),intent(inout) :: irank,ista,iend

      integer(4) :: iwork1, iwork2

      iwork1=(n2-n1+1)/npart
      iwork2=MOD(n2-n1+1,npart)
      ista=irank*iwork1+n1+MIN(irank,iwork2)
      iend=ista+iwork1-1
      if(iwork2.gt.irank) iend=iend+1

      return
      END SUBROUTINE para_range

      ! ----------------------------------------------------------------

      integer(4) function get_pid_from_globalijk(ig,jg,kg,bctype)
      implicit none
      integer(4),intent(in) :: ig,jg,kg
      character(len=*),intent(in) :: bctype
      
      integer(4) :: igf, jgf, kgf
      integer(4) :: pidx, pidy, pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif

      call clamp_globalijk(igf,jgf,kgf,ig,jg,kg,bctype)

      pidx = int((igf-1)/bc_iml)
      pidy = int((jgf-1)/bc_jml)
      pidz = int((kgf-1)/bc_kml)

      get_pid_from_globalijk = get_pid(pidx, pidy, pidz)
      return
      end function get_pid_from_globalijk

      ! ----------------------------------------------------------------

      subroutine get_localijk_from_globalijk(il,jl,kl,ig,jg,kg,bctype)
      implicit none
      integer(4),intent(inout) :: il,jl,kl
      integer(4),intent(in) :: ig,jg,kg
      character(len=*) :: bctype

      integer(4) :: igf, jgf, kgf
      integer(4) :: pid, pidx, pidy, pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      pid = get_pid_from_globalijk(ig,jg,kg,bctype)
      call get_pidxyz(pidx, pidy, pidz, pid)

      call clamp_globalijk(igf,jgf,kgf,ig,jg,kg,bctype)

      il=igf-pidx*bc_iml
      jl=jgf-pidy*bc_jml
      kl=kgf-pidz*bc_kml

      return
      end subroutine get_localijk_from_globalijk

      ! ----------------------------------------------------------------

      subroutine clamp_globalijk(igo,jgo,kgo,ig,jg,kg,bctype)
      integer(4),intent(inout) :: igo, jgo, kgo
      integer(4),intent(in) :: ig,jg,kg
      character(len=*) :: bctype

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif

      igo=ig
      jgo=jg
      kgo=kg

      select case(trim(bctype))
      case('PERIODIC')
            igo=mod(ig+min(ig-1,0)*(-bc_imori)-1,bc_imori)+1
            jgo=mod(jg+min(jg-1,0)*(-bc_jmori)-1,bc_jmori)+1
            kgo=mod(kg+min(kg-1,0)*(-bc_kmori)-1,bc_kmori)+1
      case('ADIABATIC')
            igo=min(max(ig,1),bc_imori)
            jgo=min(max(jg,1),bc_jmori)
            kgo=min(max(kg,1),bc_kmori)
      case('ADIABATIC_X')
            igo=min(max(ig,1),bc_imori)
            jgo=mod(jg+min(jg-1,0)*(-bc_jmori)-1,bc_jmori)+1
            kgo=mod(kg+min(kg-1,0)*(-bc_kmori)-1,bc_kmori)+1
      case('ADIABATIC_Y')
            igo=mod(ig+min(ig-1,0)*(-bc_imori)-1,bc_imori)+1
            jgo=min(max(jg,1),bc_jmori)
            kgo=mod(kg+min(kg-1,0)*(-bc_kmori)-1,bc_kmori)+1
      case('ADIABATIC_Z')
            igo=mod(ig+min(ig-1,0)*(-bc_imori)-1,bc_imori)+1
            jgo=mod(jg+min(jg-1,0)*(-bc_jmori)-1,bc_jmori)+1
            kgo=min(max(kg,1),bc_kmori)
      case('ADIABATIC_XY')
            igo=min(max(ig,1),bc_imori)
            jgo=min(max(jg,1),bc_jmori)
            kgo=mod(kg+min(kg-1,0)*(-bc_kmori)-1,bc_kmori)+1
      case('ADIABATIC_YZ')
            igo=mod(ig+min(ig-1,0)*(-bc_imori)-1,bc_imori)+1
            jgo=min(max(jg,1),bc_jmori)
            kgo=min(max(kg,1),bc_kmori)
      case('ADIABATIC_ZX')
            igo=min(max(ig,1),bc_imori)
            jgo=mod(jg+min(jg-1,0)*(-bc_jmori)-1,bc_jmori)+1
            kgo=min(max(kg,1),bc_kmori)
      end select

      return
      end subroutine clamp_globalijk


      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_x_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_x_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_x_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_y_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_y_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_y_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_z_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_z_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_z_bc_r8

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_r8(a,ngx,ngy,ngz)
      implicit none
      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif
      
      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      end subroutine apply_periodic_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_periodic_x_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_x_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_x_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_periodic_y_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_y_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_y_bc_r8

      ! ----------------------------------------------------------------

      subroutine apply_periodic_z_bc_r8(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(8),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_z_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_z_bc_r8

      ! ----------------------------------------------------------------

      subroutine boundary_condition_r8(bctype,a,ng)
      implicit none

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a( 1-ng*bc_dimx:bc_iml+ng*bc_dimx,
     &                            1-ng*bc_dimy:bc_jml+ng*bc_dimy,
     &                            1-ng*bc_dimz:bc_kml+ng*bc_dimz )

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      type(MPI_Request) :: isend01, isend02
      type(MPI_Request) :: isend03, isend04
      type(MPI_Request) :: isend05, isend06
      type(MPI_Request) :: irecv01, irecv02
      type(MPI_Request) :: irecv03, irecv04
      type(MPI_Request) :: irecv05, irecv06

      type(MPI_Status) :: istatus
      integer(4) :: ierr
      integer(4) :: basetag
      integer(4) :: ngx, ngy, ngz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_r8_modern:'
        write(*,*)'boundary condition uninitialized. Aborting.'
        call abort
      endif

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      ngx = bc_dimx*ng
      ngy = bc_dimy*ng
      ngz = bc_dimz*ng

      il = 1-ngx
      iu = bc_iml+ngx
      jl = 1-ngy
      ju = bc_jml+ngy
      kl = 1-ngz
      ku = bc_kml+ngz

      ! Internal transfers
      ! i planes
      if(bc_npx.gt.1.and.ngx.gt.0)then
        pidto01=bc_pidg_xm ! il_s
        call MPI_ISEND(a(1,jl,kl),       1, bc_xslab_r8(ngx), pidto01,
     &                                basetag+1, bc_comm, isend01, ierr)

        pidfrom02=bc_pidg_xp ! iu_r
        call MPI_IRECV(a(bc_iml+1,jl,kl),1, bc_xslab_r8(ngx), pidfrom02,
     &                                basetag+1, bc_comm, irecv02, ierr)

        pidto02=bc_pidg_xp ! iu_s
        call MPI_ISEND(a(bc_iml-ngx+1,jl,kl),1,bc_xslab_r8(ngx),pidto02,
     &                                basetag+2, bc_comm, isend02, ierr)      

        pidfrom01=bc_pidg_xm ! il_r
        call MPI_IRECV(a(il,jl,kl),      1, bc_xslab_r8(ngx), pidfrom01,
     &                                basetag+2, bc_comm, irecv01, ierr)
            
        call MPI_WAIT(isend01,istatus,ierr)
        call MPI_WAIT(isend02,istatus,ierr)
        call MPI_WAIT(irecv01,istatus,ierr)
        call MPI_WAIT(irecv02,istatus,ierr)
      endif !if(bc_npx.gt.1.and.ngx.gt.0)then
      ! Done only for x direction - not complete yet

      ! j planes
      if(bc_npy.gt.1.and.ngy.gt.0)then
        pidto03=bc_pidg_ym ! jl_s
        call MPI_ISEND(a(il,1,kl), 1,        bc_yslab_r8(ngy), pidto03,
     &                                basetag+3, bc_comm, isend03, ierr)

        pidfrom04=bc_pidg_yp ! ju_r
        call MPI_IRECV(a(il,bc_jml+1,kl), 1,bc_yslab_r8(ngy), pidfrom04, 
     &                                basetag+3, bc_comm, irecv04, ierr)

        pidto04=bc_pidg_yp ! ju_s
        call MPI_ISEND(a(il,bc_jml-ngy+1,kl),1,bc_yslab_r8(ngy),pidto04,
     &                                basetag+4, bc_comm, isend04, ierr)

        pidfrom03=bc_pidg_ym ! jl_r      
        call MPI_IRECV(a(il,jl,kl), 1,     bc_yslab_r8(ngy), pidfrom03,
     &                                basetag+4, bc_comm, irecv03, ierr)

        call MPI_WAIT(isend03,istatus,ierr)
        call MPI_WAIT(isend04,istatus,ierr)
        call MPI_WAIT(irecv03,istatus,ierr)
        call MPI_WAIT(irecv04,istatus,ierr)
      endif !if(bc_npy.gt.1.and.ngy.gt.0)then
      ! Done for x and y direction - not completed yet

      ! k planes
      if(bc_npz.gt.1.and.ngz.gt.0)then
        pidto05=bc_pidg_zm ! kl_s
        call MPI_ISEND(a(il,jl,1), 1,        bc_zslab_r8(ngz), pidto05,
     &                                basetag+5, bc_comm, isend05, ierr)

        pidfrom06=bc_pidg_zp ! ku_r
        call MPI_IRECV(a(il,jl,bc_kml+1), 1,bc_zslab_r8(ngz), pidfrom06,
     &                                basetag+5, bc_comm, irecv06, ierr)
            
        pidto06=bc_pidg_zp ! ku_s
        call MPI_ISEND(a(il,jl,bc_kml-ngz+1),1,bc_zslab_r8(ngz),pidto06,
     &                                basetag+6, bc_comm, isend06, ierr)

        pidfrom05=bc_pidg_zm ! kl_r
        call MPI_IRECV(a(il,jl,kl), 1,      bc_zslab_r8(ngz), pidfrom05,
     &                                basetag+6, bc_comm, irecv05, ierr)

        call MPI_WAIT(isend05,istatus,ierr)
        call MPI_WAIT(isend06,istatus,ierr)
        call MPI_WAIT(irecv05,istatus,ierr)
        call MPI_WAIT(irecv06,istatus,ierr)
      endif !if(bc_npz.gt.1.and.ngz.gt.0)then
      ! Internal tranfoer done for all directions - complete

      ! Applying the actual boundary conditions
      select case(trim(bctype))
      case('PERIODIC')
        call apply_periodic_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC')
        call apply_adiabatic_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_XY')
        call apply_periodic_z_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_YZ')
        call apply_periodic_x_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_ZX')
        call apply_periodic_y_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_X')
        call apply_periodic_y_bc_r8(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_Y')
        call apply_periodic_x_bc_r8(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r8(a, ngx, ngy, ngz)
      case('ADIABATIC_Z')
        call apply_periodic_x_bc_r8(a, ngx, ngy, ngz)
        call apply_periodic_y_bc_r8(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r8(a, ngx, ngy, ngz)
      case default 
        write(*,*)'boundary_conditions module: fatal error:'
        write(*,*)'  bctype = ',trim(bctype)
        write(*,*)'  unknown boundary condition type. Aborting.'
        call abort
      end select !select case(trim(bctype))

      return 
      end subroutine boundary_condition_r8

      ! ----------------------------------------------------------------


      subroutine apply_adiabatic_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_x_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_x_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_x_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_y_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_y_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_y_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_z_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_z_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_z_bc_r4

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_r4(a,ngx,ngy,ngz)
      implicit none
      integer(4),intent(in) :: ngx,ngy,ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif
      
      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      end subroutine apply_periodic_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_x_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_x_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_x_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_y_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_y_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_y_bc_r4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_z_bc_r4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      real(4),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_z_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_z_bc_r4

      ! ----------------------------------------------------------------

      subroutine boundary_condition_r4(bctype,a,ng)
      implicit none

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a( 1-ng*bc_dimx:bc_iml+ng*bc_dimx,
     &                            1-ng*bc_dimy:bc_jml+ng*bc_dimy,
     &                            1-ng*bc_dimz:bc_kml+ng*bc_dimz )

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      type(MPI_Request) :: isend01, isend02
      type(MPI_Request) :: isend03, isend04
      type(MPI_Request) :: isend05, isend06
      type(MPI_Request) :: irecv01, irecv02
      type(MPI_Request) :: irecv03, irecv04
      type(MPI_Request) :: irecv05, irecv06

      type(MPI_Status) :: istatus
      integer(4) :: ierr
      integer(4) :: basetag
      integer(4) :: ngx, ngy, ngz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_r8_modern:'
        write(*,*)'boundary condition uninitialized. Aborting.'
        call abort
      endif

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      ngx = bc_dimx*ng
      ngy = bc_dimy*ng
      ngz = bc_dimz*ng

      il = 1-ngx
      iu = bc_iml+ngx
      jl = 1-ngy
      ju = bc_jml+ngy
      kl = 1-ngz
      ku = bc_kml+ngz

      ! Internal transfers
      ! i planes
      if(bc_npx.gt.1.and.ngx.gt.0)then
        pidto01=bc_pidg_xm ! il_s
        call MPI_ISEND(a(1,jl,kl),       1, bc_xslab_r4(ngx), pidto01,
     &                                basetag+1, bc_comm, isend01, ierr)

        pidfrom02=bc_pidg_xp ! iu_r
        call MPI_IRECV(a(bc_iml+1,jl,kl),1, bc_xslab_r4(ngx), pidfrom02,
     &                                basetag+1, bc_comm, irecv02, ierr)

        pidto02=bc_pidg_xp ! iu_s
        call MPI_ISEND(a(bc_iml-ngx+1,jl,kl),1,bc_xslab_r4(ngx),pidto02,
     &                                basetag+2, bc_comm, isend02, ierr)      

        pidfrom01=bc_pidg_xm ! il_r
        call MPI_IRECV(a(il,jl,kl),      1, bc_xslab_r4(ngx), pidfrom01,
     &                                basetag+2, bc_comm, irecv01, ierr)
            
        call MPI_WAIT(isend01,istatus,ierr)
        call MPI_WAIT(isend02,istatus,ierr)
        call MPI_WAIT(irecv01,istatus,ierr)
        call MPI_WAIT(irecv02,istatus,ierr)
      endif !if(bc_npx.gt.1.and.ngx.gt.0)then
      ! Done only for x direction - not complete yet

      ! j planes
      if(bc_npy.gt.1.and.ngy.gt.0)then
        pidto03=bc_pidg_ym ! jl_s
        call MPI_ISEND(a(il,1,kl), 1,        bc_yslab_r4(ngy), pidto03,
     &                                basetag+3, bc_comm, isend03, ierr)

        pidfrom04=bc_pidg_yp ! ju_r
        call MPI_IRECV(a(il,bc_jml+1,kl), 1,bc_yslab_r4(ngy), pidfrom04, 
     &                                basetag+3, bc_comm, irecv04, ierr)

        pidto04=bc_pidg_yp ! ju_s
        call MPI_ISEND(a(il,bc_jml-ngy+1,kl),1,bc_yslab_r4(ngy),pidto04,
     &                                basetag+4, bc_comm, isend04, ierr)

        pidfrom03=bc_pidg_ym ! jl_r      
        call MPI_IRECV(a(il,jl,kl), 1,     bc_yslab_r4(ngy), pidfrom03,
     &                                basetag+4, bc_comm, irecv03, ierr)

        call MPI_WAIT(isend03,istatus,ierr)
        call MPI_WAIT(isend04,istatus,ierr)
        call MPI_WAIT(irecv03,istatus,ierr)
        call MPI_WAIT(irecv04,istatus,ierr)
      endif !if(bc_npy.gt.1.and.ngy.gt.0)then
      ! Done for x and y direction - not completed yet

      ! k planes
      if(bc_npz.gt.1.and.ngz.gt.0)then
        pidto05=bc_pidg_zm ! kl_s
        call MPI_ISEND(a(il,jl,1), 1,        bc_zslab_r4(ngz), pidto05,
     &                                basetag+5, bc_comm, isend05, ierr)

        pidfrom06=bc_pidg_zp ! ku_r
        call MPI_IRECV(a(il,jl,bc_kml+1), 1,bc_zslab_r4(ngz), pidfrom06,
     &                                basetag+5, bc_comm, irecv06, ierr)
            
        pidto06=bc_pidg_zp ! ku_s
        call MPI_ISEND(a(il,jl,bc_kml-ngz+1),1,bc_zslab_r4(ngz),pidto06,
     &                                basetag+6, bc_comm, isend06, ierr)

        pidfrom05=bc_pidg_zm ! kl_r
        call MPI_IRECV(a(il,jl,kl), 1,      bc_zslab_r4(ngz), pidfrom05,
     &                                basetag+6, bc_comm, irecv05, ierr)

        call MPI_WAIT(isend05,istatus,ierr)
        call MPI_WAIT(isend06,istatus,ierr)
        call MPI_WAIT(irecv05,istatus,ierr)
        call MPI_WAIT(irecv06,istatus,ierr)
      endif !if(bc_npz.gt.1.and.ngz.gt.0)then
      ! Internal tranfoer done for all directions - complete

      ! Applying the actual boundary conditions
      select case(trim(bctype))
      case('PERIODIC')
        call apply_periodic_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC')
        call apply_adiabatic_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_XY')
        call apply_periodic_z_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_YZ')
        call apply_periodic_x_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_ZX')
        call apply_periodic_y_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_X')
        call apply_periodic_y_bc_r4(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_Y')
        call apply_periodic_x_bc_r4(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_r4(a, ngx, ngy, ngz)
      case('ADIABATIC_Z')
        call apply_periodic_x_bc_r4(a, ngx, ngy, ngz)
        call apply_periodic_y_bc_r4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_r4(a, ngx, ngy, ngz)
      case default 
        write(*,*)'boundary_conditions module: fatal error:'
        write(*,*)'  bctype = ',trim(bctype)
        write(*,*)'  unknown boundary condition type. Aborting.'
        call abort
      end select !select case(trim(bctype))

      return 
      end subroutine boundary_condition_r4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_x_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_x_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidx .eq. 0 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(1-lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(pidx .eq. bc_npx-1 .and. ngx .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(bc_iml+lg,j,k) = a(bc_iml+1-lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_x_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_y_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_y_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidy .eq. 0 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,1-lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      if(pidy .eq. bc_npy-1 .and. ngy .gt. 0) then
        do k=1-ngz,bc_kml+ngz
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngy
            a(i,bc_jml+lg,k) = a(i,bc_jml+1-lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_y_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_adiabatic_z_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_z_bc_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      call get_pidxyz(pidx, pidy, pidz, bc_pidg)

      if(pidz .eq. 0 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,1-lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      if(pidz .eq. bc_npz-1 .and. ngz .gt. 0) then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,bc_kml+lg) = a(i,j,bc_kml+1-lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_adiabatic_z_bc_i4

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_i4(a,ngx,ngy,ngz)
      implicit none
      integer(4),intent(in) :: ngx,ngy,ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif
      
      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      end subroutine apply_periodic_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_x_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_x_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npx.eq.1.and.ngx.gt.0)then
        do k=1-ngz,bc_kml+ngz
        do j=1-ngy,bc_jml+ngy
          do lg=1,ngx
            a(-ngx+lg,j,k) = a(bc_iml-ngx+lg,j,k)
            a(bc_iml+lg,j,k) = a(lg,j,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_x_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_y_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) :: 
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_y_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npy.eq.1.and.ngy.gt.0)then
        do i=1-ngx,bc_iml+ngx
        do k=1-ngz,bc_kml+ngz
          do lg=1,ngy
            a(i,-ngy+lg,k) = a(i,bc_jml-ngy+lg,k)
            a(i,bc_jml+lg,k) = a(i,lg,k)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_y_bc_i4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_z_bc_i4(a, ngx, ngy, ngz)
      implicit none
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(inout) ::
     &             a(1-ngx:bc_iml+ngx,1-ngy:bc_jml+ngy,1-ngz:bc_kml+ngz)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_z_r8:'
        write(*,*)'  boundary condition uninitialized. Aborting.'
        call abort()
      endif

      if(bc_npz.eq.1.and.ngz.gt.0)then
        do j=1-ngy,bc_jml+ngy
        do i=1-ngx,bc_iml+ngx
          do lg=1,ngz
            a(i,j,-ngz+lg) = a(i,j,bc_kml-ngz+lg)
            a(i,j,bc_kml+lg) = a(i,j,lg)
          enddo
        enddo
        enddo
      endif

      return 
      end subroutine apply_periodic_z_bc_i4

      ! ----------------------------------------------------------------

      subroutine boundary_condition_i4(bctype,a,ng)
      implicit none

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a( 1-ng*bc_dimx:bc_iml+ng*bc_dimx,
     &                               1-ng*bc_dimy:bc_jml+ng*bc_dimy,
     &                               1-ng*bc_dimz:bc_kml+ng*bc_dimz )

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      type(MPI_Request) :: isend01, isend02
      type(MPI_Request) :: isend03, isend04
      type(MPI_Request) :: isend05, isend06
      type(MPI_Request) :: irecv01, irecv02
      type(MPI_Request) :: irecv03, irecv04
      type(MPI_Request) :: irecv05, irecv06

      type(MPI_Status) :: istatus
      integer(4) :: ierr
      integer(4) :: basetag
      integer(4) :: ngx, ngy, ngz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_r8_modern:'
        write(*,*)'boundary condition uninitialized. Aborting.'
        call abort
      endif

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      ngx = bc_dimx*ng
      ngy = bc_dimy*ng
      ngz = bc_dimz*ng

      il = 1-ngx
      iu = bc_iml+ngx
      jl = 1-ngy
      ju = bc_jml+ngy
      kl = 1-ngz
      ku = bc_kml+ngz

      ! Internal transfers
      ! i planes
      if(bc_npx.gt.1.and.ngx.gt.0)then
        pidto01=bc_pidg_xm ! il_s
        call MPI_ISEND(a(1,jl,kl),        1, bc_xslab_i4(ngx), pidto01,
     &                                basetag+1, bc_comm, isend01, ierr)

        pidfrom02=bc_pidg_xp ! iu_r
        call MPI_IRECV(a(bc_iml+1,jl,kl),1, bc_xslab_i4(ngx), pidfrom02,
     &                                basetag+1, bc_comm, irecv02, ierr)

        pidto02=bc_pidg_xp ! iu_s
        call MPI_ISEND(a(bc_iml-ngx+1,jl,kl),1,bc_xslab_i4(ngx),pidto02,
     &                                basetag+2, bc_comm, isend02, ierr)      

        pidfrom01=bc_pidg_xm ! il_r
        call MPI_IRECV(a(il,jl,kl),      1, bc_xslab_i4(ngx), pidfrom01,
     &                                basetag+2, bc_comm, irecv01, ierr)
            
        call MPI_WAIT(isend01,istatus,ierr)
        call MPI_WAIT(isend02,istatus,ierr)
        call MPI_WAIT(irecv01,istatus,ierr)
        call MPI_WAIT(irecv02,istatus,ierr)
      endif !if(bc_npx.gt.1.and.ngx.gt.0)then
      ! Done only for x direction - not complete yet

      ! j planes
      if(bc_npy.gt.1.and.ngy.gt.0)then
        pidto03=bc_pidg_ym ! jl_s
        call MPI_ISEND(a(il,1,kl), 1,        bc_yslab_i4(ngy), pidto03,
     &                                basetag+3, bc_comm, isend03, ierr)

        pidfrom04=bc_pidg_yp ! ju_r
        call MPI_IRECV(a(il,bc_jml+1,kl), 1,bc_yslab_i4(ngy), pidfrom04, 
     &                                basetag+3, bc_comm, irecv04, ierr)

        pidto04=bc_pidg_yp ! ju_s
        call MPI_ISEND(a(il,bc_jml-ngy+1,kl),1,bc_yslab_i4(ngy),pidto04,
     &                                basetag+4, bc_comm, isend04, ierr)

        pidfrom03=bc_pidg_ym ! jl_r      
        call MPI_IRECV(a(il,jl,kl), 1,       bc_yslab_i4(ngy),pidfrom03,
     &                                basetag+4, bc_comm, irecv03, ierr)

        call MPI_WAIT(isend03,istatus,ierr)
        call MPI_WAIT(isend04,istatus,ierr)
        call MPI_WAIT(irecv03,istatus,ierr)
        call MPI_WAIT(irecv04,istatus,ierr)
      endif !if(bc_npy.gt.1.and.ngy.gt.0)then
      ! Done for x and y direction - not completed yet

      ! k planes
      if(bc_npz.gt.1.and.ngz.gt.0)then
        pidto05=bc_pidg_zm ! kl_s
        call MPI_ISEND(a(il,jl,1), 1,        bc_zslab_i4(ngz), pidto05,
     &                                basetag+5, bc_comm, isend05, ierr)

        pidfrom06=bc_pidg_zp ! ku_r
        call MPI_IRECV(a(il,jl,bc_kml+1), 1, bc_zslab_i4(ngz),pidfrom06,
     &                                basetag+5, bc_comm, irecv06, ierr)
            
        pidto06=bc_pidg_zp ! ku_s
        call MPI_ISEND(a(il,jl,bc_kml-ngz+1),1,bc_zslab_i4(ngz),pidto06,
     &                                basetag+6, bc_comm, isend06, ierr)

        pidfrom05=bc_pidg_zm ! kl_r
        call MPI_IRECV(a(il,jl,kl), 1,      bc_zslab_i4(ngz), pidfrom05,
     &                                basetag+6, bc_comm, irecv05, ierr)

        call MPI_WAIT(isend05,istatus,ierr)
        call MPI_WAIT(isend06,istatus,ierr)
        call MPI_WAIT(irecv05,istatus,ierr)
        call MPI_WAIT(irecv06,istatus,ierr)
      endif !if(bc_npz.gt.1.and.ngz.gt.0)then
      ! Internal tranfoer done for all directions - complete

      ! Applying the actual boundary conditions
      select case(trim(bctype))
      case('PERIODIC')
        call apply_periodic_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC')
        call apply_adiabatic_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_XY')
        call apply_periodic_z_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_YZ')
        call apply_periodic_x_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_ZX')
        call apply_periodic_y_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_X')
        call apply_periodic_y_bc_i4(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_x_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_Y')
        call apply_periodic_x_bc_i4(a, ngx, ngy, ngz)
        call apply_periodic_z_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_y_bc_i4(a, ngx, ngy, ngz)
      case('ADIABATIC_Z')
        call apply_periodic_x_bc_i4(a, ngx, ngy, ngz)
        call apply_periodic_y_bc_i4(a, ngx, ngy, ngz)
        call apply_adiabatic_z_bc_i4(a, ngx, ngy, ngz)
      case default 
        write(*,*)'boundary_conditions module: fatal error:'
        write(*,*)'  bctype = ',trim(bctype)
        write(*,*)'  unknown boundary condition type. Aborting.'
        call abort
      end select !select case(trim(bctype))

      return 
      end subroutine boundary_condition_i4

      ! ----------------------------------------------------------------

      end module boundary_conditions