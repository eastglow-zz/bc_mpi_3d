      module boundary_conditions
      use mpi
      implicit none

      integer(4) :: bc_call_count
      integer(4) :: bc_max_mpi_tag = 0
      integer(4),parameter :: bc_mpitag_jumpstep = 13

      integer(4) :: bc_npx=1, bc_npy=1, bc_npz=1 !Number of partitioning along x,y, and z axis
      integer(4) :: iml,jml,kml
      integer(4) :: bc_imori, bc_jmori, bc_kmori
      integer(4) :: bc_pid

      logical :: bc_initialized = .false.
      

      contains
      subroutine init_bc(inpx,inpy,inpz,i_all,j_all,k_all,irank)
      implicit none
      integer(4),intent(in) :: inpx,inpy,inpz !saved in bc_npx, bc_npy, bc_npz
      integer(4),intent(in) :: i_all, j_all, k_all !saved in bc_imori,bc_jmori,bc_kmori
      integer(4),intent(in) :: irank !saved in bc_pid

      integer(4) :: ilb,iub,jlb,jub,klb,kub

      bc_npx = inpx
      bc_npy = inpy
      bc_npz = inpz

      bc_imori = i_all
      bc_jmori = j_all
      bc_kmori = k_all
      
      bc_pid = irank

      call get_ijkrange(ilb,iub,jlb,jub,klb,kub,bc_pid)
      iml=iub-ilb+1
      jml=jub-jlb+1
      kml=kub-klb+1

      bc_initialized = .true.

      call init_bc_call_count(0)

      end subroutine init_bc


      subroutine init_bc_call_count(val)
      implicit none
      integer(4),intent(in) :: val

       bc_call_count = val

      end subroutine init_bc_call_count


      subroutine verbose_bc_call_count()
      implicit none
      
       write(*,*) 'bc_call_count:', bc_call_count

      end subroutine verbose_bc_call_count


      subroutine verbose_bc_max_mpi_tag()
      implicit none

      write(*,*) 'bc_max_mpi_tag:', bc_max_mpi_tag

      end subroutine verbose_bc_max_mpi_tag


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
      if(i.gt.iml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(iU)"
      endif
      if(j.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jL)"
      endif
      if(j.gt.jml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jU)"
      endif
      if(k.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kL)"
      endif
      if(k.gt.kml+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kU)"
      endif

      if(go_nogo.eqv. .false.) then
            write(*,*)"Tried to access the out of array bound:"
            write(*,*)"flagid=",trim(flagid),", i,j,k = ",i,j,k
            write(*,*)"err_sign: ",trim(err_sign)," ,rank=",bc_pid
            call abort
      endif

      end subroutine bound_check


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


      subroutine get_pidxyz(pidx, pidy, pidz, pid)
      implicit none
      integer(4),intent(inout) :: pidx, pidy, pidz
      integer(4),intent(in) :: pid

      pidz = int(pid/bc_npx/bc_npy)
      pidy = int((pid - pidz*bc_npx*bc_npy)/bc_npx)
      pidx = pid - pidz*bc_npx*bc_npy - pidy*bc_npx

      end subroutine get_pidxyz


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

      end subroutine get_ijkrange


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

      end subroutine get_global_ijk


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

      END SUBROUTINE para_range


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

      pidx = int((igf-1)/iml)
      pidy = int((jgf-1)/jml)
      pidz = int((kgf-1)/kml)

      get_pid_from_globalijk = get_pid(pidx, pidy, pidz)
      return
      end function get_pid_from_globalijk


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

      il=igf-pidx*iml
      jl=jgf-pidy*jml
      kl=kgf-pidz*kml

      end subroutine get_localijk_from_globalijk


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

      end subroutine clamp_globalijk

      !-----------------------------------------------------------------

      subroutine boundary_condition_2d_i2(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv05, irecv06

      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 1
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER2, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER2, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_2d_i2(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_2d_i2(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_2d_i2(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_2d_i2(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_i2

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_i2(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      integer(4) :: i,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(iml-ng+lg,1,k)
                        a(iml+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,kml-ng+lg)
                        a(i,1,kml+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_i2

      !-----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_i2(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, lh
      logical :: gconx, gconz
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gconz=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,1,kml)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,1,l1)=a(iml+1-lg,1,l1)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,1,kml+lg)=a(l1,1,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
                  a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
            enddo
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,1,kml+1+lg,ng)
                  a(1-ng+l1,1,kml+1+lg)=a(1+lg,1,kml+1-ng+l1)
            enddo
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_i2


      !-----------------------------------------------------------------

      subroutine boundary_condition_2d_i4(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv05, irecv06

      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 1
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_2d_i4(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_2d_i4(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_2d_i4(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_2d_i4(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_i4

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_i4(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      integer(4) :: i,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(iml-ng+lg,1,k)
                        a(iml+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,kml-ng+lg)
                        a(i,1,kml+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_i4

      !-----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_i4(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, lh
      logical :: gconx, gconz
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gconz=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,1,kml)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,1,l1)=a(iml+1-lg,1,l1)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,1,kml+lg)=a(l1,1,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
                  a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
            enddo
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,1,kml+1+lg,ng)
                  a(1-ng+l1,1,kml+1+lg)=a(1+lg,1,kml+1-ng+l1)
            enddo
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_i4

      !-----------------------------------------------------------------

      subroutine boundary_condition_i2(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend03, isend04
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv03, irecv04
      integer(4) :: irecv05, irecv06


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 2*ng + jml
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1-ng
      ju = jml+ng
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER2, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER2, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER2, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_i2(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_i2(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_i2


      subroutine apply_periodic_bc_i2(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(iml-ng+lg,j,k)
                        a(iml+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,jml-ng+lg,k)
                        a(i,jml+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,kml-ng+lg)
                        a(i,j,kml+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_i2

      
      subroutine apply_nonperiodic_bc_i2(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(2),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, l2, lh
      logical :: gconx, gcony, gconz
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gcony=.false.; gconz=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gcony=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gcony=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gcony=.false.; gconz=.false.
      case('ADIABATIC_Y')
            gconx=.false.; gcony=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gcony=.false.; gconz=.true.
      case('ADIABATIC_XY')
            gconx=.true.; gcony=.true.; gconz=.false.
      case('ADIABATIC_YZ')
            gconx=.false.; gcony=.true.; gconz=.true.
      case('ADIABATIC_XZ')
            gconx=.true.; gcony=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconjl=lconjl.and.gcony
      lconju=lconju.and.gcony
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,jml,kml)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.jml+ng.and.l2.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,l1,l2)=a(iml+1-lg,l1,l2)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
      if(l1.le.kml+ng.and.l2.le.iml+ng)then  
            do lg=1,ng
                  ! jl
                  if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
                  ! ju
                  if(lconju.eqv..true.)a(l2,jml+lg,l1)=a(l2,jml+1-lg,l1)
            enddo
      endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng.and.l2.le.jml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,l2,kml+lg)=a(l1,l2,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                        a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
                  enddo
                  endif
            endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jukl",l2,jml+1+lg,1-ng+l1,ng)
                        a(l2,jml+1+lg,1-ng+l1)=a(l2,jml+1-ng+l1,1+lg)
                  enddo
                  endif
            endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlku",l2,1-ng+lg,kml+1+l1,ng)
                        a(l2,1-ng+lg,kml+1+l1)=a(l2,1+l1,kml+1-ng+lg)
                  enddo
                  endif
            endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("juku",l2,jml+1+lg,kml+1+l1,ng)
                        a(l2,jml+1+lg,kml+1+l1)=a(l2,jml-l1,kml-lg)
                  enddo
                  endif
            endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                        a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
                  enddo
                  endif
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuil",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kliu",1-ng+l1,l2,kml+1+lg,ng)
                        a(1-ng+l1,l2,kml+1+lg)=a(1+lg,l2,kml+1-ng+l1)
                  enddo
                  endif
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuiu",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                        a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iujl",iml+1+lg,1-ng+l1,l2,ng)
                        a(iml+1+lg,1-ng+l1,l2)=a(iml+1-ng+l1,1+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("ilju",1-ng+lg,jml+1+l1,l2,ng)
                        a(1-ng+lg,jml+1+l1,l2)=a(1+l1,jml+1-ng+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iuju",iml+1+lg,jml+1+l1,l2,ng)
                        a(iml+1+lg,jml+1+l1,l2)=a(iml-l1,jml-lg,l2)
                  enddo
                  endif
            endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,1-ng+lg)=a(iml-ng+1-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,1-ng+lg)=a(ng-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,1-ng+lg)=
     &                                  a(iml-ng+1-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,kml+1+lg)=a(ng-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,kml+1+lg)=
     &                                  a(iml-ng+1-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,kml+1+lg)=
     &                                  a(ng-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,kml+1+lg)=
     &                            a(iml-ng+1-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_i2

      !-----------------------------------------------------------------

      subroutine boundary_condition_i4(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06
c      integer(4) :: 
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06
c      integer(4) :: 


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 2*ng + jml
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1-ng
      ju = jml+ng
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_i4(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_i4(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_i4

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_i4(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(iml-ng+lg,j,k)
                        a(iml+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,jml-ng+lg,k)
                        a(i,jml+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,kml-ng+lg)
                        a(i,j,kml+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_i4

      !-----------------------------------------------------------------
      
      subroutine apply_nonperiodic_bc_i4(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, l2, lh
      logical :: gconx, gcony, gconz
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gcony=.false.; gconz=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gcony=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gcony=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gcony=.false.; gconz=.false.
      case('ADIABATIC_Y')
            gconx=.false.; gcony=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gcony=.false.; gconz=.true.
      case('ADIABATIC_XY')
            gconx=.true.; gcony=.true.; gconz=.false.
      case('ADIABATIC_YZ')
            gconx=.false.; gcony=.true.; gconz=.true.
      case('ADIABATIC_XZ')
            gconx=.true.; gcony=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconjl=lconjl.and.gcony
      lconju=lconju.and.gcony
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,jml,kml)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.jml+ng.and.l2.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,l1,l2)=a(iml+1-lg,l1,l2)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
      if(l1.le.kml+ng.and.l2.le.iml+ng)then  
            do lg=1,ng
                  ! jl
                  if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
                  ! ju
                  if(lconju.eqv..true.)a(l2,jml+lg,l1)=a(l2,jml+1-lg,l1)
            enddo
      endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng.and.l2.le.jml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,l2,kml+lg)=a(l1,l2,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                        a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
                  enddo
                  endif
            endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jukl",l2,jml+1+lg,1-ng+l1,ng)
                        a(l2,jml+1+lg,1-ng+l1)=a(l2,jml+1-ng+l1,1+lg)
                  enddo
                  endif
            endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlku",l2,1-ng+lg,kml+1+l1,ng)
                        a(l2,1-ng+lg,kml+1+l1)=a(l2,1+l1,kml+1-ng+lg)
                  enddo
                  endif
            endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("juku",l2,jml+1+lg,kml+1+l1,ng)
                        a(l2,jml+1+lg,kml+1+l1)=a(l2,jml-l1,kml-lg)
                  enddo
                  endif
            endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                        a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
                  enddo
                  endif
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuil",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kliu",1-ng+l1,l2,kml+1+lg,ng)
                        a(1-ng+l1,l2,kml+1+lg)=a(1+lg,l2,kml+1-ng+l1)
                  enddo
                  endif
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuiu",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                        a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iujl",iml+1+lg,1-ng+l1,l2,ng)
                        a(iml+1+lg,1-ng+l1,l2)=a(iml+1-ng+l1,1+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("ilju",1-ng+lg,jml+1+l1,l2,ng)
                        a(1-ng+lg,jml+1+l1,l2)=a(1+l1,jml+1-ng+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iuju",iml+1+lg,jml+1+l1,l2,ng)
                        a(iml+1+lg,jml+1+l1,l2)=a(iml-l1,jml-lg,l2)
                  enddo
                  endif
            endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,1-ng+lg)=a(iml-ng+1-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,1-ng+lg)=a(ng-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,1-ng+lg)=
     &                                  a(iml-ng+1-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,kml+1+lg)=a(ng-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,kml+1+lg)=
     &                                  a(iml-ng+1-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,kml+1+lg)=
     &                                  a(ng-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,kml+1+lg)=
     &                            a(iml-ng+1-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_i4

      !-----------------------------------------------------------------


      subroutine boundary_condition_2d_r4(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv05, irecv06

      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      arrsize(1) = 2*ng + iml
      arrsize(2) = 1
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_REAL, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_REAL, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_2d_r4(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_2d_r4(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_2d_r4(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_2d_r4(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_r4

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_r4(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      integer(4) :: i,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(iml-ng+lg,1,k)
                        a(iml+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,kml-ng+lg)
                        a(i,1,kml+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_r4

      !-----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_r4(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, lh
      logical :: gconx, gconz
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gconz=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,1,kml)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,1,l1)=a(iml+1-lg,1,l1)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,1,kml+lg)=a(l1,1,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
                  a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
            enddo
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,1,kml+1+lg,ng)
                  a(1-ng+l1,1,kml+1+lg)=a(1+lg,1,kml+1-ng+l1)
            enddo
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_r4

      !-----------------------------------------------------------------

      subroutine boundary_condition_2d_r8(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv05, irecv06


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      arrsize(1) = 2*ng + iml
      arrsize(2) = 1
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_2d_r8(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_2d_r8(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_2d_r8(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_2d_r8(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_r8

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_r8(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      integer(4) :: i,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(iml-ng+lg,1,k)
                        a(iml+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,kml-ng+lg)
                        a(i,1,kml+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_r8

      !-----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_r8(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, lh
      logical :: gconx, gconz
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gconz=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,1,kml)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,1,l1)=a(iml+1-lg,1,l1)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,1,kml+lg)=a(l1,1,kml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
                  a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
            enddo
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,1,kml+1+lg,ng)
                  a(1-ng+l1,1,kml+1+lg)=a(1+lg,1,kml+1-ng+l1)
            enddo
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",iml+1+l1,1,kml+1+lg,ng)
                  a(iml+1+l1,1,kml+1+lg)=a(iml-lg,1,kml-l1)
            enddo
            endif
      endif
      
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_r8
      
      !-----------------------------------------------------------------

      subroutine boundary_condition_r4(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend03, isend04
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv03, irecv04
      integer(4) :: irecv05, irecv06


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 2*ng + jml
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1-ng
      ju = jml+ng
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_REAL, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_REAL, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_REAL, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_r4(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_r4(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_r4


      subroutine apply_periodic_bc_r4(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(iml-ng+lg,j,k)
                        a(iml+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,jml-ng+lg,k)
                        a(i,jml+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,kml-ng+lg)
                        a(i,j,kml+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_r4

      
      subroutine apply_nonperiodic_bc_r4(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, l2, lh
      logical :: gconx, gcony, gconz
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gcony=.false.; gconz=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gcony=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gcony=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gcony=.false.; gconz=.false.
      case('ADIABATIC_Y')
            gconx=.false.; gcony=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gcony=.false.; gconz=.true.
      case('ADIABATIC_XY')
            gconx=.true.; gcony=.true.; gconz=.false.
      case('ADIABATIC_YZ')
            gconx=.false.; gcony=.true.; gconz=.true.
      case('ADIABATIC_XZ')
            gconx=.true.; gcony=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconjl=lconjl.and.gcony
      lconju=lconju.and.gcony
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,jml,kml)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.jml+ng.and.l2.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,l1,l2)=a(iml+1-lg,l1,l2)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
      if(l1.le.kml+ng.and.l2.le.iml+ng)then  
            do lg=1,ng
                  ! jl
                  if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
                  ! ju
                  if(lconju.eqv..true.)a(l2,jml+lg,l1)=a(l2,jml+1-lg,l1)
            enddo
      endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng.and.l2.le.jml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,l2,jml+lg)=a(l1,l2,jml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                        a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
                  enddo
                  endif
            endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jukl",l2,jml+1+lg,1-ng+l1,ng)
                        a(l2,jml+1+lg,1-ng+l1)=a(l2,jml+1-ng+l1,1+lg)
                  enddo
                  endif
            endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlku",l2,1-ng+lg,kml+1+l1,ng)
                        a(l2,1-ng+lg,kml+1+l1)=a(l2,1+l1,kml+1-ng+lg)
                  enddo
                  endif
            endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("juku",l2,jml+1+lg,kml+1+l1,ng)
                        a(l2,jml+1+lg,kml+1+l1)=a(l2,jml-l1,kml-lg)
                  enddo
                  endif
            endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                        a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
                  enddo
                  endif
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuil",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kliu",1-ng+l1,l2,kml+1+lg,ng)
                        a(1-ng+l1,l2,kml+1+lg)=a(1+lg,l2,kml+1-ng+l1)
                  enddo
                  endif
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuiu",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                        a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iujl",iml+1+lg,1-ng+l1,l2,ng)
                        a(iml+1+lg,1-ng+l1,l2)=a(iml+1-ng+l1,1+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("ilju",1-ng+lg,jml+1+l1,l2,ng)
                        a(1-ng+lg,jml+1+l1,l2)=a(1+l1,jml+1-ng+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iuju",iml+1+lg,jml+1+l1,l2,ng)
                        a(iml+1+lg,jml+1+l1,l2)=a(iml-l1,jml-lg,l2)
                  enddo
                  endif
            endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,1-ng+lg)=a(iml-ng+1-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,1-ng+lg)=a(ng-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,1-ng+lg)=
     &                                  a(iml-ng+1-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,kml+1+lg)=a(ng-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,kml+1+lg)=
     &                                  a(iml-ng+1-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,kml+1+lg)=
     &                                  a(ng-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,kml+1+lg)=
     &                            a(iml-ng+1-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_r4

      !-----------------------------------------------------------------

      subroutine boundary_condition_r8(bctype,a,ng)
      implicit none
      !include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02
      integer(4) :: pidto03, pidto04
      integer(4) :: pidto05, pidto06
      integer(4) :: pidfrom01, pidfrom02
      integer(4) :: pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06

      integer(4) :: isend01, isend02
      integer(4) :: isend03, isend04
      integer(4) :: isend05, isend06
      integer(4) :: irecv01, irecv02
      integer(4) :: irecv03, irecv04
      integer(4) :: irecv05, irecv06

      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + iml
      arrsize(2) = 2*ng + jml
      arrsize(3) = 2*ng + kml

      il = 1-ng
      iu = iml+ng
      jl = 1-ng
      ju = jml+ng
      kl = 1-ng
      ku = kml+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bctype))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_r8(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(iml+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(iml-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,jml+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,jml-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,kml+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,kml-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_r8(bctype, a, ng)

      end select !select case(trim(bctype))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_r8

      !-----------------------------------------------------------------

      subroutine apply_periodic_bc_r8(a,ng)
      implicit none
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(iml-ng+lg,j,k)
                        a(iml+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,jml-ng+lg,k)
                        a(i,jml+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,kml-ng+lg)
                        a(i,j,kml+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_r8

      !-----------------------------------------------------------------
      
      subroutine apply_nonperiodic_bc_r8(bctype, a, ng)
      implicit none
      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)

      !Local variables
      integer(4) :: lg
      integer(4) :: l1, l2, lh
      logical :: gconx, gcony, gconz
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'boundary_conditions module: fatal error:'
            write(*,*)'boundary condition uninitialized. Aborting.'
            call abort
      endif


      gconx=.false.; gcony=.false.; gconz=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      select case(trim(bctype))
      case('PERIODIC')
            gconx=.false.; gcony=.false.; gconz=.false.
      case('ADIABATIC')
            gconx=.true.; gcony=.true.; gconz=.true.
      case('ADIABATIC_X')
            gconx=.true.; gcony=.false.; gconz=.false.
      case('ADIABATIC_Y')
            gconx=.false.; gcony=.true.; gconz=.false.
      case('ADIABATIC_Z')
            gconx=.false.; gcony=.false.; gconz=.true.
      case('ADIABATIC_XY')
            gconx=.true.; gcony=.true.; gconz=.false.
      case('ADIABATIC_YZ')
            gconx=.false.; gcony=.true.; gconz=.true.
      case('ADIABATIC_XZ')
            gconx=.true.; gcony=.false.; gconz=.true.
      end select

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconx
      lconiu=lconiu.and.gconx
      lconjl=lconjl.and.gcony
      lconju=lconju.and.gcony
      lconkl=lconkl.and.gconz
      lconku=lconku.and.gconz

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(iml,jml,kml)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.jml+ng.and.l2.le.kml+ng)then  
            do lg=1,ng
                  ! il
                  if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
                  ! iu
                  if(lconiu.eqv..true.)a(iml+lg,l1,l2)=a(iml+1-lg,l1,l2)
            enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
      if(l1.le.kml+ng.and.l2.le.iml+ng)then  
            do lg=1,ng
                  ! jl
                  if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
                  ! ju
                  if(lconju.eqv..true.)a(l2,jml+lg,l1)=a(l2,jml+1-lg,l1)
            enddo
      endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.iml+ng.and.l2.le.jml+ng)then  
            do lg=1,ng
                  ! kl
                  if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
                  ! ku
                  if(lconku.eqv..true.)a(l1,l2,jml+lg)=a(l1,l2,jml+1-lg)
            enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                        a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
                  enddo
                  endif
            endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jukl",l2,jml+1+lg,1-ng+l1,ng)
                        a(l2,jml+1+lg,1-ng+l1)=a(l2,jml+1-ng+l1,1+lg)
                  enddo
                  endif
            endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("jlku",l2,1-ng+lg,kml+1+l1,ng)
                        a(l2,1-ng+lg,kml+1+l1)=a(l2,1+l1,kml+1-ng+lg)
                  enddo
                  endif
            endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
            if(l2.ge.1.and.l2.le.iml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("juku",l2,jml+1+lg,kml+1+l1,ng)
                        a(l2,jml+1+lg,kml+1+l1)=a(l2,jml-l1,kml-lg)
                  enddo
                  endif
            endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                        a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
                  enddo
                  endif
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuil",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kliu",1-ng+l1,l2,kml+1+lg,ng)
                        a(1-ng+l1,l2,kml+1+lg)=a(1+lg,l2,kml+1-ng+l1)
                  enddo
                  endif
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l2.ge.1.and.l2.le.jml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("kuiu",iml+1+l1,l2,kml+1+lg,ng)
                        a(iml+1+l1,l2,kml+1+lg)=a(iml-lg,l2,kml-l1)
                  enddo
                  endif
            endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                        a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iujl",iml+1+lg,1-ng+l1,l2,ng)
                        a(iml+1+lg,1-ng+l1,l2)=a(iml+1-ng+l1,1+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("ilju",1-ng+lg,jml+1+l1,l2,ng)
                        a(1-ng+lg,jml+1+l1,l2)=a(1+l1,jml+1-ng+lg,l2)
                  enddo
                  endif
            endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
            if(l2.ge.1.and.l2.le.kml) then
                  if(l1.ge.0.and.l1.le.ng-1) then
                  do lg=0,ng-1
                        call bound_check("iuju",iml+1+lg,jml+1+l1,l2,ng)
                        a(iml+1+lg,jml+1+l1,l2)=a(iml-l1,jml-lg,l2)
                  enddo
                  endif
            endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,1-ng+lg)=a(iml-ng+1-lg,ng-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,1-ng+lg)=a(ng-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,1-ng+lg)=
     &                                  a(iml-ng+1-lg,jml-ng+1-l2,ng-l1)
            enddo
            endif
            endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,1-ng+l2,kml+1+lg)=a(ng-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,1-ng+l2,kml+1+lg)=
     &                                  a(iml-ng+1-lg,ng-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(1-ng+l1,jml+1+l2,kml+1+lg)=
     &                                  a(ng-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
            if(l2.ge.0.and.l2.le.ng-1) then
            if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
            a(iml+1+l1,jml+1+l2,kml+1+lg)=
     &                            a(iml-ng+1-lg,jml-ng+1-l2,kml-ng+1-l1)
            enddo
            endif
            endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_r8

      !-----------------------------------------------------------------

      end module boundary_conditions