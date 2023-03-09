      module boundary_conditions
      use input
      implicit none

      integer(4) :: bc_call_count
      integer(4) :: bc_max_mpi_tag = 0
      integer(4),parameter :: bc_mpitag_jumpstep = 13
      ! bc_call_count increases at every single call of either boundary_condition_i4_par() or boundary_condition_r8_par().
      ! Purpose of this variable is not to overlap the ranges of MPI_TAG used during MPI_ISEND/IRECV from a calls of bc to another bc calls to prevent mix-up of the data buffer.
      ! The allowed maximum value of MPI_TAG is (most conservatively) known as 32767, but possibly higher library by library or version by version.
      ! Highly recommend reinitialize bc_call_count at every time step via the subroutine init_bc_call_count(val).
      ! how bc_call_count work:
      !  1. bc_call_count increases by 1 at every call of bc.
      !  2. The base MPI_TAG value is computed by bc_call_count*bc_mpitag_jumpstep
      !     where, bc_mpitag_jumpstep is a parameter limiting the maximum independent mpi_isend/irecv in a single bc call.
      !     For now, minimum value of bc_mpitag_jumpstep is 13 (now, MPI_TAG of base+1~12 is used in a single bc call).
      !  3. Lower bc_mpitag_jumpstep, higher bc_call_count possible.

      ! Layer indexing and naming: ij, jk, ki, 6 planes: il, iu, jl, ju, kl, ku
      !1: il plane: jk indexing near i=1-ng, stack of 2D planes with the thickness of the # of ghost layers
      !2: iu plane: jk indexing near i=im+ng, stack of 2D planes with the thickness of the # of ghost layers
      !3: jl plane: ki indexing near j=1-ng, stack of 2D planes with the thickness of the # of ghost layers
      !4: ju plane: ki indexing near j=jm+ng, stack of 2D planes with the thickness of the # of ghost layers
      !5: kl plane: ij indexing near k=1-ng, stack of 2D planes with the thickness of the # of ghost layers
      !6: ku plane: ij indexing near k=km+ng, stack of 2D planes with the thickness of the # of ghost layers

      ! Edge naming
      !iljl edge: intersecting il and jl plane
      !ilju edge: intersecting il and ju plane
      !iujl edge: intersecting iu and jl plane
      !iuju edge: intersecting iu and ju plane
      !jlkl edge: intersecting jl and kl plane
      !jlku edge: intersecting jl and ku plane
      !jukl edge: intersecting ju and kl plane
      !juku edge: intersecting ju and ku plane
      !klil edge: intersecting kl and il plane
      !kliu edge: intersecting kl and iu plane
      !kuil edge: intersecting ku and il plane
      !kuiu edge: intersecting ku and iu plane

      ! Corner naming
      !lll: (i,j,k)= nearby( 1, 1, 1)
      !ull: (i,j,k)= nearby(im, 1, 1)
      !lul: (i,j,k)= nearby( 1,jm, 1)
      !uul: (i,j,k)= nearby(im,jm, 1)
      !llu: (i,j,k)= nearby( 1, 1,km)
      !ulu: (i,j,k)= nearby(im, 1,km)
      !luu: (i,j,k)= nearby( 1,jm,km)
      !uuu: (i,j,k)= nearby(im,jm,km)

      ! Index permutation is being kept during the copy of the data between differently structured arrays, i.e. 3D array to 2D array
      ! a(i,j,*) --> aa(i,j,$)
      ! a(i,*,k) --> aa(k,i,$)
      ! a(*,j,k) --> aa(j,k,$)

      ! Process ID (=myrank) management
      ! Raw process ID from MPI library: PID = 0~NP-1
      ! Auxiliary process ID coordinate in partitioning coordinate: PIDX=0~NPX, PIDY=0~NPY, PIDZ=0~NPZ

      ! PIDX, PIDY, PIDZ v.s. PID Relation  --> can be obrained through 'get_pidxyz(pidx,pidy,pidz,pid)' subroutine
      ! pidz = int(pid/NPX/NPY) 
      ! pidy = int((pid - pidz*NPX*NPY)/NPX) 
      ! pidx = pid - pidz*NPX*NPY - pidy*NPX 

      !Converting arbitrary (PIDX, PIDY, PIDZ) into PID_temporary
      ! PID = pidxw + pidyw*NPX + pidzw*NPX*NPY --> can be converted through 'get_pid(pidx,pidy,pidz,bcmode)' function

      ! BC mode names
      ! periodic: periodic BC along all directions
      ! adiabatic: adiabatic BC for all boundaries
      ! adiabatic_x: adiabatic BC along x direction (for il & iu), but periodic for the rest boundaries
      ! adiabatic_y: adiabatic BC along y direction (for jl & ju), but periodic for the rest boundaries
      ! adiabatic_z: adiabatic BC along z direction (for kl & ku), but periodic for the rest boundaries
      ! adiabatic_xy: adiabatic BC along x and y directions (il,iu,jl,ju), but periodic for z direction (kl,ku)
      ! adiabatic_yz: adiabatic BC along y and z directions (jl,ju,kl,ku), but periodic for x direction (il,iu)
      ! adiabatic_zx: adiabatic BC along z and x directions (kl,ku,il,iu), but periodic for y direction (jl,ju)
      ! Further BC modes will be coming up...


      contains

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
      include 'mpif.h'
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


      integer(4) function get_pid(pidx, pidy, pidz)
      use input
      implicit none
      integer(4),intent(in) :: pidx, pidy, pidz

      integer(4) :: pidxw, pidyw, pidzw

      pidxw = pidx
      if(pidxw.gt.NPX-1)pidxw=0
      if(pidxw.lt.0)pidxw=NPX-1

      pidyw = pidy
      if(pidyw.gt.NPY-1)pidyw=0
      if(pidyw.lt.0)pidyw=NPY-1

      pidzw = pidz
      if(pidzw.gt.NPZ-1)pidzw=0
      if(pidzw.lt.0)pidzw=NPZ-1

      get_pid = pidxw + pidyw*NPX + pidzw*NPX*NPY
      return
      end function get_pid


      subroutine get_pidxyz(pidx, pidy, pidz, pid)
      use input
      implicit none
      integer(4),intent(inout) :: pidx, pidy, pidz
      integer(4),intent(in) :: pid

      pidz = int(pid/NPX/NPY)
      pidy = int((pid - pidz*NPX*NPY)/NPX)
      pidx = pid - pidz*NPX*NPY - pidy*NPX

      end subroutine get_pidxyz


      subroutine get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      use input
      implicit none
      integer(4),intent(inout) :: il,iu,jl,ju,kl,ku
      integer(4),intent(in) :: pid

      integer(4) :: pidx, pidy, pidz

      call get_pidxyz(pidx, pidy, pidz, pid)

      ! x bounds
      call para_range(1,im,NPX,pidx,il,iu)
      ! y bounds
      call para_range(1,jm,NPY,pidy,jl,ju)
      ! z bounds
      call para_range(1,km,NPZ,pidz,kl,ku)

      end subroutine get_ijkrange


      subroutine get_global_ijk(ig,jg,kg, ilocal,jlocal,klocal, pid)
      use input
      implicit none
      integer(4),intent(inout) :: ig, jg, kg !global ijk index
      integer(4),intent(in) :: ilocal, jlocal, klocal !local ijk index
      integer(4),intent(in) :: pid

      integer(4) :: il,iu, jl,ju, kl,ku
      integer(4) :: pidx, pidy, pidz

      call get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      ig = il + (ilocal-1)  !for spatial index starts from 1, not 0
      jg = jl + (jlocal-1)  !for spatial index starts from 1, not 0
      kg = kl + (klocal-1)  !for spatial index starts from 1, not 0

      end subroutine get_global_ijk


      SUBROUTINE para_range(n1,n2,npart,irank,ista,iend)
      ! Deviding system size KMORI into KM = KMORI/(porcess number)
      ! nprocs: total number of porcesses (nprocs=48 in my computer)
      ! myrank: current process name (myrank=1~47)
      ! ista and iend: first and last grids (J-direction) on each process in the total grid system
      ! use: CALL para_range(1,KMORI,nprocs,myrank,ista,iend)
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

      igo=ig
      jgo=jg
      kgo=kg

      select case(trim(bctype))
      case('periodic')
            igo=mod(ig+min(ig-1,0)*(-IM)-1,IM)+1
            jgo=mod(jg+min(jg-1,0)*(-JM)-1,JM)+1
            kgo=mod(kg+min(kg-1,0)*(-KM)-1,KM)+1
      case('adiabatic')
            igo=min(max(ig,1),IM)
            jgo=min(max(jg,1),JM)
            kgo=min(max(kg,1),KM)
      case('adiabatic_x')
            igo=min(max(ig,1),IM)
            jgo=mod(jg+min(jg-1,0)*(-JM)-1,JM)+1
            kgo=mod(kg+min(kg-1,0)*(-KM)-1,KM)+1
      case('adiabatic_y')
            igo=mod(ig+min(ig-1,0)*(-IM)-1,IM)+1
            jgo=min(max(jg,1),JM)
            kgo=mod(kg+min(kg-1,0)*(-KM)-1,KM)+1
      case('adiabatic_z')
            igo=mod(ig+min(ig-1,0)*(-IM)-1,IM)+1
            jgo=mod(jg+min(jg-1,0)*(-JM)-1,JM)+1
            kgo=min(max(kg,1),KM)
      case('adiabatic_xy')
            igo=min(max(ig,1),IM)
            jgo=min(max(jg,1),JM)
            kgo=mod(kg+min(kg-1,0)*(-KM)-1,KM)+1
      case('adiabatic_yz')
            igo=mod(ig+min(ig-1,0)*(-IM)-1,IM)+1
            jgo=min(max(jg,1),JM)
            kgo=min(max(kg,1),KM)
      case('adiabatic_zx')
            igo=min(max(ig,1),IM)
            jgo=mod(jg+min(jg-1,0)*(-JM)-1,JM)+1
            kgo=min(max(kg,1),KM)
      end select

      end subroutine clamp_globalijk

      subroutine boundary_condition_i4_par(bctype,a,ng)
      use input
      implicit none
      include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      !Local variables
      !Contiguous arrays for MPI communications, s: send, r: recv
      !planes for send
      integer(4) :: a_il_s(1-ng:jml+ng,1-ng:kml+ng,ng)
      integer(4) :: a_iu_s(1-ng:jml+ng,1-ng:kml+ng,ng)
      integer(4) :: a_jl_s(1-ng:kml+ng,1-ng:iml+ng,ng)
      integer(4) :: a_ju_s(1-ng:kml+ng,1-ng:iml+ng,ng)
      integer(4) :: a_kl_s(1-ng:iml+ng,1-ng:jml+ng,ng)
      integer(4) :: a_ku_s(1-ng:iml+ng,1-ng:jml+ng,ng)
      !edges for send
      integer(4) :: a_iljl_s(1-ng:kml+ng,ng,ng) !edge intersecting il and jl plane
      integer(4) :: a_ilju_s(1-ng:kml+ng,ng,ng) !edge intersecting il and ju plane
      integer(4) :: a_iujl_s(1-ng:kml+ng,ng,ng) !edge intersecting iu and jl plane
      integer(4) :: a_iuju_s(1-ng:kml+ng,ng,ng) !edge intersecting iu and ju plane
      integer(4) :: a_jlkl_s(1-ng:iml+ng,ng,ng) !edge intersecting jl and kl plane
      integer(4) :: a_jlku_s(1-ng:iml+ng,ng,ng) !edge intersecting jl and ku plane
      integer(4) :: a_jukl_s(1-ng:iml+ng,ng,ng) !edge intersecting ju and kl plane
      integer(4) :: a_juku_s(1-ng:iml+ng,ng,ng) !edge intersecting ju and ku plane
      integer(4) :: a_klil_s(1-ng:jml+ng,ng,ng) !edge intersecting kl and il plane
      integer(4) :: a_kliu_s(1-ng:jml+ng,ng,ng) !edge intersecting kl and iu plane
      integer(4) :: a_kuil_s(1-ng:jml+ng,ng,ng) !edge intersecting ku and il plane
      integer(4) :: a_kuiu_s(1-ng:jml+ng,ng,ng) !edge intersecting ku and iu plane
      !corners for send
      integer(4) :: a_lll_s(ng,ng,ng) !(i,j,k)= nearby( 1, 1, 1)
      integer(4) :: a_ull_s(ng,ng,ng) !(i,j,k)= nearby(im, 1, 1)
      integer(4) :: a_lul_s(ng,ng,ng) !(i,j,k)= nearby( 1,jm, 1)
      integer(4) :: a_uul_s(ng,ng,ng) !(i,j,k)= nearby(im,jm, 1)
      integer(4) :: a_llu_s(ng,ng,ng) !(i,j,k)= nearby( 1, 1,km)
      integer(4) :: a_ulu_s(ng,ng,ng) !(i,j,k)= nearby(im, 1,km)
      integer(4) :: a_luu_s(ng,ng,ng) !(i,j,k)= nearby( 1,jm,km)
      integer(4) :: a_uuu_s(ng,ng,ng) !(i,j,k)= nearby(im,jm,km)

      !planes for recv
      integer(4) :: a_il_r(1-ng:jml+ng,1-ng:kml+ng,ng)
      integer(4) :: a_iu_r(1-ng:jml+ng,1-ng:kml+ng,ng)
      integer(4) :: a_jl_r(1-ng:kml+ng,1-ng:iml+ng,ng)
      integer(4) :: a_ju_r(1-ng:kml+ng,1-ng:iml+ng,ng)
      integer(4) :: a_kl_r(1-ng:iml+ng,1-ng:jml+ng,ng)
      integer(4) :: a_ku_r(1-ng:iml+ng,1-ng:jml+ng,ng)
      !edges for recv
      integer(4) :: a_iljl_r(1-ng:kml+ng,ng,ng) !edge intersecting il and jl plane
      integer(4) :: a_ilju_r(1-ng:kml+ng,ng,ng) !edge intersecting il and ju plane
      integer(4) :: a_iujl_r(1-ng:kml+ng,ng,ng) !edge intersecting iu and jl plane
      integer(4) :: a_iuju_r(1-ng:kml+ng,ng,ng) !edge intersecting iu and ju plane
      integer(4) :: a_jlkl_r(1-ng:iml+ng,ng,ng) !edge intersecting jl and kl plane
      integer(4) :: a_jlku_r(1-ng:iml+ng,ng,ng) !edge intersecting jl and ku plane
      integer(4) :: a_jukl_r(1-ng:iml+ng,ng,ng) !edge intersecting ju and kl plane
      integer(4) :: a_juku_r(1-ng:iml+ng,ng,ng) !edge intersecting ju and ku plane
      integer(4) :: a_klil_r(1-ng:jml+ng,ng,ng) !edge intersecting kl and il plane
      integer(4) :: a_kliu_r(1-ng:jml+ng,ng,ng) !edge intersecting kl and iu plane
      integer(4) :: a_kuil_r(1-ng:jml+ng,ng,ng) !edge intersecting ku and il plane
      integer(4) :: a_kuiu_r(1-ng:jml+ng,ng,ng) !edge intersecting ku and iu plane
      !corners for recv
      integer(4) :: a_lll_r(ng,ng,ng) !(i,j,k)= nearby( 1, 1, 1)
      integer(4) :: a_ull_r(ng,ng,ng) !(i,j,k)= nearby(im, 1, 1)
      integer(4) :: a_lul_r(ng,ng,ng) !(i,j,k)= nearby( 1,jm, 1)
      integer(4) :: a_uul_r(ng,ng,ng) !(i,j,k)= nearby(im,jm, 1)
      integer(4) :: a_llu_r(ng,ng,ng) !(i,j,k)= nearby( 1, 1,km)
      integer(4) :: a_ulu_r(ng,ng,ng) !(i,j,k)= nearby(im, 1,km)
      integer(4) :: a_luu_r(ng,ng,ng) !(i,j,k)= nearby( 1,jm,km)
      integer(4) :: a_uuu_r(ng,ng,ng) !(i,j,k)= nearby(im,jm,km)

      integer(4) :: i,j,k,n1,n2

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      !Process coordinate
      call get_pidxyz(pidx, pidy, pidz, myrank)

      !Preparation for send buffers
      ! Planes
      call extract_plane_i4(a_il_s, 'il', a, 1,jml,1,kml,ng)
      call extract_plane_i4(a_iu_s, 'iu', a, 1,jml,1,kml,ng)
      call extract_plane_i4(a_jl_s, 'jl', a, 1,kml,1,iml,ng)
      call extract_plane_i4(a_ju_s, 'ju', a, 1,kml,1,iml,ng)
      call extract_plane_i4(a_kl_s, 'kl', a, 1,iml,1,jml,ng)
      call extract_plane_i4(a_ku_s, 'ku', a, 1,iml,1,jml,ng)

      !edges
      call extract_edge_i4(a_iljl_s, 'iljl', a, 1,kml,ng)
      call extract_edge_i4(a_ilju_s, 'ilju', a, 1,kml,ng)
      call extract_edge_i4(a_iujl_s, 'iujl', a, 1,kml,ng)
      call extract_edge_i4(a_iuju_s, 'iuju', a, 1,kml,ng)
      call extract_edge_i4(a_jlkl_s, 'jlkl', a, 1,iml,ng)
      call extract_edge_i4(a_jlku_s, 'jlku', a, 1,iml,ng)
      call extract_edge_i4(a_jukl_s, 'jukl', a, 1,iml,ng)
      call extract_edge_i4(a_juku_s, 'juku', a, 1,iml,ng)
      call extract_edge_i4(a_klil_s, 'klil', a, 1,jml,ng)
      call extract_edge_i4(a_kliu_s, 'kliu', a, 1,jml,ng)
      call extract_edge_i4(a_kuil_s, 'kuil', a, 1,jml,ng)
      call extract_edge_i4(a_kuiu_s, 'kuiu', a, 1,jml,ng)

      !corners
      call extract_corner_i4(a_lll_s, 'lll', a, ng)
      call extract_corner_i4(a_ull_s, 'ull', a, ng)
      call extract_corner_i4(a_lul_s, 'lul', a, ng)
      call extract_corner_i4(a_uul_s, 'uul', a, ng)
      call extract_corner_i4(a_llu_s, 'llu', a, ng)
      call extract_corner_i4(a_ulu_s, 'ulu', a, ng)
      call extract_corner_i4(a_luu_s, 'luu', a, ng)
      call extract_corner_i4(a_uuu_s, 'uuu', a, ng)

      ! Layer send/recv
      ! i planes
      pidto01=get_pid(pidx-1,pidy,pidz) ! a_il_s
      pidto02=get_pid(pidx+1,pidy,pidz) ! a_iu_s
      pidto03=get_pid(pidx,pidy-1,pidz) ! a_jl_s
      pidto04=get_pid(pidx,pidy+1,pidz) ! a_ju_s
      pidto05=get_pid(pidx,pidy,pidz-1) ! a_kl_s
      pidto06=get_pid(pidx,pidy,pidz+1) ! a_ku_s
      call MPI_ISEND(a_il_s(1-ng,1-ng,1),size(a_il_s),MPI_INTEGER,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_iu_s(1-ng,1-ng,1),size(a_iu_s),MPI_INTEGER,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_jl_s(1-ng,1-ng,1),size(a_jl_s),MPI_INTEGER,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_ju_s(1-ng,1-ng,1),size(a_ju_s),MPI_INTEGER,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_kl_s(1-ng,1-ng,1),size(a_kl_s),MPI_INTEGER,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_ku_s(1-ng,1-ng,1),size(a_ku_s),MPI_INTEGER,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

      pidfrom01=get_pid(pidx-1,pidy,pidz) ! a_il_r
      pidfrom02=get_pid(pidx+1,pidy,pidz) ! a_iu_r
      pidfrom03=get_pid(pidx,pidy-1,pidz) ! a_jl_r
      pidfrom04=get_pid(pidx,pidy+1,pidz) ! a_ju_r
      pidfrom05=get_pid(pidx,pidy,pidz-1) ! a_kl_r
      pidfrom06=get_pid(pidx,pidy,pidz+1) ! a_ku_r
      call MPI_IRECV(a_il_r(1-ng,1-ng,1),size(a_il_r),MPI_INTEGER,
     & pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_iu_r(1-ng,1-ng,1),size(a_iu_r),MPI_INTEGER,
     & pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_jl_r(1-ng,1-ng,1),size(a_jl_r),MPI_INTEGER,
     & pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_ju_r(1-ng,1-ng,1),size(a_ju_r),MPI_INTEGER,
     & pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_kl_r(1-ng,1-ng,1),size(a_kl_r),MPI_INTEGER,
     & pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_ku_r(1-ng,1-ng,1),size(a_ku_r),MPI_INTEGER,
     & pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)

      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)


      ! Edge send/recv
      pidto01=get_pid(pidx-1,pidy-1,pidz) ! a_iljl_s
      pidto02=get_pid(pidx-1,pidy+1,pidz) ! a_ilju_s
      pidto03=get_pid(pidx+1,pidy-1,pidz) ! a_iujl_s
      pidto04=get_pid(pidx+1,pidy+1,pidz) ! a_iuju_s
      pidto05=get_pid(pidx,pidy-1,pidz-1) ! a_jlkl_s
      pidto06=get_pid(pidx,pidy-1,pidz+1) ! a_jlku_s
      pidto07=get_pid(pidx,pidy+1,pidz-1) ! a_jukl_s
      pidto08=get_pid(pidx,pidy+1,pidz+1) ! a_juku_s
      pidto09=get_pid(pidx-1,pidy,pidz-1) ! a_klil_s
      pidto10=get_pid(pidx+1,pidy,pidz-1) ! a_kliu_s
      pidto11=get_pid(pidx-1,pidy,pidz+1) ! a_kuil_s
      pidto12=get_pid(pidx+1,pidy,pidz+1) ! a_kuiu_s
      call MPI_ISEND(a_iljl_s(1-ng,1,1),size(a_iljl_s),MPI_INTEGER,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_ilju_s(1-ng,1,1),size(a_ilju_s),MPI_INTEGER,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_iujl_s(1-ng,1,1),size(a_iujl_s),MPI_INTEGER,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_iuju_s(1-ng,1,1),size(a_iuju_s),MPI_INTEGER,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_jlkl_s(1-ng,1,1),size(a_jlkl_s),MPI_INTEGER,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_jlku_s(1-ng,1,1),size(a_jlku_s),MPI_INTEGER,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
      call MPI_ISEND(a_jukl_s(1-ng,1,1),size(a_jukl_s),MPI_INTEGER,
     & pidto07,basetag+7,MPI_COMM_WORLD,isend07,ierr)
      call MPI_ISEND(a_juku_s(1-ng,1,1),size(a_juku_s),MPI_INTEGER,
     & pidto08,basetag+8,MPI_COMM_WORLD,isend08,ierr)
      call MPI_ISEND(a_klil_s(1-ng,1,1),size(a_klil_s),MPI_INTEGER,
     & pidto09,basetag+9,MPI_COMM_WORLD,isend09,ierr)
      call MPI_ISEND(a_kliu_s(1-ng,1,1),size(a_kliu_s),MPI_INTEGER,
     & pidto10,basetag+10,MPI_COMM_WORLD,isend10,ierr)
      call MPI_ISEND(a_kuil_s(1-ng,1,1),size(a_kuil_s),MPI_INTEGER,
     & pidto11,basetag+11,MPI_COMM_WORLD,isend11,ierr)
      call MPI_ISEND(a_kuiu_s(1-ng,1,1),size(a_kuiu_s),MPI_INTEGER,
     & pidto12,basetag+12,MPI_COMM_WORLD,isend12,ierr)

      pidfrom01=get_pid(pidx-1,pidy-1,pidz) ! a_iljl_r
      pidfrom02=get_pid(pidx-1,pidy+1,pidz) ! a_ilju_r
      pidfrom03=get_pid(pidx+1,pidy-1,pidz) ! a_iujl_r
      pidfrom04=get_pid(pidx+1,pidy+1,pidz) ! a_iuju_r
      pidfrom05=get_pid(pidx,pidy-1,pidz-1) ! a_jlkl_r
      pidfrom06=get_pid(pidx,pidy-1,pidz+1) ! a_jlku_r
      pidfrom07=get_pid(pidx,pidy+1,pidz-1) ! a_jukl_r
      pidfrom08=get_pid(pidx,pidy+1,pidz+1) ! a_juku_r
      pidfrom09=get_pid(pidx-1,pidy,pidz-1) ! a_klil_r
      pidfrom10=get_pid(pidx+1,pidy,pidz-1) ! a_kliu_r
      pidfrom11=get_pid(pidx-1,pidy,pidz+1) ! a_kuil_r
      pidfrom12=get_pid(pidx+1,pidy,pidz+1) ! a_kuiu_r
      call MPI_IRECV(a_iljl_r(1-ng,1,1),size(a_iljl_r),MPI_INTEGER,
     & pidfrom01,basetag+4,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_ilju_r(1-ng,1,1),size(a_ilju_r),MPI_INTEGER,
     & pidfrom02,basetag+3,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_iujl_r(1-ng,1,1),size(a_iujl_r),MPI_INTEGER,
     & pidfrom03,basetag+2,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_iuju_r(1-ng,1,1),size(a_iuju_r),MPI_INTEGER,
     & pidfrom04,basetag+1,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_jlkl_r(1-ng,1,1),size(a_jlkl_r),MPI_INTEGER,
     & pidfrom05,basetag+8,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_jlku_r(1-ng,1,1),size(a_jlku_r),MPI_INTEGER,
     & pidfrom06,basetag+7,MPI_COMM_WORLD,irecv06,ierr)
      call MPI_IRECV(a_jukl_r(1-ng,1,1),size(a_jukl_r),MPI_INTEGER,
     & pidfrom07,basetag+6,MPI_COMM_WORLD,irecv07,ierr)
      call MPI_IRECV(a_juku_r(1-ng,1,1),size(a_juku_r),MPI_INTEGER,
     & pidfrom08,basetag+5,MPI_COMM_WORLD,irecv08,ierr)
      call MPI_IRECV(a_klil_r(1-ng,1,1),size(a_klil_r),MPI_INTEGER,
     & pidfrom09,basetag+12,MPI_COMM_WORLD,irecv09,ierr)
      call MPI_IRECV(a_kliu_r(1-ng,1,1),size(a_kliu_r),MPI_INTEGER,
     & pidfrom10,basetag+11,MPI_COMM_WORLD,irecv10,ierr)
      call MPI_IRECV(a_kuil_r(1-ng,1,1),size(a_kuil_r),MPI_INTEGER,
     & pidfrom11,basetag+10,MPI_COMM_WORLD,irecv11,ierr)
      call MPI_IRECV(a_kuiu_r(1-ng,1,1),size(a_kuiu_r),MPI_INTEGER,
     & pidfrom12,basetag+9,MPI_COMM_WORLD,irecv12,ierr)


      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(isend07,istatus,ierr)
      call MPI_WAIT(isend08,istatus,ierr)
      call MPI_WAIT(isend09,istatus,ierr)
      call MPI_WAIT(isend10,istatus,ierr)
      call MPI_WAIT(isend11,istatus,ierr)
      call MPI_WAIT(isend12,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)
      call MPI_WAIT(irecv07,istatus,ierr)
      call MPI_WAIT(irecv08,istatus,ierr)
      call MPI_WAIT(irecv09,istatus,ierr)
      call MPI_WAIT(irecv10,istatus,ierr)
      call MPI_WAIT(irecv11,istatus,ierr)
      call MPI_WAIT(irecv12,istatus,ierr)


      ! Corner send/recv
      pidto01=get_pid(pidx-1,pidy-1,pidz-1) !a_lll_s 
      pidto02=get_pid(pidx+1,pidy-1,pidz-1) !a_ull_s 
      pidto03=get_pid(pidx-1,pidy+1,pidz-1) !a_lul_s 
      pidto04=get_pid(pidx+1,pidy+1,pidz-1) !a_uul_s 
      pidto05=get_pid(pidx-1,pidy-1,pidz+1) !a_llu_s 
      pidto06=get_pid(pidx+1,pidy-1,pidz+1) !a_ulu_s 
      pidto07=get_pid(pidx-1,pidy+1,pidz+1) !a_luu_s 
      pidto08=get_pid(pidx+1,pidy+1,pidz+1) !a_uuu_s 
      call MPI_ISEND(a_lll_s(1,1,1),size(a_lll_s),MPI_INTEGER,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_ull_s(1,1,1),size(a_ull_s),MPI_INTEGER,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_lul_s(1,1,1),size(a_lul_s),MPI_INTEGER,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_uul_s(1,1,1),size(a_uul_s),MPI_INTEGER,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_llu_s(1,1,1),size(a_llu_s),MPI_INTEGER,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_ulu_s(1,1,1),size(a_ulu_s),MPI_INTEGER,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
      call MPI_ISEND(a_luu_s(1,1,1),size(a_luu_s),MPI_INTEGER,
     & pidto07,basetag+7,MPI_COMM_WORLD,isend07,ierr)
      call MPI_ISEND(a_uuu_s(1,1,1),size(a_uuu_s),MPI_INTEGER,
     & pidto08,basetag+8,MPI_COMM_WORLD,isend08,ierr)

      pidfrom01=get_pid(pidx-1,pidy-1,pidz-1) !a_lll_r 
      pidfrom02=get_pid(pidx+1,pidy-1,pidz-1) !a_ull_r 
      pidfrom03=get_pid(pidx-1,pidy+1,pidz-1) !a_lul_r 
      pidfrom04=get_pid(pidx+1,pidy+1,pidz-1) !a_uul_r 
      pidfrom05=get_pid(pidx-1,pidy-1,pidz+1) !a_llu_r 
      pidfrom06=get_pid(pidx+1,pidy-1,pidz+1) !a_ulu_r 
      pidfrom07=get_pid(pidx-1,pidy+1,pidz+1) !a_luu_r 
      pidfrom08=get_pid(pidx+1,pidy+1,pidz+1) !a_uuu_r 
      call MPI_IRECV(a_lll_r(1,1,1),size(a_lll_r),MPI_INTEGER,
     & pidfrom01,basetag+8,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_ull_r(1,1,1),size(a_ull_r),MPI_INTEGER,
     & pidfrom02,basetag+7,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_lul_r(1,1,1),size(a_lul_r),MPI_INTEGER,
     & pidfrom03,basetag+6,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_uul_r(1,1,1),size(a_uul_r),MPI_INTEGER,
     & pidfrom04,basetag+5,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_llu_r(1,1,1),size(a_llu_r),MPI_INTEGER,
     & pidfrom05,basetag+4,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_ulu_r(1,1,1),size(a_ulu_r),MPI_INTEGER,
     & pidfrom06,basetag+3,MPI_COMM_WORLD,irecv06,ierr)
      call MPI_IRECV(a_luu_r(1,1,1),size(a_luu_r),MPI_INTEGER,
     & pidfrom07,basetag+2,MPI_COMM_WORLD,irecv07,ierr)
      call MPI_IRECV(a_uuu_r(1,1,1),size(a_uuu_r),MPI_INTEGER,
     & pidfrom08,basetag+1,MPI_COMM_WORLD,irecv08,ierr)

      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(isend07,istatus,ierr)
      call MPI_WAIT(isend08,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)
      call MPI_WAIT(irecv07,istatus,ierr)
      call MPI_WAIT(irecv08,istatus,ierr)

      !Update the ghost layers with the data from neighboring processes
      !update planes
      call update_plane_i4(a, 'il', a_il_r, 1,jml,1,kml,ng,'regular')
      call update_plane_i4(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'regular')
      call update_plane_i4(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'regular')
      call update_plane_i4(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'regular')
      call update_plane_i4(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'regular')
      call update_plane_i4(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'regular')

      !update edges
      call update_edge_i4(a, 'iljl', a_iljl_r, 1,kml,ng,'regular')
      call update_edge_i4(a, 'ilju', a_ilju_r, 1,kml,ng,'regular')
      call update_edge_i4(a, 'iujl', a_iujl_r, 1,kml,ng,'regular')
      call update_edge_i4(a, 'iuju', a_iuju_r, 1,kml,ng,'regular')
      call update_edge_i4(a, 'jlkl', a_jlkl_r, 1,iml,ng,'regular')
      call update_edge_i4(a, 'jlku', a_jlku_r, 1,iml,ng,'regular')
      call update_edge_i4(a, 'jukl', a_jukl_r, 1,iml,ng,'regular')
      call update_edge_i4(a, 'juku', a_juku_r, 1,iml,ng,'regular')
      call update_edge_i4(a, 'klil', a_klil_r, 1,jml,ng,'regular')
      call update_edge_i4(a, 'kliu', a_kliu_r, 1,jml,ng,'regular')
      call update_edge_i4(a, 'kuil', a_kuil_r, 1,jml,ng,'regular')
      call update_edge_i4(a, 'kuiu', a_kuiu_r, 1,jml,ng,'regular')

      !corners
      call update_corner_i4(a, 'lll', a_lll_r, ng,'regular')
      call update_corner_i4(a, 'ull', a_ull_r, ng,'regular')
      call update_corner_i4(a, 'lul', a_lul_r, ng,'regular')
      call update_corner_i4(a, 'uul', a_uul_r, ng,'regular')
      call update_corner_i4(a, 'llu', a_llu_r, ng,'regular')
      call update_corner_i4(a, 'ulu', a_ulu_r, ng,'regular')
      call update_corner_i4(a, 'luu', a_luu_r, ng,'regular')
      call update_corner_i4(a, 'uuu', a_uuu_r, ng,'regular')

      ! Actual boundary condition is applied here.
      ! BC mode names(=bctype)
      ! periodic: periodic BC along all directions(default)
      ! adiabatic: adiabatic BC for all boundaries
      ! adiabatic_x: adiabatic BC along x direction (for il & iu), but periodic for the rest boundaries
      ! adiabatic_y: adiabatic BC along y direction (for jl & ju), but periodic for the rest boundaries
      ! adiabatic_z: adiabatic BC along z direction (for kl & ku), but periodic for the rest boundaries
      ! adiabatic_xy: adiabatic BC along x and y directions (il,iu,jl,ju), but periodic for z direction (kl,ku)
      ! adiabatic_yz: adiabatic BC along y and z directions (jl,ju,kl,ku), but periodic for x direction (il,iu)
      ! adiabatic_zx: adiabatic BC along z and x directions (kl,ku,il,iu), but periodic for y direction (jl,ju)
      ! Further BC modes will be coming up...
      select case(trim(bctype))
      case('adiabatic')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_i4(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_i4(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_i4(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_i4(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_i4(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_i4(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      ! Global iljl edge
      if(pidx.eq.0.and.pidy.eq.0)then
      call extract_edge_i4(a_iljl_s,'iljl', a, 1,kml,ng)
      call update_edge_i4(a, 'iljl', a_iljl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0)then
      ! Global ilju edge
      if(pidx.eq.0.and.pidy.eq.NPY-1)then
      call extract_edge_i4(a_ilju_s,'ilju', a, 1,kml,ng)
      call update_edge_i4(a, 'ilju', a_ilju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1)then 
      ! Global iujl edge
      if(pidx.eq.NPX-1.and.pidy.eq.0)then
      call extract_edge_i4(a_iujl_s,'iujl', a, 1,kml,ng)
      call update_edge_i4(a, 'iujl', a_iujl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0)then
      ! Global iuju edge
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      call extract_edge_i4(a_iuju_s,'iuju', a, 1,kml,ng)
      call update_edge_i4(a, 'iuju', a_iuju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      ! Global jlkl edge
      if(pidy.eq.0.and.pidz.eq.0)then
      call extract_edge_i4(a_jlkl_s,'jlkl', a, 1,iml,ng)
      call update_edge_i4(a, 'jlkl', a_jlkl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.0)then
      ! Global jlku edge
      if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_edge_i4(a_jlku_s,'jlku', a, 1,iml,ng)
      call update_edge_i4(a, 'jlku', a_jlku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global jukl edge
      if(pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_edge_i4(a_jukl_s,'jukl', a, 1,iml,ng)
      call update_edge_i4(a, 'jukl', a_jukl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global juku edge
      if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_edge_i4(a_juku_s,'juku', a, 1,iml,ng)
      call update_edge_i4(a, 'juku', a_juku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      ! Global klil edge
      if(pidz.eq.0.and.pidx.eq.0)then
      call extract_edge_i4(a_klil_s,'klil', a, 1,jml,ng)
      call update_edge_i4(a, 'klil', a_klil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.0)then
      ! Global kliu edge
      if(pidz.eq.0.and.pidx.eq.NPX-1)then
      call extract_edge_i4(a_kliu_s,'kliu', a, 1,jml,ng)
      call update_edge_i4(a, 'kliu', a_kliu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.NPX-1)then
      ! Global kuil edge
      if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      call extract_edge_i4(a_kuil_s,'kuil', a, 1,jml,ng)
      call update_edge_i4(a, 'kuil', a_kuil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      ! Global kuiu edge
      if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      call extract_edge_i4(a_kuiu_s,'kuiu', a, 1,jml,ng)
      call update_edge_i4(a, 'kuiu', a_kuiu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then

      ! Global lll corner
      if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.0)then
      call extract_corner_i4(a_lll_s, 'lll', a, ng)
      call update_corner_i4(a, 'lll', a_lll_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.0)then
      ! Global ull corner
      if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.0)then
      call extract_corner_i4(a_ull_s, 'ull', a, ng)
      call update_corner_i4(a, 'ull', a_ull_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.0)then
      ! Global lul corner
      if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_corner_i4(a_lul_s, 'lul', a, ng)
      call update_corner_i4(a, 'lul', a_lul_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global uul corner
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_corner_i4(a_uul_s, 'uul', a, ng)
      call update_corner_i4(a, 'uul', a_uul_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global llu corner
      if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_corner_i4(a_llu_s, 'llu', a, ng)
      call update_corner_i4(a, 'llu', a_llu_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global ulu corner
      if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_corner_i4(a_ulu_s, 'ulu', a, ng)
      call update_corner_i4(a, 'ulu', a_ulu_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global luu corner
      if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_corner_i4(a_luu_s, 'luu', a, ng)
      call update_corner_i4(a, 'luu', a_luu_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      ! Global uuu corner
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_corner_i4(a_uuu_s, 'uuu', a, ng)
      call update_corner_i4(a, 'uuu', a_uuu_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then

      case('adiabatic_x')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_i4(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_i4(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then

      case('adiabatic_y')
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_i4(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_i4(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then

      case('adiabatic_z')
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_i4(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_i4(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      case('adiabatic_xy')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_i4(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_i4(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_i4(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_i4(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then

      ! Global iljl edge
      if(pidx.eq.0.and.pidy.eq.0)then
      call extract_edge_i4(a_iljl_s,'iljl', a, 1,kml,ng)
      call update_edge_i4(a, 'iljl', a_iljl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0)then
      ! Global ilju edge
      if(pidx.eq.0.and.pidy.eq.NPY-1)then
      call extract_edge_i4(a_ilju_s,'ilju', a, 1,kml,ng)
      call update_edge_i4(a, 'ilju', a_ilju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1)then 
      ! Global iujl edge
      if(pidx.eq.NPX-1.and.pidy.eq.0)then
      call extract_edge_i4(a_iujl_s,'iujl', a, 1,kml,ng)
      call update_edge_i4(a, 'iujl', a_iujl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0)then
      ! Global iuju edge
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      call extract_edge_i4(a_iuju_s,'iuju', a, 1,kml,ng)
      call update_edge_i4(a, 'iuju', a_iuju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then

      case('adiabatic_yz')
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_i4(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_i4(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_i4(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_i4(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_i4(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      ! Global jlkl edge
      if(pidy.eq.0.and.pidz.eq.0)then
      call extract_edge_i4(a_jlkl_s,'jlkl', a, 1,iml,ng)
      call update_edge_i4(a, 'jlkl', a_jlkl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.0)then
      ! Global jlku edge
      if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_edge_i4(a_jlku_s,'jlku', a, 1,iml,ng)
      call update_edge_i4(a, 'jlku', a_jlku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global jukl edge
      if(pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_edge_i4(a_jukl_s,'jukl', a, 1,iml,ng)
      call update_edge_i4(a, 'jukl', a_jukl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global juku edge
      if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_edge_i4(a_juku_s,'juku', a, 1,iml,ng)
      call update_edge_i4(a, 'juku', a_juku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then

      case('adiabatic_zx')
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_i4(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_i4(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_i4(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_i4(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_i4(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_i4(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then

      ! Global klil edge
      if(pidz.eq.0.and.pidx.eq.0)then
      call extract_edge_i4(a_klil_s,'klil', a, 1,jml,ng)
      call update_edge_i4(a, 'klil', a_klil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.0)then
      ! Global kliu edge
      if(pidz.eq.0.and.pidx.eq.NPX-1)then
      call extract_edge_i4(a_kliu_s,'kliu', a, 1,jml,ng)
      call update_edge_i4(a, 'kliu', a_kliu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.NPX-1)then
      ! Global kuil edge
      if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      call extract_edge_i4(a_kuil_s,'kuil', a, 1,jml,ng)
      call update_edge_i4(a, 'kuil', a_kuil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      ! Global kuiu edge
      if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      call extract_edge_i4(a_kuiu_s,'kuiu', a, 1,jml,ng)
      call update_edge_i4(a, 'kuiu', a_kuiu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      end select

      end subroutine boundary_condition_i4_par



      subroutine print_ghost_layers_i4(a,ng)
      use input
      implicit none

      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      integer(4) :: a_il(1-ng:jml+ng,1-ng:kml+ng)
      integer(4) :: a_iu(1-ng:jml+ng,1-ng:kml+ng)
      integer(4) :: a_jl(1-ng:kml+ng,1-ng:iml+ng)
      integer(4) :: a_ju(1-ng:kml+ng,1-ng:iml+ng)
      integer(4) :: a_kl(1-ng:iml+ng,1-ng:jml+ng)
      integer(4) :: a_ku(1-ng:iml+ng,1-ng:jml+ng)

      integer(4) :: i,j,k

      write(*,*) "Subroutine print_ghost_layers() called"

      !il and iu
      write(*,*) 'il and iu'
      do k=1-ng,kml+ng
      do j=1-ng,jml+ng
            a_il(j,k)=a(1-ng,j,k)
            a_iu(j,k)=a(iml+ng,j,k)
            write(*,*)j,k,a_il(j,k),a_iu(j,k)
      enddo
      enddo

      !jl and ju
      write(*,*) 'jl and ju'
      do i=1-ng,iml+ng
      do k=1-ng,kml+ng
            a_jl(k,i)=a(i,1-ng,k)
            a_ju(k,i)=a(i,jml+ng,k)
            write(*,*)k,i,a_jl(k,i),a_ju(k,i)
      enddo
      enddo

      !kl and ku
      write(*,*) 'kl and ku'
      do j=1-ng,jml+ng
      do i=1-ng,iml+ng
            a_kl(i,j)=a(i,j,1-ng)
            a_ku(i,j)=a(i,j,kml+ng)
            write(*,*)i,j,a_kl(i,j),a_ku(i,j)
      enddo
      enddo

      end subroutine print_ghost_layers_i4


      subroutine print_array_i4(a,ng,il,iu,jl,ju,kl,ku)
      use input
      implicit none
      integer(4),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: il,iu,jl,ju,kl,ku

      integer(4) :: i,j,k,n1

      write(*,*)'print_array() called from rank = ',myrank

      do k=kl,ku
            write(*,*)'k=',k,", i: horizontal, j: vertical"
            do j=jl,ju
                  write(*,*) (a(i,j,k), i=il,iu)
            enddo
      enddo

      end subroutine print_array_i4


      subroutine extract_plane_i4(aout, planetype, ain, il,iu,jl,ju,ng)
      use input
      implicit none

      integer(4),intent(inout) :: aout(il-ng:iu+ng,jl-ng:ju+ng,ng)
      character(len=*),intent(in) :: planetype
      integer(4),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: il,iu,jl,ju,ng

      integer(4) :: i,j,k,n1

      select case(trim(planetype))
      case('il')
            do n1=1,ng
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  aout(j,k,n1)=ain(n1,j,k)
            enddo
            enddo
            enddo
      case('iu')
            do n1=1,ng
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  aout(j,k,n1)=ain(iml-ng+n1,j,k)
            enddo
            enddo
            enddo
      case('jl')
            do n1=1,ng
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  aout(k,i,n1)=ain(i,n1,k)
            enddo
            enddo
            enddo
      case('ju')
            do n1=1,ng
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  aout(k,i,n1)=ain(i,jml-ng+n1,k)
            enddo
            enddo
            enddo
      case('kl')
            do n1=1,ng
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  aout(i,j,n1)=ain(i,j,n1)
            enddo
            enddo
            enddo
      case('ku')
            do n1=1,ng
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  aout(i,j,n1)=ain(i,j,kml-ng+n1)
            enddo
            enddo
            enddo
      end select
      end subroutine extract_plane_i4


      subroutine extract_edge_i4(aout, edgetype, ain, il,iu,ng)
      use input
      implicit none
      integer(4),intent(inout) :: aout(il-ng:iu+ng,ng,ng)
      character(len=*),intent(in) :: edgetype
      integer(4),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: il,iu,ng

      integer(4) :: ll,n1,n2

      select case(trim(edgetype))
      case('iljl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n1,n2,ll)
                  enddo
            enddo
            enddo
      case('ilju')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n1,jml-ng+n2,ll)
                  enddo
            enddo
            enddo
      case('iujl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n1,n2,ll)
                  enddo
            enddo
            enddo
      case('iuju')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n1,jml-ng+n2,ll)
                  enddo
            enddo
            enddo
      case('jlkl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,n1,n2)
                  enddo
            enddo
            enddo
      case('jlku')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,n1,kml-ng+n2)
                  enddo
            enddo
            enddo
      case('jukl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,jml-ng+n1,n2)
                  enddo
            enddo
            enddo
      case('juku')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,jml-ng+n1,kml-ng+n2)
                  enddo
            enddo
            enddo
      case('klil')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n2,ll,n1)
                  enddo
            enddo
            enddo
      case('kliu')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n2,ll,n1)
                  enddo
            enddo
            enddo
      case('kuil')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n2,ll,kml-ng+n1)
                  enddo
            enddo
            enddo
      case('kuiu')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n2,ll,kml-ng+n1)
                  enddo
            enddo
            enddo
      end select

      end subroutine extract_edge_i4


      subroutine extract_corner_i4(aout, cornertype, ain, ng)
      use input
      implicit none
      integer(4),intent(inout) :: aout(ng,ng,ng)
      character(len=*),intent(in) :: cornertype
      integer(4),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      integer(4) :: n1,n2,n3

      select case(trim(cornertype))
      case('lll')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,n2,n3)
            enddo
            enddo
            enddo
      case('ull')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,n2,n3)
            enddo
            enddo
            enddo
      case('lul')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,jml-ng+n2,n3)
            enddo
            enddo
            enddo
      case('uul')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,jml-ng+n2,n3)
            enddo
            enddo
            enddo
      case('llu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('ulu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('luu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,jml-ng+n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('uuu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,jml-ng+n2,kml-ng+n3)
            enddo
            enddo
            enddo
      end select

      end subroutine extract_corner_i4

      subroutine update_plane_i4(aout, planetype, ain 
     &                                     , il,iu,jl,ju,ng,copymode)
      use input
      implicit none
      integer(4),intent(inout) :: 
     &                   aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: planetype  !plane name
      integer(4),intent(inout) :: ain(il-ng:iu+ng,jl-ng:ju+ng,ng)
      integer(4),intent(in) :: il,iu,jl,ju,ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: i,j,k,n1

      select case(trim(copymode))
      case('regular')
            select case(trim(planetype))
            case('il')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(1-ng+n1,j,k)=ain(j,k,1+n1)
                  enddo
                  enddo
                  enddo
            case('iu')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(iml+1+n1,j,k)=ain(j,k,1+n1)
                  enddo
                  enddo
                  enddo
            case('jl')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,1-ng+n1,k)=ain(k,i,1+n1)
                  enddo
                  enddo
                  enddo
            case('ju')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,jml+1+n1,k)=ain(k,i,1+n1)
                  enddo
                  enddo
                  enddo
            case('kl')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,1-ng+n1)=ain(i,j,1+n1)
                  enddo
                  enddo
                  enddo
            case('ku')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,kml+1+n1)=ain(i,j,1+n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(planetype))
      case('mirror')
            select case(trim(planetype))
            case('il')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(1-ng+n1,j,k)=ain(j,k,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iu')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(iml+1+n1,j,k)=ain(j,k,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jl')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,1-ng+n1,k)=ain(k,i,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ju')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,jml+1+n1,k)=ain(k,i,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kl')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,1-ng+n1)=ain(i,j,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ku')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,kml+1+n1)=ain(i,j,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(planetype))
      end select !select case(trim(copymode))

      end subroutine update_plane_i4


      subroutine update_edge_i4(aout, edgetype, ain, il,iu,ng,copymode)
      use input
      implicit none

      integer(4),intent(inout):: 
     &                   aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: edgetype  !edge name
      integer(4),intent(inout) :: ain(il-ng:iu+ng,ng,ng)
      integer(4),intent(in) :: il,iu,ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: ll,n1,n2

      select case(trim(copymode))
      case('regular')
            select case(trim(edgetype))
            case('iljl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,1-ng+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('ilju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,jml+1+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('iujl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,1-ng+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('iuju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,jml+1+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jlkl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,1-ng+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jlku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,kml+1+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jukl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,1-ng+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('juku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,kml+1+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('klil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,1-ng+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kliu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,1-ng+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kuil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,kml+1+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kuiu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,kml+1+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(edgetype))
      case('mirror')
            select case(trim(edgetype))
            case('iljl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,1-ng+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ilju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,jml+1+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iujl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,1-ng+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iuju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,jml+1+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jlkl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,1-ng+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jlku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,kml+1+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jukl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,1-ng+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('juku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,kml+1+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('klil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,1-ng+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kliu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,1-ng+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kuil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,kml+1+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kuiu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,kml+1+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(edgetype))
      end select !select case(trim(copymode))

      end subroutine update_edge_i4


      subroutine update_corner_i4(aout, cornertype, ain, ng,copymode)
      use input
      implicit none
      integer(4),intent(inout) :: 
     &                  aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(inout) :: ain(ng,ng,ng)
      character(len=*),intent(in) :: cornertype !corner name
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: n1,n2,n3

      select case(trim(copymode))
      case('regular')
            select case(trim(cornertype))
            case('lll')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('ull')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('lul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('uul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('llu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('ulu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('luu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('uuu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(cornertype))
      case('mirror')
            select case(trim(cornertype))
            case('lll')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ull')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('lul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('uul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('llu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ulu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('luu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('uuu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,kml+1+n3)= 
     &                                     ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(cornertype))
      end select !select case(trim(copymode))

      end subroutine update_corner_i4



      subroutine boundary_condition_r8_par(bctype,a,ng)
      use input
      implicit none
      include 'mpif.h'

      !Subroutine arguments
      character(len=*),intent(in) :: bctype
      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      !Local variables
      !Contiguous arrays for MPI communications, s: send, r: recv
      !planes for send
      real(8) :: a_il_s(1-ng:jml+ng,1-ng:kml+ng,ng)
      real(8) :: a_iu_s(1-ng:jml+ng,1-ng:kml+ng,ng)
      real(8) :: a_jl_s(1-ng:kml+ng,1-ng:iml+ng,ng)
      real(8) :: a_ju_s(1-ng:kml+ng,1-ng:iml+ng,ng)
      real(8) :: a_kl_s(1-ng:iml+ng,1-ng:jml+ng,ng)
      real(8) :: a_ku_s(1-ng:iml+ng,1-ng:jml+ng,ng)
      !edges for send
      real(8) :: a_iljl_s(1-ng:kml+ng,ng,ng) !edge intersecting il and jl plane
      real(8) :: a_ilju_s(1-ng:kml+ng,ng,ng) !edge intersecting il and ju plane
      real(8) :: a_iujl_s(1-ng:kml+ng,ng,ng) !edge intersecting iu and jl plane
      real(8) :: a_iuju_s(1-ng:kml+ng,ng,ng) !edge intersecting iu and ju plane
      real(8) :: a_jlkl_s(1-ng:iml+ng,ng,ng) !edge intersecting jl and kl plane
      real(8) :: a_jlku_s(1-ng:iml+ng,ng,ng) !edge intersecting jl and ku plane
      real(8) :: a_jukl_s(1-ng:iml+ng,ng,ng) !edge intersecting ju and kl plane
      real(8) :: a_juku_s(1-ng:iml+ng,ng,ng) !edge intersecting ju and ku plane
      real(8) :: a_klil_s(1-ng:jml+ng,ng,ng) !edge intersecting kl and il plane
      real(8) :: a_kliu_s(1-ng:jml+ng,ng,ng) !edge intersecting kl and iu plane
      real(8) :: a_kuil_s(1-ng:jml+ng,ng,ng) !edge intersecting ku and il plane
      real(8) :: a_kuiu_s(1-ng:jml+ng,ng,ng) !edge intersecting ku and iu plane
      !corners for send
      real(8) :: a_lll_s(ng,ng,ng) !(i,j,k)= nearby( 1, 1, 1)
      real(8) :: a_ull_s(ng,ng,ng) !(i,j,k)= nearby(im, 1, 1)
      real(8) :: a_lul_s(ng,ng,ng) !(i,j,k)= nearby( 1,jm, 1)
      real(8) :: a_uul_s(ng,ng,ng) !(i,j,k)= nearby(im,jm, 1)
      real(8) :: a_llu_s(ng,ng,ng) !(i,j,k)= nearby( 1, 1,km)
      real(8) :: a_ulu_s(ng,ng,ng) !(i,j,k)= nearby(im, 1,km)
      real(8) :: a_luu_s(ng,ng,ng) !(i,j,k)= nearby( 1,jm,km)
      real(8) :: a_uuu_s(ng,ng,ng) !(i,j,k)= nearby(im,jm,km)

      !planes for recv
      real(8) :: a_il_r(1-ng:jml+ng,1-ng:kml+ng,ng)
      real(8) :: a_iu_r(1-ng:jml+ng,1-ng:kml+ng,ng)
      real(8) :: a_jl_r(1-ng:kml+ng,1-ng:iml+ng,ng)
      real(8) :: a_ju_r(1-ng:kml+ng,1-ng:iml+ng,ng)
      real(8) :: a_kl_r(1-ng:iml+ng,1-ng:jml+ng,ng)
      real(8) :: a_ku_r(1-ng:iml+ng,1-ng:jml+ng,ng)
      !edges for recv
      real(8) :: a_iljl_r(1-ng:kml+ng,ng,ng) !edge intersecting il and jl plane
      real(8) :: a_ilju_r(1-ng:kml+ng,ng,ng) !edge intersecting il and ju plane
      real(8) :: a_iujl_r(1-ng:kml+ng,ng,ng) !edge intersecting iu and jl plane
      real(8) :: a_iuju_r(1-ng:kml+ng,ng,ng) !edge intersecting iu and ju plane
      real(8) :: a_jlkl_r(1-ng:iml+ng,ng,ng) !edge intersecting jl and kl plane
      real(8) :: a_jlku_r(1-ng:iml+ng,ng,ng) !edge intersecting jl and ku plane
      real(8) :: a_jukl_r(1-ng:iml+ng,ng,ng) !edge intersecting ju and kl plane
      real(8) :: a_juku_r(1-ng:iml+ng,ng,ng) !edge intersecting ju and ku plane
      real(8) :: a_klil_r(1-ng:jml+ng,ng,ng) !edge intersecting kl and il plane
      real(8) :: a_kliu_r(1-ng:jml+ng,ng,ng) !edge intersecting kl and iu plane
      real(8) :: a_kuil_r(1-ng:jml+ng,ng,ng) !edge intersecting ku and il plane
      real(8) :: a_kuiu_r(1-ng:jml+ng,ng,ng) !edge intersecting ku and iu plane
      !corners for recv
      real(8) :: a_lll_r(ng,ng,ng) !(i,j,k)= nearby( 1, 1, 1)
      real(8) :: a_ull_r(ng,ng,ng) !(i,j,k)= nearby(im, 1, 1)
      real(8) :: a_lul_r(ng,ng,ng) !(i,j,k)= nearby( 1,jm, 1)
      real(8) :: a_uul_r(ng,ng,ng) !(i,j,k)= nearby(im,jm, 1)
      real(8) :: a_llu_r(ng,ng,ng) !(i,j,k)= nearby( 1, 1,km)
      real(8) :: a_ulu_r(ng,ng,ng) !(i,j,k)= nearby(im, 1,km)
      real(8) :: a_luu_r(ng,ng,ng) !(i,j,k)= nearby( 1,jm,km)
      real(8) :: a_uuu_r(ng,ng,ng) !(i,j,k)= nearby(im,jm,km)

      integer(4) :: i,j,k,n1,n2

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      !Process coordinate
      call get_pidxyz(pidx, pidy, pidz, myrank)

      !Preparation for send buffers
      ! Planes
      call extract_plane_r8(a_il_s, 'il', a, 1,jml,1,kml,ng)
      call extract_plane_r8(a_iu_s, 'iu', a, 1,jml,1,kml,ng)
      call extract_plane_r8(a_jl_s, 'jl', a, 1,kml,1,iml,ng)
      call extract_plane_r8(a_ju_s, 'ju', a, 1,kml,1,iml,ng)
      call extract_plane_r8(a_kl_s, 'kl', a, 1,iml,1,jml,ng)
      call extract_plane_r8(a_ku_s, 'ku', a, 1,iml,1,jml,ng)

      !edges
      call extract_edge_r8(a_iljl_s, 'iljl', a, 1,kml,ng)
      call extract_edge_r8(a_ilju_s, 'ilju', a, 1,kml,ng)
      call extract_edge_r8(a_iujl_s, 'iujl', a, 1,kml,ng)
      call extract_edge_r8(a_iuju_s, 'iuju', a, 1,kml,ng)
      call extract_edge_r8(a_jlkl_s, 'jlkl', a, 1,iml,ng)
      call extract_edge_r8(a_jlku_s, 'jlku', a, 1,iml,ng)
      call extract_edge_r8(a_jukl_s, 'jukl', a, 1,iml,ng)
      call extract_edge_r8(a_juku_s, 'juku', a, 1,iml,ng)
      call extract_edge_r8(a_klil_s, 'klil', a, 1,jml,ng)
      call extract_edge_r8(a_kliu_s, 'kliu', a, 1,jml,ng)
      call extract_edge_r8(a_kuil_s, 'kuil', a, 1,jml,ng)
      call extract_edge_r8(a_kuiu_s, 'kuiu', a, 1,jml,ng)

      !corners
      call extract_corner_r8(a_lll_s, 'lll', a, ng)
      call extract_corner_r8(a_ull_s, 'ull', a, ng)
      call extract_corner_r8(a_lul_s, 'lul', a, ng)
      call extract_corner_r8(a_uul_s, 'uul', a, ng)
      call extract_corner_r8(a_llu_s, 'llu', a, ng)
      call extract_corner_r8(a_ulu_s, 'ulu', a, ng)
      call extract_corner_r8(a_luu_s, 'luu', a, ng)
      call extract_corner_r8(a_uuu_s, 'uuu', a, ng)

      ! Layer send/recv
      ! i planes
      pidto01=get_pid(pidx-1,pidy,pidz) ! a_il_s
      pidto02=get_pid(pidx+1,pidy,pidz) ! a_iu_s
      pidto03=get_pid(pidx,pidy-1,pidz) ! a_jl_s
      pidto04=get_pid(pidx,pidy+1,pidz) ! a_ju_s
      pidto05=get_pid(pidx,pidy,pidz-1) ! a_kl_s
      pidto06=get_pid(pidx,pidy,pidz+1) ! a_ku_s
      call MPI_ISEND(a_il_s(1-ng,1-ng,1),size(a_il_s),
     & MPI_DOUBLE_PRECISION,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_iu_s(1-ng,1-ng,1),size(a_iu_s),
     & MPI_DOUBLE_PRECISION,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_jl_s(1-ng,1-ng,1),size(a_jl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_ju_s(1-ng,1-ng,1),size(a_ju_s),
     & MPI_DOUBLE_PRECISION,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_kl_s(1-ng,1-ng,1),size(a_kl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_ku_s(1-ng,1-ng,1),size(a_ku_s),
     & MPI_DOUBLE_PRECISION,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

      pidfrom01=get_pid(pidx-1,pidy,pidz) ! a_il_r
      pidfrom02=get_pid(pidx+1,pidy,pidz) ! a_iu_r
      pidfrom03=get_pid(pidx,pidy-1,pidz) ! a_jl_r
      pidfrom04=get_pid(pidx,pidy+1,pidz) ! a_ju_r
      pidfrom05=get_pid(pidx,pidy,pidz-1) ! a_kl_r
      pidfrom06=get_pid(pidx,pidy,pidz+1) ! a_ku_r
      call MPI_IRECV(a_il_r(1-ng,1-ng,1),size(a_il_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_iu_r(1-ng,1-ng,1),size(a_iu_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_jl_r(1-ng,1-ng,1),size(a_jl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_ju_r(1-ng,1-ng,1),size(a_ju_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_kl_r(1-ng,1-ng,1),size(a_kl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_ku_r(1-ng,1-ng,1),size(a_ku_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)

      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)


      ! Edge send/recv
      pidto01=get_pid(pidx-1,pidy-1,pidz) ! a_iljl_s
      pidto02=get_pid(pidx-1,pidy+1,pidz) ! a_ilju_s
      pidto03=get_pid(pidx+1,pidy-1,pidz) ! a_iujl_s
      pidto04=get_pid(pidx+1,pidy+1,pidz) ! a_iuju_s
      pidto05=get_pid(pidx,pidy-1,pidz-1) ! a_jlkl_s
      pidto06=get_pid(pidx,pidy-1,pidz+1) ! a_jlku_s
      pidto07=get_pid(pidx,pidy+1,pidz-1) ! a_jukl_s
      pidto08=get_pid(pidx,pidy+1,pidz+1) ! a_juku_s
      pidto09=get_pid(pidx-1,pidy,pidz-1) ! a_klil_s
      pidto10=get_pid(pidx+1,pidy,pidz-1) ! a_kliu_s
      pidto11=get_pid(pidx-1,pidy,pidz+1) ! a_kuil_s
      pidto12=get_pid(pidx+1,pidy,pidz+1) ! a_kuiu_s
      call MPI_ISEND(a_iljl_s(1-ng,1,1),size(a_iljl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_ilju_s(1-ng,1,1),size(a_ilju_s),
     & MPI_DOUBLE_PRECISION,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_iujl_s(1-ng,1,1),size(a_iujl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_iuju_s(1-ng,1,1),size(a_iuju_s),
     & MPI_DOUBLE_PRECISION,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_jlkl_s(1-ng,1,1),size(a_jlkl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_jlku_s(1-ng,1,1),size(a_jlku_s),
     & MPI_DOUBLE_PRECISION,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
      call MPI_ISEND(a_jukl_s(1-ng,1,1),size(a_jukl_s),
     & MPI_DOUBLE_PRECISION,
     & pidto07,basetag+7,MPI_COMM_WORLD,isend07,ierr)
      call MPI_ISEND(a_juku_s(1-ng,1,1),size(a_juku_s),
     & MPI_DOUBLE_PRECISION,
     & pidto08,basetag+8,MPI_COMM_WORLD,isend08,ierr)
      call MPI_ISEND(a_klil_s(1-ng,1,1),size(a_klil_s),
     & MPI_DOUBLE_PRECISION,
     & pidto09,basetag+9,MPI_COMM_WORLD,isend09,ierr)
      call MPI_ISEND(a_kliu_s(1-ng,1,1),size(a_kliu_s),
     & MPI_DOUBLE_PRECISION,
     & pidto10,basetag+10,MPI_COMM_WORLD,isend10,ierr)
      call MPI_ISEND(a_kuil_s(1-ng,1,1),size(a_kuil_s),
     & MPI_DOUBLE_PRECISION,
     & pidto11,basetag+11,MPI_COMM_WORLD,isend11,ierr)
      call MPI_ISEND(a_kuiu_s(1-ng,1,1),size(a_kuiu_s),
     & MPI_DOUBLE_PRECISION,
     & pidto12,basetag+12,MPI_COMM_WORLD,isend12,ierr)

      pidfrom01=get_pid(pidx-1,pidy-1,pidz) ! a_iljl_r
      pidfrom02=get_pid(pidx-1,pidy+1,pidz) ! a_ilju_r
      pidfrom03=get_pid(pidx+1,pidy-1,pidz) ! a_iujl_r
      pidfrom04=get_pid(pidx+1,pidy+1,pidz) ! a_iuju_r
      pidfrom05=get_pid(pidx,pidy-1,pidz-1) ! a_jlkl_r
      pidfrom06=get_pid(pidx,pidy-1,pidz+1) ! a_jlku_r
      pidfrom07=get_pid(pidx,pidy+1,pidz-1) ! a_jukl_r
      pidfrom08=get_pid(pidx,pidy+1,pidz+1) ! a_juku_r
      pidfrom09=get_pid(pidx-1,pidy,pidz-1) ! a_klil_r
      pidfrom10=get_pid(pidx+1,pidy,pidz-1) ! a_kliu_r
      pidfrom11=get_pid(pidx-1,pidy,pidz+1) ! a_kuil_r
      pidfrom12=get_pid(pidx+1,pidy,pidz+1) ! a_kuiu_r
      call MPI_IRECV(a_iljl_r(1-ng,1,1),size(a_iljl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom01,basetag+4,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_ilju_r(1-ng,1,1),size(a_ilju_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom02,basetag+3,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_iujl_r(1-ng,1,1),size(a_iujl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom03,basetag+2,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_iuju_r(1-ng,1,1),size(a_iuju_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom04,basetag+1,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_jlkl_r(1-ng,1,1),size(a_jlkl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom05,basetag+8,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_jlku_r(1-ng,1,1),size(a_jlku_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom06,basetag+7,MPI_COMM_WORLD,irecv06,ierr)
      call MPI_IRECV(a_jukl_r(1-ng,1,1),size(a_jukl_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom07,basetag+6,MPI_COMM_WORLD,irecv07,ierr)
      call MPI_IRECV(a_juku_r(1-ng,1,1),size(a_juku_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom08,basetag+5,MPI_COMM_WORLD,irecv08,ierr)
      call MPI_IRECV(a_klil_r(1-ng,1,1),size(a_klil_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom09,basetag+12,MPI_COMM_WORLD,irecv09,ierr)
      call MPI_IRECV(a_kliu_r(1-ng,1,1),size(a_kliu_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom10,basetag+11,MPI_COMM_WORLD,irecv10,ierr)
      call MPI_IRECV(a_kuil_r(1-ng,1,1),size(a_kuil_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom11,basetag+10,MPI_COMM_WORLD,irecv11,ierr)
      call MPI_IRECV(a_kuiu_r(1-ng,1,1),size(a_kuiu_r),
     & MPI_DOUBLE_PRECISION,
     & pidfrom12,basetag+9,MPI_COMM_WORLD,irecv12,ierr)


      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(isend07,istatus,ierr)
      call MPI_WAIT(isend08,istatus,ierr)
      call MPI_WAIT(isend09,istatus,ierr)
      call MPI_WAIT(isend10,istatus,ierr)
      call MPI_WAIT(isend11,istatus,ierr)
      call MPI_WAIT(isend12,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)
      call MPI_WAIT(irecv07,istatus,ierr)
      call MPI_WAIT(irecv08,istatus,ierr)
      call MPI_WAIT(irecv09,istatus,ierr)
      call MPI_WAIT(irecv10,istatus,ierr)
      call MPI_WAIT(irecv11,istatus,ierr)
      call MPI_WAIT(irecv12,istatus,ierr)


      ! Corner send/recv
      pidto01=get_pid(pidx-1,pidy-1,pidz-1) !a_lll_s 
      pidto02=get_pid(pidx+1,pidy-1,pidz-1) !a_ull_s 
      pidto03=get_pid(pidx-1,pidy+1,pidz-1) !a_lul_s 
      pidto04=get_pid(pidx+1,pidy+1,pidz-1) !a_uul_s 
      pidto05=get_pid(pidx-1,pidy-1,pidz+1) !a_llu_s 
      pidto06=get_pid(pidx+1,pidy-1,pidz+1) !a_ulu_s 
      pidto07=get_pid(pidx-1,pidy+1,pidz+1) !a_luu_s 
      pidto08=get_pid(pidx+1,pidy+1,pidz+1) !a_uuu_s 
      call MPI_ISEND(a_lll_s(1,1,1),size(a_lll_s),MPI_DOUBLE_PRECISION,
     & pidto01,basetag+1,MPI_COMM_WORLD,isend01,ierr)
      call MPI_ISEND(a_ull_s(1,1,1),size(a_ull_s),MPI_DOUBLE_PRECISION,
     & pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)
      call MPI_ISEND(a_lul_s(1,1,1),size(a_lul_s),MPI_DOUBLE_PRECISION,
     & pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
      call MPI_ISEND(a_uul_s(1,1,1),size(a_uul_s),MPI_DOUBLE_PRECISION,
     & pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
      call MPI_ISEND(a_llu_s(1,1,1),size(a_llu_s),MPI_DOUBLE_PRECISION,
     & pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
      call MPI_ISEND(a_ulu_s(1,1,1),size(a_ulu_s),MPI_DOUBLE_PRECISION,
     & pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
      call MPI_ISEND(a_luu_s(1,1,1),size(a_luu_s),MPI_DOUBLE_PRECISION,
     & pidto07,basetag+7,MPI_COMM_WORLD,isend07,ierr)
      call MPI_ISEND(a_uuu_s(1,1,1),size(a_uuu_s),MPI_DOUBLE_PRECISION,
     & pidto08,basetag+8,MPI_COMM_WORLD,isend08,ierr)

      pidfrom01=get_pid(pidx-1,pidy-1,pidz-1) !a_lll_r 
      pidfrom02=get_pid(pidx+1,pidy-1,pidz-1) !a_ull_r 
      pidfrom03=get_pid(pidx-1,pidy+1,pidz-1) !a_lul_r 
      pidfrom04=get_pid(pidx+1,pidy+1,pidz-1) !a_uul_r 
      pidfrom05=get_pid(pidx-1,pidy-1,pidz+1) !a_llu_r 
      pidfrom06=get_pid(pidx+1,pidy-1,pidz+1) !a_ulu_r 
      pidfrom07=get_pid(pidx-1,pidy+1,pidz+1) !a_luu_r 
      pidfrom08=get_pid(pidx+1,pidy+1,pidz+1) !a_uuu_r 
      call MPI_IRECV(a_lll_r(1,1,1),size(a_lll_r),MPI_DOUBLE_PRECISION,
     & pidfrom01,basetag+8,MPI_COMM_WORLD,irecv01,ierr)
      call MPI_IRECV(a_ull_r(1,1,1),size(a_ull_r),MPI_DOUBLE_PRECISION,
     & pidfrom02,basetag+7,MPI_COMM_WORLD,irecv02,ierr)
      call MPI_IRECV(a_lul_r(1,1,1),size(a_lul_r),MPI_DOUBLE_PRECISION,
     & pidfrom03,basetag+6,MPI_COMM_WORLD,irecv03,ierr)
      call MPI_IRECV(a_uul_r(1,1,1),size(a_uul_r),MPI_DOUBLE_PRECISION,
     & pidfrom04,basetag+5,MPI_COMM_WORLD,irecv04,ierr)
      call MPI_IRECV(a_llu_r(1,1,1),size(a_llu_r),MPI_DOUBLE_PRECISION,
     & pidfrom05,basetag+4,MPI_COMM_WORLD,irecv05,ierr)
      call MPI_IRECV(a_ulu_r(1,1,1),size(a_ulu_r),MPI_DOUBLE_PRECISION,
     & pidfrom06,basetag+3,MPI_COMM_WORLD,irecv06,ierr)
      call MPI_IRECV(a_luu_r(1,1,1),size(a_luu_r),MPI_DOUBLE_PRECISION,
     & pidfrom07,basetag+2,MPI_COMM_WORLD,irecv07,ierr)
      call MPI_IRECV(a_uuu_r(1,1,1),size(a_uuu_r),MPI_DOUBLE_PRECISION,
     & pidfrom08,basetag+1,MPI_COMM_WORLD,irecv08,ierr)

      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(isend07,istatus,ierr)
      call MPI_WAIT(isend08,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)
      call MPI_WAIT(irecv07,istatus,ierr)
      call MPI_WAIT(irecv08,istatus,ierr)

      !Update the ghost layers with the data from neighboring processes
      !update planes
      call update_plane_r8(a, 'il', a_il_r, 1,jml,1,kml,ng,'regular')
      call update_plane_r8(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'regular')
      call update_plane_r8(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'regular')
      call update_plane_r8(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'regular')
      call update_plane_r8(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'regular')
      call update_plane_r8(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'regular')

      !update edges
      call update_edge_r8(a, 'iljl', a_iljl_r, 1,kml,ng,'regular')
      call update_edge_r8(a, 'ilju', a_ilju_r, 1,kml,ng,'regular')
      call update_edge_r8(a, 'iujl', a_iujl_r, 1,kml,ng,'regular')
      call update_edge_r8(a, 'iuju', a_iuju_r, 1,kml,ng,'regular')
      call update_edge_r8(a, 'jlkl', a_jlkl_r, 1,iml,ng,'regular')
      call update_edge_r8(a, 'jlku', a_jlku_r, 1,iml,ng,'regular')
      call update_edge_r8(a, 'jukl', a_jukl_r, 1,iml,ng,'regular')
      call update_edge_r8(a, 'juku', a_juku_r, 1,iml,ng,'regular')
      call update_edge_r8(a, 'klil', a_klil_r, 1,jml,ng,'regular')
      call update_edge_r8(a, 'kliu', a_kliu_r, 1,jml,ng,'regular')
      call update_edge_r8(a, 'kuil', a_kuil_r, 1,jml,ng,'regular')
      call update_edge_r8(a, 'kuiu', a_kuiu_r, 1,jml,ng,'regular')

      !corners
      call update_corner_r8(a, 'lll', a_lll_r, ng,'regular')
      call update_corner_r8(a, 'ull', a_ull_r, ng,'regular')
      call update_corner_r8(a, 'lul', a_lul_r, ng,'regular')
      call update_corner_r8(a, 'uul', a_uul_r, ng,'regular')
      call update_corner_r8(a, 'llu', a_llu_r, ng,'regular')
      call update_corner_r8(a, 'ulu', a_ulu_r, ng,'regular')
      call update_corner_r8(a, 'luu', a_luu_r, ng,'regular')
      call update_corner_r8(a, 'uuu', a_uuu_r, ng,'regular')

      ! Actual boundary condition is applied here.
      ! BC mode names(=bctype)
      ! periodic: periodic BC along all directions(default)
      ! adiabatic: adiabatic BC for all boundaries
      ! adiabatic_x: adiabatic BC along x direction (for il & iu), but periodic for the rest boundaries
      ! adiabatic_y: adiabatic BC along y direction (for jl & ju), but periodic for the rest boundaries
      ! adiabatic_z: adiabatic BC along z direction (for kl & ku), but periodic for the rest boundaries
      ! adiabatic_xy: adiabatic BC along x and y directions (il,iu,jl,ju), but periodic for z direction (kl,ku)
      ! adiabatic_yz: adiabatic BC along y and z directions (jl,ju,kl,ku), but periodic for x direction (il,iu)
      ! adiabatic_zx: adiabatic BC along z and x directions (kl,ku,il,iu), but periodic for y direction (jl,ju)
      ! Further BC modes will be coming up...
      select case(trim(bctype))
      case('adiabatic')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_r8(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_r8(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_r8(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_r8(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_r8(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_r8(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      ! Global iljl edge
      if(pidx.eq.0.and.pidy.eq.0)then
      call extract_edge_r8(a_iljl_s,'iljl', a, 1,kml,ng)
      call update_edge_r8(a, 'iljl', a_iljl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0)then
      ! Global ilju edge
      if(pidx.eq.0.and.pidy.eq.NPY-1)then
      call extract_edge_r8(a_ilju_s,'ilju', a, 1,kml,ng)
      call update_edge_r8(a, 'ilju', a_ilju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1)then 
      ! Global iujl edge
      if(pidx.eq.NPX-1.and.pidy.eq.0)then
      call extract_edge_r8(a_iujl_s,'iujl', a, 1,kml,ng)
      call update_edge_r8(a, 'iujl', a_iujl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0)then
      ! Global iuju edge
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      call extract_edge_r8(a_iuju_s,'iuju', a, 1,kml,ng)
      call update_edge_r8(a, 'iuju', a_iuju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      ! Global jlkl edge
      if(pidy.eq.0.and.pidz.eq.0)then
      call extract_edge_r8(a_jlkl_s,'jlkl', a, 1,iml,ng)
      call update_edge_r8(a, 'jlkl', a_jlkl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.0)then
      ! Global jlku edge
      if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_edge_r8(a_jlku_s,'jlku', a, 1,iml,ng)
      call update_edge_r8(a, 'jlku', a_jlku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global jukl edge
      if(pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_edge_r8(a_jukl_s,'jukl', a, 1,iml,ng)
      call update_edge_r8(a, 'jukl', a_jukl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global juku edge
      if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_edge_r8(a_juku_s,'juku', a, 1,iml,ng)
      call update_edge_r8(a, 'juku', a_juku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      ! Global klil edge
      if(pidz.eq.0.and.pidx.eq.0)then
      call extract_edge_r8(a_klil_s,'klil', a, 1,jml,ng)
      call update_edge_r8(a, 'klil', a_klil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.0)then
      ! Global kliu edge
      if(pidz.eq.0.and.pidx.eq.NPX-1)then
      call extract_edge_r8(a_kliu_s,'kliu', a, 1,jml,ng)
      call update_edge_r8(a, 'kliu', a_kliu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.NPX-1)then
      ! Global kuil edge
      if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      call extract_edge_r8(a_kuil_s,'kuil', a, 1,jml,ng)
      call update_edge_r8(a, 'kuil', a_kuil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      ! Global kuiu edge
      if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      call extract_edge_r8(a_kuiu_s,'kuiu', a, 1,jml,ng)
      call update_edge_r8(a, 'kuiu', a_kuiu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then

      ! Global lll corner
      if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.0)then
      call extract_corner_r8(a_lll_s, 'lll', a, ng)
      call update_corner_r8(a, 'lll', a_lll_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.0)then
      ! Global ull corner
      if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.0)then
      call extract_corner_r8(a_ull_s, 'ull', a, ng)
      call update_corner_r8(a, 'ull', a_ull_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.0)then
      ! Global lul corner
      if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_corner_r8(a_lul_s, 'lul', a, ng)
      call update_corner_r8(a, 'lul', a_lul_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global uul corner
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_corner_r8(a_uul_s, 'uul', a, ng)
      call update_corner_r8(a, 'uul', a_uul_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global llu corner
      if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_corner_r8(a_llu_s, 'llu', a, ng)
      call update_corner_r8(a, 'llu', a_llu_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global ulu corner
      if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_corner_r8(a_ulu_s, 'ulu', a, ng)
      call update_corner_r8(a, 'ulu', a_ulu_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global luu corner
      if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_corner_r8(a_luu_s, 'luu', a, ng)
      call update_corner_r8(a, 'luu', a_luu_s, ng, 'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      ! Global uuu corner
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_corner_r8(a_uuu_s, 'uuu', a, ng)
      call update_corner_r8(a, 'uuu', a_uuu_s, ng, 'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then

      case('adiabatic_x')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_r8(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_r8(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then

      case('adiabatic_y')
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_r8(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_r8(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then

      case('adiabatic_z')
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_r8(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_r8(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      case('adiabatic_xy')
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_r8(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_r8(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_r8(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_r8(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then

      ! Global iljl edge
      if(pidx.eq.0.and.pidy.eq.0)then
      call extract_edge_r8(a_iljl_s,'iljl', a, 1,kml,ng)
      call update_edge_r8(a, 'iljl', a_iljl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.0)then
      ! Global ilju edge
      if(pidx.eq.0.and.pidy.eq.NPY-1)then
      call extract_edge_r8(a_ilju_s,'ilju', a, 1,kml,ng)
      call update_edge_r8(a, 'ilju', a_ilju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.0.and.pidy.eq.NPY-1)then 
      ! Global iujl edge
      if(pidx.eq.NPX-1.and.pidy.eq.0)then
      call extract_edge_r8(a_iujl_s,'iujl', a, 1,kml,ng)
      call update_edge_r8(a, 'iujl', a_iujl_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.0)then
      ! Global iuju edge
      if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then
      call extract_edge_r8(a_iuju_s,'iuju', a, 1,kml,ng)
      call update_edge_r8(a, 'iuju', a_iuju_s, 1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1.and.pidy.eq.NPY-1)then

      case('adiabatic_yz')
      ! Global jl plane
      if(pidy.eq.0)then
      call extract_plane_r8(a_jl_r,'jl',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'jl', a_jl_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.0)then
      ! Global ju plane
      if(pidy.eq.NPY-1)then
      call extract_plane_r8(a_ju_r,'ju',a, 1,kml,1,iml,ng)
      call update_plane_r8(a, 'ju', a_ju_r, 1,kml,1,iml,ng,'mirror') 
      endif !if(pidy.eq.NPY-1)then
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_r8(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_r8(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then

      ! Global jlkl edge
      if(pidy.eq.0.and.pidz.eq.0)then
      call extract_edge_r8(a_jlkl_s,'jlkl', a, 1,iml,ng)
      call update_edge_r8(a, 'jlkl', a_jlkl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.0)then
      ! Global jlku edge
      if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      call extract_edge_r8(a_jlku_s,'jlku', a, 1,iml,ng)
      call update_edge_r8(a, 'jlku', a_jlku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.0.and.pidz.eq.NPZ-1)then
      ! Global jukl edge
      if(pidy.eq.NPY-1.and.pidz.eq.0)then
      call extract_edge_r8(a_jukl_s,'jukl', a, 1,iml,ng)
      call update_edge_r8(a, 'jukl', a_jukl_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.0)then
      ! Global juku edge
      if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then
      call extract_edge_r8(a_juku_s,'juku', a, 1,iml,ng)
      call update_edge_r8(a, 'juku', a_juku_s, 1,iml,ng,'mirror')
      endif !if(pidy.eq.NPY-1.and.pidz.eq.NPZ-1)then

      case('adiabatic_zx')
      ! Global kl plane
      if(pidz.eq.0)then
      call extract_plane_r8(a_kl_r,'kl',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'kl', a_kl_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.0)then
      ! Global ku plane
      if(pidz.eq.NPZ-1)then
      call extract_plane_r8(a_ku_r,'ku',a, 1,iml,1,jml,ng)
      call update_plane_r8(a, 'ku', a_ku_r, 1,iml,1,jml,ng,'mirror') 
      endif !if(pidz.eq.NPZ-1)then
      ! Global il plane
      if(pidx.eq.0)then
      call extract_plane_r8(a_il_r,'il',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'il', a_il_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.0)then
      ! Global iu plane
      if(pidx.eq.NPX-1)then
      call extract_plane_r8(a_iu_r,'iu',a, 1,jml,1,kml,ng)
      call update_plane_r8(a, 'iu', a_iu_r, 1,jml,1,kml,ng,'mirror')
      endif !if(pidx.eq.NPX-1)then

      ! Global klil edge
      if(pidz.eq.0.and.pidx.eq.0)then
      call extract_edge_r8(a_klil_s,'klil', a, 1,jml,ng)
      call update_edge_r8(a, 'klil', a_klil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.0)then
      ! Global kliu edge
      if(pidz.eq.0.and.pidx.eq.NPX-1)then
      call extract_edge_r8(a_kliu_s,'kliu', a, 1,jml,ng)
      call update_edge_r8(a, 'kliu', a_kliu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.0.and.pidx.eq.NPX-1)then
      ! Global kuil edge
      if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      call extract_edge_r8(a_kuil_s,'kuil', a, 1,jml,ng)
      call update_edge_r8(a, 'kuil', a_kuil_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.0)then
      ! Global kuiu edge
      if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      call extract_edge_r8(a_kuiu_s,'kuiu', a, 1,jml,ng)
      call update_edge_r8(a, 'kuiu', a_kuiu_s, 1,jml,ng,'mirror')
      endif !if(pidz.eq.NPZ-1.and.pidx.eq.NPX-1)then
      end select

      end subroutine boundary_condition_r8_par



      subroutine print_ghost_layers_r8(a,ng)
      use input
      implicit none

      real(8),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      real(8) :: a_il(1-ng:jml+ng,1-ng:kml+ng)
      real(8) :: a_iu(1-ng:jml+ng,1-ng:kml+ng)
      real(8) :: a_jl(1-ng:kml+ng,1-ng:iml+ng)
      real(8) :: a_ju(1-ng:kml+ng,1-ng:iml+ng)
      real(8) :: a_kl(1-ng:iml+ng,1-ng:jml+ng)
      real(8) :: a_ku(1-ng:iml+ng,1-ng:jml+ng)

      integer(4) :: i,j,k

      write(*,*) "Subroutine print_ghost_layers() called"

      !il and iu
      write(*,*) 'il and iu'
      do k=1-ng,kml+ng
      do j=1-ng,jml+ng
            a_il(j,k)=a(1-ng,j,k)
            a_iu(j,k)=a(iml+ng,j,k)
            write(*,*)j,k,a_il(j,k),a_iu(j,k)
      enddo
      enddo

      !jl and ju
      write(*,*) 'jl and ju'
      do i=1-ng,iml+ng
      do k=1-ng,kml+ng
            a_jl(k,i)=a(i,1-ng,k)
            a_ju(k,i)=a(i,jml+ng,k)
            write(*,*)k,i,a_jl(k,i),a_ju(k,i)
      enddo
      enddo

      !kl and ku
      write(*,*) 'kl and ku'
      do j=1-ng,jml+ng
      do i=1-ng,iml+ng
            a_kl(i,j)=a(i,j,1-ng)
            a_ku(i,j)=a(i,j,kml+ng)
            write(*,*)i,j,a_kl(i,j),a_ku(i,j)
      enddo
      enddo

      end subroutine print_ghost_layers_r8


      subroutine print_array_r8(a,ng,il,iu,jl,ju,kl,ku)
      use input
      implicit none
      real(8),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: il,iu,jl,ju,kl,ku

      integer(4) :: i,j,k,n1

      write(*,*)'print_array() called from rank = ',myrank

      do k=kl,ku
            write(*,*)'k=',k,", i: horizontal, j: vertical"
            do j=jl,ju
                  write(*,*) (a(i,j,k), i=il,iu)
            enddo
      enddo

      end subroutine print_array_r8


      subroutine extract_plane_r8(aout, planetype, ain, il,iu,jl,ju,ng)
      use input
      implicit none

      real(8),intent(inout) :: aout(il-ng:iu+ng,jl-ng:ju+ng,ng)
      character(len=*),intent(in) :: planetype
      real(8),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: il,iu,jl,ju,ng

      integer(4) :: i,j,k,n1

      select case(trim(planetype))
      case('il')
            do n1=1,ng
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  aout(j,k,n1)=ain(n1,j,k)
            enddo
            enddo
            enddo
      case('iu')
            do n1=1,ng
            do k=1-ng,kml+ng
            do j=1-ng,jml+ng
                  aout(j,k,n1)=ain(iml-ng+n1,j,k)
            enddo
            enddo
            enddo
      case('jl')
            do n1=1,ng
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  aout(k,i,n1)=ain(i,n1,k)
            enddo
            enddo
            enddo
      case('ju')
            do n1=1,ng
            do i=1-ng,iml+ng
            do k=1-ng,kml+ng
                  aout(k,i,n1)=ain(i,jml-ng+n1,k)
            enddo
            enddo
            enddo
      case('kl')
            do n1=1,ng
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  aout(i,j,n1)=ain(i,j,n1)
            enddo
            enddo
            enddo
      case('ku')
            do n1=1,ng
            do j=1-ng,jml+ng
            do i=1-ng,iml+ng
                  aout(i,j,n1)=ain(i,j,kml-ng+n1)
            enddo
            enddo
            enddo
      end select
      end subroutine extract_plane_r8


      subroutine extract_edge_r8(aout, edgetype, ain, il,iu,ng)
      use input
      implicit none
      real(8),intent(inout) :: aout(il-ng:iu+ng,ng,ng)
      character(len=*),intent(in) :: edgetype
      real(8),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: il,iu,ng

      integer(4) :: ll,n1,n2

      select case(trim(edgetype))
      case('iljl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n1,n2,ll)
                  enddo
            enddo
            enddo
      case('ilju')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n1,jml-ng+n2,ll)
                  enddo
            enddo
            enddo
      case('iujl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n1,n2,ll)
                  enddo
            enddo
            enddo
      case('iuju')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n1,jml-ng+n2,ll)
                  enddo
            enddo
            enddo
      case('jlkl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,n1,n2)
                  enddo
            enddo
            enddo
      case('jlku')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,n1,kml-ng+n2)
                  enddo
            enddo
            enddo
      case('jukl')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,jml-ng+n1,n2)
                  enddo
            enddo
            enddo
      case('juku')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(ll,jml-ng+n1,kml-ng+n2)
                  enddo
            enddo
            enddo
      case('klil')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n2,ll,n1)
                  enddo
            enddo
            enddo
      case('kliu')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n2,ll,n1)
                  enddo
            enddo
            enddo
      case('kuil')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(n2,ll,kml-ng+n1)
                  enddo
            enddo
            enddo
      case('kuiu')
            do n2=1,ng
            do n1=1,ng
                  do ll=il-ng,iu+ng
                        aout(ll,n1,n2)=ain(iml-ng+n2,ll,kml-ng+n1)
                  enddo
            enddo
            enddo
      end select

      end subroutine extract_edge_r8


      subroutine extract_corner_r8(aout, cornertype, ain, ng)
      use input
      implicit none
      real(8),intent(inout) :: aout(ng,ng,ng)
      character(len=*),intent(in) :: cornertype
      real(8),intent(inout)::ain(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ng

      integer(4) :: n1,n2,n3

      select case(trim(cornertype))
      case('lll')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,n2,n3)
            enddo
            enddo
            enddo
      case('ull')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,n2,n3)
            enddo
            enddo
            enddo
      case('lul')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,jml-ng+n2,n3)
            enddo
            enddo
            enddo
      case('uul')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,jml-ng+n2,n3)
            enddo
            enddo
            enddo
      case('llu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('ulu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('luu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(n1,jml-ng+n2,kml-ng+n3)
            enddo
            enddo
            enddo
      case('uuu')
            do n3=1,ng
            do n2=1,ng
            do n1=1,ng
                  aout(n1,n2,n3)=ain(iml-ng+n1,jml-ng+n2,kml-ng+n3)
            enddo
            enddo
            enddo
      end select

      end subroutine extract_corner_r8

      subroutine update_plane_r8(aout, planetype, ain 
     &                                     , il,iu,jl,ju,ng,copymode)
      use input
      implicit none
      real(8),intent(inout) :: 
     &                   aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: planetype  !plane name
      real(8),intent(inout) :: ain(il-ng:iu+ng,jl-ng:ju+ng,ng)
      integer(4),intent(in) :: il,iu,jl,ju,ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: i,j,k,n1

      select case(trim(copymode))
      case('regular')
            select case(trim(planetype))
            case('il')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(1-ng+n1,j,k)=ain(j,k,1+n1)
                  enddo
                  enddo
                  enddo
            case('iu')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(iml+1+n1,j,k)=ain(j,k,1+n1)
                  enddo
                  enddo
                  enddo
            case('jl')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,1-ng+n1,k)=ain(k,i,1+n1)
                  enddo
                  enddo
                  enddo
            case('ju')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,jml+1+n1,k)=ain(k,i,1+n1)
                  enddo
                  enddo
                  enddo
            case('kl')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,1-ng+n1)=ain(i,j,1+n1)
                  enddo
                  enddo
                  enddo
            case('ku')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,kml+1+n1)=ain(i,j,1+n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(planetype))
      case('mirror')
            select case(trim(planetype))
            case('il')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(1-ng+n1,j,k)=ain(j,k,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iu')
                  do n1=0,ng-1
                  do k=1-ng,kml+ng
                  do j=1-ng,jml+ng
                        aout(iml+1+n1,j,k)=ain(j,k,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jl')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,1-ng+n1,k)=ain(k,i,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ju')
                  do n1=0,ng-1
                  do i=1-ng,iml+ng
                  do k=1-ng,kml+ng
                        aout(i,jml+1+n1,k)=ain(k,i,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kl')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,1-ng+n1)=ain(i,j,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ku')
                  do n1=0,ng-1
                  do j=1-ng,jml+ng
                  do i=1-ng,iml+ng
                        aout(i,j,kml+1+n1)=ain(i,j,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(planetype))
      end select !select case(trim(copymode))

      end subroutine update_plane_r8


      subroutine update_edge_r8(aout, edgetype, ain, il,iu,ng,copymode)
      use input
      implicit none

      real(8),intent(inout):: 
     &                   aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: edgetype  !edge name
      real(8),intent(inout) :: ain(il-ng:iu+ng,ng,ng)
      integer(4),intent(in) :: il,iu,ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: ll,n1,n2

      select case(trim(copymode))
      case('regular')
            select case(trim(edgetype))
            case('iljl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,1-ng+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('ilju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,jml+1+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('iujl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,1-ng+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('iuju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,jml+1+n2,ll)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jlkl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,1-ng+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jlku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,kml+1+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('jukl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,1-ng+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('juku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,kml+1+n2)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('klil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,1-ng+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kliu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,1-ng+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kuil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,kml+1+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            case('kuiu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,kml+1+n1)=ain(ll,1+n1,1+n2)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(edgetype))
      case('mirror')
            select case(trim(edgetype))
            case('iljl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,1-ng+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ilju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n1,jml+1+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iujl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,1-ng+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('iuju')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n1,jml+1+n2,ll)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jlkl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,1-ng+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jlku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,1-ng+n1,kml+1+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('jukl')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,1-ng+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('juku')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(ll,jml+1+n1,kml+1+n2)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('klil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,1-ng+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kliu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,1-ng+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kuil')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(1-ng+n2,ll,kml+1+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('kuiu')
                  do n2=0,ng-1
                  do n1=0,ng-1
                  do ll=il-ng,iu+ng
                        aout(iml+1+n2,ll,kml+1+n1)=ain(ll,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(edgetype))
      end select !select case(trim(copymode))

      end subroutine update_edge_r8


      subroutine update_corner_r8(aout, cornertype, ain, ng,copymode)
      use input
      implicit none
      real(8),intent(inout) :: 
     &                  aout(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      real(8),intent(inout) :: ain(ng,ng,ng)
      character(len=*),intent(in) :: cornertype !corner name
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: copymode  !data copy mode: regular or mirror

      integer(4) :: n1,n2,n3

      select case(trim(copymode))
      case('regular')
            select case(trim(cornertype))
            case('lll')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('ull')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('lul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('uul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,1-ng+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('llu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('ulu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('luu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            case('uuu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,kml+1+n3)=ain(1+n1,1+n2,1+n3)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(cornertype))
      case('mirror')
            select case(trim(cornertype))
            case('lll')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ull')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('lul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('uul')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,1-ng+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('llu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,1-ng+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('ulu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,1-ng+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('luu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(1-ng+n1,jml+1+n2,kml+1+n3)=ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            case('uuu')
                  do n3=0,ng-1
                  do n2=0,ng-1
                  do n1=0,ng-1
                  aout(iml+1+n1,jml+1+n2,kml+1+n3)= 
     &                                     ain(ng-n3,ng-n2,ng-n1)
                  enddo
                  enddo
                  enddo
            end select !select case(trim(cornertype))
      end select !select case(trim(copymode))

      end subroutine update_corner_r8


      end module boundary_conditions