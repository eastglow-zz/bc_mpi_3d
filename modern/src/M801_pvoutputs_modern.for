      module pvoutputs

        use mpi_f08
        implicit none
        character(len=300) :: PVOPATH='./outputs/'
        integer(4),parameter :: pvdheaderlines = 3
        integer(4),parameter :: pvdfooterlines = 2

        integer(4) :: pv_dimen
        integer(4) :: pv_dimx, pv_dimy, pv_dimz

        real(8) :: pv_x0_o, pv_y0_o, pv_z0_o  ! origin shift
        real(8) :: pv_dx_o, pv_dy_o, pv_dz_o  ! grid sizes
        integer(4) :: pv_nprocs_o ! # of processors

        integer(4) :: pv_img, pv_jmg, pv_kmg
        integer(4) :: pv_iml, pv_jml, pv_kml

        integer(4) :: pv_gpid
        integer(4),allocatable :: pv_pidx(:), pv_pidy(:), pv_pidz(:)
        integer(4),allocatable :: pv_il(:), pv_iu(:) 
        integer(4),allocatable :: pv_jl(:), pv_ju(:) 
        integer(4),allocatable :: pv_kl(:), pv_ku(:) 

        type(MPI_Comm) :: pv_comm

        logical :: is_initialized = .false.


        !---------------------------------------------------------------
        ! !Included subroutine list
        !  init_pvoutputs 
        !  finalize_pvoutputs
        !  pv_get_dimension
        !  pv_get_file_lines
        !  output_i4_raw
        !  output_pvd
        !  output_i4_pvtr
        !  output_r8_raw
        !  output_r8_pvtr
        !  checkpoint_i4_pvtr
        !  checkpoint_r8_pvtr
      contains
      subroutine init_pvoutputs(imglo,jmglo,kmglo,imloc, jmloc, kmloc, 
     &                                           il, iu, jl, ju, kl, ku,
     &                                           xo, yo, zo, dx, dy, dz, 
     &                           numprocs, pidx, pidy, pidz, gpid, comm)
      
      ! This routine assumes rectilinear coordinate system
      ! 3 directional parallel partitioning is also assumed. This can be achieved via MPI_COMM_CART() or any use defined subroutines.
      ! imglo: The number of grid along x-axis, before parallel partitioning (global # of grid, x)
      ! jmglo: The number of grid along y-axis, before parallel partitioning (global # of grid, y)
      ! kmglo: The nubmer of grid along z-axis, before parallel partitioning (global # of grid, z)
      ! imloc: The number of grid along x-axis of in a local process (local # of grid, x)
      ! jmloc: The number of grid along y-axis of in a local process (local # of grid, y)
      ! kmloc: The number of grid along z-axis of in a local process (local # of grid, z)
      ! il: Lower bound of global grid index along x-axis 
      ! iu: Upper bound of global grid index along x-axis
      ! jl: Lower bound of global grid index along y-axis
      ! ju: Upper bound of global grid index along y-axis
      ! kl: Lower bound of global grid index along z-axis 
      ! ku: Upper bound of global grid index along z-axis
      ! xo: X coordinate of the origin of the grid in meter unit
      ! yo: Y coordinate of the origin of the grid in meter unit
      ! zo: Z coordinate of the origin of the grid in meter unit
      ! dx: Grid size along x-axis in meter unit 
      ! dy: Grid size along y-axis in meter unit 
      ! dz: Grid size along z-axis in meter unit 
      ! numprocs: Total # of MPI process 
      ! pidx: Process coordinate in 3 directional parallel partitioning, along x-axis 
      ! pidy: Process coordinate in 3 directional parallel partitioning, along y-axis 
      ! pidz: Process coordinate in 3 directional parallel partitioning, along z-axis 
      ! gpid: Global process coordinate which may identical to MPI rank of a process
      ! comm: MPI communicator handle

      implicit none
      integer(4),intent(in) :: imglo, jmglo, kmglo
      integer(4),intent(in) :: imloc, jmloc, kmloc 
      integer(4),intent(in) :: il,iu,jl,ju,kl,ku 
      real(8),intent(in) :: xo, yo, zo
      real(8),intent(in) :: dx, dy, dz
      integer(4), intent(in) :: numprocs
      integer(4), intent(in) :: pidx, pidy, pidz, gpid 
      type(MPI_Comm),intent(in) :: comm

      integer(4) :: ierr 
      integer(4) :: nn 
      integer(4) :: tmp 

      allocate(pv_pidx(0:numprocs-1))
      allocate(pv_pidy(0:numprocs-1))
      allocate(pv_pidz(0:numprocs-1))
      allocate(pv_il(0:numprocs-1))
      allocate(pv_iu(0:numprocs-1))
      allocate(pv_jl(0:numprocs-1))
      allocate(pv_ju(0:numprocs-1))
      allocate(pv_kl(0:numprocs-1))
      allocate(pv_ku(0:numprocs-1))

      call mpi_allgather(
     &        pidx, 1, MPI_INTEGER, pv_pidx, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &        pidy, 1, MPI_INTEGER, pv_pidy, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &        pidz, 1, MPI_INTEGER, pv_pidz, 1, MPI_INTEGER, comm, ierr)

      call mpi_allgather(
     &            il, 1, MPI_INTEGER, pv_il, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &            iu, 1, MPI_INTEGER, pv_iu, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &            jl, 1, MPI_INTEGER, pv_jl, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &            ju, 1, MPI_INTEGER, pv_ju, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &            kl, 1, MPI_INTEGER, pv_kl, 1, MPI_INTEGER, comm, ierr)
      call mpi_allgather(
     &            ku, 1, MPI_INTEGER, pv_ku, 1, MPI_INTEGER, comm, ierr)

      call pv_get_dimension(imglo, jmglo, kmglo)

      pv_img = imglo 
      pv_jmg = jmglo 
      pv_kmg = kmglo 
      pv_iml = imloc 
      pv_jml = jmloc
      pv_kml = kmloc 
      pv_x0_o = xo
      pv_y0_o = yo
      pv_z0_o = zo
      pv_dx_o = dx
      pv_dy_o = dy
      pv_dz_o = dz
      pv_nprocs_o = numprocs
      pv_gpid = gpid

      pv_comm = comm

      if(gpid .eq. 0) then 
        write(*,*)'=== Initialization of pvoutputs ==='
        write(*,*)'===    pvoutputs parameters     ==='
        write(*,*)'pv_img = ', pv_img 
        write(*,*)'pv_jmg = ', pv_jmg
        write(*,*)'pv_kmg = ', pv_kmg 
        write(*,*)'pv_iml = ', pv_iml 
        write(*,*)'pv_jml = ', pv_jml
        write(*,*)'pv_kml = ', pv_kml 
        write(*,*)'pv_x0_o = ', pv_x0_o 
        write(*,*)'pv_y0_o = ', pv_y0_o 
        write(*,*)'pv_z0_o = ', pv_z0_o 
        write(*,*)'pv_dx_o = ', pv_dx_o
        write(*,*)'pv_dy_o = ', pv_dy_o
        write(*,*)'pv_dz_o = ', pv_dz_o
        write(*,*)'pv_nprocs_o = ', pv_nprocs_o

        write(*,*)'  pid, pv_pidx, pv_pidy, pv_pidz'
        do nn=0,numprocs-1
          write(*,*) nn, pv_pidx(nn), pv_pidy(nn), pv_pidz(nn)
        enddo

        write(*,*)'  pid, pv_il, pv_iu, pv_jl, pv_ju, pv_kl, pv_ku'
        do nn=0,numprocs-1
          write(*,*)nn,
     &       pv_il(nn),pv_iu(nn),pv_jl(nn),pv_ju(nn),pv_kl(nn),pv_ku(nn)
        enddo
        write(*,*)'==================================='
      endif

      return
      end subroutine init_pvoutputs

      !-----------------------------------------------------------------

      subroutine finalize_pvoutputs()
      implicit none 

      deallocate(pv_pidx)
      deallocate(pv_pidy)
      deallocate(pv_pidz)
      deallocate(pv_il)
      deallocate(pv_iu)
      deallocate(pv_jl)
      deallocate(pv_ju)
      deallocate(pv_kl)
      deallocate(pv_ku)

      return
      end subroutine finalize_pvoutputs

      !-----------------------------------------------------------------

      subroutine pv_get_dimension(im, jm, km)
      implicit none 

      integer(4),intent(in) :: im, jm, km

      pv_dimx = 1 
      pv_dimy = 1
      pv_dimz = 1
      if(im .le. 1) pv_dimx = 0
      if(jm .le. 1) pv_dimy = 0
      if(km .le. 1) pv_dimz = 0 
      pv_dimen = pv_dimx + pv_dimy + pv_dimz
      
      if(pv_dimen .eq. 0) then
        write(*,*)'Error, pv_get_dimension() in M801_pvoutputs.for'
        write(*,*)'  No dimension specified'
        write(*,*)'  Aborting...'
        call abort()
      end if

      return 
      end subroutine pv_get_dimension

      !-----------------------------------------------------------------

      subroutine pv_get_file_lines(nlines, fname)
      implicit none
      integer(4),intent(inout) :: nlines
      character(len=*),intent(in) :: fname

      integer(4):: io

      nlines = 0
      open(100, FILE=fname, status='old',action='read')
      DO
        READ(100,*,iostat=io)
        IF (io/=0) exit
        nlines = nlines + 1
      END DO
      close(100)
      return
      end subroutine pv_get_file_lines

      !-----------------------------------------------------------------

      subroutine pv_get_pidxyz(pidx, pidy, pidz, gpid)
      implicit none
      integer(4),intent(inout) :: pidx, pidy, pidz 
      integer(4),intent(in) :: gpid

      pidx = pv_pidx(gpid)
      pidy = pv_pidy(gpid)
      pidz = pv_pidz(gpid)
      
      return 
      end subroutine pv_get_pidxyz

      !-----------------------------------------------------------------

      subroutine pv_get_ijkrange(il,iu,jl,ju,kl,ku, gpid)
      implicit none 
      integer(4),intent(inout) :: il,iu,jl,ju,kl,ku
      integer(4),intent(in) :: gpid 

      il = pv_il(gpid) 
      iu = pv_iu(gpid) 
      jl = pv_jl(gpid)
      ju = pv_ju(gpid)
      kl = pv_kl(gpid)
      ku = pv_ku(gpid)

      return 
      end subroutine pv_get_ijkrange

      !-----------------------------------------------------------------

      subroutine output_i4_raw(a, ng, fnamebase, timestep, 
     & gpid, description)
      implicit none
      !!include 'mpif.h'
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                           1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                           1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep, gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k
      character(len=300) :: fnameout
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku

      call pv_get_pidxyz(pidx, pidy, pidz, gpid)
      call pv_get_ijkrange(il,iu,jl,ju,kl,ku, gpid)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      fnameout = trim(PVOPATH)//'/'//trim(fnamebase)//'_'//trim(ts2s)//
     & '_'//trim(pid2s)//'.raw'

      open(1000,file=trim(fnameout),status='unknown',action='write')
      ! Header output
      write(1000,*)description
      write(1000,*)
     & 'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
      write(1000,*)
     & gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

      ! Data output - fortran order
      do k=1,pv_kml
      do j=1,pv_jml
      do i=1,pv_iml
            write(1000,*) a(i,j,k)
      enddo
      enddo
      enddo

      close(1000)

      end subroutine output_i4_raw

      !-----------------------------------------------------------------

      subroutine output_pvd(fnamebase,subdir,timestep,time,description)
      implicit none
      character(len=*),intent(in) :: fnamebase
      character(len=*),intent(in) :: subdir
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: time
      character(len=*),intent(in) :: description

      character(len=300) :: fnamepvtr, fnamepvd
      character(len=8) :: ts2s
      character(len=20) :: time_str
      logical :: pvdexist
      integer(4) :: nlines, nbody, line
      write(ts2s,"(i8.8)") int(timestep)
      write(time_str,"(E20.15)") time
      fnamepvtr = trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'
      fnamepvd=trim(PVOPATH)//'/'//trim(fnamebase)//'.pvd'

      ! Check if the file exists or not
      inquire(file=trim(fnamepvd),exist=pvdexist)
      if (pvdexist .eqv. .false.)then  ! .pvd doesn't exist. Make a new one.

        open(unit=900,file=trim(fnamepvd),
     &      status='unknown', action='write')

        ! Write .pvd header
        write(900,'(A)')'<?xml version="1.0"?>'
        write(900,'(A)')'<VTKFile type="Collection" version="0.1" '//
     &    'byte_order="LittleEndian">'
        write(900,'(A)')'<Collection>'

        ! Write the pvd data
        write(900,'(A)')'<DataSet timestep="'//trim(time_str)//
     &    '" group="" part="0" file="'//trim(fnamepvtr)//'" name="'//
     &    trim(description)//'"/>'

        ! Write .pvd footer
        write(900,'(A)')'</Collection>'
        write(900,'(A)')'</VTKFile>'
        close(900)

      else ! .pvd already exists. Update the lines.
        ! Find the # of lines in the pvd file
        call pv_get_file_lines(nlines, trim(fnamepvd))
        nbody = nlines - pvdheaderlines - pvdfooterlines

        open(unit=900,file=trim(fnamepvd),
     &      status='old', action='readwrite')

        ! Skip the lines over prior data
        do line = 1, pvdheaderlines + nbody
          read(900,*)
        enddo

        ! Write the pvd data
        write(900,'(A)')'<DataSet timestep="'//trim(time_str)//
     &    '" group="" part="0" file="'//trim(fnamepvtr)//'" name="'//
     &    trim(description)//'"/>'

        ! Write .pvd footer
        write(900,'(A)')'</Collection>'
        write(900,'(A)')'</VTKFile>'
        close(900)

      endif ! if (pdvexist .eqv. .false.)then ... else

      return
      end subroutine output_pvd

      !-----------------------------------------------------------------

      subroutine output_i4_pvtr(a,ng,fnamebase,timestep,
     &                                         timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                           1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                           1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(fnamebase,subdir,timestep,timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)
       
      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_pvtr

      !-----------------------------------------------------------------

      subroutine output_r8_raw(a, ng, fnamebase, timestep, 
     & gpid, description)

      implicit none
      !include 'mpif.h'
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                        1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                        1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep, gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k
      character(len=300) :: fnameout, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: ierr

      call pv_get_pidxyz(pidx, pidy, pidz, gpid)
      call pv_get_ijkrange(il,iu,jl,ju,kl,ku, gpid)


      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      fnameout = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.raw'

      open(1000,file=trim(fnameout),status='unknown',action='write')
      ! Header output
      write(1000,*)description
      write(1000,*)
     & 'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
      write(1000,*)
     & gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

      ! Data output - fortran order
      do k=1,pv_kml
      do j=1,pv_jml
      do i=1,pv_iml
            write(1000,*) a(i,j,k)
      enddo
      enddo
      enddo

      close(1000)

      end subroutine output_r8_raw

      !-----------------------------------------------------------------

      subroutine output_r4_pvtr(a,ng,fnamebase,timestep,timeval,gpid,
     & description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      real(4),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                        1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                        1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr


      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(fnamebase,subdir,timestep,timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid) 
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(E21.15)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r4_pvtr

      !-----------------------------------------------------------------
      subroutine output_r8_pvtr(a,ng,fnamebase,timestep,
     &                                         timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                        1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                        1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(fnamebase,subdir,timestep,timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)
  

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(E21.15)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_pvtr

      !-----------------------------------------------------------------

      subroutine checkpoint_i4_pvtr(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                           1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                           1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr


      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/')
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/'//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd('checkpoint/'//fnamebase,subdir,timestep,
     &    timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_i4_pvtr

      !-----------------------------------------------------------------

      subroutine checkpoint_r4_pvtr(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      real(4),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                        1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                        1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr


      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/')
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/'//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd('checkpoint/'//fnamebase,subdir,timestep,
     &    timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)   
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(E21.15)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_r4_pvtr

      !-----------------------------------------------------------------

      subroutine checkpoint_r8_pvtr(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng*pv_dimx:pv_iml+ng*pv_dimx,
     &                        1-ng*pv_dimy:pv_jml+ng*pv_dimy,
     &                        1-ng*pv_dimz:pv_kml+ng*pv_dimz )
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr


      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/')
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/'//trim(subdir))
      endif
      call mpi_barrier(pv_comm, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd('checkpoint/'//fnamebase,subdir,timestep,
     &    timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,pv_nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call pv_get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call pv_get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,pv_nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call pv_get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call pv_get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,'(E21.15)') i*pv_dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,'(E21.15)') j*pv_dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,'(E21.15)') k*pv_dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(E21.15)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_r8_pvtr

      !-----------------------------------------------------------------

      end module pvoutputs  