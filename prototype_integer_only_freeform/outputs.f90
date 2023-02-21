module outputs
implicit none

contains
subroutine output_i4_raw(a, ng, fnamebase, timestep, gpid, description)
use input
use boundary_conditions
implicit none
include 'mpif.h'
integer(4),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
integer(4),intent(in) :: ng
character(len=*),intent(in) :: fnamebase
integer(4),intent(in) :: timestep, gpid
character(len=*),intent(in) :: description

integer(4) :: i,j,k,nn
character(len=300) :: fnameout
character(len=8) :: ts2s, pid2s
integer(4) :: pidx, pidy, pidz
integer(4) :: il, iu, jl, ju, kl, ku

call get_pidxyz(pidx, pidy, pidz, gpid)
call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

write(ts2s,"(i8.8)") int(timestep)
write(pid2s,"(i4.4)") int(gpid)
fnameout = trim(opath)//'/'//trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.raw'

open(1000,file=trim(fnameout),status='unknown',action='write')
! Header output
write(1000,*)description
write(1000,*)'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
write(1000,*)gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

! Data output - fortran order
do k=1,kml
do j=1,jml
do i=1,iml
      write(1000,*) a(i,j,k)
enddo
enddo
enddo

close(1000)

end subroutine output_i4_raw



subroutine output_i4_pvtr(a, ng, fnamebase, timestep, gpid, description)
! Output of pvtr requires at least one ghost layer (ng>=1)
use input
use boundary_conditions
implicit none
include 'mpif.h'
integer(4),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
integer(4),intent(in) :: ng
character(len=*),intent(in) :: fnamebase
integer(4),intent(in) :: timestep, gpid
character(len=*),intent(in) :: description

integer(4) :: i,j,k,nn
character(len=300) :: fnameg, fnamel
character(len=8) :: ts2s, pid2s
integer(4) :: pidx, pidy, pidz
integer(4) :: il, iu, jl, ju, kl, ku
integer(4) :: itmp, jtmp, ktmp
character(len=300) :: str_extent_g, str_extent_l

call get_pidxyz(pidx, pidy, pidz, gpid)
call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

write(ts2s,"(i8.8)") int(timestep)
write(pid2s,"(i4.4)") int(gpid)
write(str_extent_g,"(6(i4,1x))")0,im-1,0,jm-1,0,km-1
fnameg = trim(opath)//'/'//trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

! Output of pvtr file (meta file to combine each local data(vtr))
if(gpid.eq.0)then
      open(1001,FILE=fnameg,status='unknown',action='write')
      write(1001,'(A)')'<?xml version="1.0"?>'
      write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//'version="0.1" byte_order="LittleEndian">'
      write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//trim(adjustl(str_extent_g))//'">'

      write(1001,'(A)')'<PCoordinates>'
      write(1001,'(A)')'<DataArray Name="x" type="Float32" format="ascii">'
      write(1001,'(A)')'</DataArray>'
      write(1001,'(A)')'<DataArray Name="y" type="Float32" format="ascii">'
      write(1001,'(A)')'</DataArray>'
      write(1001,'(A)')'<DataArray Name="z" type="Float32" format="ascii">'
      write(1001,'(A)')'</DataArray>'
      write(1001,'(A)')'</PCoordinates>'

      write(1001,'(A)')'<PPointData Scalars="'//trim(fnamebase)//'">'
      write(1001,'(A)')'<DataArray type="Int32" Name="'//trim(fnamebase)//'" format="ascii">'
      write(1001,'(A)')'</DataArray>'
      write(1001,'(A)')'</PPointData>'

      do nn=0,nprocs-1
            write(pid2s,"(i4.4)") int(nn)
            fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
            call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
            call get_pidxyz(pidx, pidy, pidz, nn)

            itmp=il-1
            jtmp=jl-1
            ktmp=kl-1
            if(pidx.ge.1)itmp=il-2
            if(pidy.ge.1)jtmp=jl-2
            if(pidz.ge.1)ktmp=kl-2
            write(str_extent_l,"(6(i4,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1

            write(1001,'(A)')'<Piece Extent="'//trim(adjustl(str_extent_l))//'" Source="'//trim(fnamel)//'">'
            write(1001,'(A)')'</Piece>'
      enddo !do nn=0,nprocs-1

      write(1001,'(A)')'</PRectilinearGrid>'
      write(1001,'(A)')'</VTKFile>'
      close(1001)      
endif !if(gpid.eq.0)then


! Output of vtr file (actual local data)
write(pid2s,"(i4.4)") int(gpid)
call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
call get_pidxyz(pidx, pidy, pidz, gpid)

itmp=il-1
jtmp=jl-1
ktmp=kl-1
if(pidx.ge.1)itmp=il-2
if(pidy.ge.1)jtmp=jl-2
if(pidz.ge.1)ktmp=kl-2
write(str_extent_l,"(6(i4,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
fnamel=trim(opath)//'/'//trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
open(1002,FILE=fnamel,status='unknown',action='write')
write(1002,'(A)')'<?xml version="1.0"?>'
write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//'version="0.1" byte_order="LittleEndian">'
write(1002,'(A)')'<RectilinearGrid WholeExtent="'//trim(adjustl(str_extent_g))//'">'
write(1002,'(A)')'<Piece Extent="'//trim(adjustl(str_extent_l))//'">'
write(1002,'(A)')'<Coordinates>'
write(1002,'(A)')'<DataArray Name="x" type="Float32" format="ascii">'
do i=itmp,iu-1
      write(1002,'(F0.15)') i*DXL
enddo
write(1002,'(A)')'</DataArray>'
write(1002,'(A)')'<DataArray Name="y" type="Float32" format="ascii">'
do j=jtmp,ju-1
      write(1002,'(F0.15)') j*DYL
enddo
write(1002,'(A)')'</DataArray>'
write(1002,'(A)')'<DataArray Name="z" type="Float32" format="ascii">'
do k=ktmp,ku-1
      write(1002,'(F0.15)') k*DZL
enddo
write(1002,'(A)')'</DataArray>'
write(1002,'(A)')'</Coordinates>'

write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
write(1002,'(A)')'<DataArray type="Int32" Name="'//trim(fnamebase)//'" format="ascii">'
itmp=1-1
jtmp=1-1
ktmp=1-1
if(pidx.ge.1)itmp=1-2
if(pidy.ge.1)jtmp=1-2
if(pidz.ge.1)ktmp=1-2
do k=ktmp,kml-1
do j=jtmp,jml-1
do i=itmp,iml-1
      write(1002,'(I0)')a(i+1,j+1,k+1)
enddo
enddo
enddo

write(1002,'(A)')'</DataArray>'
write(1002,'(A)')'</PointData>'
write(1002,'(A)')'</Piece>'
write(1002,'(A)')'</RectilinearGrid>'
write(1002,'(A)')'</VTKFile>'

end subroutine output_i4_pvtr


subroutine merge_i4_raw()
end subroutine merge_i4_raw
end module outputs