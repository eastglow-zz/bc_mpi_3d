      module draw
      use input
      use boundary_conditions
      implicit none

      contains

      subroutine put_sphere_i4(a, ig, jg, kg, r,val, ng, bctype,pid)
      implicit none
      integer(4),intent(inout) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      integer(4),intent(in) :: ig,jg,kg,r
      integer(4),intent(in) :: val
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: bctype
      integer(4),intent(in) :: pid

      integer(4) :: i,j,k
      integer(4) :: il,jl,kl
      integer(4) :: diff,kkrank
      
      do k = kg-2*r, kg+2*r
      do j = jg-2*r, jg+2*r
      do i = ig-2*r, ig+2*r
            diff = (i-ig)**2 + (j-jg)**2 + (k-kg)**2 - r*r
            if(diff.le.0)then
                  KKRANK = get_pid_from_globalijk(I,J,K,bctype)
                  call get_localijk_from_globalijk(il,jl,kl, I,J,K
     &                                                   ,bctype)
                  if(kkrank.eq.pid)then
                        a(il,jl,kl)=val 
                  endif
            endif
      enddo
      enddo
      enddo

      end subroutine put_sphere_i4

      end module draw