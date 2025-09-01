      module diffusion_equation

      use mpi_f08

      use input
      use boundary_conditions

      implicit none 

      real(8), parameter :: Mobility = 10.0d0

      contains 

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

      subroutine diffusion_sourced(a,ao,ngx,ngy,ngz,timestep,rate, 
     &                                                       upperlimit)
      implicit none 

      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      real(8),intent(inout) :: 
     &                     ao(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: rate  ! Source strength
      real(8),intent(in) :: upperlimit  ! Upper limit for concentration

      real(8) :: M, dx, dt
      real(8) :: jxp, jxn 
      real(8) :: jyp, jyn
      real(8) :: jzp, jzn
      real(8) :: div_x, div_y, div_z
      real(8) :: cnew_try, dcdt 
      real(8) :: jxp_block, jxn_block 
      real(8) :: jyp_block, jyn_block
      real(8) :: jzp_block, jzn_block
      real(8) :: Meff_xp, Meff_xn 
      real(8) :: Meff_yp, Meff_yn
      real(8) :: Meff_zp, Meff_zn
      real(8) :: source_term 

      integer(4) :: i,j,k
      integer(4) :: ig,jg,kg

      M = Mobility
      dt=0.8*0.15*DXL*DXL/M
      dttime = dt

      do k=1,kml
      do j=1,jml
      do i=1,iml
        ! Get global coorinate 
        call get_global_ijk(ig,jg,kg,i,j,k,myrank)

        if(ig.eq.max(IMORI/2,1).and.  jg.eq.max(JMORI/2,1).and.
     &                                        kg.eq.max(KMORI/2,1)) then
          !write(*,*)'Source term applied at (',ig,',',jg,',',kg,')'
          source_term = rate
        else
          source_term = 0.0d0
        end if

        ! Initialization with uniform mobility
        Meff_xp = M
        Meff_xn = M
        Meff_yp = M
        Meff_yn = M
        Meff_zp = M
        Meff_zn = M

        ! Fluxes in x direction
        if(ngx.gt.0) then
          jxn = - Meff_xn * (ao(i,j,k) - ao(i-1,j,k)) / DXL
          jxp = - Meff_xp * (ao(i+1,j,k) - ao(i,j,k)) / DXL
        else
          jxn = 0.0d0
          jxp = 0.0d0
        end if

        ! Fluxes in y direction
        if(ngy.gt.0) then
          jyn = - Meff_yn * (ao(i,j,k) - ao(i,j-1,k)) / DYL
          jyp = - Meff_yp * (ao(i,j+1,k) - ao(i,j,k)) / DYL
        else
          jyn = 0.0d0
          jyp = 0.0d0
        end if

        ! Fluxes in z direction
        if(ngz.gt.0) then
          jzn = - Meff_zn * (ao(i,j,k) - ao(i,j,k-1)) / DZL
          jzp = - Meff_zp * (ao(i,j,k+1) - ao(i,j,k)) / DZL
        else
          jzn = 0.0d0
          jzp = 0.0d0
        end if

        div_x = (jxp - jxn) / DXL 
        div_y = (jyp - jyn) / DYL
        div_z = (jzp - jzn) / DZL

        dcdt = - (div_x + div_y + div_z) + source_term

        cnew_try = ao(i,j,k) + dcdt*dt

        ! Put zero mobility at this grid point if cnew_try > upperlimit.
        if(cnew_try .gt. upperlimit) then
          write(*,*)
     &     'Concentration upper limit reached at (',ig,',',jg,',',kg,')'
          write(*,*)'  mobility set zero here'
          Meff_xn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
          Meff_xp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_yp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_yn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
          Meff_zp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_zn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
          source_term = 0.0d0
        endif

        ! Fluxes in x direction
        if(ngx.gt.0) then
          jxn = - Meff_xn * (ao(i,j,k) - ao(i-1,j,k)) / DXL
          jxp = - Meff_xp * (ao(i+1,j,k) - ao(i,j,k)) / DXL
        else
          jxn = 0.0d0
          jxp = 0.0d0
        end if

        ! Fluxes in y direction
        if(ngy.gt.0) then
          jyn = - Meff_yn * (ao(i,j,k) - ao(i,j-1,k)) / DYL
          jyp = - Meff_yp * (ao(i,j+1,k) - ao(i,j,k)) / DYL
        else
          jyn = 0.0d0
          jyp = 0.0d0
        end if

        ! Fluxes in z direction
        if(ngz.gt.0) then
          jzn = - Meff_zn * (ao(i,j,k) - ao(i,j,k-1)) / DZL
          jzp = - Meff_zp * (ao(i,j,k+1) - ao(i,j,k)) / DZL
        else
          jzp = 0.0d0
          jzn = 0.0d0
        end if

        div_x = (jxp - jxn) / DXL 
        div_y = (jyp - jyn) / DYL
        div_z = (jzp - jzn) / DZL

        dcdt = - (div_x + div_y + div_z) + source_term

        a(i,j,k)=ao(i,j,k) + dcdt*dt

      enddo
      enddo
      enddo

      return 
      end subroutine diffusion_sourced

      ! ----------------------------------------------------------------

      subroutine diffusion_biphase_box(a,ao,ngx,ngy,ngz,timestep, 
     &      ilb, iub, jlb, jub, klb, kub, upperlimit_in, upperlimit_out)
      implicit none 

      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      real(8),intent(inout) :: 
     &                     ao(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: timestep
      integer(4),intent(in) :: ilb, iub, jlb, jub, klb, kub
      real(8),intent(in) :: upperlimit_in, upperlimit_out

      real(8) :: M, dx, dt
      real(8) :: jxp, jxn 
      real(8) :: jyp, jyn
      real(8) :: jzp, jzn
      real(8) :: div_x, div_y, div_z
      real(8) :: cnew_try, dcdt 
      real(8) :: jxp_block, jxn_block 
      real(8) :: jyp_block, jyn_block
      real(8) :: jzp_block, jzn_block
      real(8) :: Meff_xp, Meff_xn 
      real(8) :: Meff_yp, Meff_yn
      real(8) :: Meff_zp, Meff_zn
      real(8) :: upperlimit
      

      integer(4) :: i,j,k
      integer(4) :: ig,jg,kg

      M = Mobility
      dt=0.8*0.15*DXL*DXL/M
      dttime = dt

      do k=1,kml
      do j=1,jml
      do i=1,iml
        ! Get global coorinate 
        call get_global_ijk(ig,jg,kg,i,j,k,myrank)

        ! Initialization with uniform mobility
        Meff_xp = M
        Meff_xn = M
        Meff_yp = M
        Meff_yn = M
        Meff_zp = M
        Meff_zn = M

        ! Fluxes in x direction
        if(ngx.gt.0) then
          jxn = - Meff_xn * (ao(i,j,k) - ao(i-1,j,k)) / DXL
          jxp = - Meff_xp * (ao(i+1,j,k) - ao(i,j,k)) / DXL
        else
          jxn = 0.0d0
          jxp = 0.0d0
        end if

        ! Fluxes in y direction
        if(ngy.gt.0) then
          jyn = - Meff_yn * (ao(i,j,k) - ao(i,j-1,k)) / DYL
          jyp = - Meff_yp * (ao(i,j+1,k) - ao(i,j,k)) / DYL
        else
          jyn = 0.0d0
          jyp = 0.0d0
        end if

        ! Fluxes in z direction
        if(ngz.gt.0) then
          jzn = - Meff_zn * (ao(i,j,k) - ao(i,j,k-1)) / DZL
          jzp = - Meff_zp * (ao(i,j,k+1) - ao(i,j,k)) / DZL
        else
          jzn = 0.0d0
          jzp = 0.0d0
        end if

        div_x = (jxp - jxn) / DXL 
        div_y = (jyp - jyn) / DYL
        div_z = (jzp - jzn) / DZL

        dcdt = - (div_x + div_y + div_z)

        cnew_try = ao(i,j,k) + dcdt*dt

        ! Identify the upper limit at this grid point 
        if (ig.ge.ilb .and. ig.le.iub .and. 
     &                                   jg.ge.jlb .and. jg.le.jub .and. 
     &                                   kg.ge.klb .and. kg.le.kub) then
          upperlimit = upperlimit_in
        else
          upperlimit = upperlimit_out
        end if

        ! Put zero mobility at this grid point if cnew_try > upperlimit.
        if(cnew_try .gt. upperlimit) then
    !       write(*,*)
    !  &     'Concentration upper limit reached at (',ig,',',jg,',',kg,')'
    !       write(*,*)'  mobility set zero here'
          Meff_xn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
          Meff_xp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_yp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_yn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
          Meff_zp = 2.0*0.0d0*M/(0.0d0+M)
          Meff_zn = 2.0*M*0.0d0/(M+0.0d0) ! Harmonic mean
        endif

        ! Fluxes in x direction
        if(ngx.gt.0) then
          jxn = - Meff_xn * (ao(i,j,k) - ao(i-1,j,k)) / DXL
          jxp = - Meff_xp * (ao(i+1,j,k) - ao(i,j,k)) / DXL
        else
          jxn = 0.0d0
          jxp = 0.0d0
        end if

        ! Fluxes in y direction
        if(ngy.gt.0) then
          jyn = - Meff_yn * (ao(i,j,k) - ao(i,j-1,k)) / DYL
          jyp = - Meff_yp * (ao(i,j+1,k) - ao(i,j,k)) / DYL
        else
          jyn = 0.0d0
          jyp = 0.0d0
        end if

        ! Fluxes in z direction
        if(ngz.gt.0) then
          jzn = - Meff_zn * (ao(i,j,k) - ao(i,j,k-1)) / DZL
          jzp = - Meff_zp * (ao(i,j,k+1) - ao(i,j,k)) / DZL
        else
          jzp = 0.0d0
          jzn = 0.0d0
        end if

        div_x = (jxp - jxn) / DXL 
        div_y = (jyp - jyn) / DYL
        div_z = (jzp - jzn) / DZL

        dcdt = - (div_x + div_y + div_z)

        a(i,j,k)=ao(i,j,k) + dcdt*dt

      enddo
      enddo
      enddo

      return 
      end subroutine diffusion_biphase_box

      ! ----------------------------------------------------------------

      subroutine diffu_ic_stepfunction_x(a,ngx,ngy,ngz,ilow,ihigh,val)
      implicit none
      real(8),intent(inout) :: 
     &                      a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: ngx,ngy,ngz
      integer(4),intent(in) :: ilow,ihigh
      real(8),intent(in) :: val

      integer(4) :: i,j,k
      integer(4) :: ig,jg,kg

      do k=1,kml 
      do j=1,jml 
      do i=1,iml
        call get_global_ijk(ig,jg,kg,i,j,k,myrank)
        if (ig.ge.ilow .and. ig.le.ihigh) then
          write(*,*)
     &             'Setting initial condition at (',ig,',',jg,',',kg,')'
          write(*,*)'  Setting value to: ', val
          a(i,j,k) = val
        else
          a(i,j,k) = 0.0d0
        end if
      enddo
      enddo
      enddo

      return 
      end subroutine diffu_ic_stepfunction_x

      end module diffusion_equation