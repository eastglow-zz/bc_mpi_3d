      module diffusion_equation

      use mpi_f08

      use input
      use boundary_conditions

      implicit none 

      real(8), parameter :: Mobility = 10.0d0
      real(8), parameter :: Mob_p1 = 10.0d0
      real(8), parameter :: Mob_p2 = 1.0d0
      real(8), parameter :: Mob_p3 = 10.0d0

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
      !dx=1.0
      dx = 1.0
      dt=0.4*0.15*DXL*DXL/D
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

      subroutine diffusion_by_flux(a,ao,mob,pf,ngx,ngy,ngz,timestep, dt)
      implicit none
      real(8),intent(inout) :: a(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8),intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8),intent(inout) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4),intent(in) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4),intent(in) :: ngx,ngy,ngz
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: dt

      real(8) :: dcdt 
      real(8) :: lowerlimit, upperlimit

      integer(4) :: i,j,k

      do k=1,kml
      do j=1,jml  
      do i=1,iml
        
        ! call calc_dcdt(dcdt, i,j,k, ao, mob, ngx, ngy, ngz)
        call calc_dcdt_arithmetic(dcdt, i,j,k, ao, mob, ngx, ngy, ngz)

        a(i,j,k)=ao(i,j,k) + dcdt*dt

        ! Sanity check for valid concentration value 
        ! Identify concentration limits based on phase
        call comp_range_by_phase(pf(i,j,k), lowerlimit, upperlimit)
        if(a(i,j,k).lt.lowerlimit) then
          error_flag = .true.
          write(error_log,*)'Error, diffusion_by_flux(), ',
     &       'Concentration lower limit reached at (',i,',',j,',',k,')',
     &          ', phase=',pf(i,j,k), ', ao=',ao(i,j,k),', a=',a(i,j,k),
     &                                               ', mob=',mob(i,j,k)         
        else if(a(i,j,k).gt.upperlimit) then
          error_flag = .true.
          write(error_log,*)'Error, diffusion_by_flux(), ',
     &       'Concentration upper limit reached at (',i,',',j,',',k,')',
     &          ', phase=',pf(i,j,k), ', ao=',ao(i,j,k),', a=',a(i,j,k),
     &                                               ', mob=',mob(i,j,k)
        end if

      enddo
      enddo
      enddo

      return 
      end subroutine diffusion_by_flux

      ! ----------------------------------------------------------------

      subroutine calc_diffusion_dt(dt, mob, dx)
      implicit none
      real(8), intent(inout) :: dt
      real(8), intent(in) :: mob(:,:,:)
      real(8), intent(in) :: dx

      dt = 0.8 * 0.15 * dx * dx / maxval(mob)

      return 
      end subroutine calc_diffusion_dt

      ! ----------------------------------------------------------------

      subroutine calc_dcdt(dcdt, i,j,k, ao, mob, ngx, ngy, ngz)
      implicit none
      real(8), intent(inout) :: dcdt
      integer(4), intent(in) :: i,j,k
      real(8), intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: ngx, ngy, ngz

      real(8) :: jxp, jxn 
      real(8) :: jyp, jyn
      real(8) :: jzp, jzn
      real(8) :: div_x, div_y, div_z
      real(8) :: Meff_xp, Meff_xn 
      real(8) :: Meff_yp, Meff_yn
      real(8) :: Meff_zp, Meff_zn
      real(8),parameter :: smallnumber = 1.0d-10

      if(ngx.gt.0) then 
        Meff_xp=2.0*mob(i+1,j,k)*mob(i,j,k)/
     &                             (mob(i+1,j,k)+mob(i,j,k)+smallnumber)
        Meff_xn=2.0*mob(i,j,k)*mob(i-1,j,k)/
     &                             (mob(i,j,k)+mob(i-1,j,k)+smallnumber)
      else
        Meff_xp = 0.0d0
        Meff_xn = 0.0d0
      end if
      if(ngy.gt.0) then 
        Meff_yp=2.0*mob(i,j+1,k)*mob(i,j,k)/
     &                             (mob(i,j+1,k)+mob(i,j,k)+smallnumber)
        Meff_yn=2.0*mob(i,j,k)*mob(i,j-1,k)/
     &                             (mob(i,j,k)+mob(i,j-1,k)+smallnumber)
      else
        Meff_yp = 0.0d0
        Meff_yn = 0.0d0
      end if
      if(ngz.gt.0) then
        Meff_zp=2.0*mob(i,j,k+1)*mob(i,j,k)/
     &                             (mob(i,j,k+1)+mob(i,j,k)+smallnumber)
        Meff_zn=2.0*mob(i,j,k)*mob(i,j,k-1)/
     &                             (mob(i,j,k)+mob(i,j,k-1)+smallnumber)
      else
        Meff_zp = 0.0d0
        Meff_zn = 0.0d0
      end if

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

      return 
      end subroutine calc_dcdt

      ! ----------------------------------------------------------------


      subroutine calc_dcdt_arithmetic(dcdt, i,j,k, ao, mob, ngx,ngy,ngz)
      implicit none
      real(8), intent(inout) :: dcdt
      integer(4), intent(in) :: i,j,k
      real(8), intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: ngx, ngy, ngz

      real(8) :: jxp, jxn 
      real(8) :: jyp, jyn
      real(8) :: jzp, jzn
      real(8) :: div_x, div_y, div_z
      real(8) :: Meff_xp, Meff_xn 
      real(8) :: Meff_yp, Meff_yn
      real(8) :: Meff_zp, Meff_zn

      if(ngx.gt.0) then 
        Meff_xp=(mob(i+1,j,k)+mob(i,j,k))/2.0
        Meff_xn=(mob(i,j,k)+mob(i-1,j,k))/2.0
      else
        Meff_xp = 0.0d0
        Meff_xn = 0.0d0
      end if
      if(ngy.gt.0) then 
        Meff_yp=(mob(i,j+1,k)+mob(i,j,k))/2.0
        Meff_yn=(mob(i,j,k)+mob(i,j-1,k))/2.0
      else
        Meff_yp = 0.0d0
        Meff_yn = 0.0d0
      end if
      if(ngz.gt.0) then
        Meff_zp=(mob(i,j,k+1)+mob(i,j,k))/2.0
        Meff_zn=(mob(i,j,k)+mob(i,j,k-1))/2.0
      else
        Meff_zp = 0.0d0
        Meff_zn = 0.0d0
      end if

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

      return 
      end subroutine calc_dcdt_arithmetic

      ! ----------------------------------------------------------------

      subroutine calc_diffusion_mobility(mob_op,flow_factor,mob,ao,pf,
     &                                                   ngx,ngy,ngz,dt)
      implicit none
      real(8), intent(inout) :: mob_op(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(inout) :: flow_factor(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: ngx, ngy, ngz
      real(8), intent(in) :: dt 

      integer(4) :: i,j,k

      do k=1,kml
      do j=1,jml
      do i=1,iml
    !     call calc_flow_factor(flow_factor,i,j,k,ao,pf,mob,ngx,ngy,ngz,
    !  &                                                               dt)
        call calc_flow_factor_arithmetic(flow_factor,i,j,k,ao,pf,mob,
     &                                                   ngx,ngy,ngz,dt)
        !mob_op(i,j,k) = flow_factor(i,j,k) * mob(i,j,k)
        if (pf(i,j,k) .eq. 1) then
          mob_op(i,j,k) = Mob_p1 * flow_factor(i,j,k)
        else if( pf(i,j,k) .eq. 2) then
          mob_op(i,j,k) = Mob_p2 * flow_factor(i,j,k)
        else if ( pf(i,j,k) .eq. 3) then
          mob_op(i,j,k) = Mob_p3 * flow_factor(i,j,k)
        else
          mob_op(i,j,k) = flow_factor(i,j,k) * mob(i,j,k)
        end if
      enddo
      enddo
      enddo

      return 

      end subroutine calc_diffusion_mobility

      ! ----------------------------------------------------------------

      subroutine calc_flow_factor(flow_factor,i,j,k,ao,pf,mob,
     &                                                   ngx,ngy,ngz,dt)
      implicit none
      real(8), intent(inout) :: flow_factor(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: i,j,k
      real(8), intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: ngx, ngy, ngz
      real(8), intent(in) :: dt

      real(8) :: dcmax_over_dc ! Ratio of max allowable concentration change to predicted concentration change
      real(8) :: dcdt 
      real(8) :: cnew_try, delta_c
      real(8) :: lowerlimit, upperlimit
      real(8) :: af_xp, af_xn
      real(8) :: af_yp, af_yn
      real(8), parameter :: smallnumber = 1.0d-10

      flow_factor(i,j,k) = 1.0d0 ! Initialization with full flow capacity

      ! Calculate dcdt at this grid point
      call calc_dcdt(dcdt, i,j,k, ao, mob, ngx, ngy, ngz)

      cnew_try = ao(i,j,k) + dcdt*dt
      delta_c = cnew_try - ao(i,j,k)
      ! Identify concentration limits based on phase
      call comp_range_by_phase(pf(i,j,k), lowerlimit, upperlimit)

      if(cnew_try <= lowerlimit) then 
        dcmax_over_dc = (lowerlimit - ao(i,j,k)) / delta_c
      else if(cnew_try >= upperlimit) then
        dcmax_over_dc = (upperlimit - ao(i,j,k)) / delta_c
      else
        return ! No adjustment needed
      end if

      af_xn = dcmax_over_dc*mob(i-1,j,k)/
     &       ((1.0d0-dcmax_over_dc)*mob(i,j,k)+mob(i-1,j,k)+smallnumber)
      af_xp = dcmax_over_dc*mob(i+1,j,k)/
     &       ((1.0d0-dcmax_over_dc)*mob(i,j,k)+mob(i+1,j,k)+smallnumber)
      af_yn = dcmax_over_dc*mob(i,j-1,k)/
     &       ((1.0d0-dcmax_over_dc)*mob(i,j,k)+mob(i,j-1,k)+smallnumber)
      af_yp = dcmax_over_dc*mob(i,j+1,k)/
     &       ((1.0d0-dcmax_over_dc)*mob(i,j,k)+mob(i,j+1,k)+smallnumber)

      flow_factor(i,j,k) = min(af_xn, af_xp, af_yn, af_yp)

      if(flow_factor(i,j,k) .lt. 0.0d0) then 
        error_flag = .true.
        write(error_log,*) 
     &                'Error, calc_flow_factor(): Negative flow factor',
     &               ' at (', i,',',j,',',k,') pid=', myrank, 
     &               'pf=', pf(i,j,k),' lowerlimit=', lowerlimit,
     &               ' upperlimit=', upperlimit,
     &               ' cnew_try=', cnew_try, ' ao=', ao(i,j,k),
     &               ' dcdt=', dcdt, ' delta_c=', delta_c,
     &               ' dcmax_over_dc=', dcmax_over_dc,
     &               ' af_xn=', af_xn, ' af_xp=', af_xp,
     &               ' af_yn=', af_yn, ' af_yp=', af_yp
      end if


      return 
      end subroutine calc_flow_factor

      ! ----------------------------------------------------------------

      subroutine calc_flow_factor_arithmetic(flow_factor,i,j,k,ao,pf,
     &                                               mob,ngx,ngy,ngz,dt)
      implicit none
      real(8), intent(inout) :: flow_factor(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: i,j,k
      real(8), intent(in) :: ao(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      real(8), intent(in) :: mob(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4), intent(in) :: ngx, ngy, ngz
      real(8), intent(in) :: dt

      real(8) :: dcmax ! Ratio of max allowable concentration change to predicted concentration change
      real(8) :: dcdt 
      real(8) :: cnew_try, delta_c
      real(8) :: lowerlimit, upperlimit
      real(8) :: af_xp, af_xn
      real(8) :: af_yp, af_yn
      real(8), parameter :: smallnumber = 1.0d-10
      real(8) :: denominator

      flow_factor(i,j,k) = 1.0d0 ! Initialization with full flow capacity

      ! Calculate dcdt at this grid point
      call calc_dcdt_arithmetic(dcdt, i,j,k, ao, mob, ngx, ngy, ngz)

      cnew_try = ao(i,j,k) + dcdt*dt
      ! Identify concentration limits based on phase
      call comp_range_by_phase(pf(i,j,k), lowerlimit, upperlimit)

      if(cnew_try <= lowerlimit) then 
        dcmax = cnew_try - lowerlimit
        denominator = 0.5d0*mob(i,j,k)*dt*( smallnumber +
     &             (ao(i+1,j,k)-2.0d0*ao(i,j,k)+ao(i-1,j,k))/(DXL*DXL)+
     &             (ao(i,j+1,k)-2.0d0*ao(i,j,k)+ao(i,j-1,k))/(DYL*DYL) )
        flow_factor(i,j,k) = 1.0d0 - dcmax / denominator
      else if(cnew_try >= upperlimit) then
        dcmax = cnew_try - upperlimit
        denominator = 0.5d0*mob(i,j,k)*dt*( smallnumber +
     &             (ao(i+1,j,k)-2.0d0*ao(i,j,k)+ao(i-1,j,k))/(DXL*DXL)+
     &             (ao(i,j+1,k)-2.0d0*ao(i,j,k)+ao(i,j-1,k))/(DYL*DYL) )
        flow_factor(i,j,k) = 1.0d0 - dcmax / denominator
      else
        return ! No adjustment needed
      end if

      ! Negative flow_factor allowed here since arithmetic mean is used.

      return 
      end subroutine calc_flow_factor_arithmetic

      ! ----------------------------------------------------------------

      subroutine comp_range_by_phase(phaseid, comp_min, comp_max)
      implicit none
      integer(4),intent(in) :: phaseid
      real(8),intent(inout) :: comp_min, comp_max

      select case(phaseid)
      case(1)  ! Compressible background phase such as air
        comp_min = 0.0d0
        comp_max = 1000.0d0
      case(2)  ! Imcompressible liquid phase 
        comp_min = 0.0d0
        comp_max = 1.0d0
      case(3)  ! Solid solution with limited solubility
        comp_min = 0.20d0
        comp_max = 0.50d0
      case default 
        write(*,*)'comp_range_by_phase: Unknown phase ID'
        write(*,*)'  phaseid = ', phaseid
        write(*,*)'  Aborting...'
        call abort()
      end select

      return 
      end subroutine comp_range_by_phase

      ! ----------------------------------------------------------------

      subroutine phaseid_ic_box(pf, ngx, ngy, ngz, ilb, iub, jlb, jub, 
     &                                                klb, kub, phaseid)
      implicit none
      integer(4),intent(inout) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4),intent(in) :: ngx, ngy, ngz
      integer(4),intent(in) :: ilb, iub, jlb, jub, klb, kub
      integer(4),intent(in) :: phaseid

      integer(4) :: i, j, k
      integer(4) :: ig, jg, kg

      do k=1,kml
      do j=1,jml
      do i=1,iml
        call get_global_ijk(ig,jg,kg,i,j,k,myrank)
        if (ig.ge.ilb .and. ig.le.iub .and. jg.ge.jlb .and. jg.le.jub 
     &                             .and. kg.ge.klb .and. kg.le.kub) then
          pf(i,j,k) = phaseid
        end if
      end do
      end do
      end do

      return
      end subroutine phaseid_ic_box

      ! ----------------------------------------------------------------

      subroutine diffu_ic_by_phase(a,pf,ngx,ngy,ngz,
     &                                       c_phase1,c_phase2,c_phase3)
      implicit none
      real(8),intent(inout) :: a(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4),intent(in) :: pf(1-ngx:iml+ngx,1-ngy:jml+ngy,
     &                                                    1-ngz:kml+ngz)
      integer(4),intent(in) :: ngx,ngy,ngz
      real(8),intent(in) :: c_phase1,c_phase2,c_phase3

      integer(4) :: i,j,k

      do k=1,kml
      do j=1,jml
      do i=1,iml
        select case (pf(i,j,k))
        case (1)
          a(i,j,k) = c_phase1
          ! write(*,*) 'Setting phase 1 at (', i, ',', j, ',', k, ')'
        case (2)
          ! write(*,*) 'Setting phase 2 at (', i, ',', j, ',', k, ')'
          a(i,j,k) = c_phase2
        case (3)
          ! write(*,*) 'Setting phase 3 at (', i, ',', j, ',', k, ')'
          a(i,j,k) = c_phase3
        case default
          write(*,*) 'Unknown phase ID at (', i, ',', j, ',', k, ')',
     &                                                  'pid = ', myrank
          write(*,*) '  Check phase ID mapping.'
          write(*,*) '  Aborting...'
          call abort()
        end select
      end do
      end do
      end do

      return 
      end subroutine diffu_ic_by_phase

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


      ! ----------------------------------------------------------------

      real(8) function total_amount_c(a,ngx,ngy,ngz)
      implicit none
      real(8),intent(in) :: a(1-ngx:iml+ngx,1-ngy:jml+ngy,1-ngz:kml+ngz)
      integer(4),intent(in) :: ngx,ngy,ngz
      real(8) :: total

      integer(4) :: i,j,k

      total = 0.0d0
      do k=1,kml
      do j=1,jml
      do i=1,iml
        total = total + a(i,j,k)
      end do
      end do
      end do

      call MPI_Allreduce(MPI_IN_PLACE, total, 1, MPI_DOUBLE_PRECISION, 
     &                                          MPI_SUM, MPI_COMM_WORLD)

      total_amount_c = total 
      
      end function total_amount_c

      end module diffusion_equation