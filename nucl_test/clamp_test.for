      program clamp_test
      use input
      use boundary_conditions
      implicit none

      integer(4) :: ig,jg,kg
      integer(4) :: igf,jgf,kgf

      im=100
      jm=100
      km=100

      write(*,*) 'IM, KM, KM = ',im,jm,km
      write(*,*) 'Enter ig, jg, kg:'
      read(*,*) ig, jg, kg

      call clamp_globalijk(igf,jgf,kgf,ig,jg,kg,'periodic')

      write(*,*) 'After clamping'
      write(*,*) igf, jgf, kgf

      end program clamp_test