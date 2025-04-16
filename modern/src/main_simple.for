      program modern_mpi_test
      use mpi_f08
      implicit none

      type(MPI_Comm) :: comm
      type(MPI_Status) :: status
      integer :: rank, size, ierr, next, prev, tag
      integer :: send_data, recv_data

      call MPI_Init(ierr)
      comm = MPI_COMM_WORLD
      call MPI_Comm_rank(comm, rank, ierr)
      call MPI_Comm_size(comm, size, ierr)

      next = mod(rank + 1, size)
      prev = mod(rank - 1 + size, size)
      tag = 0

      send_data = rank
      call MPI_Send(send_data, 1, MPI_INTEGER, next, tag, comm, ierr)
      call MPI_Recv(recv_data, 1, MPI_INTEGER, prev, tag, comm, status, 
     &                                                             ierr)

      print *, 
     & 'Process ',rank,' received data ',recv_data,' from process ',prev

      call MPI_Finalize(ierr)
      end program modern_mpi_test

