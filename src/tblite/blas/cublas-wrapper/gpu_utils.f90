module tblite_blas_gpu_utils
   use iso_c_binding
   
   implicit none
   private
   public :: cuda_runtime, cuda_dmatrix, dp, get_cublas_error, get_cuda_error, cublasCreate
   public :: cublasDestroy, cudaDeviceSynchronize, sp, cuda_smatrix
   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: sp = selected_real_kind(6)
   type cuda_runtime
      integer(kind=c_intptr_t) :: cuda_device
      
      logical :: double = .false.
      type(c_ptr) :: handle
   contains
      procedure :: setup
      procedure :: create_handle
      procedure :: destroy_handle
   end type

   type cuda_dmatrix
      type(c_ptr), allocatable :: ptr
      real(c_double), pointer :: mat(:, :)
      integer(kind=c_int) :: rows, cols
      integer(kind=c_int) :: elemSize = 8
   contains 
      generic :: copy_2device => copy_2device_dp
      generic :: copy_2host => copy_2host_dp
      generic :: free_device => free_device_dp
      procedure :: copy_2device_dp
      procedure :: copy_2host_dp
      procedure :: free_device_dp
   end type

   type cuda_smatrix
      type(c_ptr), allocatable :: ptr
      real(c_float), pointer :: mat(:, :)
      integer(kind=c_int) :: rows, cols
      integer(kind=c_int) :: elemSize = 4
   contains 
      generic :: copy_2device => copy_2device_sp
      generic :: copy_2host => copy_2host_sp
      generic :: free_device => free_device_sp
      procedure :: copy_2device_sp
      procedure :: copy_2host_sp
      procedure :: free_device_sp
   end type

   type cublasHandle
      type(c_ptr) :: handle
   end type

   interface
   integer(c_int) function cublasGetVersion(handle, version) bind(C, name="cublasGetVersion_v2")
      use iso_c_binding
      type(c_ptr), value :: handle
      integer(c_int) :: version
   end function
   integer(c_int) function cublasCreate(handle) bind(C, name="cublasCreate_v2")
      use iso_c_binding
      type(c_ptr) :: handle
   end function cublasCreate
   integer(c_int) function cublasDestroy(handle) bind(C, name="cublasDestroy_v2")
      use iso_c_binding
      type(c_ptr), value :: handle
   end function cublasDestroy
   integer  (c_int) function cudaMalloc(ptr, size) bind(C,name='cudaMalloc')
      use iso_c_binding
      integer(c_int), value:: size
      type(c_ptr) :: ptr
   end function cudaMalloc
   integer  (c_int) function cudaFree(ptr) bind(C,name='cudaFree')
      use iso_c_binding
      type(c_ptr),value :: ptr
   end function cudaFree
   integer  (c_int) function cublasSetMatrix(rows, cols, elemSize, A, lda, devA, lda2) bind(C,name='cublasSetMatrix')
      use iso_c_binding
      integer(c_int), value :: rows, cols, elemSize, lda, lda2
      type(c_ptr), value :: A
      type(c_ptr), value :: devA
   end function cublasSetMatrix
   integer  (c_int) function cublasGetMatrix(rows, cols, elemSize, devA, lda, B, ldb) bind(C,name='cublasGetMatrix')
      use iso_c_binding
      integer(c_int), value:: rows,cols, elemSize,lda,ldb
      type(c_ptr), value :: devA
      type(c_ptr), value :: B
   end function cublasGetMatrix
   integer (c_int) function cudaMemGetInfo(fre, tot) bind(C, name="cudaMemGetInfo")
      use iso_c_binding
      implicit none
      integer(kind=c_intptr_t) :: fre
      integer(kind=c_intptr_t) :: tot
    end function cudaMemGetInfo
   integer (c_int) function cudaSetDevice(dev) bind(C, name="cudaSetDevice")
      use iso_c_binding
      implicit none
      integer(c_intptr_t), value :: dev
   end function cudaSetDevice
   integer (c_int) function cudaGetDeviceCount(n_dev) bind(C, name="cudaGetDeviceCount")
      use iso_c_binding
      implicit none
      integer(c_int) :: n_dev
   end function cudaGetDeviceCount
   type(c_ptr) function cudaGetErrorName(error) bind(C, name="cudaGetErrorName")
      use iso_c_binding
      implicit none
      integer(c_int), value :: error
   end function
   type(c_ptr) function cudaGetErrorString(error) bind(C, name="cudaGetErrorString")
      use iso_c_binding
      implicit none
      integer(c_int), value :: error
   end function
   type(c_ptr) function cublasGetStatusName(error) bind(C, name="cublasGetStatusName")
      use iso_c_binding
      implicit none
      integer(c_int), value :: error
   end function
   type(c_ptr) function cublasGetStatusString(error) bind(C, name="cublasGetStatusString")
      use iso_c_binding
      implicit none
      integer(c_int), value :: error
   end function
   subroutine cudaDeviceSynchronize() BIND(C, NAME='cudaDeviceSynchronize')
   end subroutine
   end interface

contains

subroutine create_handle(self)
   class(cuda_runtime) :: self
   integer(c_int) :: stat = 0 
   integer(c_int) :: version = 0
   stat = cudaSetDevice(self%cuda_device)
   self%handle = c_null_ptr
   call get_cuda_error(stat)
   stat = cublasCreate(self%handle)
   call get_cublas_error(stat)
   stat = cublasGetVersion(self%handle, version)
   call get_cublas_error(stat)
   write(*,*) version 
end subroutine

subroutine destroy_handle(self)
   class(cuda_runtime) :: self
   integer(c_int) :: stat = 0 
   stat = cublasDestroy(self%handle)
   call get_cublas_error(stat)
end subroutine


subroutine copy_2device_dp(self, mat, cudart)
   class(cuda_dmatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=8), allocatable :: mat(:,:)
   integer(c_int) :: stat = 0
   real(kind=c_double) :: size_
   self%cols = size(mat, dim=2)
   self%rows = size(mat, dim=1)
   self%elemSize = c_sizeof(size_)

   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)
   
   allocate(self%ptr)
   stat = cudaMalloc(self%ptr, self%cols*self%rows*self%elemSize)
   call get_cuda_error(stat)
   allocate(self%mat(self%rows, self%cols))
   self%mat = mat
   stat = cublasSetMatrix(self%rows, self%cols, self%elemSize, c_loc(self%mat), self%rows, self%ptr, self%rows)
   call get_cublas_error(stat)
   deallocate(mat)
end subroutine

subroutine free_device_dp(self, cudart)
   class(cuda_dmatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=8), allocatable :: mat(:,:)
   integer(c_int) :: stat = 0
   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)
   if (.not. allocated(self%ptr)) return

   stat = cudaFree(self%ptr)
   call get_cuda_error(stat)
   deallocate(self%ptr)
end subroutine

subroutine copy_2host_dp(self, mat, cudart)
   class(cuda_dmatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=8), allocatable :: mat(:,:)
   
   integer(c_int) :: stat = 0 
   if (.not.(allocated(self%ptr))) return
   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)
 
   stat = cublasGetMatrix(self%rows, self%cols, self%elemSize, self%ptr, self%rows, c_loc(self%mat), self%rows)
   call get_cublas_error(stat)
   stat = cudaFree(self%ptr)
   call get_cuda_error(stat)
   deallocate(self%ptr)
   mat = self%mat
   deallocate(self%mat)
end subroutine

subroutine copy_2device_sp(self, mat, cudart)
   class(cuda_smatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=4), allocatable :: mat(:,:)
   integer(c_int) :: stat = 0
   real(kind=c_float) :: size_
   self%cols = size(mat, dim=2)
   self%rows = size(mat, dim=1)
   self%elemSize = c_sizeof(size_)
   

   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)
   
   allocate(self%ptr)
   stat = cudaMalloc(self%ptr, self%cols*self%rows*self%elemSize)
   call get_cuda_error(stat)
   allocate(self%mat(self%rows, self%cols))
   self%mat = mat 
   stat = cublasSetMatrix(self%rows, self%cols, self%elemSize, c_loc(self%mat), self%rows, self%ptr, self%rows)
   call get_cublas_error(stat)
   deallocate(mat)

end subroutine

subroutine free_device_sp(self, cudart)
   class(cuda_smatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=8), allocatable :: mat(:,:)
   integer(c_int) :: stat = 0
   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)
   if (.not. allocated(self%ptr)) return

   stat = cudaFree(self%ptr)
   call get_cuda_error(stat)
   deallocate(self%ptr)
end subroutine

subroutine copy_2host_sp(self, mat, cudart)
   class(cuda_smatrix) :: self
   type(cuda_runtime) :: cudart
   real(kind=4), allocatable :: mat(:,:)
   real(kind=c_float), allocatable :: tmp(:, :)
   integer(c_int) :: stat = 0 
   if (.not.(allocated(self%ptr))) return
   stat = cudaSetDevice(cudart%cuda_device)
   call get_cuda_error(stat)

   allocate(tmp(self%rows,self%cols))
   
   stat = cublasGetMatrix(self%rows, self%cols, self%elemSize, self%ptr, self%rows, c_loc(self%mat), self%rows)
   call get_cublas_error(stat)
   stat = cudaFree(self%ptr)
   call get_cuda_error(stat)
   deallocate(self%ptr)
   mat = self%mat
   deallocate(self%mat)
end subroutine

subroutine get_cuda_error(stat)
   character(kind=c_char), pointer, dimension(:) :: errorstr
   integer(kind=c_size_t) :: stringlen = 255
   character(len=:), allocatable :: fchar, totchar
   integer(kind=c_int) :: stat
   if(stat /= 0) then
      call c_f_pointer(cudaGetErrorName(stat), errorstr, [ stringlen ])
      call c_f_character(errorstr, fchar)
      totchar = fchar
      call c_f_pointer(cudaGetErrorString(stat), errorstr, [ stringlen ])
      call c_f_character(errorstr, fchar)
      totchar = totchar//": "//fchar
      write(*,*) totchar
   end if
end subroutine

subroutine get_cublas_error(stat)
   character(kind=c_char), pointer, dimension(:) :: errorstr
   integer(kind=c_size_t) :: stringlen = 255
   character(len=:), allocatable :: fchar, totchar
   integer(kind=c_int) :: stat
   if(stat /= 0) then
      call c_f_pointer(cublasGetStatusName(stat), errorstr, [ stringlen ])
      call c_f_character(errorstr, fchar)
      totchar = fchar
      call c_f_pointer(cublasGetStatusString(stat), errorstr, [ stringlen ])
      call c_f_character(errorstr, fchar)
      totchar = totchar//": "//fchar
      write(*,*) totchar
   end if
end subroutine

subroutine setup(self, max_dim, n_mat)
   class(cuda_runtime) :: self
   !> maximum size of the matrices that need to be allocated
   integer, intent(in) :: max_dim
   !> number of matrcies that need to be allocated
   integer, intent(in) :: n_mat
   integer(kind=c_intptr_t) :: freemem, totmem, availmem
   integer(kind=c_int) :: stat, ndev
   integer(kind=c_size_t) :: request
   real(c_double) :: array(max_dim,max_dim)
   real(c_float) :: array_s(max_dim,max_dim)
   integer(c_intptr_t) :: i
   
   freemem = 0
   totmem = 0
   stat = cudaGetDeviceCount(ndev)
   call get_cuda_error(stat)
   availmem = 0
   do i = 0, ndev -1
   ! to be done at a later point setup_cuda(self%max_dim, self%n_mat)
      stat = cudaSetDevice(i)
      stat = cudaMemGetInfo(freemem, totmem)
      if (stat .ne. 0 ) then
         call get_cuda_error(stat)
         self%cuda_device = -1
      end if
      if (freemem > availmem) then 
         availmem = freemem
         self%cuda_device = i
      end if
      
   end do
  
   request = n_mat*c_sizeof(array)
   if (request > availmem) self%double = .true.
   request = n_mat*c_sizeof(array_s)
   if (request > availmem) self%cuda_device = -1
end subroutine

subroutine c_f_character(rhs, lhs)
   character(kind=c_char), intent(in) :: rhs(*)
   character(len=:), allocatable, intent(out) :: lhs

   integer :: ii

   do ii = 1, huge(ii) - 1
      if (rhs(ii) == c_null_char) then
         exit
      end if
   end do
   allocate(character(len=ii-1) :: lhs)
   lhs = transfer(rhs(1:ii-1), lhs)

end subroutine c_f_character
end module




