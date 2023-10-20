module tblite_cublas_level3
   use iso_c_binding
   use tblite_blas_gpu_utils

   implicit none
   private
   public :: gemm
   
   !> Possible solvers provided by LAPACK
   type :: enum_cublas_operation
      integer(c_int) :: CUBLAS_OP_N = 0
      integer(c_int) :: CUBLAS_OP_T = 1
      integer(c_int) :: CUBLAS_OP_C = 2
   end type enum_cublas_operation

   !> Actual enumerator of possible solvers
   type(enum_cublas_operation), parameter :: cublasOperation_t = enum_cublas_operation()

   interface gemm
      module procedure :: wrap_dgemm
      module procedure :: wrap_sgemm
   end interface

   interface
      integer(c_int) function cublasDgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         & bind(C, name="cublasDgemm_v2")
         use iso_c_binding
         type(c_ptr), value :: handle
         integer(c_int), value :: transa, transb
         integer(c_int), value :: m, n, k, lda, ldb, ldc
         real(c_double) :: alpha, beta
         type(c_ptr), value :: A, B
         type(c_ptr), value :: C
      end function
   integer(c_int) function cublasSgemm(handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         & bind(C, name="cublasSgemm_v2")
         use iso_c_binding
         type(c_ptr), value :: handle
         integer(c_int), value :: transa, transb
         integer(c_int), value :: m, n, k, lda, ldb, ldc
         real(c_float) :: alpha, beta
         type(c_ptr), value :: A, B
         type(c_ptr), value :: C
      end function
   end interface
contains
   function get_cublasOperation(char) result(n)
      character(len=1) :: char
      integer(c_int) :: n
      if (char == "n" .or. char == "N") n = cublasOperation_t%CUBLAS_OP_N
      if (char == "t" .or. char == "T") n = cublasOperation_t%CUBLAS_OP_T
      if (char == "c" .or. char == "C") n = cublasOperation_t%CUBLAS_OP_C
   end function

   subroutine wrap_dgemm(amat, bmat, cmat, cudart, transa, transb, alpha, beta)
      type(cuda_dmatrix) :: amat, bmat, cmat
      type(cuda_runtime) :: cudart
      character(len=1), intent(in), optional :: transa
      character(len=1), intent(in), optional :: transb
      real(dp), intent(in), optional :: alpha
      real(dp), intent(in), optional :: beta
      real(c_double) :: a, b
      integer(c_int) :: tra, trb
      integer(c_int) :: lda, ldb, ldc, m, n, k
      integer(c_int) :: stat
      if (present(alpha)) then
         a = alpha
      else
         a = 1.0
      end if
      if (present(beta)) then
         b = beta
      else
         b = 0.0
      end if
      if (present(transa)) then
         tra = get_cublasOperation(transa)
      else
         tra = cublasOperation_t%CUBLAS_OP_N
      end if
      if (present(transb)) then
         trb = get_cublasOperation(transb)
      else
         trb = cublasOperation_t%CUBLAS_OP_N
      end if
      if ((tra.eq.cublasOperation_t%CUBLAS_OP_N)) then
         k = amat%cols
      else
         k = amat%rows
      end if
      lda = max(1, amat%rows)
      ldb = max(1, bmat%rows)
      ldc = max(1, cmat%rows)
      m = cmat%rows
      n = cmat%cols
      
      call cudaDeviceSynchronize()
      stat = cublasDgemm(cudart%handle, 0, 0, m, n, k, a, amat%ptr, lda, bmat%ptr, ldb, b, cmat%ptr, ldc)
      write(*,*) stat
      call get_cublas_error(stat)
      
      call cudaDeviceSynchronize()
      call get_cublas_error(stat)
      
   end subroutine

   subroutine wrap_sgemm(amat, bmat, cmat, cudart, transa, transb, alpha, beta)
      type(cuda_smatrix) :: amat, bmat, cmat
      type(cuda_runtime) :: cudart
      character(len=1), intent(in), optional :: transa
      character(len=1), intent(in), optional :: transb
      real(sp), intent(in), optional :: alpha
      real(sp), intent(in), optional :: beta
      real(c_float) :: a, b
      integer(c_int) :: tra, trb
      integer(c_int) :: lda, ldb, ldc, m, n, k
      integer(c_int) :: stat
      if (present(alpha)) then
         a = alpha
      else
         a = 1.0
      end if
      if (present(beta)) then
         b = beta
      else
         b = 0.0
      end if
      if (present(transa)) then
         tra = get_cublasOperation(transa)
      else
         tra = cublasOperation_t%CUBLAS_OP_N
      end if
      if (present(transb)) then
         trb = get_cublasOperation(transb)
      else
         trb = cublasOperation_t%CUBLAS_OP_N
      end if
      if ((tra.eq.cublasOperation_t%CUBLAS_OP_N)) then
         k = amat%cols
      else
         k = amat%rows
      end if
      lda = max(1, amat%rows)
      ldb = max(1, bmat%rows)
      ldc = max(1, cmat%rows)
      m = cmat%rows
      n = cmat%cols
      
      call cudaDeviceSynchronize()
      stat = cublasSgemm(cudart%handle, 0, 0, m, n, k, a, amat%ptr, lda, bmat%ptr, ldb, b, cmat%ptr, ldc)
      write(*,*) stat
      call get_cublas_error(stat)
      
      call cudaDeviceSynchronize()
      call get_cublas_error(stat)
      
   end subroutine
end module