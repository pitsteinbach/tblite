module tblite_cublas_level1
    use iso_c_binding
    use tblite_blas_gpu_utils, only : cuda_dmatrix, cuda_smatrix, cuda_runtime, cudaDeviceSynchronize, sp, dp, get_cublas_error
    implicit none
    private

    public :: Hadamard, axpy

    interface Hadamard
        module procedure :: wrap_dHadamard
        module procedure :: wrap_sHadamard
    end interface
    interface axpy
        module procedure :: wrap_daxpy
        module procedure :: wrap_saxpy
    end interface
    interface
        subroutine sHadamard(vecinout, vecin, ndim2) bind(C, name="sHadamard")
            use iso_c_binding
            type(c_ptr), value :: vecinout
            type(c_ptr), value :: vecin
            integer(c_size_t), value :: ndim2
        end subroutine
        subroutine dHadamard(vecinout, vecin, ndim2) bind(C, name="dHadamard")
            use iso_c_binding
            type(c_ptr), value :: vecinout
            type(c_ptr), value :: vecin
            integer(c_size_t), value :: ndim2
        end subroutine
        integer(c_int) function cublasSaxpy(handle, n, alpha, x, incx, y, incy) bind(C, name="cublasSaxpy_v2")
            use iso_c_binding
            type(c_ptr), value :: handle
            integer(c_int), value :: n, incx ,incy
            real(c_float), value :: alpha
            type(c_ptr), value :: x, y
        end function
        integer(c_int) function cublasDaxpy(handle, n, alpha, x, incx, y, incy) bind(C, name="cublasDaxpy_v2")
            use iso_c_binding
            type(c_ptr), value :: handle
            integer(c_int), value :: n, incx ,incy
            real(c_double), value :: alpha
            type(c_ptr), value :: x, y
        end function
    end interface
contains

subroutine wrap_sHadamard(vecinout, vecin, cudart)
    type(cuda_smatrix) :: vecinout, vecin
    type(cuda_runtime) :: cudart
    integer(c_size_t) :: vec1, vec2
    vec1 = vecinout%cols*vecinout%rows
    vec2 = vecin%cols*vecin%rows
    if (vec1 /= vec2) return
    call cudaDeviceSynchronize()
    call sHadamard(vecinout%ptr, vecin%ptr, vec1)
    call cudaDeviceSynchronize()
end subroutine

subroutine wrap_dHadamard(vecinout, vecin, cudart)
    type(cuda_dmatrix) :: vecinout, vecin
    type(cuda_runtime) :: cudart
    integer(c_size_t) :: vec1, vec2
    vec1 = vecinout%cols*vecinout%rows
    vec2 = vecin%cols*vecin%rows
    if (vec1 /= vec2) return
    call cudaDeviceSynchronize()
    call dHadamard(vecinout%ptr, vecin%ptr, vec1)
    call cudaDeviceSynchronize()
end subroutine

subroutine wrap_saxpy(a, x, y, cudart, incx, incy)
    type(cuda_smatrix) :: x, y
    real(kind=sp) :: a
    type(cuda_runtime) :: cudart
    integer, optional :: incx, incy
    integer(c_int) :: inx, iny, n
    integer(c_int) :: stat
    if (present(incx)) then
        inx = incx
    else
        inx = 1
    end if
    if (present(incy)) then
        iny = incy
    else
        iny = 1
    end if
    n = (x%cols*x%rows)/abs(inx)
    call cudaDeviceSynchronize()
    stat = cublasSaxpy(cudart%handle, n, a, x%ptr, inx, y%ptr, iny)
    call get_cublas_error(stat)
    call cudaDeviceSynchronize()
end subroutine

subroutine wrap_daxpy(a, x, y, cudart, incx, incy)
    type(cuda_dmatrix) :: x, y
    real(kind=dp) :: a
    type(cuda_runtime) :: cudart
    integer, optional :: incx, incy
    integer(c_int) :: inx, iny, n
    integer(c_int) :: stat
    if (present(incx)) then
        inx = incx
    else
        inx = 1
    end if
    if (present(incy)) then
        iny = incy
    else
        iny = 1
    end if
    n = (x%cols*x%rows)/abs(inx)
    call cudaDeviceSynchronize()
    stat = cublasDaxpy(cudart%handle, n, a, x%ptr, inx, y%ptr, iny)
    call get_cublas_error(stat)
    call cudaDeviceSynchronize()
end subroutine

end module