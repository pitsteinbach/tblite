module tblite_xtbml_feature_dictionary
    use mctc_env_accuracy, only : wp, i8
    implicit none
    private
 
    public :: timer_type, format_time
 
    type :: feature_record
       character(len=:), allocatable :: label
       real(wp), allocatable :: array1(:) 
       real(wp), allocatable :: array2(:, :) 
       real(wp), allocatable :: array3(:, :, :) 
    end type feature_record
 
    type :: feature_dictionary_type
       integer :: n = 0
       integer, allocatable :: nat
       character(len=:), allocatable :: last
       type(feature_record), allocatable :: record(:)
    contains
       procedure :: add_entry !use an interface to distinguish if a vector or 
       procedure :: update_entry ! label and array pair
       procedure :: get_entry !check
       procedure :: get_label ! return label label
       procedure, private :: push
    end type feature_dictionary_type


    interface add_entry
        subroutine add_entry_label(self,label)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
        end subroutine
        subroutine add_entry_label_array1(self, label, ndim1)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            integer :: ndim1
        end subroutine
        subroutine add_entry_label_array2(self, label, ndim1, ndim2)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            integer :: ndim1, ndim2
        end subroutine
        subroutine add_entry_label_array3(self, label, ndim1, ndim2, ndim3)
            import : feature_dictionary_type, wp 
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            integer :: ndim1, ndim2, ndim3
        end subroutine
        subroutine add_entry_label_input_array1(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:) 
        end subroutine
        subroutine add_entry_label_input_array2(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:,:) 
        end subroutine
        subroutine add_entry_label_input_array3(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:, :, :) 
        end subroutine 
    end interface

    interface get_entry
        subroutine get_entry_label_input_array1(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:) 
        end subroutine
        subroutine get_entry_label_input_array2(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:,:) 
        end subroutine
        subroutine get_entry_label_input_array3(self, label, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            character(len=*) :: label
            real(wp) :: array(:, :, :) 
        end subroutine 
        subroutine get_entry_index_input_array1(self, index, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            integer :: index
            real(wp) :: array(:) 
        end subroutine
        subroutine get_entry_index_input_array2(self, index, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            integer :: index
            real(wp) :: array(:,:) 
        end subroutine
        subroutine get_entry_index_input_array3(self, index, array)
            import : feature_dictionary_type, wp
            class(feature_dictionary_type) :: self
            integer :: index
            real(wp) :: array(:, :, :) 
        end subroutine 
    end interface
 
 contains

subroutine get_label(self, index, label)
    class(feature_dictionary_type) :: self
    integer :: index
    character(len=:), allocatable :: label
    if (index > self%n) return
    allocate(label(len(self%record(index)%label)))
    label = self%record(index)%label
end subroutine

subroutine get_entry_label_input_array1(self, label, array)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:) 
    integer : it

    call self%get_entry(return_label_index(self, label), array)

end subroutine

subroutine get_entry_label_input_array2(self, label, array)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:,:)

    call self%get_entry(return_label_index(self, label), array)

end subroutine

subroutine get_entry_label_input_array3(self, label, array)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:, :, :) 

    call self%get_entry(return_label_index(self, label), array)

end subroutine 

subroutine get_entry_index_input_array1(self, index, array)
    class(feature_dictionary_type) :: self
    integer :: index
    real(wp) :: array(:)

    if (index > self%n) return
    array = self%record(index)%array1

end subroutine

subroutine get_entry_index_input_array2(self, index, array)
    class(feature_dictionary_type) :: self
    integer :: index
    real(wp) :: array(:,:) 

    if (index > self%n) return
    array = self%record(index)%array2
end subroutine

subroutine get_entry_index_input_array3(self, index, array)
    class(feature_dictionary_type) :: self
    integer :: index
    real(wp) :: array(:, :, :) 

    if (index > self%n) return
    array = self%record(index)%array3
end subroutine 
 
subroutine push(self, label, it)
    class(timer_type), intent(inout) :: self
    character(len=*), intent(in) :: label
 
    integer, intent(out) :: it
 
    if (.not.allocated(self%record)) call resize(self%record)
    it = find(self%record(:self%n), label)
 
    if (it == 0) then
       if (self%n >= size(self%record)) then
          call resize(self%record)
       end if
 
       self%n = self%n + 1
       it = self%n
       self%record(it) = time_record(label)
    end if
 end subroutine push
 
subroutine add_entry_label(self,label)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    integer :: it

    call self%push(label, it)
end subroutine

subroutine add_entry_label_array1(self, label, ndim1)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    integer :: ndim1
    integer :: it

    call self%push(label, it)

    associate(record => self%record(it))
        allocate(record%array1(ndim1), source = 0.0_wp)
     end associate
end subroutine

subroutine add_entry_label_array2(self, label, ndim1, ndim2)
    
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    integer :: ndim1, ndim2

    integer :: it

    call self%push(label, it)

    associate(record => self%record(it))
        allocate(record%array2(ndim1, ndim2), source = 0.0_wp)
     end associate
end subroutine

subroutine add_entry_label_array3(self, label, ndim1, ndim2, ndim3)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    integer :: ndim1, ndim2, ndim3
    integer :: it

    call self%push(label, it)

    associate(record => self%record(it))
        allocate(record%array3(ndim1, ndim2, ndim3), source = 0.0_wp)
     end associate
end subroutine

subroutine add_entry_label_input_array1(self, label, array)    
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:) 
    integer :: it 

    call self%push(label, it)

    associate(record => self%record(it))
        call move_alloc(array, record%array2)
    end associate
end subroutine

subroutine add_entry_label_input_array2(self, label, array)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:,:) 
    integer :: it 

    call self%push(label, it)

    associate(record => self%record(it))
        call move_alloc(array, record%array1)
    end associate

end subroutine

subroutine add_entry_label_input_array3(self, label, array)
    class(feature_dictionary_type) :: self
    character(len=*) :: label
    real(wp) :: array(:, :, :) 
    integer :: it 

    call self%push(label, it)

    associate(record => self%record(it))
        call move_alloc(array, record%array3)
    end associate
end subroutine 
  
 
function return_label_index(self, label) result(it)
    class(feature_dictionary_type), intent(in) :: self
    character(len=*), intent(in) :: label
    real(wp) :: time
 
    integer :: it
 
    if (self%n <= 0) return
    it = find(self%record(:self%n), label)
    if (it == 0) return
 
end function return_label_index
 
 
pure function find(record, label) result(pos)
    type(time_record), intent(in) :: record(:)
    character(len=*), intent(in), optional :: label
    integer :: pos
 
    integer :: i
 
    pos = 0
    if (present(label)) then
       do i = size(record), 1, -1
          if (allocated(record(i)%label)) then
             if (label == record(i)%label) then
                pos = i
                exit
             end if
          end if
       end do
    else
       do i = size(record), 1, -1
          if (record(i)%running) then
             pos = i
             exit
          end if
       end do
    end if
 end function find
 
 
 function format_time(time) result(label)
    real(wp), intent(in) :: time
    character(len=:), allocatable :: label
 
    real(wp) :: secs
    integer :: mins, hours, days
 
    secs = time
    days = int(secs/86400.0_wp)
    secs = secs - days*86400.0_wp
    hours = int(secs/3600.0_wp)
    secs = secs - hours*3600.0_wp
    mins = int(secs/60.0_wp)
    secs = time - mins*60.0_wp
 
    if (days > 0) then
       label = format_label(days, '(i0, " d,")')
    else
       label = repeat(" ", 4)
    end if
    if (hours > 0) then
       label = label // format_label(hours, '(1x, i2, " h,")')
    else
       label = label // repeat(" ", 6)
    end if
    if (mins > 0) then
       label = label // format_label(mins, '(1x, i2, " min,")')
    else
       label = label // repeat(" ", 8)
    end if
    label = label // format_label(secs, '(f6.3)')//" sec"
 end function format_time
 
 
 function timing() result(time)
    real(wp) :: time
 
    integer(i8) :: time_count, time_rate, time_max
    call system_clock(time_count, time_rate, time_max)
    time = real(time_count, wp)/real(time_rate, wp)
 end function timing
 
 
 !> Reallocate list of timing records
 pure subroutine resize(var, n)
    !> Instance of the array to be resized
    type(time_record), allocatable, intent(inout) :: var(:)
    !> Dimension of the final array size
    integer, intent(in), optional :: n
 
    type(time_record), allocatable :: tmp(:)
    integer :: this_size, new_size
    integer, parameter :: initial_size = 20
 
    if (allocated(var)) then
       this_size = size(var, 1)
       call move_alloc(var, tmp)
    else
       this_size = initial_size
    end if
 
    if (present(n)) then
       new_size = n
    else
       new_size = this_size + this_size/2 + 1
    end if
 
    allocate(var(new_size))
 
    if (allocated(tmp)) then
       this_size = min(size(tmp, 1), size(var, 1))
       var(:this_size) = tmp(:this_size)
       deallocate(tmp)
    end if
 
 end subroutine resize
 
 end module tblite_timer