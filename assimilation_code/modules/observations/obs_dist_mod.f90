module obs_dist_mod
    use mpi
    use iso_c_binding
    ! use ifport, only : random_number
    use     mpi_utilities_mod, only : task_count, my_task_id, send_to, receive_from
    use        types_mod, only : r8, i8, MISSING_R8, metadatalength

    use      obs_def_mod, only : obs_def_type, get_obs_def_time, read_obs_def, &
                                 write_obs_def, destroy_obs_def, copy_obs_def, &
                                 interactive_obs_def, get_obs_def_location, &
                                 get_obs_def_type_of_obs, get_obs_def_key,  &
                                 get_obs_def_error_variance, &
                                 set_obs_def_key, set_obs_def_time, set_obs_def_location, &
                                 set_obs_def_type_of_obs, set_obs_def_error_variance
    ! use obs_sequence_mod, only : obs_type
    use     location_mod, only : location_type, is_location_in_region, get_location, set_location, query_location
    use time_manager_mod, only : time_type, set_time, print_time, print_date, get_time, &
                                 operator(-), operator(+), &
                                 operator(>), operator(<), &
                                 operator(>=), operator(/=), operator(==)

    ! Place this here for now
    ! not ideal since type is no longer private but will have to do for now
    type obs_type
       ! private
    ! The key is needed to indicate the element number in the storage for the obs_sequence
    ! Do I want to enforce the identity of the particular obs_sequence?
       integer :: key
       integer(i8) :: time_actual
       type(obs_def_type) :: def
       real(r8), allocatable :: values(:)
       real(r8), allocatable :: qc(:)
       ! real(r8), allocatable :: values(:)
       ! real(r8), allocatable :: qc(:)      
       ! Put sort indices directly into the data structure
       integer :: prev_time, next_time
       integer :: cov_group
    end type obs_type

    type obs_values_qc_type
       integer :: time_order
       integer(i8) :: time_actual
       real(r8) :: val
       real(r8) :: qc
    end type obs_values_qc_type

    type obs_type_send
    ! The key is needed to indicate the element number in the storage for the obs_sequence
    ! Do I want to enforce the identity of the particular obs_sequence?
    ! Declare a simplified data structure for sending and receiving using MPI
       integer :: key
       ! Remove obs_def_type when sending obs to reduce complexity
       ! type(obs_def_type) :: def
       ! Define what is in obs_def_type here
       ! Location
       !real(r8) :: vloc
       !type(location_type) :: location
       integer  :: which_vert
       !type(time_type)     :: time
       integer  :: kind 
       ! time_type
       integer  :: seconds
       integer  :: days
       ! Other obs_def_type values
       integer  :: obs_def_key
       ! real(r8), allocatable    :: external_FO(:)
       !integer  :: ens_size
       ! cannot send pointers as component of derived types
       ! real(r8), pointer :: values(:)  => NULL()
       ! real(r8), pointer :: qc(:)      => NULL()
       ! Put sort indices directly into the data structure
       ! todo: We can't do that when we're parallelizing; think of something better
       integer :: prev_time, next_time
       integer :: cov_group
       integer :: time_order
       integer(i8) :: time_actual
       real(r8) :: lon
       real(r8) :: lat
       real(r8) :: vloc
       real(r8) :: error_variance
    end type obs_type_send

    type obs_dist_type
        ! use these for our unpacked observations and values
        ! will be accessed whenever a process attempts a one-sided get
        !
        ! Why pointers (instead of 'allocate')? I would like to switch where these are pointing 
        ! after the observations have been sorted, and would like to avoid 
        ! an accidental deallocation happening under my feet
        type(obs_values_qc_type), pointer   :: val_buf(:) => NULL()
        type(obs_type_send), pointer        :: obs_buf(:) => NULL()
        ! real(r8), pointer                       :: obs_reals(:) => NULL()
        integer                                 :: obs_mpi
        integer                                 :: val_mpi
        integer                                 :: num_obs_per_proc
        integer                                 :: our_num_obs
        integer                                 :: num_vals_per_obs
        integer                                 :: num_vals_per_proc
        integer                                 :: total_obs
        integer                                 :: rem
        integer                                 :: my_pe
        integer                                 :: nprocs   
        integer                                 :: obs_win
        integer                                 :: val_win
        integer                                 :: mpi_time
        integer                                 :: ngets
        integer                                 :: dist_type
        integer                                 :: ierror
    end type obs_dist_type
    type(obs_dist_type) :: odt
    ! integer :: obs_win, val_win

    type obs_sample_type
        private
        integer                                 :: key
        integer                                 :: selected
    end type obs_sample_type

    interface
        subroutine qsort(array, elem_count, elem_size, compare) bind(C, name="qsort")
            import
            type(c_ptr), value :: array
            integer(c_size_t), value :: elem_count
            integer(c_size_t), value :: elem_size
            type(c_funptr), value   :: compare !int(*compare)(const void *, const void *)
        end subroutine qsort

    end interface


contains

!-------------------------------------------------
subroutine get_obs_def(obs, obs_def)

type(obs_type),     intent(in)  :: obs
type(obs_def_type), intent(out) :: obs_def

! WARNING: NEED TO DEFINE A COPY ROUTINE FOR OBS_DEF !!!
call copy_obs_def(obs_def, obs%def)

end subroutine get_obs_def

!-------------------------------------------------
subroutine set_obs_def(obs, obs_def)

type(obs_type),     intent(inout) :: obs
type(obs_def_type), intent(in)    :: obs_def

call copy_obs_def(obs%def, obs_def)

end subroutine set_obs_def

!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine dist_obs_set(set, new_set, num_obs, num_values, nprocs, root, start_proc)
    type(obs_type), allocatable,            intent(inout)      :: set(:)
    type(obs_type), allocatable,           intent(inout)      :: new_set(:)
    integer,                    intent(inout)      :: num_obs
    integer,                    intent(inout)      :: num_values
    integer,                    intent(in)      :: nprocs
    integer,                    intent(in)      :: root, start_proc

    type(obs_type_send), allocatable               :: conv_set(:), new_conv_set(:)
    type(obs_type_send), allocatable               :: all_conv_set(:)
    ! real(r8), allocatable                          :: values(:)
    ! real(r8), allocatable                          :: qc(:)
    ! real(r8), allocatable                          :: all_values_qc(:)
    ! real(r8), allocatable                          :: values_qc(:)
    type(obs_values_qc_type), allocatable          :: all_values_qc(:)
    type(obs_values_qc_type), allocatable          :: values_qc(:), new_values_qc(:)
    integer                                        :: total_values, all_values, vals_per_proc, total_obs, obs_per_proc, rem, &
        gather_procs, gather_obs_per_proc, gather_vals_per_proc
    integer                                        :: obs_mpi, ierror, i, d, diff, j, l, actual_proc, gather_proc, &
        vals_mpi
    integer, allocatable                           :: disp(:), count(:), disp_vals(:), count_vals(:)


    ! calculate information about number of processes retrieving obs based on what was provided in arguments
    gather_procs = nprocs - start_proc 
    gather_proc = my_task_id() - start_proc
    gather_obs_per_proc = num_obs / gather_procs
    obs_per_proc = num_obs / nprocs
    rem = modulo(num_obs, gather_procs)
    total_values = num_values * num_obs
    if (gather_proc < rem) gather_obs_per_proc = gather_obs_per_proc + 1
    if (my_task_id() == root) then 
        allocate(all_conv_set(num_obs))
        allocate(all_values_qc(total_values))
    else
        allocate(all_conv_set(1))
        allocate(all_values_qc(1))
    endif

    gather_vals_per_proc = num_values * gather_obs_per_proc
    allocate(conv_set(gather_obs_per_proc))
    allocate(values_qc(gather_vals_per_proc))
    
    ! Convert to a struct with fully contiguous memory 
    call convert_obs_set(set, conv_set, values_qc, num_values, gather_obs_per_proc)

    call deallocate_obs_set(set, obs_per_proc + 1, num_values)

    ! Setup the dedicated data structure
    ! already performed by initialize_obs_window
    ! call setup_obs_mpi(obs_mpi, vals_mpi)
    
    ! retrieve the observations and calculate displacement
    call gather_obs_varied(all_conv_set, conv_set, all_values_qc, values_qc, odt%obs_mpi, odt%val_mpi, num_values, num_obs, & 
        nprocs)

    if (my_task_id() == root) then 
        ! sort observations in linked list traversal order using quicksort
        print *, 'after gather'
        call sort_obs_send_by_time(all_conv_set, all_values_qc, num_values, num_obs)
        print *, 'sorted by time (timestamp added)'
        ! call sort_roundrobin_inplace(all_conv_set, all_values_qc, num_obs, num_values, nprocs)
        ! print *, 'sorted with roundrobin dist'
    endif

    deallocate(conv_set)
    deallocate(values_qc)
    call mpi_barrier(MPI_COMM_WORLD, ierror)

    ! scatter observations back to all processes sorted
    rem = modulo(num_obs, nprocs)
    if (my_task_id() < rem) obs_per_proc = obs_per_proc + 1
    vals_per_proc = obs_per_proc * num_values
    allocate(new_conv_set(obs_per_proc))
    allocate(new_values_qc(vals_per_proc))
    call scatter_obs_varied(new_conv_set, all_conv_set, new_values_qc, all_values_qc, obs_mpi, vals_mpi, num_values, num_obs, &
        nprocs)

    ! call deallocate_obs_set(new_set, 1, num_values)
    call allocate_obs_set(new_set, obs_per_proc, num_values)

    call mpi_barrier(MPI_COMM_WORLD, ierror)

    ! all procs convert their sets back to packed variants
    call convert_obs_back(new_set, new_conv_set, new_values_qc, obs_per_proc, num_values)

    ! After send has completed, destroy the mpi struct datatype and deallocate memory
    ! call destroy_obs_mpi(obs_mpi, vals_mpi)
    deallocate(all_conv_set)
    deallocate(all_values_qc)

end subroutine dist_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
integer function get_obs_offset(key)
    integer,                    intent(in)      :: key
    integer                                     :: key_loc
    integer                                     :: offset_multiple
    integer                                     :: obs_excl_rem
    integer                                     :: scnd_offset
    integer                                     :: offset_in_block
    integer                                     :: orig_start
    ! if the remainder observations are scattered in roundrobin order, 
    ! our current assumptions regarding key-offset correspondence are incorrect.
    ! This function returns the correct offset.

    ! If key = -1, get outta here!
    if (key == -1) then
        get_obs_offset = -1
        return
    endif

    ! If obs is one of the remainder
    obs_excl_rem = odt%total_obs - odt%rem
    if (key > obs_excl_rem) then
        key_loc = (key - obs_excl_rem) * (odt%num_obs_per_proc + 1)
        get_obs_offset = key_loc
        return
    endif

    ! If obs is a component of one of the blocks
    scnd_offset = odt%num_obs_per_proc + 1
    offset_multiple = (key - 1) / odt%num_obs_per_proc
    offset_in_block = modulo(key, odt%num_obs_per_proc)
    orig_start = odt%num_obs_per_proc * odt%rem

    key_loc = (scnd_offset * offset_multiple) + offset_in_block
    if (offset_multiple < odt%rem) then
       key_loc = (scnd_offset * offset_multiple) + (key - (offset_multiple * odt%num_obs_per_proc))
    else
       key_loc = (scnd_offset * odt%rem) + (key - orig_start)
    endif
   get_obs_offset = key_loc

    ! get_obs_offset = key_loc
end function get_obs_offset
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine initialize_obs_window(buffer, num_obs_per_proc, num_vals_per_obs, total_obs, rem, num_alloc, dist_type, nprocs)
    ! this will initialize the obs window for one-sided communication
    ! or distributed communication generally
    type(obs_type),          intent(inout)      :: buffer(:)
    integer,                 intent(in)         :: num_obs_per_proc
    integer,                 intent(in)         :: total_obs
    integer,                 intent(in)         :: rem
    integer,                 intent(in)         :: num_alloc
    integer,                 intent(in)         :: dist_type
    integer,                 intent(in)         :: nprocs
    integer                                     :: ierror
    integer                                     :: num_vals

    ! allocate our buffers to be number of obs on this process
    num_vals = num_obs_per_proc * num_vals_per_obs
    if (dist_type == 1) then
        allocate(odt%obs_buf(num_alloc))
        allocate(odt%val_buf(num_alloc*num_vals_per_obs))
    endif

    ! set important variables related to observation distribution
    odt%num_obs_per_proc = num_obs_per_proc
    odt%num_vals_per_proc = num_vals
    odt%num_vals_per_obs = num_vals_per_obs
    odt%total_obs = total_obs
    odt%rem = rem
    odt%my_pe = my_task_id()
    odt%nprocs = nprocs
    odt%ngets = 0
    odt%mpi_time = 0.0
    odt%dist_type = dist_type

    ! set the number of obs associated with our process specifically
    odt%our_num_obs = odt%num_obs_per_proc
    if (odt%my_pe < odt%rem) then
        odt%our_num_obs = odt%our_num_obs + 1
    endif

    ! setup our datatypes
    call setup_obs_mpi(odt%obs_mpi, odt%val_mpi)

    ! convert to sendable datatype
    ! also create mpi window of observation memory
    if (dist_type == 1) then
        call convert_obs_set(buffer, odt%obs_buf, odt%val_buf, num_vals_per_obs, num_alloc)

        ! create windows
        call mpi_win_create(odt%obs_buf, num_alloc * sizeof(odt%obs_buf(1)), sizeof(odt%obs_buf(1)), MPI_INFO_NULL, MPI_COMM_WORLD, &
        odt%obs_win, &
        ierror)

        call mpi_win_create(odt%val_buf, num_alloc * num_vals_per_obs * sizeof(odt%val_buf(1)), sizeof(odt%val_buf(1)), MPI_INFO_NULL, MPI_COMM_WORLD, &
        odt%val_win, &
        ierror)
    endif

end subroutine initialize_obs_window
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_obs_loc_info(key, obs_pe, obs_offset, val_offset)
    integer,            intent(in)          :: key
    integer,            intent(inout)       :: obs_pe
    integer(KIND=MPI_ADDRESS_KIND),            intent(inout)       :: obs_offset
    integer(KIND=MPI_ADDRESS_KIND),            intent(inout)       :: val_offset

    ! determine observation location based on whether observation is component of remainder
    if (key > odt%total_obs - odt%rem) then
        ! print *, 'different calculations here'
        obs_pe = modulo(key, odt%num_obs_per_proc) - 1
        obs_offset = odt%num_obs_per_proc
        val_offset = obs_offset * odt%num_vals_per_obs
    else
        obs_pe = (key - 1) / odt%num_obs_per_proc
        obs_offset = (modulo((key - 1), odt%num_obs_per_proc))
        val_offset = obs_offset * odt%num_vals_per_obs
    endif

end subroutine get_obs_loc_info
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_obs_dist(key, obs)
    ! Used to get a single observation using one key; 
    ! requires locking before retrieving every observation; not efficient
    integer,                intent(in)              :: key
    type(obs_type),         intent(inout)           :: obs
    ! integer,                intent(inout)           :: obs
    type(obs_type_send),allocatable                 :: obs_buffer(:)
    type(obs_values_qc_type),allocatable            :: vals_buffer(:)
    integer                                         :: val_pos, obs_pos, obs_pe, rem_proc
    integer(kind=MPI_ADDRESS_KIND)                  :: obs_offset, val_offset
    integer                                         :: ierror
    real(r8)                                        :: start, end
    type(obs_type)                                  :: obs_arr(1)

    allocate(obs_buffer(1))
    allocate(vals_buffer(odt%num_vals_per_obs))
    ! todo: also determine whether the obs we are looking for is on another process
    ! if it is not, we do not need to perform a one-sided comm

    ! determine where the observation is located
    call get_obs_loc_info(key, obs_pe, obs_offset, val_offset)

    if (obs_pe == odt%my_pe) then
        ! obs_buffer(1) = odt%obs_buf(obs_offset + 1)
        obs_buffer(1) = odt%obs_buf(obs_offset + 1)
        vals_buffer(1:odt%num_vals_per_obs) = odt%val_buf(val_offset+1:val_offset+odt%num_vals_per_obs)
    else
        ! lock window of process
        call mpi_win_lock(MPI_LOCK_SHARED, obs_pe, MPI_MODE_NOCHECK, odt%obs_win, ierror)
        call mpi_win_lock(MPI_LOCK_SHARED, obs_pe, MPI_MODE_NOCHECK, odt%val_win, ierror)

        ! retrieve observation
        call mpi_get(obs_buffer, 1, odt%obs_mpi, obs_pe, obs_offset, 1, odt%obs_mpi, odt%obs_win, ierror)
        call mpi_get(vals_buffer, odt%num_vals_per_obs, odt%val_mpi, obs_pe, val_offset, odt%num_vals_per_obs, odt%val_mpi, odt%val_win, ierror)

        ! unlock window of process we're retrieving from
        call mpi_win_unlock(obs_pe, odt%obs_win, ierror)
        call mpi_win_unlock(obs_pe, odt%val_win, ierror)
    endif

    ! convert observation back to packed form
    call convert_obs_back(obs_arr, obs_buffer, vals_buffer, 1, odt%num_vals_per_obs)
    obs = obs_arr(1)
    deallocate(obs_buffer)
    deallocate(vals_buffer)

end subroutine get_obs_dist
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine dbg_print(str)
    character(len=*),   intent(in)      :: str
    if (odt%my_pe == 0) print *, str
end subroutine dbg_print
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine samplesort_obs(perc)
    integer,            intent(in)      :: perc ! what percentage of obs sequence will be our sample?

    integer                             :: all_sample_num
    integer                             :: our_sample_num
    integer                             :: sample_cnt
    integer                             :: j, i, per_proc, l, is_last, m, check
    integer(i8)                         :: k
    real                                :: start_random
    integer                             :: rand_idx
    integer                             :: our_num_obs
    integer                             :: new_obs_num
    integer                             :: new_vals_num
    integer, allocatable                :: bucket_cnt(:)
    integer, allocatable                :: bucket_disp(:)
    integer, allocatable                :: val_disp(:)
    integer, allocatable                :: new_cnt(:)
    integer, allocatable                :: new_disp(:)
    integer, allocatable                :: new_val_disp(:)
    integer, allocatable                :: val_qc_rdisp(:)
    integer, allocatable                :: val_qc_cnt(:)
    integer, allocatable                :: new_val_qc_cnt(:)
    integer, allocatable                :: val_qc_disp(:)
    integer(i8), allocatable,target     :: all_samples(:)
    integer(i8), allocatable            :: our_samples(:)
    integer, allocatable                :: is_selected(:)
    integer(i8), allocatable            :: scnd_selection(:)
    type(obs_type_send), pointer        :: new_obs_set(:) => NULL()
    type(obs_values_qc_type), pointer   :: new_val_qc(:) => NULL()
    type(obs_type), allocatable, target :: obs_set_sort(:)
    integer(i8)                         :: curr_time
    integer                             :: iden_vals, iden_rem

    our_num_obs = odt%num_obs_per_proc
    if (odt%my_pe < odt%rem) then
        our_num_obs = our_num_obs + 1
    endif

    ! Convert back to packed to make sorting obs easier
    if (odt%my_pe == 0) print *, 'Converting observations to packed set'
    call allocate_obs_set(obs_set_sort, our_num_obs, odt%num_vals_per_obs)
    call convert_obs_back(obs_set_sort, odt%obs_buf, odt%val_buf, our_num_obs, odt%num_vals_per_obs)

    ! number of samples to retrieve / send to first process
    ! note: if no parameter is specified, assume 1%
    all_sample_num = (perc * 0.01) * odt%total_obs

    ! allocate memory for all samples (only first process)

    ! calculate our number of samples and allocate an array 
    ! also: get a random number (somehow)
    our_sample_num = all_sample_num / odt%nprocs
    allocate(our_samples(our_sample_num))
    allocate(is_selected(our_num_obs))
    is_selected(1:our_num_obs) = 0

    if (odt%my_pe == 0) print *, 'allocating memory for samples'

    all_sample_num = our_sample_num * odt%nprocs
    if (odt%my_pe == 0) then
        allocate(all_samples(all_sample_num))
    else
        allocate(all_samples(1))
    endif

    if (odt%my_pe == 0) print *, 'randomly selecting obs'
    sample_cnt = 0
    do while (sample_cnt < our_sample_num)
        call random_number(start_random) ! random enough? assume yes for now (until everything explodes)
        rand_idx = floor(start_random * our_num_obs) + 1
        ! if (odt%my_pe == 0) print *, rand_idx
        if (is_selected(rand_idx) == 0) then
            is_selected(rand_idx) = 1
            sample_cnt = sample_cnt + 1
        else
            cycle
        endif
    enddo

    ! todo: also need to sort our obs seq before sending to first process
    
    ! sort observations in time order
    ! this doesn't work in parallel; for now assume that it does
    ! TODO: this **must** be fixed before testing or the world **will** end

    ! would this even work?! this is so cursed...
    if (odt%my_pe == 0) print *, 'sorting observations in time order'
    ! if (odt%my_pe == 0) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 1) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 24) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 78) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 0) print *, obs_set_sort(500)%time_actual
    call qsort(c_loc(obs_set_sort(1)), int(odt%our_num_obs, c_size_t), sizeof(obs_set_sort(1)), c_funloc(compare_time_types_alt))

    ! convert back so our unpacked obs matches the sorting of our packed observations
    call convert_obs_set(obs_set_sort, odt%obs_buf, odt%val_buf, odt%num_vals_per_obs, odt%our_num_obs)
    ! if (odt%my_pe == 0) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 1) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 24) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 78) print *, obs_set_sort(1)%time_actual
    ! if (odt%my_pe == 0) print *, obs_set_sort(500)%time_actual
    ! call sort_obs_send_by_time(odt%obs_buf, odt%val_buf, odt%num_vals_per_obs, our_num_obs)

    call dbg_print('test') ! make sure my debug print works
    ! force a leave; we need to test everything happening earlier first!
    ! return
    call dbg_print('getting samples')
    j = 1
    do i = 1, our_num_obs
        if (is_selected(i) == 1) then
            curr_time = obs_set_sort(i)%time_actual
            our_samples(j) = curr_time
            j = j + 1
            ! our_samples(j) = odt%obs_buf%
        endif
    enddo
    ! if (odt%my_pe == 0) print *, j
    ! if (odt%my_pe == 0) print *, our_sample_num
    ! if (odt%my_pe == 0) print *, 'PE: ', odt%my_pe, our_samples(1)
    ! if (odt%my_pe == 1) print *, 'PE: ', odt%my_pe, our_samples(1)
    ! if (odt%my_pe == 24) print *, 'PE: ', odt%my_pe, our_samples(1)
    ! if (odt%my_pe == 78) print *, 'PE: ', odt%my_pe, our_samples(1)
    ! if (odt%my_pe == 124) print *, 'PE: ', odt%my_pe, our_samples(1)

    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! perform gather to retrieve samples 

    call dbg_print('gathering all samples')
    call mpi_gather(our_samples, our_sample_num, MPI_INTEGER8, all_samples, our_sample_num, MPI_INTEGER8, 0, MPI_COMM_WORLD, &
    odt%ierror)

    ! if first process sort time samples using qsort and select p - 1 samples from this set
    ! allocate for all processes to store bucket vars
    call dbg_print('sorting samples and selecting new ones')
    allocate(scnd_selection(odt%nprocs - 1))
    if (odt%my_pe == 0) then
        per_proc = our_sample_num
        ! print *, all_samples(1)
        call qsort(c_loc(all_samples(1)), int(all_sample_num, c_size_t), sizeof(all_samples(1)), c_funloc(compare_large_int))
        ! print *, all_samples(1)
        ! do i = 1, all_sample_num
        !     print *, all_samples(i)
        ! enddo
        do i = 1, odt%nprocs - 1
            scnd_selection(i) = all_samples(i * per_proc)
            print *, scnd_selection(i)
        enddo
        ! print *, 'next steps...'
        ! do i = 1, odt%nprocs - 1
        !     print *, scnd_selection(i)
        ! enddo
        ! print *, scnd_selection(1)
        ! print *, scnd_selection(2)
        ! print *, scnd_selection(3)
    endif

    call dbg_print('broadcasting new samples')
    call mpi_bcast(scnd_selection, odt%nprocs - 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, odt%ierror)

    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! return
    ! Now we have buckets in which to sort elements
    ! NOTE: maybe also use these buckets as indicators as to which processes are likely to have which elements
    ! storing / duplicating such information should take a minimal amount of memory (8 bytes (time) * nprocs)
    ! Use binary search to select bucket for each element (O(nlogn))
    ! Actually...maybe not; elements are already sorted

    ! Currently assume we have no duplicate keys
    ! TODO: We'll handle that case later
    allocate(bucket_disp(odt%nprocs))
    allocate(bucket_cnt(odt%nprocs))
    check = 0
    k = odt%obs_buf(1)%time_actual
    bucket_disp(1) = 0
    l = 1
    do i = 1, odt%nprocs
        j = 0
        ! if (i < odt%nprocs) then
        !     do while (l <= our_sample_num .and. k < scnd_selection(i))
        !         j = j + 1
        !         l = l + 1
        !         if (l <= our_sample_num) k = our_samples(l)
        !     enddo
        ! else if (i == odt%nprocs) then
        !     do while (l <= our_sample_num)
        !         j = j + 1
        !         l = l + 1
        !         if (l <= our_sample_num) k = our_samples(l)
        !     enddo
        ! endif
        do while (l <= our_num_obs)
            ! do while ((i < odt%nprocs .and. k < scnd_selection(i)) .or. (i == odt%nprocs .and. k >= scnd_selection(i - 1)))
            do
                ! print *, 'l: ' , l
                ! print *, 'our_sample_num: ' , our_sample_num
                if (i < odt%nprocs) then
                    ! if (k < scnd_selection(i)) then
                    !     check = 1
                    ! else
                    !     exit
                    ! endif
                    ! print *, 'banana'
                    if (k >= scnd_selection(i)) exit
                endif
                if (i == odt%nprocs) then
                    ! if (k >= scnd_selection(i - 1)) then
                    !     check = 1
                    ! else
                    !     exit
                    ! endif
                    ! print *, 'apple'
                    if (k < scnd_selection(i - 1)) exit
                endif
                if (i == odt%nprocs .and. odt%my_pe == 1) then
                    print *, 'last scnd_selection: ', scnd_selection(i-1)
                    print *, 'l: ', l
                    print *, 'k: ', k
                endif
                j = j + 1
                l = l + 1
                ! print *, 'why????'
                if (l <= our_num_obs) k = odt%obs_buf(l)%time_actual
                if (l == our_num_obs + 1) exit
            enddo
            ! print *, 'honk' 
            exit
        enddo
        ! if (i < odt%nprocs) then
        !     do while (l <= our_sample_num .and. k < scnd_selection(i))
        !         j = j + 1
        !         l = l + 1
        !         if (l <= our_sample_num) k = our_samples(l)
        !     enddo
        !     if (l > our_sample_num) exit
        ! else
        !     ! print *, 'alt case where we are at last process'
        !     ! this is *not* correct; must fix
        !
        !     j = odt%total_obs - bucket_disp(odt%nprocs) 
        ! endif
        bucket_cnt(i) = j
        if (i /= odt%nprocs) bucket_disp(i + 1) = j + bucket_disp(i)
        ! print *, 'do something else lol'
    enddo
    if (odt%my_pe == 1) print *, 'last bucket_cnt: ', bucket_cnt(odt%nprocs)

    
    ! if (odt%my_pe == 0) then
    !     print *, 'bucket_cnt'
    !     do i = 1, odt%nprocs
    !         print *, bucket_cnt(i)
    !     enddo
    !     print *, 'bucket_disp'
    !     do i = 1, odt%nprocs
    !         print *, bucket_disp(i)
    !     enddo
    ! endif

    ! print *, odt%my_pe, 'reached barrier'
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    ! return 

    ! check for identical sets of keys
    ! If a bucket has no elements, that likely means the previous key engulfed all elements
    ! split the elements placed into that bucket across all buckets that have the same key 
    ! i = 2
    ! do while (i <= odt%nprocs)
    ! !do i = 1, odt%nprocs
    !     ! print *, 'From process', odt%my_pe, ': i = ', i
    !     m = 0
    !     is_last = 0
    !     if (bucket_cnt(i) == 0) then
    !         ! print *, 'haaah'
    !         if (i == odt%nprocs) then
    !             is_last = 0
    !         else if (scnd_selection(i) == scnd_selection(i - 1)) then
    !             is_last = 1
    !         endif
    !         if (is_last) then
    !             j = i
    !             k = 0
    !             do while (j + k < odt%nprocs .and. m == 0) 
    !                 if (bucket_cnt(j + k) == 0 .and. scnd_selection(j+k) == scnd_selection(i-1)) then
    !                     k = k + 1
    !                 else
    !                     m = 1
    !                 endif
    !             enddo
    !             k = k + 1
    !             iden_vals = bucket_cnt(i - 1) / k
    !             iden_rem = modulo(bucket_cnt(i - 1), k)
    !             l = 0
    !             do while (l < k)
    !                 bucket_cnt((i-1) + l) = iden_vals
    !                 if (l < iden_rem) then
    !                     bucket_cnt((i-1) + l) = bucket_cnt((i - 1) + l) + 1
    !                 endif
    !                 l = l + 1
    !             enddo
    !             i = j + k - 1
    !         else
    !             i = i + 1
    !         endif
    !     else
    !         i = i + 1
    !     endif
    ! enddo

    ! process i in new_cnt(i): how many elements are we receiving?
    allocate(new_cnt(odt%nprocs))

    ! displacement of the elements in the buffer for process i
    ! process(i) = where on process(i) do the elements start?
    allocate(new_disp(odt%nprocs))
    ! allocate(val_qc_disp(odt%nprocs))

    ! displacement of values on each process i 
    allocate(new_val_disp(odt%nprocs))

    ! allocate(val_qc_rdisp(odt%nprocs))
    allocate(new_val_qc_cnt(odt%nprocs))
    ! number of values to be sent to each process from *our* process
    allocate(val_qc_cnt(odt%nprocs))
    ! displacement of values/qc in buffer for *our* process
    allocate(val_qc_disp(odt%nprocs))

    ! print *, 'Process ', odt%my_pe, ' reached barrier #2'
    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! return
    ! alltoall on the displacement and count vars
    call dbg_print('attempting to run alltoall #1')

    call mpi_alltoall(bucket_cnt, 1, MPI_INTEGER, new_cnt, 1, MPI_INTEGER, MPI_COMM_WORLD, &
        odt%ierror)

    call dbg_print('attempting to run alltoall #2')

    ! call mpi_alltoall(bucket_disp, 1, MPI_INTEGER, new_disp, 1, MPI_INTEGER, MPI_COMM_WORLD, &
    !     odt%ierror)

    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! return
    ! alltoallv on the elements themselves
    ! (probably the most expensive component of this algorithm)

    ! note: need to modify count and disp for values/qc since each
    ! observation may have more than one
    val_qc_cnt(1:odt%nprocs) = bucket_cnt(1:odt%nprocs)
    val_qc_disp(1:odt%nprocs) = bucket_disp(1:odt%nprocs)

    val_qc_cnt(1:odt%nprocs) = val_qc_cnt(1:odt%nprocs) * odt%num_vals_per_obs

    ! determine displacement for bucket values
    val_qc_disp(1) = 0 
    do i = 2, odt%nprocs
        val_qc_disp(i) = val_qc_disp(i - 1) + val_qc_cnt(i - 1) 
    enddo

    call dbg_print("attempting to run alltoall #3")
    ! send new val_qc displacement to every process
    ! call mpi_alltoall(val_qc_disp, 1, MPI_INTEGER, new_val_disp, 1, MPI_INTEGER, MPI_COMM_WORLD, &
    !    odt%ierror)

    ! return

    ! determine the num of observations we will receive
    new_obs_num = 0
    do i = 1, odt%nprocs
        new_obs_num = new_obs_num + new_cnt(i)
    enddo

    ! calculate rdisp using cnt previously sent using alltoall
    ! val_qc_rdisp(1) = 0
    ! do i = 2, odt%nprocs
    !     val_qc_rdisp(i) = val_qc_rdisp(i - 1) + new_cnt(i - 1)
    ! enddo


    new_disp(1) = 0
    do i = 2, odt%nprocs
        new_disp(i) = new_disp(i - 1) + new_cnt(i - 1) ! represents sdisp
    enddo

    new_val_qc_cnt(1:odt%nprocs) = new_cnt(1:odt%nprocs) * odt%num_vals_per_obs
    new_vals_num = new_obs_num * odt%num_vals_per_obs

    new_val_disp(1) = 0
    do i = 2, odt%nprocs
        new_val_disp(i) = new_val_disp(i - 1) + new_val_qc_cnt(i - 1)
    enddo

    allocate(new_obs_set(new_obs_num))
    allocate(new_val_qc(new_vals_num))

    ! sdisp: displacement from send buffer
    ! rdisp: where in receive buffer elems will be stored
    ! if (odt%my_pe == 1) then
    !     print *, 'val_qc_cnt' 
    !     do i = 1, odt%nprocs
    !         print *, val_qc_cnt(i)
    !     enddo
    !     print *, 'val_qc_disp' 
    !     do i = 1, odt%nprocs
    !         print *, val_qc_disp(i)
    !     enddo
    !     print *, 'obs_cnt' 
    !     do i = 1, odt%nprocs
    !         print *, bucket_cnt(i)
    !     enddo
    !     print *, 'obs_disp' 
    !     do i = 1, odt%nprocs
    !         print *, bucket_disp(i)
    !     enddo
    !     print *, 'val_qc_disp' 
    !     do i = 1, odt%nprocs
    !         print *, val_qc_disp(i)
    !     enddo
    !     print *, 'new_cnt: ' 
    !     do i = 1, odt%nprocs
    !         print *, new_val_qc_cnt(i)
    !     enddo
    !     print *, 'new_cnt: ' 
    !     do i = 1, odt%nprocs
    !         print *, new_val_qc_cnt(i)
    !     enddo
    !     print *, 'new_disp: ' 
    !     do i = 1, odt%nprocs
    !         print *, new_val_disp(i)
    !     enddo
    !     print *, 'new_cnt(obs): ' 
    !     do i = 1, odt%nprocs
    !         print *, new_cnt(i)
    !     enddo
    !     print *, 'new_disp(obs): ' 
    !     do i = 1, odt%nprocs
    !         print *, new_disp(i)
    !     enddo
    ! endif
    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! if (odt%my_pe == 0) print *, 'total_obs: ', odt%total_obs
    ! call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! print *, 'pe: ', odt%my_pe, 'new_obs_num: ', new_obs_num
    ! call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! print *, 'pe: ', odt%my_pe, 'our_num_obs: ', our_num_obs
    ! call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! return 

    call dbg_print('attempting alltoallv #1')
    ! alltoallv both the observations and the values
    call mpi_alltoallv(odt%obs_buf, bucket_cnt, bucket_disp, odt%obs_mpi, new_obs_set, new_cnt, new_disp, odt%obs_mpi, &
        MPI_COMM_WORLD, odt%ierror)

    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! return

    call dbg_print('attempting alltoallv #2')

    call mpi_alltoallv(odt%val_buf, val_qc_cnt, val_qc_disp, odt%val_mpi, new_val_qc, new_val_qc_cnt, new_val_disp, odt%val_mpi, &
        MPI_COMM_WORLD, odt%ierror)

    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)

    ! return

    call dbg_print('heehaw')
    ! print *, 'pe: ', odt%my_pe, 'new_obs_num: ', new_obs_num
    ! call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! print *, 'pe: ', odt%my_pe, 'our_num_obs: ', our_num_obs
    ! if (odt%my_pe == 1) call print_obs_send(new_obs_set(112324))
    ! call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    call dbg_print("hahahahahahahahha")
    ! return
    ! Final sort of the observations on each process
    ! call convert_obs_set(obs_set_sort, new_obs_set, new_val_qc, odt%num_vals_per_obs, new_obs_num)
        ! print *, 'new_obs_num: ', new_obs_num
    call deallocate_obs_set(obs_set_sort, our_num_obs, odt%num_vals_per_obs) 
    call allocate_obs_set(obs_set_sort, new_obs_num, odt%num_vals_per_obs) 
    call convert_obs_back(obs_set_sort, new_obs_set, new_val_qc, new_obs_num, odt%num_vals_per_obs)
    call dbg_print("hahahahahahahahha")
    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! return
    call qsort(c_loc(obs_set_sort(1)), int(new_obs_num, c_size_t), sizeof(obs_set_sort(1)), c_funloc(compare_time_types_alt))
    call convert_obs_set(obs_set_sort, new_obs_set, new_val_qc, odt%num_vals_per_obs, new_obs_num)
    
    ! set our global pointers to the new sets and deallocate the old sets
    deallocate(odt%obs_buf)
    deallocate(odt%val_buf)
    odt%obs_buf => new_obs_set
    odt%val_buf => new_val_qc

    call dbg_print('made it to the end')
    call mpi_barrier(MPI_COMM_WORLD, odt%ierror)
    ! 1. select set of samples from every process's observation sequences
    !    (1% of the total observation sequence)
    ! 2. every process sends its sample set to first process
    ! 3. first process sorts samples using qsort
    ! 4. first process selects p - 1 samples from this set
    ! 5. p - 1 samples broadcasted to every process
    ! 6. elements are placed into their respective buckets
    ! 7. count and displacement sent using alltoall
    ! 8. actual elements sent using alltoallv

    ! Idea: use samplesort to sort observations into time order
    ! best case: will provide a faster, more memory efficient sort
    ! worst case: everything breaks (yay! fun!)

end subroutine samplesort_obs
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_all_obs_contiguous(obs)
    ! What happens if we attempt to retrieve all contiguous observations in obs_sequence?
    ! rather than one get for each observation?
    ! maybe would reduce startup cost?
    type(obs_type),         intent(inout)           :: obs(:)
    ! integer,                intent(inout)           :: obs
    type(obs_type_send),allocatable                 :: obs_buffer(:)
    type(obs_values_qc_type),allocatable            :: vals_buffer(:)
    integer                                         :: val_pos, obs_pos, obs_pe, rem_proc, total_values
    integer(kind=MPI_ADDRESS_KIND)                  :: obs_offset, val_offset
    integer                                         :: ierror, i, d, get_cnt, num_retrieve, val_retrieve
    integer                                         :: scnd_num_retrieve
    integer                                         :: scnd_val_retrieve
    real(r8)                                        :: stime, etime, ttime
    ! type(obs_type)                                  :: obs_arr(1)

    get_cnt = 0
    total_values = odt%num_vals_per_obs * odt%total_obs
    allocate(obs_buffer(odt%total_obs))
    allocate(vals_buffer(total_values))
    ! todo: also determine whether the obs we are looking for is on another process
    ! if it is not, we do not need to perform a one-sided comm

    ! lock window of all processes
    stime = mpi_wtime()
    call mpi_win_lock_all(MPI_MODE_NOCHECK, odt%obs_win, ierror)
    call mpi_win_lock_all(MPI_MODE_NOCHECK, odt%val_win, ierror)


    num_retrieve = 0
    val_retrieve = 0

    obs_offset = 0
    val_offset = 0

    print *, 'odt%num_obs_per_proc: ', odt%num_obs_per_proc
    print *, 'odt%nprocs: ', odt%nprocs
    scnd_num_retrieve = odt%num_obs_per_proc
    obs_pos = 1
    val_pos = 1
    do i = 0, odt%nprocs - 1
        ! print *, 'i: ', i
        ! d = 3
        scnd_num_retrieve = odt%num_obs_per_proc
        scnd_val_retrieve = odt%num_obs_per_proc * odt%num_vals_per_obs
        ! scnd_val_retrieve = 1

        if (i < odt%rem) then
            scnd_num_retrieve = scnd_num_retrieve + 1
            scnd_val_retrieve = scnd_val_retrieve + odt%num_vals_per_obs
        endif

        if (i == odt%my_pe) then
            obs_buffer(obs_pos:obs_pos+scnd_num_retrieve-1) = odt%obs_buf(1:scnd_num_retrieve)
            vals_buffer(val_pos:val_pos+scnd_val_retrieve-1) = odt%val_buf(1:scnd_val_retrieve)
        else
            call mpi_get(obs_buffer(obs_pos:obs_pos+scnd_num_retrieve-1), scnd_num_retrieve, odt%obs_mpi, i, obs_offset, scnd_num_retrieve, odt%obs_mpi, odt%obs_win, ierror)
            ! call mpi_get(vals_buffer(val_pos:val_pos+scnd_val_retrieve-1), scnd_val_retrieve, odt%val_mpi, i, 0, scnd_val_retrieve, &
            !     odt%val_mpi, odt%val_win, ierror)
            call mpi_get(vals_buffer(val_pos:val_pos+scnd_val_retrieve-1), scnd_val_retrieve, odt%val_mpi, i, val_offset, &
            scnd_val_retrieve, odt%val_mpi, odt%val_win, ierror)
        endif

        obs_pos = obs_pos + scnd_num_retrieve
        val_pos = val_pos + scnd_val_retrieve
    enddo

    ! unlock window of all processes
    call mpi_win_unlock_all(odt%obs_win, ierror)
    call mpi_win_unlock_all(odt%val_win, ierror)
    etime = mpi_wtime()

    print *, 'finished get'
    print *, 'time to get: ', etime - stime
    print *, 'converting back'
    ! convert observation(s) back to packed form
    call convert_obs_back(obs, obs_buffer, vals_buffer, odt%total_obs, odt%num_vals_per_obs)

    ! print *, 'Average time to get: ', ttime / get_cnt
    deallocate(obs_buffer)
    deallocate(vals_buffer)

end subroutine get_all_obs_contiguous

!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine get_obs_set(keys, obs, num_keys)
    ! Retrieve a set of observations simultaneously, rather than one observation at a time
    ! far faster than retrieving observations one key at a time
    ! however, speed varies depending on whether obs are on or off-node
    ! mpi_get takes longer to startup when retrieving off-node
    integer,                intent(in)              :: keys(:)
    integer,                intent(in)              :: num_keys
    type(obs_type),         intent(inout)           :: obs(:)
    ! integer,                intent(inout)           :: obs
    type(obs_type_send),allocatable                 :: obs_buffer(:)
    type(obs_values_qc_type),allocatable            :: vals_buffer(:)
    integer                                         :: val_pos, obs_pos, obs_pe, rem_proc, total_values
    integer(kind=MPI_ADDRESS_KIND)                  :: obs_offset, val_offset
    integer                                         :: ierror, i, d, get_cnt
    real(r8)                                        :: stime, etime, ttime
    ! type(obs_type)                                  :: obs_arr(1)

    get_cnt = 0
    total_values = odt%num_vals_per_obs * num_keys
    allocate(obs_buffer(num_keys))
    allocate(vals_buffer(total_values))
    ! todo: also determine whether the obs we are looking for is on another process
    ! if it is not, we do not need to perform a one-sided comm

    ! lock all windows in shared mode; this should be fine since obs sequence is not modified by any process
    stime = mpi_wtime()
    call mpi_win_lock_all(MPI_MODE_NOCHECK, odt%obs_win, ierror)
    call mpi_win_lock_all(MPI_MODE_NOCHECK, odt%val_win, ierror)

    ! Loop through all of the keys we would like to retrieve
    do i = 1, num_keys
        if (modulo(i, 1000000) == 0) print *, 'i: ', i
        d = (i-1) * odt%num_vals_per_obs
        call get_obs_loc_info(keys(i), obs_pe, obs_offset, val_offset)
        if (obs_pe == odt%my_pe) then
            obs_buffer(i) = odt%obs_buf(obs_offset + 1)
            vals_buffer(d+1:d+odt%num_vals_per_obs) = odt%val_buf(val_offset+1:val_offset+odt%num_vals_per_obs)
        else

            ! todo: maybe abstract timers away in separate type?
            ! also: provide option to enable/disable timers
            ! stime = MPI_WTime()
            call mpi_get(obs_buffer(i), 1, odt%obs_mpi, obs_pe, obs_offset, 1, odt%obs_mpi, odt%obs_win, ierror)
            ! etime = MPI_WTime()
            ! ttime = ttime + (etime - stime)

            ! stime = MPI_WTime()
            call mpi_get(vals_buffer(d+1:d+odt%num_vals_per_obs), odt%num_vals_per_obs, odt%val_mpi, obs_pe, val_offset, odt%num_vals_per_obs, odt%val_mpi, odt%val_win, ierror)
            ! etime = MPI_WTime()
            ! ttime = ttime + (etime - stime)

            get_cnt = get_cnt + 2
        endif
    enddo

    ! unlock all windows after retrieval has finished
    call mpi_win_unlock_all(odt%obs_win, ierror)
    call mpi_win_unlock_all(odt%val_win, ierror)
    etime = mpi_wtime()
    print *, 'time to get (get_obs_set): ', etime - stime

    ! convert observation(s) back to packed form
    call convert_obs_back(obs, obs_buffer, vals_buffer, num_keys, odt%num_vals_per_obs)

    ! print *, 'Average time to get: ', ttime / get_cnt
    ! obs = obs_arr(1)
    deallocate(obs_buffer)
    deallocate(vals_buffer)

end subroutine get_obs_set

!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine destroy_obs_window()
    ! this will destroy the obs window (it's in the name, chief)
    integer :: ierror
    call destroy_obs_mpi(odt%obs_mpi, odt%val_mpi)
    if (odt%dist_type == 1) then
        call mpi_win_free(odt%obs_win, ierror)
        call mpi_win_free(odt%val_win, ierror)
        deallocate(odt%val_buf)
        deallocate(odt%obs_buf)
    endif
end subroutine destroy_obs_window

!------------------------------------------------------------------

subroutine destroy_obs_mpi(mpi_obstype, mpi_valstype)
    integer,                    intent(inout)       :: mpi_obstype, mpi_valstype
    integer                                      :: ierror

    call mpi_type_free(mpi_obstype, ierror)
    call mpi_type_free(mpi_valstype, ierror)

end subroutine destroy_obs_mpi
!------------------------------------------------------------------
subroutine setup_obs_mpi(mpi_obstype, mpi_valstype)
    integer,                    intent(inout)    :: mpi_obstype
    integer,                    intent(inout)    :: mpi_valstype
    integer :: rank, nprocs, ierror, i, num_ints, num_vars, num_doubles, num_logicals, num_longs
    integer :: num_vars_vals, num_ints_vals, num_doubles_vals
    integer(MPI_ADDRESS_KIND) :: offsets(15), offsets_vals(4)
    integer(MPI_ADDRESS_KIND) :: address(15), address_vals(4)
    integer :: oldtypes(15), oldtypes_vals(4)
    integer :: bl_var(15), bl_var_vals(4)
    type(obs_type_send) :: initial
    type(obs_values_qc_type) :: init_val

    ! for obs_type_send
    num_ints = 10 
    num_doubles = 4
    num_logicals = 0
    num_longs = 1
    num_vars = 15
    bl_var(1:num_vars) = 1 

    ! for  obs_values_qc_type
    num_ints_vals = 1
    num_doubles_vals = 2
    num_longs_vals = 1
    num_vars_vals = 4
    bl_var_vals(1:num_vars_vals) = 1

    ! get the addresses of every variable; will let us calculate the displacement
    call mpi_get_address(initial, address(1), ierror)
    call mpi_get_address(initial%which_vert, address(2), ierror)
    call mpi_get_address(initial%kind, address(3), ierror)
    call mpi_get_address(initial%seconds, address(4), ierror)
    call mpi_get_address(initial%days, address(5), ierror)
    call mpi_get_address(initial%obs_def_key, address(6), ierror)
    call mpi_get_address(initial%prev_time, address(7), ierror)
    call mpi_get_address(initial%next_time, address(8), ierror)
    call mpi_get_address(initial%cov_group, address(9), ierror)
    call mpi_get_address(initial%time_order, address(10), ierror)
    call mpi_get_address(initial%time_actual, address(11), ierror)
    call mpi_get_address(initial%lon, address(12), ierror)
    call mpi_get_address(initial%lat, address(13), ierror)
    call mpi_get_address(initial%vloc, address(14), ierror)
    call mpi_get_address(initial%error_variance, address(15), ierror)

    ! also do this for values_qc derived type
    call mpi_get_address(init_val, address_vals(1), ierror)
    call mpi_get_address(init_val%time_actual, address_vals(2), ierror)
    call mpi_get_address(init_val%val, address_vals(3), ierror)
    call mpi_get_address(init_val%qc, address_vals(4), ierror)

    ! define types in struct in terms of base MPI datatypes
    oldtypes(1:num_ints) = MPI_INTEGER
    oldtypes(11) = MPI_INTEGER8
    oldtypes(num_ints+num_longs+1:num_ints+num_longs+num_doubles) = MPI_REAL8

    ! same with values_qc
    oldtypes_vals(1:num_ints_vals) = MPI_INTEGER
    oldtypes_vals(2) = MPI_INTEGER8
    oldtypes_vals(num_ints_vals+num_longs+1:num_ints_vals+num_longs+num_doubles_vals) = MPI_REAL8

    ! set the offsets of each of the variables from the first
    offsets(1) = 0
    do i = 2, num_vars
        offsets(i) = address(i) - address(1)
    enddo 

    ! you know the drill...
    offsets_vals(1) = 0
    do i = 2, num_vars_vals
        offsets_vals(i) = address_vals(i) - address_vals(1)
    enddo 

    ! Commit the new data types
    call mpi_type_create_struct(num_vars, bl_var, offsets, oldtypes, mpi_obstype, ierror)
    call mpi_type_commit(mpi_obstype, ierror)

    call mpi_type_create_struct(num_vars_vals, bl_var_vals, offsets_vals, oldtypes_vals, mpi_valstype, ierror)
    call mpi_type_commit(mpi_valstype, ierror)

end subroutine setup_obs_mpi
!------------------------------------------------------------------
!------------------------------------------------------------------
function compare_obs(obs1, obs2) Bind(C)
    type(obs_type_send),       intent(in)       :: obs1, obs2 
    integer(c_int)                              :: compare_obs

    if (obs1%time_order < obs2%time_order) then
        compare_obs = -1
    else if (obs1%time_order > obs2%time_order) then
        compare_obs = 1
    else
        compare_obs = 0
    endif

end function compare_obs
!------------------------------------------------------------------
!------------------------------------------------------------------
function compare_vals(val1, val2) Bind(C)
    type(obs_values_qc_type),       intent(in)      :: val1, val2 
    integer(c_int)                                  :: compare_vals

    if (val1%time_order < val2%time_order) then
        compare_vals = -1
    else if (val1%time_order > val2%time_order) then
        compare_vals = 1
    else
        compare_vals = 0
    endif

end function compare_vals
!------------------------------------------------------------------
!------------------------------------------------------------------
function compare_time_types(val1, val2) Bind(C)
    type(obs_values_qc_type),       intent(in)      :: val1, val2 
    integer(c_int)                                  :: compare_time_types

    if (val1%time_actual < val2%time_actual) then
        compare_time_types = -1
    else if (val1%time_actual > val2%time_actual) then
        compare_time_types = 1
    else
        compare_time_types = 0
    endif

end function compare_time_types
!------------------------------------------------------------------
!------------------------------------------------------------------
function compare_time_types_alt(val1, val2) Bind(C)
    type(obs_type),       intent(in)      :: val1, val2 
    integer(c_int)                                  :: compare_time_types_alt

    if (val1%time_actual < val2%time_actual) then
        compare_time_types_alt = -1
    else if (val1%time_actual > val2%time_actual) then
        compare_time_types_alt = 1
    else
        compare_time_types_alt = 0
    endif

end function compare_time_types_alt
!------------------------------------------------------------------
!------------------------------------------------------------------
function compare_large_int(val1, val2) Bind(C)
    integer(i8),                    intent(in)      :: val1, val2 
    integer(c_int)                                  :: compare_large_int

    if (val1 < val2) then
        compare_large_int = -1
    else if (val1 > val2) then
        compare_large_int = 1
    else
        compare_large_int = 0
    endif

end function compare_large_int
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine sort_obs_send_by_time(obs_set, values_qc, num_values, num_obs)
    ! todo: make this sort in-place; too much memory used currently
    ! need to allocate two massive arrays
    type(obs_type_send), target,        intent(inout)       :: obs_set(num_obs)
    type(obs_values_qc_type), target,   intent(inout)       :: values_qc(num_obs * num_values)
    integer,                            intent(in)          :: num_values, num_obs
    integer                                                 :: i, next_time, j, k, val_idx, x, total_values, l
    integer(C_SIZE_T)                                       :: num_obs_c, total_vals_c, sizeof_val, sizeof_obs
    ! integer                                                 :: test(4), test_2(4), 

    print *, 'odt%total_obs: ', odt%total_obs
    print *, 'odt%rem: ', odt%rem
    i = 1
    j = 1
    k = 1
    do while (i /= -1)
        ! if (modulo(i, 10000) == 0) print *, 'i = ', i
        l = get_obs_offset(i)

        obs_set(l)%time_order = j 
        ! obs_set(i)%time_order = num_obs - j

        ! need to deal with values_qc as well
        val_idx = ((l - 1) * num_values) + 1
        do x = val_idx, (val_idx + num_values) - 1
            values_qc(x)%time_order = k
            ! values_qc(x)%time_order = (num_obs * num_values) - k 
            k = k + 1
        enddo 

        ! move indices
        i = obs_set(l)%next_time
        ! print *, 'i = ', i
        ! if (i == 16777217) print *, 'i (l = 16777217): ', i
        j = j + 1
    enddo
    ! print *, 'Made it to qsort!'

    total_values = num_obs * num_values
    total_vals_c = total_values
    num_obs_c = num_obs
    sizeof_val = sizeof(values_qc(1))
    sizeof_obs = sizeof(obs_set(1))

    print *, 'sort 1'
    call qsort(c_loc(obs_set(1)), int(num_obs, c_size_t), sizeof_obs, c_funloc(compare_obs))
    print *, 'sort 2'
    call qsort(c_loc(values_qc(1)), int(num_obs, c_size_t), sizeof_val, c_funloc(compare_vals))
    


end subroutine sort_obs_send_by_time
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine sort_obs_by_time_alt(obs_set, values_qc, num_values, num_obs)
    ! todo: make this sort in-place; too much memory used currently
    ! need to allocate two massive arrays
    type(obs_type_send), target,        intent(inout)       :: obs_set(num_obs)
    type(obs_values_qc_type), target,   intent(inout)       :: values_qc(num_obs * num_values)
    integer,                            intent(in)          :: num_values, num_obs
    integer                                                 :: i, next_time, j, k, val_idx, x, total_values, l
    integer(C_SIZE_T)                                       :: num_obs_c, total_vals_c, sizeof_val, sizeof_obs
    ! integer                                                 :: test(4), test_2(4), 

    print *, 'odt%total_obs: ', odt%total_obs
    print *, 'odt%rem: ', odt%rem
    i = 1
    j = 1
    k = 1
    do while (i /= -1)
        ! if (modulo(i, 10000) == 0) print *, 'i = ', i
        l = get_obs_offset(i)

        obs_set(l)%time_order = j 
        ! obs_set(i)%time_order = num_obs - j

        ! need to deal with values_qc as well
        val_idx = ((l - 1) * num_values) + 1
        do x = val_idx, (val_idx + num_values) - 1
            values_qc(x)%time_order = k
            ! values_qc(x)%time_order = (num_obs * num_values) - k 
            k = k + 1
        enddo 

        ! move indices
        i = obs_set(l)%next_time
        ! print *, 'i = ', i
        ! if (i == 16777217) print *, 'i (l = 16777217): ', i
        j = j + 1
    enddo
    ! print *, 'Made it to qsort!'

    total_values = num_obs * num_values
    total_vals_c = total_values
    num_obs_c = num_obs
    sizeof_val = sizeof(values_qc(1))
    sizeof_obs = sizeof(obs_set(1))

    print *, 'sort 1'
    call qsort(c_loc(obs_set(1)), int(num_obs, c_size_t), sizeof_obs, c_funloc(compare_obs))
    print *, 'sort 2'
    call qsort(c_loc(values_qc(1)), int(num_obs, c_size_t), sizeof_val, c_funloc(compare_vals))
    

end subroutine sort_obs_by_time_alt
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine send_obs_set(set, proc, num_obs, num_values)
    type(obs_type),             intent(inout)      :: set(:)
    integer,                    intent(in)      :: proc
    integer,                    intent(inout)      :: num_obs
    integer,                    intent(inout)      :: num_values

    type(obs_type_send), allocatable                :: conv_set(:)
    type(obs_values_qc_type), allocatable            :: values_qc(:)
    integer                                        :: total_values
    integer                                        :: obs_mpi, ierror, i, d, diff, val_mpi

    total_values = num_values * num_obs
    allocate(conv_set(num_obs))
    print *, 'num_obs: ', num_obs
    allocate(values_qc(total_values * 10))
    
    ! Convert to a struct with fully contiguous memory 
    call convert_obs_set(set, conv_set, values_qc, num_values, num_obs)

    ! if (proc == 1) then
    !     print *, 'values(1) (before recv): ', values(1)
    !     print *, 'set(1)%values(1): ', set(1)%values(1)
    ! endif 

    ! Setup the dedicated data structure
    call setup_obs_mpi(obs_mpi, val_mpi)

    ! Send the full set to processor specified
    ! also send values and qc
    ! KY this isn't ideal; takes way too long to communicate
    ! print *, 'total_values: ', total_values
    call mpi_ssend(values_qc, total_values, val_mpi, proc, 0, MPI_COMM_WORLD, ierror)
    call mpi_ssend(conv_set, num_obs, obs_mpi, proc, 0, MPI_COMM_WORLD, ierror)
    
    ! Let's try using MPI_Scatter instead! 
    ! Or perhaps use MPI_Scatterv? (fewer scatter calls)
    ! Test simplest possible option to determine how fast it will be (not 100% accurate distribution, but still)
    ! call mpi_scatter(values, total_values, MPI_REAL8, total_values, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    ! call mpi_scatter(conv_set, num_obs, obs_mpi, num_obs, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)
    


    ! After send has completed, destroy the mpi struct and deallocate memory
    call destroy_obs_mpi(obs_mpi, val_mpi)
    deallocate(conv_set)
    deallocate(values_qc)


end subroutine send_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine recv_obs_set(set, proc, num_obs, num_values)
    type(obs_type),             intent(inout)      :: set(:)
    integer,                    intent(in)         :: proc
    integer,                    intent(inout)      :: num_obs
    integer,                    intent(inout)      :: num_values

    type(obs_type_send), allocatable               :: conv_set(:)
    type(obs_values_qc_type), allocatable          :: values_qc(:)
    integer                                        :: total_values
    integer                                        :: obs_mpi, ierror, i, diff, d, val_mpi

    total_values = num_values * num_obs
    allocate(conv_set(num_obs))
    allocate(values_qc(total_values * 10))
    
    ! print *, 'set(1)%values(1): ', set(1)%values(1)
    ! Convert to a struct with fully contiguous memory 
    ! call convert_obs_set(set, conv_set, values, qc, num_values, num_obs)

    ! Setup the dedicated data structure
    call setup_obs_mpi(obs_mpi, val_mpi)
    ! Retrieve the full set from proc specified
    ! also retrieve values and qc
    ! print *, 'trying to receive observations' 
    ! print *, 'total_values :', total_values
    call mpi_recv(values_qc, total_values, val_mpi, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    ! print *, '1'
    ! print *, '2'
    call mpi_recv(conv_set, num_obs, obs_mpi, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    print *, 'made past mpi_recv'
    ! print *, '3'

    ! print *, 'set(1)%values(1) (before convert_obs_back): ', set(1)%values(1)
    call convert_obs_back(set, conv_set, values_qc, num_obs, num_values)
    ! print *, 'set(1)%key (before setting values): ', set(1)%key

    ! if (allocated(set(1)%values)) then
    ! call print_obs_send(conv_set(1))
        ! print *, 'set(1)%values(1) (after convert_obs_back): ', set(1)%values(1)
    ! endif
    ! print *, 'hello from process ', my_task_id() 

    deallocate(conv_set)
    deallocate(values_qc)
    ! After send has completed, destroy the mpi struct
    print *, 'before destroy_obs_mpi call'
    call destroy_obs_mpi(obs_mpi, val_mpi)
    print *, 'after destroy_obs_mpi call'


end subroutine recv_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine scatter_obs_set(set, new_set, num_obs_per_proc, num_values, nprocs, root)
    type(obs_type),             intent(inout)      :: set(:)
    type(obs_type),             intent(inout)      :: new_set(:)
    integer,                    intent(inout)      :: num_obs_per_proc
    integer,                    intent(inout)      :: num_values
    integer,                    intent(in)      :: nprocs
    integer,                    intent(in)      :: root

    type(obs_type_send), allocatable               :: conv_set(:)
    type(obs_type_send), allocatable               :: all_conv_set(:)
    ! real(r8), allocatable                          :: values(:)
    ! real(r8), allocatable                          :: qc(:)
    type(obs_values_qc_type), allocatable           :: all_values_qc(:)
    type(obs_values_qc_type), allocatable           :: values_qc(:)
    integer                                        :: total_values, all_values, vals_per_proc, total_obs, val_mpi
    integer                                        :: obs_mpi, ierror, i, d, diff, j, l

    total_obs = num_obs_per_proc * nprocs
    total_values = num_values * total_obs
    if (my_task_id() == root) then
        allocate(all_conv_set(total_obs))
    else
        allocate(all_conv_set(1))
    endif
    allocate(conv_set(num_obs_per_proc))
    ! allocate(values(total_values * 10))

    all_values = total_values * 2 ! values + qc = 2
    vals_per_proc = all_values / nprocs
    if (my_task_id() == root) then
        allocate(all_values_qc(all_values))
    else
        allocate(all_values_qc(1))
    endif
    allocate(values_qc(vals_per_proc))
    ! allocate(qc(total_values * 10))
    
    ! Convert to a struct with fully contiguous memory 
    if (my_task_id() == root) then
        call convert_obs_set(set, all_conv_set, all_values_qc, num_values, total_obs)
    endif

    ! Setup the dedicated data structure
    call setup_obs_mpi(obs_mpi, val_mpi)

    ! Let's try using MPI_Scatter instead! 
    ! Or perhaps use MPI_Scatterv? (fewer scatter calls)
    ! Test simplest possible option to determine how fast it will be (not 100% accurate distribution, but still)
    call mpi_scatter(all_conv_set, num_obs_per_proc, obs_mpi, conv_set, num_obs_per_proc, obs_mpi, 0, MPI_COMM_WORLD, ierror)
    call mpi_scatter(all_values_qc, vals_per_proc, val_mpi, values_qc, vals_per_proc, val_mpi, 0, MPI_COMM_WORLD, ierror)

    ! Repack obs sequence
    call convert_obs_back(new_set, conv_set, values_qc, num_obs_per_proc, num_values)
    ! print *, 'set(1)%key (before setting values): ', set(1)%key
    ! do values_qc conversion, but in reverse
    
    
    ! After send has completed, destroy the mpi struct datatype and deallocate memory
    call destroy_obs_mpi(obs_mpi, val_mpi)
    if (my_task_id() == 0) then
        deallocate(all_conv_set)
        deallocate(all_values_qc)
    endif

    deallocate(conv_set)
    deallocate(values_qc)

end subroutine scatter_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine gather_obs_varied(dest, src, dest_val, src_val, obs_mpi, val_mpi, num_values, num_obs, num_procs)
    integer,                     intent(in)         :: obs_mpi, val_mpi, num_values, num_procs, num_obs
    type(obs_type_send),          intent(inout)      :: dest(:), src(:)
    type(obs_values_qc_type),          intent(inout)      :: dest_val(:), src_val(:)
    integer                                         :: obs_per_proc, rem, i, j, k, ierror, vals_per_proc
    integer , allocatable                                        :: disp(:), disp_vals(:), count(:), count_vals(:)

    ! calculate observations and values per process
    obs_per_proc = num_obs / num_procs
    vals_per_proc = obs_per_proc * num_values
    rem = modulo(num_obs, num_procs)

    ! allocate temporary buffers
    allocate(disp(num_procs))
    allocate(disp_vals(num_procs))
    allocate(count(num_procs))
    allocate(count_vals(num_procs))

    ! determine count and displacement for gatherv
    do i = 1, num_procs
        count(i) = obs_per_proc 
        count_vals(i) = (obs_per_proc) * num_values
        if (i <= rem) then
            count(i) = count(i) + 1
            count_vals(i) = count_vals(i) + num_values
        endif
    enddo

    disp(1) = 0
    disp_vals(1) = 0
    do i = 2, num_procs
        disp(i) = disp(i - 1) + count(i - 1)
        disp_vals(i) = disp_vals(i - 1) + count_vals(i - 1)
    enddo

    ! gather observations and values
    call mpi_gatherv(src, count(my_task_id() + 1), odt%obs_mpi, dest, count, disp, odt%obs_mpi, 0, MPI_COMM_WORLD, ierror)
    call mpi_gatherv(src_val, count_vals(my_task_id() + 1), odt%val_mpi, dest_val, count_vals, disp_vals, odt%val_mpi, 0, MPI_COMM_WORLD, &
    ierror)

    ! deallocate used memory
    deallocate(disp)
    deallocate(disp_vals)
    deallocate(count)
    deallocate(count_vals)


end subroutine gather_obs_varied
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine scatter_obs_varied(dest, src, dest_val, src_val, obs_mpi, val_mpi, num_values, num_obs, num_procs)
    integer,                     intent(in)         :: obs_mpi, val_mpi, num_values, num_procs, num_obs
    type(obs_type_send),          intent(inout)      :: dest(:), src(:)
    type(obs_values_qc_type),          intent(inout)      :: dest_val(:), src_val(:)
    integer                                         :: obs_per_proc, rem, i, j, k, ierror, vals_per_proc
    integer , allocatable                                        :: disp(:), disp_vals(:), count(:), count_vals(:)

    ! calculate observations and values per process
    obs_per_proc = num_obs / num_procs
    vals_per_proc = obs_per_proc * num_values
    rem = modulo(num_obs, num_procs)

    ! allocate temporary buffers
    allocate(disp(num_procs))
    allocate(disp_vals(num_procs))
    allocate(count(num_procs))
    allocate(count_vals(num_procs))

    ! determine count and displacement for scatterv
    do i = 1, num_procs
        count(i) = obs_per_proc 
        count_vals(i) = (obs_per_proc) * num_values
        if (i <= rem) then
            count(i) = count(i) + 1
            count_vals(i) = count_vals(i) + num_values
        endif
    enddo

    disp(1) = 0
    disp_vals(1) = 0
    do i = 2, num_procs
        disp(i) = disp(i - 1) + count(i - 1)
        disp_vals(i) = disp_vals(i - 1) + count_vals(i - 1)
    enddo

    ! scatter observations and values
    call mpi_scatterv(src, count, disp, odt%obs_mpi, dest, count(my_task_id() + 1), odt%obs_mpi, 0, MPI_COMM_WORLD, ierror)
    call mpi_scatterv(src_val, count_vals, disp_vals, odt%val_mpi, dest_val, count_vals(my_task_id() + 1), odt%val_mpi, 0, MPI_COMM_WORLD, &
    ierror)

    ! deallocate used memory
    deallocate(disp)
    deallocate(disp_vals)
    deallocate(count)
    deallocate(count_vals)


end subroutine scatter_obs_varied
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine gather_obs_set(set, new_set, num_obs_per_proc, num_values, nprocs, root)
    type(obs_type),             intent(inout)      :: set(:)
    type(obs_type),             intent(inout)      :: new_set(:)
    integer,                    intent(inout)      :: num_obs_per_proc
    integer,                    intent(inout)      :: num_values
    integer,                    intent(in)      :: nprocs
    integer,                    intent(in)      :: root 

    type(obs_type_send), allocatable               :: conv_set(:)
    type(obs_type_send), allocatable               :: all_conv_set(:)
    ! real(r8), allocatable                          :: values(:)
    ! real(r8), allocatable                          :: qc(:)
    type(obs_values_qc_type), allocatable           :: all_values_qc(:)
    type(obs_values_qc_type), allocatable           :: values_qc(:)
    integer(i8)                                    :: total_values, all_values, vals_per_proc
    integer                                        :: obs_mpi, ierror, i, d, diff, j, l, first, val_mpi, total_obs

    total_obs = num_obs_per_proc * nprocs
    total_values = num_values * total_obs ! values and qc
    allocate(conv_set(num_obs_per_proc))
    if (my_task_id() == root) then 
        allocate(all_conv_set(total_obs))
        allocate(all_values_qc(total_values))
    else
        allocate(all_conv_set(1))
        allocate(all_values_qc(1))
    endif
    ! allocate(values(total_values * 10))

    vals_per_proc = total_values / nprocs
    allocate(values_qc(vals_per_proc))
    ! allocate(qc(total_values * 10))
    
    ! Convert to a struct with fully contiguous memory 
    call convert_obs_set(set, conv_set, values_qc, num_values, num_obs_per_proc)

    ! Setup the dedicated data structure
    call setup_obs_mpi(obs_mpi, val_mpi)
    
    call mpi_barrier(MPI_COMM_WORLD, ierror)
    ! Send the full set to processor specified
    ! also send values and qc
    ! KY this isn't ideal; takes way too long to communicate
    ! print *, 'total_values: ', total_values
    ! call mpi_send(values, total_values, MPI_REAL8, proc, 0, MPI_COMM_WORLD, ierror)
    ! call mpi_send(qc, total_values, MPI_REAL8, proc, 0, MPI_COMM_WORLD, ierror)
    ! call mpi_send(conv_set, num_obs, obs_mpi, proc, 0, MPI_COMM_WORLD, ierror)
    
    ! Let's try using MPI_Scatter instead! 
    ! Or perhaps use MPI_Scatterv? (fewer scatter calls)
    ! Test simplest possible option to determine how fast it will be (not 100% accurate distribution, but still)
    ! will likely use padding to guarantee all obs and reals are sent 
    ! call mpi_scatter(all_conv_set, num_obs_per_proc, obs_mpi, conv_set, num_obs_per_proc, obs_mpi, 0, MPI_COMM_WORLD, ierror)
    ! call mpi_scatter(all_values_qc, vals_per_proc, MPI_REAL8, values_qc, vals_per_proc, MPI_REAL8, 0, MPI_COMM_WORLD, ierror)

    call mpi_gather(conv_set, num_obs_per_proc, obs_mpi, all_conv_set, num_obs_per_proc, obs_mpi, root, MPI_COMM_WORLD, ierror)
    call mpi_gather(values_qc, vals_per_proc, val_mpi, all_values_qc, vals_per_proc, val_mpi, root, MPI_COMM_WORLD, ierror)


    if (my_task_id() == root) then
        ! Test code here
        ! do j = num_obs_per_proc * 128 + 1, total_obs
        !     if (all_conv_set(j)%key == 0) then
        !         ! print *, 'key is zero at ', j - (num_obs_per_proc * 128)
        !     endif
        ! enddo
        
        call convert_obs_back(new_set, all_conv_set, all_values_qc, total_obs, num_values)
        ! print *, 'set(1)%key (before setting values): ', set(1)%key
        ! do values_qc conversion, but in reverse
        
    endif

    call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    ! After send has completed, destroy the mpi struct datatype and deallocate memory
    call destroy_obs_mpi(obs_mpi, val_mpi)
    if (my_task_id() == root) then
        deallocate(all_conv_set)
        deallocate(all_values_qc)
    endif

    deallocate(conv_set)
    deallocate(values_qc)

end subroutine gather_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine convert_obs_set(orig, out_obs, out_vals, num_values, num_obs)
    type(obs_type),             intent(inout)  :: orig(:)
    type(obs_type_send),        intent(inout) :: out_obs(:)
    type(obs_values_qc_type),   intent(inout) :: out_vals(:)
    integer,          intent(in)  :: num_values
    integer,          intent(in)  :: num_obs

    integer :: i, d, j, seconds, days, diff
    real(r8)    :: location(3)
    integer     :: which_vert
    type(location_type)  :: orig_location
    type(obs_def_type)   :: obs_def
    type(time_type)      :: obs_time


    do i = 1, num_obs

        ! get obs_def
        call get_obs_def(orig(i), obs_def)

        ! get location
        orig_location = get_obs_def_location(obs_def)
        location = get_location(orig_location)
        out_obs(i)%lon = location(1)
        out_obs(i)%lat = location(2)
        out_obs(i)%vloc = location(3)
        out_obs(i)%which_vert = nint(query_location(orig_location))

        ! get time
        obs_time = get_obs_def_time(obs_def)
        call get_time(obs_time, seconds, days)
        out_obs(i)%seconds = seconds
        out_obs(i)%days = days
        out_obs(i)%time_actual = (out_obs(i)%days * 24 * 60 * 60) + out_obs(i)%seconds

        ! get other values
        out_obs(i)%kind = get_obs_def_type_of_obs(obs_def)
        out_obs(i)%key = orig(i)%key
        out_obs(i)%prev_time = orig(i)%prev_time
        out_obs(i)%next_time = orig(i)%next_time
        out_obs(i)%cov_group = orig(i)%cov_group
        out_obs(i)%error_variance = get_obs_def_error_variance(obs_def)
        out_obs(i)%obs_def_key = get_obs_def_key(obs_def)

    enddo

    d = 1
    diff = num_values - 1 
    do i = 1, num_obs
        ! values
        ! values_qc(d:d+diff) = set(i)%values(1:num_values)
        do j = 1, num_values
            out_vals(j+d-1)%val = orig(i)%values(j)
            out_vals(j+d-1)%qc = orig(i)%qc(j)
            out_vals(j+d-1)%time_actual = out_obs(i)%time_actual
        enddo
        ! values_qc(d:d+diff)%qc = set(i)%qc(1:num_values)

        ! qc
        ! l = d + num_values
        ! values_qc(l:l+diff) = set(i)%qc(1:num_values)
        
        ! increment our offset
        d = d + (num_values)
    enddo 

end subroutine convert_obs_set
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine print_obs_send(obs_send)
    ! for debugging purposes
    type(obs_type_send),            intent(in)      :: obs_send
    print *,'obs_send%key: ', obs_send%key
    print *,'obs_send%which_vert: ',obs_send%which_vert
    print *,'obs_send%kind: ',obs_send%kind
    print *,'obs_send%seconds: ',obs_send%seconds
    print *,'obs_send%days: ',obs_send%days
    print *,'obs_send%obs_def_key: ',obs_send%obs_def_key
    print *,'obs_send%prev_time: ',obs_send%prev_time
    print *,'obs_send%next_time: ',obs_send%next_time
    print *,'obs_send%cov_group: ',obs_send%cov_group
    print *,'obs_send%lon: ',obs_send%lon
    print *,'obs_send%lat: ',obs_send%lat
    print *,'obs_send%vloc: ',obs_send%vloc
    print *,'obs_send%error_variance: ',obs_send%error_variance
end subroutine print_obs_send
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine convert_obs_back(recv, simple_obs, simple_val_qc, num_obs, num_values)
    type(obs_type),             intent(inout)  :: recv(:)
    type(obs_type_send),        intent(inout) :: simple_obs(:)
    type(obs_values_qc_type),   intent(inout) :: simple_val_qc(:)
    ! integer,          intent(in)  :: num_values
    integer,          intent(in)  :: num_obs

    integer :: i, d, j, seconds, days
    real(r8)    :: location(3)
    type(location_type)  :: location_def
    type(obs_def_type)   :: obs_def


    ! print *, 'start of loop 1'
    do i = 1, num_obs
        ! if (odt%my_pe == 1 .and. i > 110000) then
        !     print *, 'i = ', i
        ! endif

        ! print *, 'i: ', i
        ! set location
        call set_obs_def_location(recv(i)%def, set_location(simple_obs(i)%lon, simple_obs(i)%lat, simple_obs(i)%vloc, &
            simple_obs(i)%which_vert))

        ! get time
        call set_obs_def_time(recv(i)%def, set_time(simple_obs(i)%seconds, simple_obs(i)%days))
        recv(i)%time_actual = (simple_obs(i)%days*24*60*60)+simple_obs(i)%seconds

        ! get other values
        if (i == 1) then
            ! call print_obs_send(simple_obs(1))
        endif
        call set_obs_def_type_of_obs(recv(i)%def, simple_obs(i)%kind)
        recv(i)%prev_time = simple_obs(i)%prev_time
        recv(i)%next_time = simple_obs(i)%next_time
        recv(i)%cov_group = simple_obs(i)%cov_group
        recv(i)%key = simple_obs(i)%key
        call set_obs_def_error_variance(recv(i)%def, simple_obs(i)%error_variance)
        call set_obs_def_key(recv(i)%def, simple_obs(i)%obs_def_key)

    enddo
    ! print *, 'end of loop 1'


    ! allocatable components returned home
    d = 1
    diff = num_values - 1 
    do i = 1, num_obs
        ! print *, 'i: ', i
        ! values
        ! values_qc(d:d+diff) = set(i)%values(1:num_values)
        do j = 1, num_values
            recv(i)%values(j) = simple_val_qc(j+d-1)%val
            recv(i)%qc(j) = simple_val_qc(j+d-1)%qc
        enddo
        ! values_qc(d:d+diff)%qc = set(i)%qc(1:num_values)

        ! qc
        ! l = d + num_values
        ! values_qc(l:l+diff) = set(i)%qc(1:num_values)
        
        ! increment our offset
        d = d + (num_values)
    enddo 

end subroutine convert_obs_back

!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine get_job_info(nranks, nnodes)
    integer,            intent(out)  :: nranks
    integer,            intent(out)  :: nnodes
    character(len=512)               :: env_value
    character(len=64), dimension(13)  :: envs_sep
    
    integer                          :: i

    ! retrieve environment variable
    call get_environment_variable('PBS_SELECT', env_value)

    ! replace separators with comma to make read easier
    do i = 1, len_trim(env_value)
        if (env_value(i:i) == ':' .or. env_value(i:i) == '=') env_value(i:i) = ','
    enddo

    ! read each value
    read(env_value, *) envs_sep(1:13)
    read(envs_sep(1), *) nnodes
    read(envs_sep(3), *) nranks


end subroutine get_job_info
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine allocate_obs_set(buffer, size, num_values)
    type(obs_type), allocatable,    intent(inout)       :: buffer(:)
    integer,                        intent(in)          :: size
    integer,                        intent(in)          :: num_values
    integer                                             :: i, m
    real(r8)                                            :: k, j, l

    m = -1
    k = 0.0
    j = 0.0
    l = 0.0
    allocate(buffer(size))
    do i = 1, size
        allocate(buffer(i)%values(num_values))
        allocate(buffer(i)%qc(num_values))
        call set_obs_def_location(buffer(i)%def, set_location(k, j, l, m))
    enddo

end subroutine allocate_obs_set
!------------------------------------------------------------------

!------------------------------------------------------------------
subroutine deallocate_obs_set(buffer, size, num_values)
    type(obs_type), allocatable,            intent(inout)       :: buffer(:)
    integer,                                intent(in)          :: size
    integer,                                intent(in)          :: num_values
    integer                                                     :: i

    do i = 1, size
        deallocate(buffer(i)%values)
        deallocate(buffer(i)%qc)
    enddo
    deallocate(buffer)

end subroutine deallocate_obs_set
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine sort_roundrobin(buffer, new_buffer, num_obs, num_procs)
    type(obs_type),             intent(in)       :: buffer(:)
    integer,                    intent(in)       :: num_procs, num_obs
    type(obs_type),             intent(out)      :: new_buffer(:)
    integer                                      :: i, j, obs_per_proc, rem, k
    
    obs_per_proc = num_obs / num_procs
    rem = modulo(num_obs, num_procs)

    k = 0
    print *, 'num_procs: ', num_procs
    print *, 'obs_per_proc: ', obs_per_proc
    do i = 1, num_procs
        do j = i, obs_per_proc
            ! print *, '(', i, ',', j, ')'
            if (buffer(((i-1)*obs_per_proc) + j)%key == 0 .and. k == 0) then
                print *, ' 0 key was found at obs ', ((i - 1) * obs_per_proc) + j
                k = 1
            endif
            new_buffer(i*j) = buffer(((i-1)*obs_per_proc) + j)
        enddo
    enddo 

    if (rem > 0) then
        do i = 1, rem
            new_buffer((obs_per_proc * (i + 1)) - 1) = buffer((obs_per_proc * (i + 1)) - 1) 
        enddo 
    endif

end subroutine sort_roundrobin
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine sort_roundrobin_inplace(obs_set, values_qc, num_obs, num_values, num_procs)
    type(obs_type_send), target,             intent(inout)       :: obs_set(:)
    integer,                    intent(in)       :: num_procs, num_obs, num_values
    type(obs_values_qc_type), target,            intent(inout)      :: values_qc(:)
    integer                                      :: i, j, obs_per_proc, rem, k, obs_per_proc_rem, total_values, these_obs, d
    integer(C_SIZE_T)                            :: num_obs_c, total_vals_c, sizeof_val, sizeof_obs
    integer, allocatable                                      :: procs(:)


    obs_per_proc = num_obs / num_procs
    rem = modulo(num_obs, num_procs)

    obs_per_proc_rem = obs_per_proc + 1
    allocate(procs(num_procs))
    do i = 1, num_procs
        procs(i) = obs_per_proc
        if (i <= rem) procs(i) = procs(i) + 1
    enddo

    
    k = 1
    do i = 1, obs_per_proc
        do j = 1, num_procs
            if (j <= rem) then
                d = obs_per_proc_rem * (j - 1)
            else 
                d = (obs_per_proc_rem * rem) + (obs_per_proc * ((j - rem) - 1))
            endif
            obs_set(d + i)%time_order = k
            values_qc(d + i)%time_order = k
            k = k + 1
        enddo
    enddo
    
    do i = 1, rem
        d = obs_per_proc_rem * i
        obs_set(d)%time_order = d
        values_qc(d)%time_order = d 
    enddo

    ! k = 1 
    ! d = 1
    ! do i = 1, num_procs
    !     if (i <= rem) then
    !         d = obs_per_proc_rem * (i - 1)
    !     else 
    !         d = (obs_per_proc_rem * rem) + (obs_per_proc * ((i - rem) - 1))
    !     endif
    !     do j = 1, procs(i)
    !         if (j <= obs_per_proc) then
    !             obs_set(d + j)%time_order = (obs_per_proc * (i - 1)) + j
    !             values_qc(d + j)%time_order = (obs_per_proc * (i - 1)) + j
    !         else
    !             obs_set(d + j)%time_order = d + j
    !             values_qc(d + j)%time_order = d + j
    !         endif
    !     enddo
    ! enddo

    ! sort both lists
    total_values = num_obs * num_values
    total_vals_c = total_values
    num_obs_c = num_obs
    sizeof_val = sizeof(values_qc(1))
    sizeof_obs = sizeof(obs_set(1))

    call qsort(c_loc(obs_set(1)), int(num_obs, c_size_t), sizeof_obs, c_funloc(compare_obs))
    call qsort(c_loc(values_qc(1)), int(num_obs, c_size_t), sizeof_val, c_funloc(compare_vals))
    
end subroutine sort_roundrobin_inplace
!------------------------------------------------------------------
end module obs_dist_mod