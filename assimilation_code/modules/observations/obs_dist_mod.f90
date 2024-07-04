module obs_dist_mod
    use mpi
    use iso_c_binding
    use     mpi_utilities_mod, only : task_count, my_task_id, send_to, receive_from
    use        types_mod, only : r8, i8, MISSING_R8, metadatalength

    use      obs_def_mod, only : obs_def_type, get_obs_def_time, read_obs_def, &
                                 write_obs_def, destroy_obs_def, copy_obs_def, &
                                 interactive_obs_def, get_obs_def_location, &
                                 get_obs_def_type_of_obs, get_obs_def_key,  &
                                 get_obs_def_error_variance, &
                                 set_obs_def_key, set_obs_def_time, set_obs_def_location, &
                                 set_obs_def_type_of_obs, set_obs_def_error_variance
    use obs_sequence_mod, only : obs_type

    type obs_values_qc_type
       private
       integer :: time_order
       real(r8) :: val
       real(r8) :: qc
    end type obs_values_qc_type

    type obs_type_send
    ! The key is needed to indicate the element number in the storage for the obs_sequence
    ! Do I want to enforce the identity of the particular obs_sequence?
    ! Declare a simplified data structure for sending and receiving using MPI
       private 
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
       real(r8) :: lon
       real(r8) :: lat
       real(r8) :: vloc
       real(r8) :: error_variance
    end type obs_type_send

    interface
        subroutine qsort(array, elem_count, elem_size, compare) bind(C, name="qsort")
            import
            type(c_ptr), value :: array
            integer(c_size_t), value :: elem_count
            integer(c_size_t), value :: elem_size
            type(c_funptr), value   :: compare !int(*compare)(const void *, const void *)
        end subroutine qsort

    end interface

    ! use these for our unpacked observations and values
    ! will be accessed whenever a process attempts a one-sided get
    type(obs_values_qc_type), allocatable       :: val_buf(:)
    type(obs_type_send), allocatable            :: obs_buf(:)

contains

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
    ! allocate(values(total_values * 10))

    gather_vals_per_proc = num_values * gather_obs_per_proc
    ! if (gather_proc < rem) gather_vals_per_proc = gather_vals_per_proc + (num_values)
    allocate(conv_set(gather_obs_per_proc))
    allocate(values_qc(gather_vals_per_proc))
    ! allocate(qc(total_values * 10))
    
    ! Convert to a struct with fully contiguous memory 
    call convert_obs_set(set, conv_set, num_values, gather_obs_per_proc)

    ! copy values from allocatable arrays
    d = 1
    diff = num_values - 1 
    do i = 1, gather_obs_per_proc
        ! values
        ! values_qc(d:d+diff) = set(i)%values(1:num_values)
        do j = 1, num_values
            values_qc(j+d-1)%val = set(i)%values(j)
            values_qc(j+d-1)%qc = set(i)%qc(j)
        enddo
        ! values_qc(d:d+diff)%qc = set(i)%qc(1:num_values)

        ! qc
        ! l = d + num_values
        ! values_qc(l:l+diff) = set(i)%qc(1:num_values)
        
        ! increment our offset
        d = d + (num_values)
    enddo 
    ! if (proc == 1) then
    !     print *, 'values(1) (before recv): ', values(1)
    !     print *, 'set(1)%values(1): ', set(1)%values(1)
    ! deallocate(set)

    ! Setup the dedicated data structure
    call setup_obs_mpi(obs_mpi, vals_mpi)
    
    ! prepare the gatherv; allocate displacement and count arrays
    allocate(disp(nprocs))
    allocate(count(nprocs))
    allocate(disp_vals(nprocs))
    allocate(count_vals(nprocs))

    ! Setup for obs mpi struct
    ! disp(1:start_proc) = 1
    count(1:start_proc) = 0
    ! disp_vals(1:start_proc) = 0
    count_vals(1:start_proc) = 0
    do i = 0, gather_procs - 1
        actual_proc = start_proc + 1 + i
        l = num_obs / gather_procs
        if (i < rem) then
            l = l + 1
        endif
        count(actual_proc) = l
    enddo

    disp(1:start_proc) = 0
    ! do i = 2, start_proc
    !     disp(i) = disp(i - 1) + count(i - 1)
    ! enddo

    disp(start_proc+1) = 0
    do i = 1, gather_procs - 1
        actual_proc = start_proc + 1 + i
        disp(actual_proc) = count(actual_proc - 1) + disp(actual_proc - 1)
    enddo

    ! Setup for values
    ! disp_vals(1:start_proc) = 0
    count_vals(1:start_proc) = 0
    do i = 0, gather_procs - 1
        actual_proc = start_proc + 1 + i
        count_vals(actual_proc) = num_values * count(actual_proc)
    enddo

    disp_vals(1:start_proc) = 0
    ! disp_vals(1) = 1
    ! do i = 2, start_proc
    !     disp_vals(i) = disp_vals(i - 1) + count_vals(i - 1)
    ! enddo

    disp_vals(start_proc + 1) = 0
    do i = 1, gather_procs - 1
        actual_proc = start_proc + 1 + i
        disp_vals(actual_proc) = count_vals(actual_proc - 1) + disp_vals(actual_proc - 1)
    enddo

    call mpi_barrier(MPI_COMM_WORLD, ierror)

    if (my_task_id() == root) print *, 'before gather'

    call mpi_gatherv(conv_set, count(my_task_id()+1), obs_mpi, all_conv_set, count, disp, &
        obs_mpi, root, MPI_COMM_WORLD, ierror)

    call mpi_gatherv(values_qc, count_vals(my_task_id()+1), vals_mpi, all_values_qc, count_vals, disp_vals, &
        vals_mpi, root, MPI_COMM_WORLD, ierror)

    if (my_task_id() == root) then 
        print *, 'after gather'
        call sort_obs_send_by_time(all_conv_set, all_values_qc, num_values, num_obs)
        print *, 'sorted by time (timestamp added)'
        ! call sort_roundrobin_inplace(all_conv_set, all_values_qc, num_obs, num_values, nprocs)
        print *, 'sorted with roundrobin dist'
    endif

    call mpi_barrier(MPI_COMM_WORLD, ierror)

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

    call convert_obs_back(new_set, new_conv_set, obs_per_proc)

    d = 1
    diff = num_values - 1 
    do i = 1, obs_per_proc
        ! values
        ! values_qc(d:d+diff) = set(i)%values(1:num_values)
        do j = 1, num_values
            new_set(i)%values(j) = new_values_qc(j+d-1)%val
            new_set(i)%qc(j) = new_values_qc(j+d-1)%qc
        enddo
        ! values_qc(d:d+diff)%qc = set(i)%qc(1:num_values)

        ! qc
        ! l = d + num_values
        ! values_qc(l:l+diff) = set(i)%qc(1:num_values)
        
        ! increment our offset
        d = d + (num_values)
    enddo 
    ! todo: sort by time, organize roundrobin, and scatter
    ! todo: maybe scatter roundrobin using point to point ops?
    ! todo: figure out how to sort in-place rather than using an additional array
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

    ! if (my_task_id() == root) then
    !     ! Repack obs sequence
    !     call print_obs_send(all_conv_set(128*num_obs_per_proc+1))
    !     call print_obs_send(all_conv_set(128*num_obs_per_proc+2))
    !     call convert_obs_back(new_set, all_conv_set, total_obs)
    !     ! print *, 'set(1)%key (before setting values): ', set(1)%key
    !     ! do values_qc conversion, but in reverse
    !     
    !     d = 1
    !     diff = num_values - 1 
    !     do i = 1, total_obs
    !         ! make sure arrays are allocated
    !         if (.not. allocated(new_set(i)%values)) then
    !             allocate(new_set(i)%values(num_values))
    !         endif
    !         if (.not. allocated(new_set(i)%qc)) then
    !             allocate(new_set(i)%qc(num_values))
    !         endif
    !
    !         ! values
    !         new_set(i)%values(1:num_values) = all_values_qc(d:d+diff)
    !
    !         ! qc
    !         l = d + num_values
    !         new_set(i)%qc(1:num_values) = all_values_qc(l:l+diff)
    !         
    !         ! increment our offset
    !         d = d + (num_values * 2)
    !     enddo 
    ! endif
    !
    ! call mpi_barrier(MPI_COMM_WORLD, ierror)
    
    ! After send has completed, destroy the mpi struct datatype and deallocate memory
    call destroy_obs_mpi(obs_mpi)
    if (my_task_id() == root) then
        deallocate(all_conv_set)
        deallocate(all_values_qc)
    endif

    deallocate(conv_set)
    deallocate(values_qc)
    deallocate(disp)
    deallocate(count)
    deallocate(disp_vals)
    deallocate(count_vals)

end subroutine dist_obs_set

subroutine initialize_obs_window()
    ! this will initialize the obs window for one-sided communication
end subroutine initialize_obs_window

subroutine destroy_obs_window()
    ! this will destroy the obs window (it's in the name, chief)
end subroutine destroy_obs_window

end module obs_dist_mod
