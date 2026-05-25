program main_draw_trajectory_on_file
    ! Generate a trajectory of a stochastic process and save it to a file
    use declarations_module
    use functions_local
    use functions
    use subroutines
    implicit none
    double precision, allocatable :: xs(:), ts(:)
    character(len=10) :: model
    double precision :: dt_sampling, x0, t0, dt_integration,tf
    integer*8 :: len,seed
    double precision, allocatable :: model_params(:)
    integer*8 :: ii
    integer :: unit
    character(len=100) :: addname,filename

    call assingments()
    seed = time()
    print *, "seed=", seed
    print *

    call dran_ini(int(seed))

    ! Parameters for the simulation
    mu = 3.5829193064160207E-002
    K = 0.44660111475710185
    D = 0.18920043956403140 
    x0 = mu ! Initial condition 
    model = "OU"                   ! Change to "EN", "CP", "DE", "GGM" as needed
    t_sampling = 0.001d0 !Sampling time step (\Delta t in the paper)
    dt_integration = dt_sampling/10.0d0      ! integration time step (only used for models without analytical solution, for example: EN)
    tf = 10.19998                 ! Final time of the trajectory
    t0 = 0.0d0                 ! Initial time of the trajectory
    ! Optional model parameters (for example: GGM)
    if (model == "GGM") then
        allocate(model_params(1))
        model_params(1) = 1.5d0  ! Set theta value for GGM model
    end if
    len = int8((tf-t0)/dt_sampling)
    print *, "Generating trajectory for model: ", trim(model)
    print *, "tf = ", tf,dt_sampling * len
    print *, "Model parameters: "
    print *, "mu = ", mu, "K = ", K, "D = ", D, "x0 = ", x0
    

    allocate(xs(len), ts(len))

    ! Call the subroutine
    call generate_trajectory_in_array(trim(model), dt_sampling, x0, t0, xs, ts, len, dt_integration, model_params)

    ! Save the trajectory to a file
    filename = "Path.dat"

    print *, "Trajectory saved to "//trim(filename)
    unit = 1001
    open(unit=unit, file=trim(filename), iostat=ios, status="unknown", action="write")
    do ii = 1, len
        write(unit, *) ts(ii), xs(ii)
    end do
    close(unit)

    print *, "Trajectory generation complete.", "xf= ", xs(len), "tf= ", ts(len), "len= ", len
    deallocate(xs, ts)
    

101 end program main_draw_trajectory_on_file

  