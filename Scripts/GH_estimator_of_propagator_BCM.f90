program Evaluate_propagator
    ! Estimate the propagator of a jump with the BCM estimator and compare it with the exact propagato (when available) or the Gaussian approximation
    use declarations_module
    use functions_local
    use functions
    use subroutines
    implicit none
    integer*8 :: n_params_target_model
    double precision, allocatable :: lk_target(:),lk_ref(:)
    character(len=3) :: target_model,bridge_model
    double precision ::  x0, t0,tf,xf,xmin,xmax,D_2,dt,dummy,dt_bridge,prop_estimator,S_err,R_err,tfmt0,prop_true,prop_Gauss
    integer*8 :: seed,N_Bessel
    double precision, allocatable :: model_params(:),params_target_model(:)
    integer*8 :: ii,N_bridges
    integer :: unit
    character(len=100) :: addname,filename
    logical :: compute_R_err

    call assingments()
    seed = time()
    print *, "seed=", seed
    print *

    call dran_ini(int(seed))

    ! Parameters for the simulation
    target_model = "DE"                       ! Change to "OU", "EN", "CP", "DE", "GGM" as needed
    if ((target_model == "OU").or.(target_model == "DE").or.(target_model == "GB")) then
        compute_R_err = .true.
    else
        compute_R_err = .false.
    endif
    bridge_model = "OU" ! Change to "OU", "WI", "GB" as needed
    mu = 1.0d0 
    K = 1.0d0 
    D = 1.0d0 
    
    x0 = mu !initial condition
    xf = x0 + 0.1 !final condition
    t0 = 0.0d0
    tf = 1.0d0

    D_2 = D**2
    N_Bessel = 20000 ! Number of terms in the truncated series expansion of the Bessel function to evaluate the propagator of the DE model.
    tfmt0 = tf-t0
    dt_bridge = 5*1.0d-2 !\Delta t_B in the paper (time step of the bridges)
    N_bridges = 1000 !N_B in the paper (number of bridges generated)
    
    if (target_model == "GB") then 
        n_params_target_model = int8(2)
        allocate(params_target_model(n_params_target_model))
        params_target_model = [K,D]
    else 
        n_params_target_model = int8(3)
        allocate(params_target_model(n_params_target_model))
        params_target_model = [mu,K,D]
    endif
    print *, "Estimating propagator for target_model: ", trim(target_model)
    print *, "from t0,x0 = ", t0,x0,"to tf,xf = ", tf,xf
    allocate(lk_target(N_bridges),lk_ref(N_bridges))
    select case (trim(target_model))
    case("OU")
        call propagator_OU(mu,K,D_2,x0,xf,tfmt0,prop_true)
    case("DE")
        call propagator_dem(mu,K,D_2,x0,xf,tfmt0,N_Bessel,prop_true)
    case("GB")
        call propagator_GB(K,D_2,x0,xf,tfmt0,prop_true)
    case default
        print*, "There is no exact propagator for requested model"
      end select
    if (trim(bridge_model)=="GB") then
        call Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,prop_estimator,S_err,R_err,compute_R_err,lk_target,lk_ref)
    else if ((trim(bridge_model)=="OU").or.(trim(bridge_model)=="WI")) then
        call Draw_lks_on_arrays_and_estimate_propagator_W_BCM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,prop_estimator,S_err,R_err,compute_R_err,lk_target,lk_ref)
    else
        print *, "Bridge model not recognized"
        print *, "location: main_evaluate_propagator_on_jump.f90"
        stop
    endif
    

    call Evaluate_propagator_jump_Gaussian_approx &
        (x0, t0, xf, tf, prop_Gauss, target_model, params_target_model, n_params_target_model)

    ! Save the propagator to a file
    filename = "propagator.dat"
    print *, "lks saved to "//trim(filename)
    unit = 1001
    open(unit=unit, file=trim(filename), iostat=ios, status="unknown", action="write")
    do ii = 1, N_bridges
        write(unit, *) lk_target(ii), lk_ref(ii)
    end do
    close(unit)

    print *, "prop_true=",prop_true,"Estimator = ", prop_estimator, "S_err = ", S_err, "R_err = ", R_err
    print *, "prop_Gauss=",prop_Gauss
    deallocate(lk_target,lk_ref) 
     
end program Evaluate_propagator