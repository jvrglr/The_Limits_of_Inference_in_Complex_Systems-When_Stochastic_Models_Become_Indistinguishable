program LKs_real_data
    use declarations_module
    use functions_local
    use functions 
    use subroutines
    implicit none
    double precision, allocatable :: ts(:), xs(:), LKs(:),ts_long(:), xs_long(:)
    integer*8 :: ii, jj, n_measures,len
    character(len=200) :: data_file, filename
    character(len=32) :: addname
    character(len=3) :: target_model, bridge_model
    double precision ::  LK, t_CPU
    real :: start, finish
    integer*8 :: seed, select_target_model, n_params_target_model
    double precision :: theta, E_D, E_D_2, E_mu, E_K, P
    double precision, dimension(3) :: params_target_model

    
    call assingments()
    ! seed = time()
    seed = 1994
    print *, "seed=", seed
    print *
    
    call dran_ini(int(seed))
    theta = 0.5d0
    n_params_target_model = 3
    P = 0.0d0
    call cpu_time(start)
    !-------------------Read data-------------------
    data_file = "Data/GH_twizzers.dat"
    print *, "  Read data from:", data_file
    call read_trajectory_unknown_nrows(trim(data_file), ts_long, xs_long, n_measures)
    print *, "Number of rows read: ", n_measures
    print *, "First time point: ", ts_long(1),"last time",ts_long(n_measures)
    print *, "First value: ", xs_long(1),"last state",xs_long(n_measures)
    print *
    !-----------------------------------------------

    len = 100
    allocate(xs(len), ts(len))
    xs = xs_long(1:len)
    ts = ts_long(1:len)


    print *
    do jj = 1,3,1
        select_target_model = jj !jj
        call select_model_name(select_target_model, target_model)


        print *, "Procesing data with "//trim(target_model)//" model ..."
        print *
        call estimate_parameters_Gaussian_and_lamperti_QV(E_D,E_D_2,E_mu,E_K,xs,ts,len,model=trim(target_model),params_model=[theta])
        print*, "Parameters estimated QV+Gauss: E_mu=",E_mu,"E_K=",E_K,"E_D=",E_D
        params_target_model = [E_mu,E_K,E_D]

        call Evaluate_log_likelihood_time_series_Gaussian_approx(xs,ts,len,LK,trim(target_model),params_target_model,n_params_target_model)
        print *, "Log-likelihood of the data under "//trim(target_model)//" model: LK=",LK
        print *
    enddo



end program LKs_real_data