program Model_distinguishability_EN_DE
    use declarations_module
    use functions_local
    use functions 
    use subroutines
    implicit none
    double precision, allocatable :: xs(:),ts(:),dts(:),Ps(:)
    integer*8, allocatable :: n_rows(:)
    integer*8 :: ii,jj,kk,ll,meas,n_short,seed,len_dts,len_n_rows,n_measures,realizs,select_model,n_params_target_model
    integer :: params,decimals
    double precision :: D2,var,dt,dt_integration,theta,x0,t0,dt_bridge,LK_EN,LK_DE,P,dummy,total_time,t_CPU
    double precision :: E_D_EN,E_D_EN_2,E_mu_EN,E_K_EN
    double precision :: E_D_DE,E_D_DE_2,E_mu_DE,E_K_DE
    double precision, dimension(3) :: params_target_model
    character(len=10) :: arg
    character(len=100) :: addname,filename
    character(len=3) :: model,model_bridge,target_model
    double precision :: dran_g, dran_gamma
    real :: start,finish

    
    call assingments()
    seed = time()
    call get_command_argument(1, arg)
    read(arg, *) select_model
    print *, "select_model=", select_model
    print *, "seed=", seed
    print *

    call dran_ini(int(seed))
    realizs = 100
    params = 2
    mu = 2.0d0 
    K = 1.0d0
    D = 0.2d0 
    D2 = D*D
    var = D2/K/2.0d0
    theta = 0.5d0
    dt_integration = 1.0d-3
    n_params_target_model = 3
    model_bridge = "WI" ! "OU" 
    ! model = "EN" !"EN" !"DE" !"OU" !"CP" !"GGM"
    if (select_model == 1) then
        model = "EN" 
    else if (select_model == 2) then
        model = "DE"
    else if (select_model == 3) then
        model = "OU"
    else if (select_model == 4) then
        model = "CP"
    else if (select_model == 5) then
        model = "GGM"
    else
        print *, "Model not recognized, fatal error,",select_model
        stop
    end if
    print *, "Model: ", trim(model)
    x0 = mu
    t0 = 0.0d0
    len_dts = 3 
    len_n_rows = 7 
    allocate(dts(len_dts),n_rows(len_n_rows),Ps(len_dts))
    call logspace(start = -3.0d0, end = -1.0d0, num = len_dts, base = 10.0d0, result = dts)
    n_rows = [3,10,30,100,300,1000,3000]


    filename = "data/P_distinguish_EN_DE_target_"//trim(model)//"_Gauss_varying_N.dat"
    print *
    print *, "Save data in: ",trim(filename)
    print *, "dts=",dts
    print *, "n_rows=",n_rows
    open(unit=3003, file=trim(filename), iostat=ios,action="write",status="unknown")
        if ( ios /= 0 ) stop "Error opening file name"
    total_time = 0.0d0
    do jj = 1, len_n_rows, 1
        n_measures = n_rows(jj)
        meas = n_measures-1
        allocate(ts(n_measures),xs(n_measures))

        do kk = 1, len_dts, 1
            dt = dts(kk)
            P = 0.0d0
            call cpu_time(start)
            do ll = 1, realizs, 1
                !-------------------Generate trajectory-------------------
                call generate_trajectory_in_array(model,dt,x0,t0,xs,ts,n_measures,dt_integration,model_params=[theta])

                call estimate_parameters_Gaussian_and_QV(E_D_EN,E_D_EN_2,E_mu_EN,E_K_EN,xs,ts,n_measures,model="EN",params_model=[theta])

                call estimate_parameters_Gaussian_and_QV(E_D_DE,E_D_DE_2,E_mu_DE,E_K_DE,xs,ts,n_measures,model="DE",params_model=[theta])

                params_target_model = [E_mu_EN,E_K_EN,E_D_EN]
                target_model = "EN"
                call Evaluate_log_likelihood_time_series_Gaussian_approx(xs,ts,n_measures,LK_EN,target_model,params_target_model,n_params_target_model)

                params_target_model = [E_mu_DE,E_K_DE,E_D_DE]
                target_model = "DE"
                call Evaluate_log_likelihood_time_series_Gaussian_approx(xs,ts,n_measures,LK_DE,target_model,params_target_model,n_params_target_model)

                if (( LK_EN>=LK_DE ).and.(model=="EN")) then
                    P = P + 1
                elseif (( LK_DE>=LK_EN ).and.(model=="DE")) then
                    P = P + 1
                end if
                
            end do
            dummy = P/dble(realizs)
            Ps(kk) = dummy
            call cpu_time(finish)
            t_CPU = finish - start
            total_time = total_time + t_CPU
            print *, "n_measures=",n_measures,"dt=",dt,"P=",dummy,"model=",model,"t_CPU=",t_CPU
            write(3003,*) n_measures,dt,dummy,model,t_CPU
            
        end do
        deallocate(xs,ts) 
    end do
          
    close(3003)  
    print * 
    
    print *, "time elapsed: ",total_time  
    print *, "End of the program",model
end program Model_distinguishability_EN_DE