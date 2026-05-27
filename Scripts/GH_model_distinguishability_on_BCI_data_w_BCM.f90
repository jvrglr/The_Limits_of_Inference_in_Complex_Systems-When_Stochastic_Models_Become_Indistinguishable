program LKs_real_data
    use declarations_module
    use functions_local
    use functions 
    use subroutines
    implicit none
    double precision, allocatable :: ts(:), xs(:), LKs(:)
    integer*8 :: ii, jj, n_measures, n_files, N_bridges
    character(len=200) :: data_file, filename
    character(len=32) :: addname
    character(len=3) :: target_model, bridge_model
    double precision :: dt_bridge, LK, t_CPU
    real :: start, finish
    integer*8 :: seed, select_target_model, n_params_target_model
    double precision :: theta, E_D, E_D_2, E_mu, E_K, P,corr_cut_off,tau
    double precision, dimension(3) :: params_target_model
    integer*8, dimension(3) :: select_bridge_models

    
    call assingments()
    ! seed = time()
    seed = 1994
    print *, "seed=", seed
    print *
    
    call dran_ini(int(seed))
    theta = 0.5d0
    tau = 5.0d0
    dt_bridge = tau/100.0d0
    N_Bridges = 1000
    n_params_target_model = 3
    select_bridge_models = [2,2,1] ! 1=WI, 2=OU, 3=GB
    corr_cut_off = 100.0d0
    P = 0.0d0
    call cpu_time(start)
    
    addname = "BCI"
    n_files = 179 !179

    allocate(Lks(n_files))
    filename = "output_file.dat"
    open(unit=3003, file=trim(filename), iostat=ios,action="write",status="unknown")
    if ( ios /= 0 ) stop "Error opening file name"
    print *, "Save data in: ",trim(filename)
    do jj = 1,3,1 !1,3,1
        select_target_model = jj
        call select_model_name(select_target_model, target_model)
        call select_bridge_model_name(select_bridge_models(jj), bridge_model)
        print *, "Procesing "//trim(target_model)//" model using "//trim(bridge_model)//" bridges ..."
        print *
        do ii = 1, n_files, 1
            if ((ii == 82)) then
                cycle
            end if
            !-------------------Read data-------------------
            data_file = "Data/GH_BCI/ts_xs_"//trim(addname)//"_"//trim(str8(ii))//".dat"
            call read_trajectory_unknown_nrows(trim(data_file), ts, xs, n_measures)
            !-------------------Estimate parameters-------------------
            call estimate_parameters_Gaussian_and_lamperti_QV(E_D,E_D_2,E_mu,E_K,xs,ts,n_measures,model=trim(target_model),params_model=[theta])
            
            params_target_model = [E_mu,E_K,E_D]
            ! print *, "parameters=",params_target_model
            call Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff(ts,xs,n_measures,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,LK,corr_cut_off)
            LKs(ii) = Lk
            deallocate(xs,ts)
            
        end do
    ! print *, "  Finished with=",target_model,"LK[-1]=",LK
    print*
    write(3003,*) Lks
    enddo
    call cpu_time(finish)
    t_CPU = finish - start
    
     
    close(3003) 
    print *, "End program. time elapsed: ",t_CPU 

end program LKs_real_data