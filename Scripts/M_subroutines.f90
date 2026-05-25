module Subroutines

    !Module with public subroutines.
    !Depends on declarations_module and functions.
  
    use declarations_module
    use functions
    use functions_local
    use, intrinsic :: ieee_arithmetic
  
  contains
    
    subroutine trapz(x, y, n, integral)
      ! Integrate y vs x using the trapezoidal rule
      implicit none
      integer(8), intent(in)    :: n
      real(8),    intent(in)    :: x(n), y(n)
      real(8),    intent(out)   :: integral
      integer(8)                :: i

      integral = 0.0d0
      do i = 1, n - 1
        integral = integral + 0.5d0 * (y(i) + y(i+1)) * (x(i+1) - x(i))
      end do
    end subroutine trapz

    subroutine compute_ACF(xs, ts, n, max_lag, num_lags, lag_window, lags, acf)
    !Compute autocorrelation function (ACF) for a discrete time series with inhomogeneous spacing.
    !num_lags = number of lags (dt bins) to compute ACF.
    !max_lag = maximum lag (dt) to compute ACF for.
    !lag_window = width of the lag bin (default is 50% of bin size).
      implicit none
      integer*8, intent(in) :: n, num_lags
      real(8), intent(in) :: xs(n), ts(n)
      real(8), intent(in), optional :: max_lag, lag_window
      real(8), allocatable, intent(out) :: lags(:), acf(:)

      real(8), allocatable :: counts(:)
      real(8) :: variance, used_max_lag, used_lag_window
      real(8) :: x_mean, dt, delta_t
      integer*8 :: i,j, k, lag_idx
      real(8) :: lag_spacing
      real(8), allocatable :: x_centered(:)
      real(8) :: abs_diff

      ! Compute default max_lag if not provided
      if (present(max_lag)) then
        used_max_lag = max_lag
      else
        used_max_lag = maxval(ts) - minval(ts)
      end if

      if (present(lag_window)) then
        used_lag_window = lag_window
      else
        used_lag_window = used_max_lag / (2.0d0 * num_lags)  ! Default lag window is 10% of max_lag
      end if

      ! Allocate output arrays
      allocate(lags(num_lags), acf(num_lags), counts(num_lags))
      lags = [(used_max_lag * dble(i - 1) / dble(num_lags - 1), i = 1, num_lags)]
      acf = 0.0d0
      counts = 0.0d0

      ! Center x
      x_mean = sum(xs) / dble(n)
      allocate(x_centered(n))
      x_centered = xs - x_mean
      ! Variance
      variance = sum(x_centered**2) / dble(n)

      ! Main loop over all pairs
      do j = 1, n
        do k = j, n
          delta_t = ts(k) - ts(j)
          if (delta_t > used_max_lag) cycle

          ! Find nearest lag bin (dt label)
          lag_idx = 1
          do while (lag_idx < num_lags)
            if (lags(lag_idx + 1) > delta_t) exit
            lag_idx = lag_idx + 1
          end do

          acf(lag_idx) = acf(lag_idx) + x_centered(j) * x_centered(k)
          counts(lag_idx) = counts(lag_idx) + 1.0d0

        end do
      end do
      
      do j = 1, num_lags ! Normalize
        if (counts(j) > 0.0d0) then
          acf(j) = acf(j) / (counts(j) * variance)
        end if
      end do
      acf(1) = 1.0d0  ! Ensure ACF[0] = 1

      deallocate(x_centered)
    end subroutine compute_ACF

    subroutine log_gamma_pdf(x, a, b, log_fx)
      implicit none
      double precision, intent(in)  :: x     ! evaluation point
      double precision, intent(in)  :: a     ! shape parameter
      double precision, intent(in)  :: b     ! rate parameter (b = 1 / scale)
      double precision, intent(out) :: log_fx! log PDF value

      double precision :: log_gamma_val
      double precision, parameter :: pi = 3.141592653589793d0
      double precision :: tmp

      ! Handle invalid inputs
      if (x <= 0.0d0 .or. a <= 0.0d0 .or. b <= 0.0d0) then
          log_fx = -1.0d300  ! return very small value instead of -Inf
          return
      end if

      ! Compute log(Gamma(a)) safely
      if (a < 171.0d0) then
          ! safe to use intrinsic GAMMA
          log_gamma_val = log(gamma(a))
      else
          ! Stirling’s approximation for log(Gamma(a))
          log_gamma_val = a * log(a) - a + 0.5d0 * log(2.0d0 * pi / a)
      end if

      ! Compute log of gamma PDF
      log_fx = a * log(b) + (a - 1.0d0) * log(x) - b * x - log_gamma_val
    end subroutine log_gamma_pdf

    subroutine gamma_pdf(x, a, b, fx)
    ! Gamma PDF with rate parameter b (i.e., f(x; a, b) where b = 1/β)
    ! f(x; a, b) = b**a * x**(a-1) * exp(-b*x) / Gamma(a)
      implicit none
      double precision, intent(in)  :: x     ! evaluation point
      double precision, intent(in)  :: a     ! shape parameter (alpha)
      double precision, intent(in)  :: b     ! rate parameter (b = 1 / beta)
      double precision, intent(out) :: fx    ! PDF value at x

      if (x <= 0.0d0 .or. a <= 0.0d0 .or. b <= 0.0d0) then
          fx = 0.0d0
          return
      end if

      fx = b**a * x**(a - 1.0d0) * exp(-b * x) / gamma(a)
    end subroutine gamma_pdf

    subroutine read_trajectory_unknown_nrows(filename, ts, xs, NN)
    ! Reads a two-column file with unknown number of rows.
    ! Allocates ts and xs arrays and returns them along with number of rows.
      implicit none
      character(len=*), intent(in) :: filename
      integer*8, intent(out) :: NN
      double precision, allocatable, intent(out) :: ts(:), xs(:)
      integer*8 :: ii
      character(len=256) :: line

      ! First pass: count number of rows
      NN = 0
      open(unit=5005, file=filename, status='old', action='read')
      do
          read(5005,'(A)', iostat=ios) line
          if (ios /= 0) exit
          NN = NN + 1
      end do
      close(5005)

      ! Allocate arrays
      allocate(ts(NN))
      allocate(xs(NN))

      ! Second pass: read data
      open(unit=5005, file=filename, status='old', action='read')
      do ii = 1, NN
          read(5005,*) ts(ii), xs(ii)
      end do
      close(5005)

    end subroutine read_trajectory_unknown_nrows

    subroutine write_trajectory(filename, ts, xs, NN)
    ! Writes ts and xs arrays to a two-column file
    implicit none
    character(len=*), intent(in) :: filename
    integer*8, intent(in) :: NN
    double precision, intent(in) :: ts(:), xs(:)
    integer*8 :: ii
    integer :: ios

    ! Basic check
    if (size(ts) /= NN .or. size(xs) /= NN) then
        print*, "Error: array sizes do not match NN"
        stop
    end if

    ! Open file (overwrite)
    open(unit=5006, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
        print*, "Error opening file for writing:", filename
        stop
    end if

    ! Write data
    do ii = 1, NN
        write(5006,*) ts(ii), xs(ii)
    end do

    close(5006)

end subroutine write_trajectory

    subroutine select_bridge_model_name(select_model, model)
    ! Given an integer code (select_model), this subroutine sets the correspondingmodel name as a string. If the code is unrecognized, the program stops with an error message.
        implicit none
        integer*8, intent(in) :: select_model
        character(len=*), intent(out) :: model
        character(len=10) :: tmp_model

        select case (select_model)
        case (1)
            tmp_model = "WI"
        case (2)
            tmp_model = "OU"
        case (3)
            tmp_model = "GB"
        case default
            print *, "Model not recognized, fatal error,", select_model
            stop
        end select

        model = tmp_model(1:len_trim(model))  ! Ensures assignment to same-length output

    end subroutine select_bridge_model_name

    subroutine select_model_name(select_model, model)
    ! Given an integer code (select_model), this subroutine sets the correspondingmodel name as a string. If the code is unrecognized, the program stops with an error message.
        implicit none
        integer*8, intent(in) :: select_model
        character(len=*), intent(out) :: model
        character(len=10) :: tmp_model

        select case (select_model)
        case (1)
            tmp_model = "OU"
        case (2)
            tmp_model = "DE"
        case (3)
            tmp_model = "EN"
        case (4)
            tmp_model = "CP"
        case (5)
            tmp_model = "GGM"
        case default
            print *, "Model not recognized, fatal error,", select_model
            stop
        end select

        model = tmp_model(1:len_trim(model))  ! Ensures assignment to same-length output

    end subroutine select_model_name

    subroutine read_trajectory(filename,ts,xs,NN)
    !Read trajectory from filename, storing times and process in ts and xs respectively. NN is the number of rows in the data file.
      implicit none
      character(len=*), intent(in) :: filename
      integer*8, intent(in) :: NN
      double precision, dimension(NN), intent(inout) :: ts,xs
      integer*8 :: ii
      
      ! print *, "open data in: ",filename
      ! print *
      open(unit=5005, file=filename, status='old', action='read')
      do ii = 1, NN, 1
          read(5005,*) ts(ii),xs(ii)
      end do
      close(5005) 

    end subroutine read_trajectory

    subroutine generate_trajectory_in_array(model, dt_sampling, x0, t0, xs, ts, len, dt_integration, model_params)
    ! Generate trajectory of specified model with given parameters.
    ! x, t, mu, K, and D inputs are global variables.
      implicit none
      character(len=*), intent(in) :: model
      double precision, intent(in) :: dt_sampling, t0, x0
      double precision, intent(in), optional :: dt_integration
      integer*8, intent(in) :: len
      double precision, intent(in), optional :: model_params(:) ! Assumed-shape array
      double precision, dimension(len), intent(out) :: xs, ts
      double precision :: next_t, theta
      integer*8 :: ii

      ! Initialize global state
      t = t0
      x = x0
      xs(1) = x0

      if (model=="EN" .or. model=="CP" .or. model=="GGM") then
          if (.not. present(dt_integration)) then
              print *, "Fata error in subroutine: generate_trajectory_in_array"
              print *, "Error: dt_integration must be provided for models EN, CP, and GGM"
              print *, "model=", model
              stop
          else if (dt_integration <= 0.0d0) then
            print *, "Fata error in subroutine: generate_trajectory_in_array"
              print *, "Error: dt_integration must be positive"
              print *, "dt_integration=", dt_integration
              stop
          else if (dt_sampling <dt_integration) then
              print *, "Fata error in subroutine: generate_trajectory_in_array"
              print *, "Error: dt_sampling must be greater than or equal to dt_integration"
              print *, "dt_sampling=", dt_sampling, "dt_integration=", dt_integration
              stop
          endif
      end if

      if (model == "GGM") then
          if (.not. present(model_params)) then
              print*, "Fata error in subroutine: generate_trajectory_in_array"
              print *, "Error: model_params must be provided for GGM model"
              stop
          else if (size(model_params) /= 1) then
            print*, "Fata error in subroutine: generate_trajectory_in_array"
              print *, "Error: model_params must be a single value for GGM model"
              stop
          endif
      end if

      if (model == "EN") then
          do ii = 2, len
              next_t = t + dt_sampling
              call Update_system_env(next_t, dt_integration)
              ts(ii) = t
              xs(ii) = x
          end do

      elseif (model == "CP") then
          do ii = 2, len
              next_t = t + dt_sampling
              call Update_system_CP(next_t, dt_integration)
              ts(ii) = t
              xs(ii) = x
          end do

      elseif (model == "OU") then
          do ii = 2, len
              next_t = t + dt_sampling
              call Update_system_OU(next_t)
              ts(ii) = t
              xs(ii) = x
          end do

      elseif (model == "DE") then
          do ii = 2, len
              next_t = t + dt_sampling
              call Update_system_dem_exact(next_t)
              ts(ii) = t
              xs(ii) = x
          end do

      elseif (model == "GGM") then
          theta = model_params(1)
          do ii = 2, len
              next_t = t + dt_sampling
              call Update_system_general_gamma_model(next_t, dt_integration, theta)
              ts(ii) = t
              xs(ii) = x
          end do

      else
          print *, "Error: model not recognized"
          stop
      end if

    end subroutine generate_trajectory_in_array

    subroutine Final_states_montecarlo_in_array(model, dt_sampling, x0, t0, xs, len, dt_integration, model_params)
    ! Generate len realizations of the process at time t0+dt_sampling conditioned to X_t0 = xo.
    ! x, t, mu, K, and D inputs are global variables.
      implicit none
      character(len=*), intent(in) :: model
      double precision, intent(in) :: dt_sampling, t0, x0
      double precision, intent(in), optional :: dt_integration
      integer*8, intent(in) :: len
      double precision, intent(in), optional :: model_params(:) ! Assumed-shape array
      double precision, dimension(len), intent(out) :: xs
      double precision :: next_t, theta
      integer*8 :: ii

      ! Initialize global state
      t = t0
      x = x0
      next_t = t + dt_sampling
      if (model == "EN") then
          do ii = 1, len
              call Update_system_env(next_t, dt_integration)
              xs(ii) = x
              t = t0
              x = x0
          end do
      elseif (model == "CP") then
          do ii = 1, len
              call Update_system_CP(next_t, dt_integration)
              xs(ii) = x
              t = t0
              x = x0
          end do

      elseif (model == "OU") then
          do ii = 1, len
              call Update_system_OU(next_t)
              xs(ii) = x
              t = t0
              x = x0
          end do

      elseif (model == "DE") then
          do ii = 1, len
              call Update_system_dem_exact(next_t)
              xs(ii) = x
              t = t0
              x = x0  
          end do
      
      elseif (model == "BE") then ! This generates trajectories of DE model where K=0 (Bessel process) (exact solution breaks)
          do ii = 1, len
              call Update_system_dem_approx(next_t, dt_integration)
              xs(ii) = x
              t = t0
              x = x0
          end do

      elseif (model == "GGM") then
          if (.not. present(model_params)) then
              print *, "Error: model_params must be provided for GGM model"
              stop
          endif
          theta = model_params(1)
          do ii = 1, len
              call Update_system_general_gamma_model(next_t, dt_integration, theta)
              xs(ii) = x
              t = t0
              x = x0
          end do

      else
          print *, "Error: model not recognized"
          stop
      end if
    end subroutine Final_states_montecarlo_in_array

    subroutine estimate_parameters_Gaussian(E_D,E_D_2,E_mu,E_K,xs,ts,len,model,params_model)
    !Estimate three parameters of model using Gaussian approximation of propagator.
      implicit none
      double precision, intent(out) :: E_mu,E_K,E_D,E_D_2
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: xs,ts
      character(len=*), intent(in) :: model
      double precision, intent(in),optional :: params_model(:) ! Assumed-shape array
      integer*8 :: meas 
      double precision, dimension(len-1) :: g_2s,xis,dxs,fs
      double precision :: dt,pa,pb,theta,dummy,tfmt0,dmeas
      double precision :: c,c1,c2,d1,d2,a,b,Z

      dt = ts(2)-ts(1)
      tfmt0 = ts(len)-ts(1)
      if (dt <= 0.0d0) then
        print *, "Error: dt must be positive","dt=", dt
        print *, "Location: estimate_parameters_Gaussian"
        stop
      end if
      if (len < 2) then
        print *, "Error: len must be greater than 1"
        print *, "Location: estimate_parameters_Gaussian"
        stop
      end if
      meas = len-1 ! number of jumps
      dmeas = dble(meas)
      xis = xs(1:meas)
      dxs = xs(2:len) - xis  !Differences of the process
      if ( model=="OU" ) then
        pa = 0.0d0
        pb = 1.0d0
        g_2s = 1.0d0
      else if ( model=="DE" ) then
        pa = 0.0d0
        pb = 1.0d0
        g_2s = xis
      else if ( model=="EN" ) then
        pa = 1.0d0
        pb = 2.0d0
        g_2s = xis**2.0d0
      else if ( model=="GGM") then
        theta = params_model(1)
        dummy = theta+1.0d0
        pa = theta
        pb = dummy
        g_2s = xis**(dummy)
      else if (model == "CP") then
        pa = 1.0d0
        pb = 2.0d0
        g_2s = xis
      end if
      c = sum(xis**(pa+pb)/g_2s)/dmeas
      c1 = sum(xis**(2*pa)/g_2s)/dmeas
      c2 = sum(xis**(2*pb)/g_2s)/dmeas
      d1 = sum(dxs*xis**(pa)/g_2s)/tfmt0
      d2 = sum(dxs*xis**(pb)/g_2s)/tfmt0
      Z = c**2.0d0 - c1*c2
      a = (c*d2-c2*d1)/Z
      b = (c*d1-c1*d2)/Z
      E_mu = -a/b
      E_K = -b
      fs = a*xis**pa + b*xis**pb
      E_D_2 = sum((fs*dt-dxs)**2.0d0/g_2s)/dmeas/dt
      ! E_D_2 = sum((dxs)**2.0d0/g_2s)/dmeas/dt
      E_D = sqrt(E_D_2)
    end subroutine estimate_parameters_Gaussian

    subroutine estimate_parameters_Gaussian_and_QV(E_D,E_D_2,E_mu,E_K,xs,ts,len,model,params_model)
    !Estimate mu and K parameters of model using Gaussian approximation of propagator, while D is estimated using quadratic variation.
      implicit none
      double precision, intent(out) :: E_mu,E_K,E_D,E_D_2
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: xs,ts
      double precision, dimension(len-1) :: dts
      character(len=*), intent(in) :: model
      double precision, intent(in),optional :: params_model(:) ! Assumed-shape array
      integer*8 :: meas,ii
      double precision, dimension(len-1) :: g_2s,xis,dxs
      double precision :: pa,pb,theta,dummy,tfmt0,dmeas
      double precision :: c,c1,c2,d1,d2,a,b,Z

      if (len < 2) then
        print *, "Error: len must be greater than 1"
        print*, "Fatal error in subroutine; estimate_parameters_Gaussian_and_QV"
        stop
      end if

      meas = len-1 ! number of jumps
      dmeas = dble(meas)
      dts = ts(2:len) - ts(1:meas)
      tfmt0 = ts(len)-ts(1)
      do ii = 1, meas, 1
        if (dts(ii) <= 0.0d0) then
        print *, "Error: dt must be positive, ii=",ii,"dt=",dts(ii)
        print*, "Fatal error in subroutine; estimate_parameters_Gaussian_and_QV"
        stop
      end if
      end do
      
      if (model == "GGM") then
          if (.not. present(params_model)) then
              print*, "Fatal error in subroutine: generate_trajectory_in_array"
              print *, "Error: params_model must be provided for GGM model"
              stop
          else if (size(params_model) < 1) then
              print*, "Fatal error in subroutine: generate_trajectory_in_array"
              print *, "Error: params_model must have at least one value for GGM model"
              stop
          endif
      end if
      
      xis = xs(1:meas)
      dxs = xs(2:len) - xis  !Differences of the process
      select case (trim(model))
      case ("OU")
          pa = 0.0d0
          pb = 1.0d0
          g_2s = 1.0d0
      case ("DE")
          pa = 0.0d0
          pb = 1.0d0
          g_2s = xis
      case ("EN")
          pa = 1.0d0
          pb = 2.0d0
          g_2s = xis**2.0d0
      case ("GGM")
          theta = params_model(1)
          dummy = theta + 1.0d0
          pa = theta
          pb = dummy
          g_2s = xis**dummy
      case ("CP")
          pa = 1.0d0
          pb = 2.0d0
          g_2s = xis
      case default
          print*, "Fatal error in subroutine; estimate_parameters_Gaussian_and_QV"
          print *, "Fatal error: model not recognized in g_2s selection"
          stop
      end select
      c  = sum(xis**(pa+pb)/g_2s)/dmeas
      c1 = sum(xis**(2*pa)/g_2s)/dmeas
      c2 = sum(xis**(2*pb)/g_2s)/dmeas
      d1 = sum(dxs*xis**(pa)/g_2s)/tfmt0
      d2 = sum(dxs*xis**(pb)/g_2s)/tfmt0
      Z = c**2.0d0 - c1*c2
      if (abs(Z) < 1.0d-14) then
          print*, "Fatal error in subroutine; estimate_parameters_Gaussian_and_QV"
          print *, "Error: Ill-conditioned system in Gaussian parameter estimation (Z too small)"
          stop
      end if
      a = (c*d2-c2*d1)/Z
      b = (c*d1-c1*d2)/Z
      if (abs(b) < 1.0d-14) then
          print*, "Fatal error in subroutine; estimate_parameters_Gaussian_and_QV"
          print*, "Fatal error: Division by zero in mu estimate (b ~ 0)"
          stop
      end if
      E_mu = -a/b
      E_K = -b
      ! E_D_2 = sum((dxs)**2.0d0/g_2s/dts)/dmeas
      E_D_2 = sum((dxs)**2.0d0/g_2s/dts)/dmeas
      E_D = sqrt(E_D_2)
    end subroutine estimate_parameters_Gaussian_and_QV

    subroutine estimate_parameters_Gaussian_and_Lamperti_QV(E_D,E_D_2,E_mu,E_K,xs,ts,len,model,params_model)
    !Estimate mu and K parameters of model using Gaussian approximation of propagator, while D is estimated using quadratic variation.
      implicit none
      double precision, intent(out) :: E_mu,E_K,E_D,E_D_2
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: xs,ts
      double precision, dimension(len-1) :: dts
      character(len=*), intent(in) :: model
      double precision, intent(in),optional :: params_model(:) ! Assumed-shape array
      integer*8 :: meas,ii
      double precision, dimension(len-1) :: g_2s,xis,dxs,dys
      double precision, dimension(len) :: ys
      double precision :: pa,pb,theta,dummy,tfmt0,dmeas,dummy2
      double precision :: c,c1,c2,d1,d2,a,b,Z

      if (len < 2) then
        print *, "Error: len must be greater than 1"
        print *, "Fatal error in subroutine; estimate_parameters_Gaussian_and_Lamperti_QV"
        stop
      end if

      meas = len-1 ! number of jumps
      dmeas = dble(meas)
      dts = ts(2:len) - ts(1:meas)
      tfmt0 = ts(len)-ts(1)
      do ii = 1, meas, 1
        if (dts(ii) <= 0.0d0) then
        print *, "Error: dt must be positive, ii=",ii,"dt=",dts(ii)
        print *, "Fatal error in subroutine; estimate_parameters_Gaussian_and_Lamperti_QV"
        stop
      end if
      end do
      
      if (model == "GGM") then
          if (.not. present(params_model)) then
              print*, "Fatal error in subroutine: generate_trajectory_in_array"
              print *, "Error: params_model must be provided for GGM model"
              stop
          else if (size(params_model) < 1) then
              print*, "Fatal error in subroutine: generate_trajectory_in_array"
              print *, "Error: params_model must have at least one value for GGM model"
              stop
          endif
      end if
      
      xis = xs(1:meas)
      dxs = xs(2:len) - xis  !Differences of the process
      select case (trim(model))
      case ("OU")
          pa = 0.0d0
          pb = 1.0d0
          g_2s = 1.0d0
          ys = xs
      case ("DE")
          pa = 0.0d0
          pb = 1.0d0
          g_2s = xis
          ys = 2.0d0*sqrt(xs)
      case ("EN")
          pa = 1.0d0
          pb = 2.0d0
          g_2s = xis**2.0d0
          ys = log(xs)
      case ("GGM")
          theta = params_model(1)
          dummy = theta + 1.0d0
          pa = theta
          pb = dummy
          g_2s = xis**dummy
          dummy2 = (1.0d0-theta)/2.0d0
          ys = (xs**dummy2)/dummy2
      case ("CP")
          pa = 1.0d0
          pb = 2.0d0
          g_2s = xis
          ys = 2.0d0*sqrt(xs)
      case default
          print *, "Fatal error in subroutine; estimate_parameters_Gaussian_and_Lamperti_QV"
          print *, "Fatal error: model not recognized in g_2s selection"
          stop
      end select
      dys = ys(2:len) - ys(1:meas)
      c  = sum(xis**(pa+pb)/g_2s)/dmeas
      c1 = sum(xis**(2*pa)/g_2s)/dmeas
      c2 = sum(xis**(2*pb)/g_2s)/dmeas
      d1 = sum(dxs*xis**(pa)/g_2s)/tfmt0
      d2 = sum(dxs*xis**(pb)/g_2s)/tfmt0
      Z = c**2.0d0 - c1*c2
      if (abs(Z) < 1.0d-300) then
          print *, "Fatal error in subroutine; estimate_parameters_Gaussian_and_Lamperti_QV"
          print *, "Error: Ill-conditioned system in Gaussian parameter estimation (Z too small), Z= ",Z
          stop
      end if
      a = (c*d2-c2*d1)/Z
      b = (c*d1-c1*d2)/Z
      if (abs(b) < 1.0d-14) then
          print *, "Fatal error in subroutine; estimate_parameters_Gaussian_and_Lamperti_QV"
          print*, "Fatal error: Division by zero in mu estimate (b ~ 0)"
          stop
      end if
      E_mu = -a/b
      E_K = -b
      E_D_2 = sum((dys)**2.0d0/dts)/dmeas
      E_D = sqrt(E_D_2)
    end subroutine estimate_parameters_Gaussian_and_Lamperti_QV

    subroutine estimate_parameters_Gaussian_and_QV_w_corr_cutoff(E_D,E_D_2,E_mu,E_K,xs,ts,len,model,params_model,cutoff_dt)
      implicit none
      double precision, intent(out) :: E_mu,E_K,E_D,E_D_2
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: xs, ts
      character(len=*), intent(in) :: model
      double precision, intent(in), optional :: params_model(:)
      double precision, intent(in) :: cutoff_dt

      integer*8 :: ii, n_valid, meas
      double precision, allocatable :: dts(:), dxs(:), xis(:), g_2s(:), mask(:)
      double precision, allocatable :: filt_dxs(:), filt_xis(:), filt_dts(:), filt_g2s(:)
      double precision :: pa, pb, theta, dummy, dmeas, tfmt0
      double precision :: c, c1, c2, d1, d2, a, b, Z

      if (len < 2) then
        print *, "Error: len must be greater than 1"
        stop
      end if

      meas = len - 1
      allocate(dts(meas), dxs(meas), xis(meas), g_2s(meas), mask(meas))
      dts = ts(2:len) - ts(1:meas)
      dxs = xs(2:len) - xs(1:meas)
      xis = xs(1:meas)

      do ii = 1, meas
        if (dts(ii) <= 0.0d0) then
          print *, "Error: dt must be positive, ii=", ii, "dt=", dts(ii)
          stop
        end if
      end do

      ! Model-dependent quantities
      select case (trim(model))
      case ("OU")
        pa = 0.0d0; pb = 1.0d0
        g_2s = 1.0d0
      case ("DE")
        pa = 0.0d0; pb = 1.0d0
        g_2s = xis
      case ("EN")
        pa = 1.0d0; pb = 2.0d0
        g_2s = xis**2.0d0
      case ("GGM")
        if (.not. present(params_model)) then
          print *, "Error: params_model must be provided for GGM model"
          stop
        end if
        if (size(params_model) < 1) then
          print *, "Error: params_model must have at least one value for GGM model"
          stop
        end if
        theta = params_model(1)
        dummy = theta + 1.0d0
        pa = theta; pb = dummy
        g_2s = xis**dummy
      case ("CP")
        pa = 1.0d0; pb = 2.0d0
        g_2s = xis
      case default
        print *, "Fatal error: model not recognized"
        stop
      end select

      ! Apply cutoff filter
      mask = 0.0d0
      n_valid = 0
      do ii = 1, meas
        if (dts(ii) < cutoff_dt) then
          n_valid = n_valid + 1
          mask(n_valid) = dble(ii)
        end if
      end do

      if (n_valid < 2) then
        print *, "Error: Not enough data points after applying cutoff_dt"
        stop
      end if

      allocate(filt_dxs(n_valid), filt_xis(n_valid), filt_dts(n_valid), filt_g2s(n_valid))
      do ii = 1, n_valid
        filt_dxs(ii) = dxs(int(mask(ii)))
        filt_xis(ii) = xis(int(mask(ii)))
        filt_dts(ii) = dts(int(mask(ii)))
        filt_g2s(ii) = g_2s(int(mask(ii)))
      end do

      dmeas = dble(n_valid)
      tfmt0 = sum(filt_dts)

      c  = sum(filt_xis**(pa+pb)/filt_g2s) / dmeas
      c1 = sum(filt_xis**(2*pa)/filt_g2s) / dmeas
      c2 = sum(filt_xis**(2*pb)/filt_g2s) / dmeas
      d1 = sum(filt_dxs * filt_xis**pa / filt_g2s) / tfmt0
      d2 = sum(filt_dxs * filt_xis**pb / filt_g2s) / tfmt0

      Z = c**2.0d0 - c1 * c2
      if (abs(Z) < 1.0d-14) then
        print *, "Error: Ill-conditioned system (Z ≈ 0)"
        stop
      end if

      a = (c*d2 - c2*d1) / Z
      b = (c*d1 - c1*d2) / Z
      if (abs(b) < 1.0d-14) then
        print *, "Error: Division by zero in mu estimate (b ≈ 0)"
        stop
      end if

      E_mu = -a / b
      E_K = -b
      E_D_2 = sum((filt_dxs)**2.0d0 / filt_g2s / filt_dts) / dmeas
      E_D = sqrt(E_D_2)

      ! Deallocate
      deallocate(dts, dxs, xis, g_2s, mask, filt_dxs, filt_xis, filt_dts, filt_g2s)

    end subroutine estimate_parameters_Gaussian_and_QV_w_corr_cutoff

    subroutine estimate_parameters_segmented(xs, ts, len, cutoff_dt, model, params_model, E_mu, E_K, E_D, E_D_2)
      implicit none
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: xs, ts
      double precision, intent(in) :: cutoff_dt
      character(len=*), intent(in) :: model
      double precision, intent(in), optional :: params_model(:)
      double precision, intent(out) :: E_mu, E_K, E_D, E_D_2

      integer*8 :: i, start_idx, end_idx, nseg, seg_len, total_meas
      double precision :: dt
      double precision :: mu_k, K_k, D_k, D2_k
      double precision :: weight, sum_wt, sum_mu, sum_K, sum_D, sum_D2
      double precision, allocatable :: seg_xs(:), seg_ts(:)

      sum_mu = 0.0d0
      sum_K  = 0.0d0
      sum_D  = 0.0d0
      sum_D2 = 0.0d0
      sum_wt = 0.0d0
      total_meas = 0
      start_idx = 1

      do i = 2, len
        dt = ts(i) - ts(i-1)

        if (dt >= cutoff_dt .or. i == len) then
          end_idx = i - 1
          if (end_idx - start_idx + 1 >= 2) then
            seg_len = end_idx - start_idx + 1
            allocate(seg_xs(seg_len), seg_ts(seg_len))
            seg_xs = xs(start_idx:end_idx)
            seg_ts = ts(start_idx:end_idx)

            ! Call the original estimator on the segment
            if (present(params_model)) then
              call estimate_parameters_Gaussian_and_QV(D_k, D2_k, mu_k, K_k, seg_xs, seg_ts, seg_len, model, params_model)
            else
              call estimate_parameters_Gaussian_and_QV(D_k, D2_k, mu_k, K_k, seg_xs, seg_ts, seg_len, model)
            end if

            weight = dble(seg_len - 1)
            sum_wt = sum_wt + weight
            total_meas = total_meas + seg_len - 1

            sum_mu  = sum_mu + weight * mu_k
            sum_K   = sum_K  + weight * K_k
            sum_D   = sum_D  + weight * D_k
            sum_D2  = sum_D2 + weight * D2_k

            deallocate(seg_xs, seg_ts)
          end if
          start_idx = i
        end if
      end do

      if (sum_wt <= 0.0d0) then
        print *, "Error: No valid segments found under cutoff_dt"
        stop
      end if

      E_mu  = sum_mu / sum_wt
      E_K   = sum_K  / sum_wt
      E_D   = sum_D  / sum_wt
      E_D_2 = sum_D2 / sum_wt
    end subroutine estimate_parameters_segmented

    subroutine Evaluate_likelihood_time_series_w_bridge_CM(ts,xs,len_time_series,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,LK_estimator)
    !Evaluate the likelihood of a time series using importance sampling with bridge measure.
      implicit none
      integer*8, intent(in) :: N_bridges,n_params_target_model,len_time_series
      double precision, intent(in) :: dt_bridge
      double precision, dimension(len_time_series), intent(in) :: ts,xs
      character(len=*), intent(in) :: target_model,bridge_model
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      double precision, intent(out) :: LK_estimator
      integer*8 :: ii
      double precision :: t0,tf,x0,xf,estimator,S_err,R_err

      
      if (len_time_series < 2) then
          print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM"
          print*, "Error: len_time_series must be at least 2"
          stop
      end if
      
      LK_estimator = 1.0d0
      t0 = ts(1)
      x0 = xs(1)
      ! print *, "   ---> t0=",t0,"x0=",x0,"params=",params_target_model,"dt_bridge=",dt_bridge,"N_bridges=",N_bridges,"bridge_model=",bridge_model
      do ii = 2,len_time_series
        tf = ts(ii)
        xf = xs(ii)
        call Evaluate_propagator_w_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err)
        LK_estimator = LK_estimator*estimator
        t0 = tf
        x0 = xf
      enddo
      
    end subroutine Evaluate_likelihood_time_series_w_bridge_CM

    subroutine Evaluate_log_likelihood_time_series_w_bridge_CM(ts,xs,len_time_series,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,LK_estimator)
    !Evaluate the normalized log-likelihood of a time series using importance sampling with bridge measure.
      implicit none
      integer*8, intent(in) :: N_bridges,n_params_target_model,len_time_series
      double precision, intent(in) :: dt_bridge
      double precision, dimension(len_time_series), intent(in) :: ts,xs
      character(len=*), intent(in) :: target_model,bridge_model
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      double precision, intent(out) :: LK_estimator
      integer*8 :: ii
      double precision :: t0,tf,x0,xf,estimator,S_err,R_err

      if (len_time_series < 2) then
          print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM"
          print*, "Error: len_time_series must be at least 2"
          stop
      end if
      
      LK_estimator = 0.0d0
      t0 = ts(1)
      x0 = xs(1)

      do ii = 2,len_time_series
        tf = ts(ii)
        xf = xs(ii)
        call Evaluate_propagator_w_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err)
        if (estimator <= 0.0d0) then
            print*
            print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM"
            print*, "Error: estimator <= 0 at ii=", ii, "estimator=", estimator
            print*, "t0=", t0, "tf=", tf, "x0=", x0, "xf=", xf
            stop
        end if
        LK_estimator = LK_estimator + log(estimator)
        t0 = tf
        x0 = xf
      enddo
      LK_estimator = LK_estimator/(tf-ts(1))
      ! print*, "-->Estimator = ",LK_estimator,"model=",target_model
    end subroutine Evaluate_log_likelihood_time_series_w_bridge_CM

    subroutine Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff(ts,xs,len_time_series,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,LK_estimator,corr_cut_off,stop_if_lt_zero)
    !Evaluate the normalized log-likelihood of a time series using importance sampling with bridge measure.
      implicit none
      integer*8, intent(in) :: N_bridges,n_params_target_model,len_time_series
      double precision, intent(in) :: dt_bridge,corr_cut_off
      double precision, dimension(len_time_series), intent(in) :: ts,xs
      character(len=*), intent(in) :: target_model,bridge_model
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      double precision, intent(out) :: LK_estimator
      logical, optional, intent(in) :: stop_if_lt_zero
      integer*8 :: ii
      double precision :: t0,tf,x0,xf,estimator,S_err,R_err

      if (len_time_series < 2) then
          print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff"
          print*, "Error: len_time_series must be at least 2"
          stop
      end if
      
      LK_estimator = 0.0d0
      t0 = ts(1)
      x0 = xs(1)

      do ii = 2,len_time_series
        tf = ts(ii)
        xf = xs(ii)
        if ((tf-t0)<= corr_cut_off) then
          if (( bridge_model == "WI" ).or.( bridge_model == "OU" )) then
            call Evaluate_propagator_w_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err)
          else if ( bridge_model == "GB" ) then
            call Evaluate_propagator_w_M1_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err)
          else
            print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff"
            print*, "Error: bridge model not recognized"
            stop
          end if
          
          if (estimator <= 0.0d0) then    !Flag an error if estimator is negative or zero
              if (present(stop_if_lt_zero)) then
                if (stop_if_lt_zero) then
                  print*
                  print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff"
                  print*, "Error: estimator <= 0 at ii=", ii, "estimator=", estimator
                  stop
                endif
              else
                print*
                print*, "Fatal error in Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff"
                print*, "Error: estimator <= 0 at ii=", ii, "estimator=", estimator
                stop
              end if
          end if
        else 
          estimator = 1.0d0 ! If the time difference is larger than corr_cut_off, we assume no correlation and set estimator to 1 (has no effect on log-likelihood)
        end if
        LK_estimator = LK_estimator + log(estimator)
        t0 = tf
        x0 = xf
      enddo
      LK_estimator = LK_estimator/(tf-ts(1))
      ! print*, "-->Estimator = ",LK_estimator,"model=",target_model
    end subroutine Evaluate_log_likelihood_time_series_w_bridge_CM_w_corr_cutoff

    subroutine Evaluate_propagator_w_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err,compute_R_err)
    !Estimate the propagator of a process using importance sampling with bridge measure.
    !Bridges are Gaussian bridges of a Lamperti-transformed process
    !Outcomes are the estimated propagator, the statistical error (S_err) and the relative error (R_err) when exact result is known.
    !Lamperti transformation is performed to compute Radon-Nikodym derivatives using Girsanov theorem.
    !Available models for bridges are: "WI", "OU", "GB"
    !Available models for target models are: "OU", "DE", "EN"
    !SHOULD COMPUTE R_err as an optional variable
      implicit none
      double precision, parameter :: EXPONENT_OVERFLOW_LIMIT = 700.0d0
      integer*8, intent(in) :: N_bridges,n_params_target_model
      double precision, intent(in) :: t0,tf,x0,xf,dt_bridge
      double precision, intent(out) :: estimator, R_err, S_err
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      character(len=*), intent(in) :: bridge_model,target_model
      integer*8 :: ii,len_bridge, N_Bessel
      logical, intent(in), optional :: compute_R_err
      double precision :: tfmt0,t0_bridge,tf_bridge,x0_bridge,xf_bridge,prop_true
      double precision, allocatable :: xs_bridge(:),ts_bridge(:)
      
      double precision :: l,prop_B,RN,RN2,Jacobian,l_ref,l_target
      double precision :: mu_OU,K_OU,D_OU,D_OU_2,dummy
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2
      double precision :: mu_bridge,K_bridge,D_bridge
      integer*8, parameter :: n_params_bridge_model = 3 
      double precision, dimension(n_params_bridge_model) :: params_bridge_model!Needed for EN model

      ! print*, "  ---> t0=",t0,"tf=",tf,"x0=",x0,"xf=",xf,"params_target_model=",params_target_model,"dt_bridge=",dt_bridge,"N_bridges=",N_bridges,"bridge_model=",bridge_model
      R_err = ieee_value(R_err, ieee_quiet_nan) ! By default, set relative error to NaN
      tfmt0 = tf-t0
      if (dt_bridge>tfmt0) then
        !"dt_bridge is larger than tf-t0--> no bridges at all"
        len_bridge = 2
      else
        len_bridge = int8(tfmt0/dt_bridge)+1
      endif
      t0_bridge = t0
      tf_bridge = tf

      if (x0 <= 0.0d0 .or. xf <= 0.0d0) then
        print *, "Fatal error in Evaluate_propagator_w_bridge_CM"
        print *, "Error: x0 and xf must be positive for EN/DE/CP/GGM models"
        stop
      end if

      allocate(xs_bridge(len_bridge),ts_bridge(len_bridge)) 

      ! ----- Parameters of the target model and Lamperti transform -----
      !proceeding in this way, I compute lamperti transform twice per point :: ENHANCE EFFICIENCY
        if ( trim(target_model) == "OU" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, OU model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_OU = params_target_model(1) 
          K_OU = params_target_model(2)  
          D_OU = params_target_model(3) 
          D_OU_2 = D_OU*D_OU
          x0_bridge = x0/D_OU ! Lamperti transformation
          xf_bridge = xf/D_OU ! Lamperti transformation
          Jacobian = 1.0d0/D_OU ! Lamperti transformation    
          if (present(compute_R_err)) then
            if (compute_R_err) then
              call propagator_OU(mu_OU,K_OU,D_OU_2,x0,xf,tfmt0,prop_true) ! Compute true value of propagator
            end if
          end if
          
        else if ( trim(target_model) == "DE" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, DE model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_DE = params_target_model(1) 
          K_DE = params_target_model(2)  
          D_DE = params_target_model(3) 
          D_DE_2 = D_DE*D_DE
          x0_bridge = 2.0d0*sqrt(x0)/D_DE ! Lamperti transformation
          xf_bridge = 2.0d0*sqrt(xf)/D_DE ! Lamperti transformation
          Jacobian = 1.0d0/sqrt(xf)/D_DE ! Lamperti transformation 
          if (present(compute_R_err)) then
            if (compute_R_err) then
               N_Bessel = 10000
              call propagator_dem(mu_DE,K_DE,D_DE_2,x0,xf,tfmt0,N_Bessel,prop_true) ! Compute true value of propagator
            end if
          end if
          
        else if ( trim(target_model) == "EN" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, EN model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_EN = params_target_model(1) 
          K_EN = params_target_model(2)  
          D_EN = params_target_model(3) 
          D_EN_2 = D_EN*D_EN
          x0_bridge = log(x0)/D_EN ! Lamperti transformation
          xf_bridge = log(xf)/D_EN ! Lamperti transformation
          Jacobian = 1.0d0 / (xf * D_EN) ! Lamperti transformation       
        else 
          print *, "Fatal error, target model not recognized"
          stop
        end if
      ! ----- _______________________________ -----
      
      ! ----- Parameters of the bridge model and RN bridge-unconstrained -----
        if ( trim(bridge_model)=="WI" ) then
          call propagator_Unbiased_Brownian(1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)  ! Evaluate RN of unconstrained process wrt its bridge
          do ii = 1, n_params_bridge_model, 1 ! Browanianian bridge does not depend on parameters
            params_bridge_model(ii) = ieee_value(params_bridge_model(ii) , ieee_quiet_nan)
          end do
        else if ( trim(bridge_model)=="GB" ) then
          k_bridge = log(xf_bridge/x0_bridge)/tfmt0 
          D_bridge = 1.0d0 
          ! mu_bridge = ieee_value(mu_bridge , ieee_quiet_nan) ! geometric Browanianian bridge only has two parameters
          params_bridge_model(1:2) = [K_bridge,D_bridge]
          call propagator_GB(k_bridge,1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)  ! Evaluate RN of unconstrained process wrt its bridge
        else if ( trim(bridge_model)=="OU" ) then
          mu_bridge = (x0_bridge+xf_bridge)/2.0d0 
          K_bridge = params_target_model(2)
          D_bridge = 1.0d0 
          params_bridge_model = [mu_bridge,K_bridge,D_bridge]
          call propagator_OU(mu_bridge,K_bridge,1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)
        else 
          print *, "Fatal error, bridge model not recognized"
          stop
        end if
      ! ----- _______________________________ -----
      
      estimator = 0.0d0
      S_err = 0.0d0
      do ii = 1,N_bridges,1
        
        call Fill_gaps_with_lamperti_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,len_bridge,xs_bridge,ts_bridge,bridge_model,n_params_bridge_model,params_bridge_model)

        ! Evaluate log-likelihood OU/Wiener USING EXACT EXPRESSIONS OF STOCHASTIC INTEGRALS
        call Evaluate_log_RN_lamperti_W_Girsanov(target_model,params_target_model,n_params_target_model,len_bridge,ts_bridge,xs_bridge,l_target)
        call Evaluate_log_RN_lamperti_W_Girsanov(bridge_model,params_bridge_model,n_params_bridge_model,len_bridge,ts_bridge,xs_bridge,l_ref   )
        l = l_target - l_ref !Chain rule of Radon-Nikodym derivatives
        dummy = tfmt0*l
        ! Compute estimator of propagator
        ! print*, "----> l=",l,"l_target=",l_target,"l_ref=",l_ref,"dummy=",dummy
        if (dummy > EXPONENT_OVERFLOW_LIMIT) then
          print*
          
          print*, "Fatal error in Evaluate_propagator_w_bridge_CM"
          print *, "Error: exp(tfmt0*l) would overflow"
          print *, "==================== DEBUG INFO ===================="
          print *, "Bridge filling and likelihood evaluation diagnostics"
          print *, "---------------------------------------------------"
          print *, "Time interval:"
          print *, "  t0_bridge  =", t0_bridge
          print *, "  tf_bridge  =", tf_bridge
          print *, "  tfmt0      =", tfmt0
          print *, "Lamperti transformed states States:"
          print *, "  x0_bridge  =", x0_bridge
          print *, "  xf_bridge  =", xf_bridge
          print *, "States:"
          print *, "  x0  =", x0
          print *, "  xf  =", xf
          print *, "Bridge parameters:"
          print *, "  len_bridge =", len_bridge
          print *, "  bridge_model =", trim(bridge_model)
          print *, "  target_model =", trim(target_model)
          print *, "  n_params_bridge_model =", n_params_bridge_model
          print *, "  params_bridge_model =", params_bridge_model
          print *, "  n_params_target_model =", n_params_target_model
          print *, "  params_target_model =", params_target_model
          print *, "Generated path (ts_bridge and xs_bridge):"
          print *, "Computed log-likelihoods:"
          print *, "  l_target =", l_target
          print *, "  l_ref    =", l_ref
          print *, "  l = l_target - l_ref =", l
          print *, "  tfmt0 * l =", dummy
          print *, "==================================================="
          stop
        endif
        RN = exp(dummy)*prop_B !Radon-Nikodym of OU wrt Wiener bridge 
        RN2 = RN*RN
        estimator = estimator + RN
        S_err = S_err + RN2  
      enddo 
        
      estimator = estimator/N_bridges*Jacobian ! Inverting Lamperti transformation
      S_err = S_err/N_bridges*(Jacobian*Jacobian) 
      S_err = sqrt(abs(S_err-estimator**2.0d0)/N_bridges) !statistical error
      R_err = abs(prop_true-estimator)/prop_true !Relative error
      deallocate(xs_bridge,ts_bridge) 
    end subroutine Evaluate_propagator_w_bridge_CM

    subroutine Evaluate_propagator_w_M1_bridge_CM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err,compute_R_err)
    !Estimate the propagator of a process using importance sampling with bridge measure.
    !Bridges are Gaussian bridges of a Lamperti-transformed process
    !Outcomes are the estimated propagator, the statistical error (S_err) and the relative error (R_err) when exact result is known.
    !Lamperti transformation is performed to compute Radon-Nikodym derivatives using Girsanov theorem.
    !Available models for bridges are: "WI", "OU", "GB"
    !Available models for target models are: "OU", "DE", "EN"
    !SHOULD COMPUTE R_err as an optional variable
      implicit none
      double precision, parameter :: EXPONENT_OVERFLOW_LIMIT = 700.0d0
      integer*8, intent(in) :: N_bridges,n_params_target_model
      double precision, intent(in) :: t0,tf,x0,xf,dt_bridge
      double precision, intent(out) :: estimator, R_err, S_err
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      character(len=*), intent(in) :: bridge_model,target_model
      integer*8 :: ii,len_bridge, N_Bessel
      logical, intent(in), optional :: compute_R_err
      double precision :: tfmt0,t0_bridge,tf_bridge,x0_bridge,xf_bridge,prop_true
      double precision, allocatable :: xs_bridge(:),ts_bridge(:)
      
      double precision :: l,prop_B,RN,RN2,Jacobian,l_ref,l_target
      double precision :: mu_OU,K_OU,D_OU,D_OU_2,dummy
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2,inv_D_EN
      double precision :: mu_GB,K_GB,D_GB,D_GB_2,inv_D_GB
      double precision :: mu_bridge,K_bridge,D_bridge
      integer*8, parameter :: n_params_bridge_model = 3 
      double precision, dimension(n_params_bridge_model) :: params_bridge_model!Needed for EN model

      R_err = ieee_value(R_err, ieee_quiet_nan) ! By default, set relative error to NaN
      tfmt0 = tf-t0
      if (dt_bridge>tfmt0) then
        !"dt_bridge is larger than tf-t0--> no bridges at all"
        len_bridge = 2
      else
        len_bridge = int8(tfmt0/dt_bridge)+1
      endif
      t0_bridge = t0
      tf_bridge = tf

      if (x0 <= 0.0d0 .or. xf <= 0.0d0) then
        print *, "Fatal error in Evaluate_propagator_w_bridge_CM"
        print *, "Error: x0 and xf must be positive for EN/DE/CP/GGM models"
        stop
      end if

      allocate(xs_bridge(len_bridge),ts_bridge(len_bridge)) 

      ! ----- Parameters of the target model and M1-transform -----
      !proceeding in this way, I compute M1-transform twice per point :: ENHANCE EFFICIENCY
       if ( trim(target_model) == "DE" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, DE model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_DE = params_target_model(1) 
          K_DE = params_target_model(2)  
          D_DE = params_target_model(3) 
          D_DE_2 = D_DE*D_DE
          x0_bridge = exp(2.0d0*sqrt(x0)/D_DE) ! M1-transformation
          xf_bridge = exp(2.0d0*sqrt(xf)/D_DE) ! M1-transformation
          Jacobian = xf_bridge/sqrt(xf)/D_DE ! M1-transformation 
          if ((present(compute_R_err)).and.(compute_R_err)) then
            N_Bessel = 10000
            call propagator_dem(mu_DE,K_DE,D_DE_2,x0,xf,tfmt0,N_Bessel,prop_true) ! Compute true value of propagator
          end if
        else if ( trim(target_model) == "GB") then
          if (n_params_target_model/=2) then
            print *, "Fatal error, GB model requires 2 parameters"
            stop
          endif
          ! Parameters of the target model 
          K_GB = params_target_model(1)  
          D_GB = params_target_model(2) 
          D_GB_2 = D_GB*D_GB
          inv_D_GB = 1.0d0/D_GB
          x0_bridge = x0**inv_D_GB ! M1-transformation
          xf_bridge = xf**inv_D_GB ! M1-transformation
          Jacobian = xf_bridge/xf*inv_D_GB ! Jacobian to invert M1-transformation     
          if ((present(compute_R_err)).and.(compute_R_err)) then
            call propagator_GB(K_GB,D_GB_2,x0,xf,tfmt0,prop_true) ! Compute true value of propagator
          end if
        else if ( trim(target_model) == "EN" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, EN model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_EN = params_target_model(1) 
          K_EN = params_target_model(2)  
          D_EN = params_target_model(3) 
          D_EN_2 = D_EN*D_EN
          inv_D_EN = 1.0d0/D_EN
          x0_bridge = x0**inv_D_EN ! M1-transformation
          xf_bridge = xf**inv_D_EN ! M1-transformation  
          Jacobian = xf_bridge/xf*inv_D_GB ! Jacobian to invert M1-transformation  
        else 
          print *, "Fatal error, target model not recognized : ",trim(target_model)
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1"
          stop
        end if
      ! ----- _______________________________ -----
      
      ! ----- Parameters of the bridge model and RN bridge-unconstrained -----
        if ( trim(bridge_model)=="GB" ) then
          k_bridge = log(xf_bridge/x0_bridge)/tfmt0 
          D_bridge = 1.0d0 
          ! mu_bridge = ieee_value(mu_bridge , ieee_quiet_nan) ! geometric Browanianian bridge only has two parameters
          params_bridge_model(1:2) = [K_bridge,D_bridge]
          call propagator_GB(k_bridge,1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)  ! Evaluate RN of unconstrained process wrt its bridge
        else if ((trim(bridge_model)=="OU").or.(trim(bridge_model)=="WI")) then
          print *, "Fatal error, bridge model cannot be OU or WI as I am using transforms to M1--> M1-Multiplicative reference measures only"
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1"
          stop
        else 
          print *, "Fatal error, bridge model not recognized"
          stop
        end if
      ! ----- _______________________________ -----
      
      estimator = 0.0d0
      S_err = 0.0d0
      do ii = 1,N_bridges,1
        
        call Fill_gaps_with_lamperti_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,len_bridge,xs_bridge,ts_bridge,bridge_model,n_params_bridge_model,params_bridge_model)

        ! Evaluate log-likelihood OU/Wiener USING EXACT EXPRESSIONS OF STOCHASTIC INTEGRALS
        call Evaluate_log_RN_transformed_W_M1(target_model,params_target_model,n_params_target_model,len_bridge,ts_bridge,xs_bridge,l_target)
        call Evaluate_log_RN_transformed_W_M1(bridge_model,params_bridge_model,n_params_bridge_model,len_bridge,ts_bridge,xs_bridge,l_ref   )
        l = l_target - l_ref !Chain rule of Radon-Nikodym derivatives
        dummy = tfmt0*l
        ! Compute estimator of propagator
        ! print*, "----> l=",l,"l_target=",l_target,"l_ref=",l_ref,"dummy=",dummy
        if (dummy > EXPONENT_OVERFLOW_LIMIT) then
          print*
          
          print*, "Fatal error in Evaluate_propagator_w_bridge_CM"
          print *, "Error: exp(tfmt0*l) would overflow"
          print *, "==================== DEBUG INFO ===================="
          print *, "Bridge filling and likelihood evaluation diagnostics"
          print *, "---------------------------------------------------"
          print *, "Time interval:"
          print *, "  t0_bridge  =", t0_bridge
          print *, "  tf_bridge  =", tf_bridge
          print *, "  tfmt0      =", tfmt0
          print *, "Lamperti transformed states States:"
          print *, "  x0_bridge  =", x0_bridge
          print *, "  xf_bridge  =", xf_bridge
          print *, "States:"
          print *, "  x0  =", x0
          print *, "  xf  =", xf
          print *, "Bridge parameters:"
          print *, "  len_bridge =", len_bridge
          print *, "  bridge_model =", trim(bridge_model)
          print *, "  target_model =", trim(target_model)
          print *, "  n_params_bridge_model =", n_params_bridge_model
          print *, "  params_bridge_model =", params_bridge_model
          print *, "  n_params_target_model =", n_params_target_model
          print *, "  params_target_model =", params_target_model
          print *, "Generated path (ts_bridge and xs_bridge):"
          print *, "Computed log-likelihoods:"
          print *, "  l_target =", l_target
          print *, "  l_ref    =", l_ref
          print *, "  l = l_target - l_ref =", l
          print *, "  tfmt0 * l =", dummy
          print *, "==================================================="
          stop
        endif
        RN = exp(dummy)*prop_B !Radon-Nikodym of OU wrt Wiener bridge 
        RN2 = RN*RN
        estimator = estimator + RN
        S_err = S_err + RN2  
      enddo 
        
      estimator = estimator/N_bridges*Jacobian ! Inverting Lamperti transformation
      S_err = S_err/N_bridges*(Jacobian*Jacobian) 
      S_err = sqrt(abs(S_err-estimator**2.0d0)/N_bridges) !statistical error
      R_err = abs(prop_true-estimator)/prop_true !Relative error
      deallocate(xs_bridge,ts_bridge) 
    end subroutine Evaluate_propagator_w_M1_bridge_CM

    subroutine Draw_lks_on_arrays_and_estimate_propagator_W_BCM(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err,compute_R_err,lk_target,lk_ref)
    !THIS SUBROUTINE SHOULD BE USED WITH ADDITIVE REFERENCE PROCESSES (OU-WI)
    !IT SHOULD BE RENAMED AS Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M0
    !Estimate the propagator of a process using importance sampling with bridge measure.
    !Bridges are Gaussian bridges of a Lamperti-transformed process
    !Outcomes are the estimated propagator, the statistical error (S_err) and the relative error (R_err) when exact result is known. Also, arrays of likelihoods for target and reference model on every bridge path are provided as output
    !Lamperti transformation is performed to compute Radon-Nikodym derivatives using Girsanov theorem.
    !Available models for bridges are: "WI", "OU"
    !Available models for target models are: "OU", "DE", "EN"
    !SHOULD COMPUTE R_err as an optional variable
      implicit none
      double precision, parameter :: EXPONENT_OVERFLOW_LIMIT = 700.0d0
      integer*8, intent(in) :: N_bridges,n_params_target_model
      double precision, intent(in) :: t0,tf,x0,xf,dt_bridge
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      character(len=*), intent(in) :: bridge_model,target_model
      logical, intent(in), optional :: compute_R_err
      double precision, intent(out) :: estimator, R_err, S_err
      double precision, intent(out) :: lk_target(N_bridges),lk_ref(N_bridges)
      double precision :: tfmt0,t0_bridge,tf_bridge,x0_bridge,xf_bridge,prop_true
      double precision, allocatable :: xs_bridge(:),ts_bridge(:)
      integer*8 :: ii,len_bridge, N_Bessel
      double precision :: l,prop_B,RN,RN2,Jacobian,l_ref,l_target
      double precision :: mu_OU,K_OU,D_OU,D_OU_2,dummy
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2
      double precision :: mu_GB,K_GB,D_GB,D_GB_2
      double precision :: mu_bridge,K_bridge,D_bridge
      integer*8, parameter :: n_params_bridge_model = 3 
      double precision, dimension(n_params_bridge_model) :: params_bridge_model!Needed for EN model

      ! print*, "  ---> t0=",t0,"tf=",tf,"x0=",x0,"xf=",xf,"params_target_model=",params_target_model,"dt_bridge=",dt_bridge,"N_bridges=",N_bridges,"bridge_model=",bridge_model
      R_err = ieee_value(R_err, ieee_quiet_nan) ! By default, set relative error to NaN
      tfmt0 = tf-t0
      if (dt_bridge>tfmt0) then
        print *, "Fatal error, dt_bridge is larger than tf-t0"
        stop
      endif

      t0_bridge = t0
      tf_bridge = tf
      len_bridge = int8(tfmt0/dt_bridge)+1
      
      allocate(xs_bridge(len_bridge),ts_bridge(len_bridge)) 

      ! ----- Parameters of the target model and Lamperti transform -----
      !proceeding in this way, I compute lamperti transform twice per point :: ENHANCE EFFICIENCY
        if ( trim(target_model) == "OU" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, OU model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_OU = params_target_model(1) 
          K_OU = params_target_model(2)  
          D_OU = params_target_model(3) 
          D_OU_2 = D_OU*D_OU
          x0_bridge = x0/D_OU ! Lamperti transformation
          xf_bridge = xf/D_OU ! Lamperti transformation
          Jacobian = 1.0d0/D_OU ! Lamperti transformation    
          if (present(compute_R_err)) then
            if (compute_R_err) then
              call propagator_OU(mu_OU,K_OU,D_OU_2,x0,xf,tfmt0,prop_true) ! Compute true value of propagator
            end if
          end if
        
        else if ( trim(target_model) == "GB" ) then
          if (n_params_target_model/=2) then
            print *, "Fatal error, GB model requires 2 parameters"
            stop
          endif
          ! Parameters of the target model 
          K_GB = params_target_model(1)  
          D_GB = params_target_model(2) 
          D_GB_2 = D_GB*D_GB
          x0_bridge = log(x0)/D_GB ! Lamperti transformation
          xf_bridge = log(xf)/D_GB ! Lamperti transformation
          Jacobian = 1.0d0 / (xf * D_GB) ! Lamperti transformation     
          if (present(compute_R_err)) then
            if (compute_R_err) then
              call propagator_GB(K_GB,D_GB_2,x0,xf,tfmt0,prop_true) ! Compute true value of propagator
            end if
          end if
          
        else if ( trim(target_model) == "DE" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, DE model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_DE = params_target_model(1) 
          K_DE = params_target_model(2)  
          D_DE = params_target_model(3) 
          D_DE_2 = D_DE*D_DE
          x0_bridge = 2.0d0*sqrt(x0)/D_DE ! Lamperti transformation
          xf_bridge = 2.0d0*sqrt(xf)/D_DE ! Lamperti transformation
          Jacobian = 1.0d0/sqrt(xf)/D_DE ! Lamperti transformation 
          if (present(compute_R_err)) then
            if (compute_R_err) then
               N_Bessel = 10000
              call propagator_dem(mu_DE,K_DE,D_DE_2,x0,xf,tfmt0,N_Bessel,prop_true) ! Compute true value of propagator
            end if
          end if
          
        else if ( trim(target_model) == "EN" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, EN model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_EN = params_target_model(1) 
          K_EN = params_target_model(2)  
          D_EN = params_target_model(3) 
          D_EN_2 = D_EN*D_EN
          x0_bridge = log(x0)/D_EN ! Lamperti transformation
          xf_bridge = log(xf)/D_EN ! Lamperti transformation
          Jacobian = 1.0d0 / (xf * D_EN) ! Lamperti transformation       
        else 
          print *, "Fatal error, target model not recognized"
          stop
        end if
      ! ----- _______________________________ -----
      
      ! ----- Parameters of the bridge model and RN bridge-unconstrained -----
        if ( trim(bridge_model)=="WI" ) then
          call propagator_Unbiased_Brownian(1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)  ! Evaluate RN of unconstrained process wrt its bridge
          do ii = 1, n_params_bridge_model, 1
            params_bridge_model(ii) = ieee_value(params_bridge_model(ii) , ieee_quiet_nan)
          end do
        else if ( trim(bridge_model)=="OU" ) then
          mu_bridge = (x0_bridge+xf_bridge)/2.0d0 
          K_bridge = params_target_model(2)
          D_bridge = 1.0d0 
          params_bridge_model = [mu_bridge,K_bridge,D_bridge]
          call propagator_OU(mu_bridge,K_bridge,1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)
        else if ( trim(bridge_model)=="GB" ) then
          print *, "Fatal error, bridge model cannot be GB as I am using lamperti transforms--> Additive reference measures only"
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM"
          stop

        else 
          print *, "Fatal error, bridge model not recognized"
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM"
          stop
        end if
      ! ----- _______________________________ -----
      
      estimator = 0.0d0
      S_err = 0.0d0
      do ii = 1,N_bridges,1
        call Fill_gaps_with_lamperti_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,len_bridge,xs_bridge,ts_bridge,bridge_model,n_params_bridge_model,params_bridge_model)

        ! Evaluate log-likelihood OU/Wiener USING EXACT EXPRESSIONS OF STOCHASTIC INTEGRALS
        call Evaluate_log_RN_lamperti_W_Girsanov(target_model,params_target_model,n_params_target_model,len_bridge,ts_bridge,xs_bridge,l_target)
        call Evaluate_log_RN_lamperti_W_Girsanov(bridge_model,params_bridge_model,n_params_bridge_model,len_bridge,ts_bridge,xs_bridge,l_ref   )
        lk_ref(ii) = l_ref
        lk_target(ii) = l_target
        l = l_target - l_ref !Chain rule of Radon-Nikodym derivatives
        dummy = tfmt0*l
        ! Compute estimator of propagator
        
        if (dummy > EXPONENT_OVERFLOW_LIMIT) then
          print *, "Error: exp(tfmt0*l_OU) would overflow"
          stop
        endif
        RN = exp(dummy)*prop_B !Radon-Nikodym of OU wrt Wiener bridge 
        RN2 = RN*RN
        estimator = estimator + RN
        S_err = S_err + RN2  
      enddo 
          ! print *, "==================== DEBUG INFO 2 ===================="
          ! print *, "Bridge filling and likelihood evaluation diagnostics"
          ! print *, "---------------------------------------------------"
          ! print *, "Time interval:"
          ! print *, "  t0_bridge  =", t0_bridge
          ! print *, "  tf_bridge  =", tf_bridge
          ! print *, "  tfmt0      =", tfmt0
          ! print *, "Lamperti transformed states States:"
          ! print *, "  x0_bridge  =", x0_bridge
          ! print *, "  xf_bridge  =", xf_bridge
          ! print *, "States:"
          ! print *, "  x0  =", x0
          ! print *, "  xf  =", xf
          ! print *, "Bridge parameters:"
          ! print *, "  len_bridge =", len_bridge
          ! print *, "  bridge_model =", trim(bridge_model)
          ! print *, "  target_model =", trim(target_model)
          ! print *, "  n_params_bridge_model =", n_params_bridge_model
          ! print *, "  params_bridge_model =", params_bridge_model
          ! print *, "  n_params_target_model =", n_params_target_model
          ! print *, "  params_target_model =", params_target_model
          ! print *, "Generated path (ts_bridge and xs_bridge):"
          ! print *, "Computed log-likelihoods:"
          ! print *, "  l_target =", l_target
          ! print *, "  l_ref    =", l_ref
          ! print *, "  l = l_target - l_ref =", l
          ! print *, "  tfmt0 * l =", dummy
          ! print *, "==================================================="
      estimator = estimator/N_bridges*Jacobian ! Inverting Lamperti transformation
      S_err = S_err/N_bridges*(Jacobian*Jacobian) 
      S_err = sqrt(abs(S_err-estimator**2.0d0)/N_bridges) !statistical error
      R_err = abs(prop_true-estimator)/prop_true !Relative error
      deallocate(xs_bridge,ts_bridge) 
    end subroutine Draw_lks_on_arrays_and_estimate_propagator_W_BCM

    subroutine Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1(t0,tf,x0,xf,target_model,params_target_model,n_params_target_model,dt_bridge,N_bridges,bridge_model,estimator,S_err,R_err,compute_R_err,lk_target,lk_ref)
    !THIS SUBROUTINE SHOULD BE USED WITH multiplicative M1 references (e.g. GB)
    !Estimate the propagator of a process using importance sampling with bridge measure.
    !Bridges are geometric Brownian bridges of a transformed process
    !Outcomes are the estimated propagator, the statistical error (S_err) and the relative error (R_err) when exact result is known. Also, arrays of likelihoods for target and reference model on every bridge path are provided as output
    !Transformation of data to get M1-noise is performed to compute Radon-Nikodym derivatives using generalized Girsanov theorem.
    !Available models for bridges are: "GM"
    !Available models for target models are: "DE"
    !SHOULD COMPUTE R_err as an optional variable
      implicit none
      double precision, parameter :: EXPONENT_OVERFLOW_LIMIT = 700.0d0
      integer*8, intent(in) :: N_bridges,n_params_target_model
      double precision, intent(in) :: t0,tf,x0,xf,dt_bridge
      double precision, dimension(n_params_target_model), intent(in) :: params_target_model
      character(len=*), intent(in) :: bridge_model,target_model
      logical, intent(in), optional :: compute_R_err
      double precision, intent(out) :: estimator, R_err, S_err
      double precision, intent(out) :: lk_target(N_bridges),lk_ref(N_bridges)
      double precision :: tfmt0,t0_bridge,tf_bridge,x0_bridge,xf_bridge,prop_true
      double precision, allocatable :: xs_bridge(:),ts_bridge(:)
      integer*8 :: ii,len_bridge, N_Bessel
      double precision :: l,prop_B,RN,RN2,Jacobian,l_ref,l_target
      double precision :: mu_OU,K_OU,D_OU,D_OU_2,dummy
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2,inv_D_EN 
      double precision :: mu_GB,K_GB,D_GB,D_GB_2,inv_D_GB
      double precision :: mu_bridge,K_bridge,D_bridge
      integer*8, parameter :: n_params_bridge_model = 3 
      double precision, dimension(n_params_bridge_model) :: params_bridge_model!Needed for EN model

      ! print*, "  ---> t0=",t0,"tf=",tf,"x0=",x0,"xf=",xf,"params_target_model=",params_target_model,"dt_bridge=",dt_bridge,"N_bridges=",N_bridges,"bridge_model=",bridge_model
      R_err = ieee_value(R_err, ieee_quiet_nan) ! By default, set relative error to NaN
      tfmt0 = tf-t0
      if (dt_bridge>tfmt0) then
        print *, "Fatal error, dt_bridge is larger than tf-t0"
        stop
      endif

      t0_bridge = t0
      tf_bridge = tf
      len_bridge = int8(tfmt0/dt_bridge)+1
      
      allocate(xs_bridge(len_bridge),ts_bridge(len_bridge)) 

      ! ----- Parameters of the target model and M1-transform -----
      !proceeding in this way, I compute M1-transform twice per point :: ENHANCE EFFICIENCY
       if ( trim(target_model) == "DE" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, DE model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_DE = params_target_model(1) 
          K_DE = params_target_model(2)  
          D_DE = params_target_model(3) 
          D_DE_2 = D_DE*D_DE
          x0_bridge = exp(2.0d0*sqrt(x0)/D_DE) ! M1-transformation
          xf_bridge = exp(2.0d0*sqrt(xf)/D_DE) ! M1-transformation
          Jacobian = xf_bridge/sqrt(xf)/D_DE ! M1-transformation 
          if (present(compute_R_err)) then
            if (compute_R_err) then
               N_Bessel = 10000
              call propagator_dem(mu_DE,K_DE,D_DE_2,x0,xf,tfmt0,N_Bessel,prop_true) ! Compute true value of propagator
            end if
          end if
        else if ( trim(target_model) == "GB") then
          if (n_params_target_model/=2) then
            print *, "Fatal error, GB model requires 2 parameters"
            stop
          endif
          ! Parameters of the target model 
          K_GB = params_target_model(1)  
          D_GB = params_target_model(2) 
          D_GB_2 = D_GB*D_GB
          inv_D_GB = 1.0d0/D_GB
          x0_bridge = x0**inv_D_GB ! M1-transformation
          xf_bridge = xf**inv_D_GB ! M1-transformation
          Jacobian = xf_bridge/xf*inv_D_GB ! Jacobian to invert M1-transformation     
          if (present(compute_R_err)) then
            if (compute_R_err) then
              call propagator_GB(K_GB,D_GB_2,x0,xf,tfmt0,prop_true) ! Compute true value of propagator
            end if
          end if
          else if ( trim(target_model) == "EN" ) then
          if (n_params_target_model/=3) then
            print *, "Fatal error, EN model requires 3 parameters"
            stop
          endif
          ! Parameters of the target model
          mu_EN = params_target_model(1) 
          K_EN = params_target_model(2)  
          D_EN = params_target_model(3) 
          D_EN_2 = D_EN*D_EN
          inv_D_EN = 1.0d0/D_EN
          x0_bridge = x0**inv_D_EN ! M1-transformation
          xf_bridge = xf**inv_D_EN ! M1-transformation  
          Jacobian = xf_bridge/xf*inv_D_GB ! Jacobian to invert M1-transformation  
        else 
          print *, "Fatal error, target model not recognized : ",trim(target_model)
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1"
          stop
        end if
      ! ----- _______________________________ -----
      
      ! ----- Parameters of the bridge model and RN bridge-unconstrained -----
        if ( trim(bridge_model)=="GB" ) then
          k_bridge = log(xf_bridge/x0_bridge)/tfmt0 
          D_bridge = 1.0d0 
          ! mu_bridge = ieee_value(mu_bridge , ieee_quiet_nan) ! geometric Browanianian bridge only has two parameters
          params_bridge_model(1:2) = [K_bridge,D_bridge]
          call propagator_GB(k_bridge,1.0d0,x0_bridge,xf_bridge,tfmt0,prop_B)  ! Evaluate RN of unconstrained process wrt its bridge
        else if ((trim(bridge_model)=="OU").or.(trim(bridge_model)=="WI")) then
          print *, "Fatal error, bridge model cannot be OU or WI as I am using transforms to M1--> M1-Multiplicative reference measures only"
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1"
          stop
        else 
          print *, "Fatal error, bridge model not recognized"
          stop
        end if
      ! ----- _______________________________ -----
      
      estimator = 0.0d0
      S_err = 0.0d0
      do ii = 1,N_bridges,1
        call Fill_gaps_with_lamperti_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,len_bridge,xs_bridge,ts_bridge,bridge_model,n_params_bridge_model,params_bridge_model)

        ! Evaluate log-likelihood OU/Wiener USING EXACT EXPRESSIONS OF STOCHASTIC INTEGRALS
        call Evaluate_log_RN_transformed_W_M1(target_model,params_target_model,n_params_target_model,len_bridge,ts_bridge,xs_bridge,l_target)
        call Evaluate_log_RN_transformed_W_M1(bridge_model,params_bridge_model,n_params_bridge_model,len_bridge,ts_bridge,xs_bridge,l_ref   )
        lk_ref(ii) = l_ref
        lk_target(ii) = l_target
        l = l_target - l_ref !Chain rule of Radon-Nikodym derivatives
        dummy = tfmt0*l
        ! Compute estimator of propagator
        
        if (dummy > EXPONENT_OVERFLOW_LIMIT) then
          print *
          print *, "Error: exp(tfmt0*l) would overflow"
          print *, "Location: Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1"
          print *
          print *, "==================== DEBUG INFO ===================="
          print *, "Bridge filling and likelihood evaluation diagnostics"
          print *, "---------------------------------------------------"
          print *, "Time interval:"
          print *, "  t0_bridge  = ", t0_bridge
          print *, "  tf_bridge  = ", tf_bridge
          print *, "  tfmt0      = ", tfmt0
          print *, "Lamperti transformed states States:"
          print *, "  x0_bridge  = ", x0_bridge
          print *, "  xf_bridge  = ", xf_bridge
          print *, "States:"
          print *, "  x0  =", x0
          print *, "  xf  =", xf
          print *, "Bridge parameters:"
          print *, "  len_bridge = ", len_bridge
          print *, "  bridge_model = ", trim(bridge_model)
          print *, "  target_model = " , trim(target_model)
          print *, "  n_params_bridge_model = ", n_params_bridge_model
          print *, "  params_bridge_model = ", params_bridge_model
          print *, "  n_params_target_model = ", n_params_target_model
          print *, "  params_target_model = ", params_target_model
          print *, "Generated path (ts_bridge and xs_bridge):"
          print *, "Computed log-likelihoods:"
          print *, "  l_target = ", l_target
          print *, "  l_ref    = ", l_ref
          print *, "  l = l_target - l_ref = ", l
          print *, "  tfmt0 * l = ", dummy
          print *, "==================================================="
          stop
        endif
        RN = exp(dummy)*prop_B !Radon-Nikodym of OU wrt Wiener bridge 
        RN2 = RN*RN
        estimator = estimator + RN
        S_err = S_err + RN2  
      enddo 
          ! print *, "==================== DEBUG INFO 2 ===================="
          ! print *, "Bridge filling and likelihood evaluation diagnostics"
          ! print *, "---------------------------------------------------"
          ! print *, "Time interval:"
          ! print *, "  t0_bridge  =", t0_bridge
          ! print *, "  tf_bridge  =", tf_bridge
          ! print *, "  tfmt0      =", tfmt0
          ! print *, "Lamperti transformed states States:"
          ! print *, "  x0_bridge  =", x0_bridge
          ! print *, "  xf_bridge  =", xf_bridge
          ! print *, "  Jacobian   =", Jacobian
          ! print *, "States:"
          ! print *, "  x0  =", x0
          ! print *, "  xf  =", xf
          ! print *, "Bridge parameters:"
          ! print *, "  len_bridge =", len_bridge
          ! print *, "  bridge_model =", trim(bridge_model)
          ! print *, "  target_model =", trim(target_model)
          ! print *, "  n_params_bridge_model =", n_params_bridge_model
          ! print *, "  params_bridge_model =", params_bridge_model
          ! print *, "  n_params_target_model =", n_params_target_model
          ! print *, "  params_target_model =", params_target_model
          ! print *, "Generated path (ts_bridge and xs_bridge):"
          ! print *, "Computed log-likelihoods:"
          ! print *, "  l_target =", l_target
          ! print *, "  l_ref    =", l_ref
          ! print *, "  l = l_target - l_ref =", l
          ! print *, "  tfmt0 * l =", dummy
          ! print *, "==================================================="
      estimator = estimator/N_bridges*Jacobian ! Inverting Lamperti transformation
      S_err = S_err/N_bridges*(Jacobian*Jacobian) 
      S_err = sqrt(abs(S_err-estimator**2.0d0)/N_bridges) !statistical error
      R_err = abs(prop_true-estimator)/prop_true !Relative error
      deallocate(xs_bridge,ts_bridge) 
    end subroutine Draw_lks_on_arrays_and_estimate_propagator_W_BCM_M1

    subroutine Evaluate_log_RN_lamperti_W_Girsanov(model,params_model,n_params_model,len_data,ts,ys,l)
    !Evaluate normalized log-RN derivative  of target model wrt Wiener process 
    ! model is assumed to be additive with constant diffusion equal to 1, Lamperti transformed
    !Input data is Lamperti transformed process yt
      implicit none
      integer*8, intent(in) :: len_data,n_params_model
      double precision, dimension(n_params_model), intent(in) :: params_model
      character(len=*), intent(in) :: model
      double precision, dimension(len_data), intent(in) :: ts,ys
      double precision, intent(out) :: l
      double precision :: a,b
      double precision :: w_a,w_b,tfmt0,m_y,m_y2,yf,y0
      double precision :: m_inv_y2,br_y,br_inv_y
      double precision :: m_E_Y,m_E_2Y,br_E_Y
      double precision :: mu_OU,K_OU,D_OU,D_OU_2
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2
      double precision :: k_GB,D_GB,v_GB
      integer*8 :: n_jumps

      n_jumps = len_data-1
      y0 = ys(1)
      yf = ys(len_data)
      tfmt0 = ts(len_data)-ts(1) 
      if (tfmt0 <= 0.0d0) then
        print *, "Error: tfmt0 must be positive"
        print *, "Location: Evaluate_log_RN_lamperti_W_Girsanov"
        stop
      end if

      if (trim(model)=="WI")then
        l = 0.0d0
      else if (trim(model)=="GB")then
        k_GB = params_model(1)
        D_GB = params_model(2)
        v_GB = k_GB/D_GB - 0.5d0*D_GB
        l = v_GB*(yf-y0)/tfmt0 - 0.5d0*v_GB**2.0d0
      else if (trim(model)=="OU") then
        mu_OU = params_model(1) 
        K_OU = params_model(2)  
        D_OU = params_model(3) 
        D_OU_2 = D_OU*D_OU
        a = K_OU*mu_OU/D_OU
        b  = K_OU
        m_y = sum(ys(1:n_jumps))/dble(n_jumps)
        m_y2 = sum(ys(1:n_jumps)**2.0d0)/dble(n_jumps)
        w_a = (yf-y0)/tfmt0
        w_b = (yf**2.0-y0**2.0)*0.5d0/tfmt0-0.5d0
        l = w_a*a-w_b*b-(a**2.0d0)*0.5d0-(b**2.0d0)*m_y2*0.5d0+a*b*m_y

      else if (trim(model)=="DE") then
        mu_DE = params_model(1) 
        K_DE = params_model(2)  
        D_DE = params_model(3) 
        D_DE_2 = D_DE*D_DE
        a = 2.0d0*K_DE*mu_DE/D_DE_2-0.5d0 
        b  = K_DE*0.5d0
        m_inv_y2 = sum(ys(1:n_jumps)**(-2.0d0))/dble(n_jumps)
        m_y2 = sum(ys(1:n_jumps)**2.0d0)/dble(n_jumps)
        br_y = (yf**2.0d0 - y0**2.0d0)*0.5d0/tfmt0 -0.5d0
        br_inv_y = log(yf/y0)/tfmt0 + 0.5d0*m_inv_y2
        l = a*br_inv_y - b*br_y - (a**2.0d0)*m_inv_y2*0.5d0-(b**2.0d0)*m_y2*0.5d0+a*b

      else if (trim(model)=="EN") then
        mu_EN = params_model(1) 
        K_EN = params_model(2)  
        D_EN = params_model(3) 
        D_EN_2 = D_EN*D_EN
        a = K_EN*mu_EN/D_EN - D_EN*0.5d0 
        b  = K_EN/D_EN
        m_E_Y = sum(exp(ys(1:n_jumps)*D_EN))/dble(n_jumps)
        m_E_2Y = sum(exp(2.0d0*ys(1:n_jumps)*D_EN))/dble(n_jumps)
        br_E_Y = (exp(yf*D_EN)-exp(y0*D_EN))/tfmt0/D_EN - D_EN*m_E_Y*0.5d0
        w_a = (yf-y0)/tfmt0
        l = w_a*a - br_E_Y*b - (a**2.0d0)*0.5d0 - (b**2.0d0)*m_E_2Y/2.0d0 + a*b*m_E_Y
     
     else
        print *, "Error: model not recognized in Evaluate_log_RN_lamperti_W_Girsanov"
        print *, "Location: Evaluate_log_RN_lamperti_W_Girsanov"
        stop
      end if

    end subroutine Evaluate_log_RN_lamperti_W_Girsanov

    subroutine Evaluate_log_RN_transformed_W_M1(model,params_model,n_params_model,len_data,ts,ys,l)
    !Evaluate normalized log-RN derivative  of target model wrt M1 process
    ! model is assumed to be multiplicative with constant diffusion equal to 1
    !Input data is transformed process yt
      implicit none
      integer*8, intent(in) :: len_data,n_params_model
      double precision, dimension(n_params_model), intent(in) :: params_model
      character(len=*), intent(in) :: model
      double precision, dimension(len_data), intent(in) :: ts,ys
      double precision, intent(out) :: l
      double precision, dimension(len_data-1):: yys,dys,dts,arg,y2s
      
      double precision :: a,b
      double precision :: w_a,w_b,tfmt0,m_y,m_y2,yf,y0
      double precision :: m_inv_y2,br_y,br_inv_y
      double precision :: m_E_Y,m_E_2Y,br_E_Y
      double precision :: br,ms
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      double precision :: mu_EN,K_EN,D_EN,D_EN_2
      double precision :: k_GB,kk_GB,D_GB,v_GB,D_GB_2
      integer*8 :: n_jumps

      n_jumps = len_data-1
      y0 = ys(1)
      yf = ys(len_data)
      tfmt0 = ts(len_data)-ts(1) 
      yys = ys(1:n_jumps)
      y2s = yys**2.0d0
      dys = ys(2:len_data)-ys(1:n_jumps)
      dts = ts(2:len_data)-ts(1:n_jumps)
      if (tfmt0 <= 0.0d0) then
        print *, "Error: tfmt0 must be positive"
        print *, "Location: Evaluate_log_RN_lamperti_W_M1"
        stop
      end if

      if (trim(model)=="DE") then
        mu_DE = params_model(1) 
        K_DE = params_model(2)  
        D_DE = params_model(3) 
        D_DE_2 = D_DE*D_DE
        ! a = 2.0d0*K_DE*mu_DE/D_DE_2-0.5d0 
        ! b  = K_DE*0.5d0
        a   = 4.0d0*K_DE*mu_DE/D_DE_2 - 1.0d0
        arg = 0.5d0*yys*(a/log(yys) - K_DE*log(yys)+ 1.0d0)
        br  = sum(arg*dys/y2s)/tfmt0
        ms  = sum(arg**2.0d0*dts/y2s)/tfmt0
        l   = br - 0.5d0*ms
      
      else if (trim(model)=="GB")then
        k_GB = params_model(1)
        D_GB = params_model(2)
        D_GB_2 = D_GB*D_GB
        kk_GB  = k_GB/D_GB + 0.5d0*(1.0d0-D_GB)
        v_GB   = 0.5d0*kk_GB*(1.0d0 - kk_GB)
        
        l = kk_GB/tfmt0 *log(yf/y0) + v_GB
      else if (trim(model)=="EN") then
        mu_EN = params_model(1) 
        K_EN = params_model(2)  
        D_EN = params_model(3) 
        a = k_EN*mu_EN/D_EN + (1.0d0-D_EN)*0.5d0
        b = k_EN/D_EN
        arg = yys*(a-b*yys**D_EN)
        br = sum(arg*dys/y2s)/tfmt0
        ms  = sum(arg**2.0d0*dts/y2s)/tfmt0
        l   = br - 0.5d0*ms
     
     else
        print *
        print *, "Error: model not recognized"
        print *, "Location: Evaluate_log_RN_lamperti_W_M1"
        print *
        stop
      end if

    end subroutine Evaluate_log_RN_transformed_W_M1

    subroutine Fill_gaps_with_lamperti_bridge(t0,x0,tf,xf,len_bridge,xs,ts,bridge_model,n_params_bridge,params_bridge)
    !Generate len-1 points as a bridge linking the points 
    !OU, WI --> The bridge of a Gaussian process is itself Gaussian.
    !GB --> lognormal process
    ! Process is a assumed to have diffusion constant = 1 (Lamperti transformed)
    !X(t0)=x0 and X(tf)=xf
    !The method uses a multiresolution algorithm
    !Store the points generated and the input extremes in the array xs.
    !Also created the array of times where the process is generated, ts.
    !The number of new points generated needs to be bigger or equal than 1
    !The grid generated by this algorithm is regular, ts(ii+1)-ts(ii) = cte for all ii.
      implicit none
      integer*8,intent(in) :: len_bridge,n_params_bridge
      double precision,intent(in) :: t0,x0,tf,xf
      double precision,dimension (len_bridge),intent(out) :: xs,ts
      double precision, dimension(n_params_bridge), intent(in) :: params_bridge
      double precision :: D_B,D_B_2
      character(len=*), intent(in) :: bridge_model
      integer*8 :: ii,extra_points
      double precision :: dt_bridge

      if ( (len_bridge < 2)) then
        print *, 'Error: len_bridge must be greater than 2.'
        return
      end if
      extra_points = len_bridge - 2 ! bridge points not present in original data
      dt_bridge = (tf-t0)/(extra_points + 1)
      do ii = 1, len_bridge, 1
        ts(ii) = t0+(ii-1)*dt_bridge
      end do

      xs(1) = x0
      do ii = 2, len_bridge-1, 1
        call Update_system_Lamperti_bridge(ts(ii) ,ts(ii-1),xs(ii-1), tf,xf,bridge_model,n_params_bridge,params_bridge)
        xs(ii) = x
      end do
      xs(len_bridge) = xf

    end subroutine Fill_gaps_with_lamperti_bridge

    subroutine Estimate_OU_propator_w_bridges(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_bridge,params_bridge,n_params_bridge,model_bridge,mu_OU,K_OU,D_OU,prop_estimator,n_bridges)
    !Estimate the propagator of the OU process using importance sampling with bridge measure.
      implicit none
      integer*8, intent(in) :: len_bridge,n_params_bridge,n_bridges
      double precision, intent(in) :: t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge
      double precision, intent(in) :: mu_OU,K_OU,D_OU
      double precision, intent(out) :: prop_estimator
      double precision, intent(in), dimension(n_params_bridge) :: params_bridge
      character(len=*), intent(in) :: model_bridge
      double precision, dimension(len_bridge) :: xs_bridge,ts_bridge
      double precision :: tfmt0_bridge,l_OU,J,prop_OU,RN
      double precision :: D_B,D_B_2,prop_B
      integer*8 :: ii
      
      tfmt0_bridge = tf_bridge-t0_bridge
      if ( model_bridge == "WI" ) then
        if (n_params_bridge/=1) then
          print *, "Fatal error, Wiener bridge only has one parameter"
          stop
        endif
        D_B = params_bridge(1)
        D_B_2 = D_B*D_B
        prop_estimator = 0.0d0
        do ii = 1, n_bridges, 1
          call Fill_gaps_with_Brownian_bridge_multiresolution(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_bridge,D_B,xs_bridge,ts_bridge)
          ! RN Wiener-Wiener bridge
          call propagator_Unbiased_Brownian(D_B_2,x0_bridge,xf_bridge,tfmt0_bridge,prop_B)
          ! RN OU-Wiener
          call Evaluate_ell_OU_EM(mu_OU/D_OU,K_OU,xs_bridge,ts_bridge,len_bridge,l_OU)
          ! Jacobian
          J = 1.0d0/D_OU
          RN = exp(tfmt0_bridge*(l_OU))*prop_B*J
          print *, "  RN=",RN,"l_OU=",l_OU,"prop_B=",prop_B
          prop_estimator = prop_estimator + RN
        end do
        prop_estimator = prop_estimator/dble(n_bridges)
      end if 
    end subroutine Estimate_OU_propator_w_bridges

    subroutine Evaluate_log_likelihood_time_series_w_bridges(xs,ts,len,l,model,params_model,n_params,model_bridge,params_bridge,n_params_bridge,len_bridge,n_bridges)
    !Evaluate the normalized log-likelihood of a time series under a target measure using bridge importance sampling
      implicit none
      integer*8, intent(in) :: len,n_params,n_params_bridge,n_bridges,len_bridge
      double precision, dimension(len), intent(in) :: xs,ts
      double precision, dimension(n_params) :: params_model
      double precision, dimension(n_params_bridge) :: params_bridge
      double precision, intent(out) :: l
      character(len=*), intent(in) :: model,model_bridge
      double precision :: t0_bridge,x0_bridge,tf_bridge,xf_bridge,prop_estimator,tfmt0,dt_bridge
      double precision :: mu_OU,K_OU,D_OU
      integer*8 :: ii,meas,len_dummy_bridge

      meas = len-1  
      t0_bridge = ts(1)
      x0_bridge = xs(1)
      l = 0.0d0
      tfmt0 = ts(len)-ts(1)
      mu_OU = params_model(1)
      K_OU = params_model(2)
      D_OU = params_model(3)
      do ii = 1, meas, 1
        tf_bridge = ts(ii +1)
        xf_bridge = xs(ii +1)
        dt_bridge = (tf_bridge-t0_bridge)/len_bridge
        len_dummy_bridge = int((tf_bridge-t0_bridge)/dt_bridge)+1
        call Estimate_OU_propator_w_bridges(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_dummy_bridge,params_bridge,n_params_bridge,model_bridge,mu_OU,K_OU,D_OU,prop_estimator,n_bridges)
        l = l + log(prop_estimator)

        x0_bridge = xf_bridge
        tf_bridge = ts(ii +1)
      end do
      l = l / tfmt0
    end subroutine Evaluate_log_likelihood_time_series_w_bridges

    subroutine Evaluate_likelihood_time_series_Gaussian_approx(xs,ts,len,Prob,Model,params_model,n_params)
    !Evaluate the likelihood of a time series using the Gaussian approximation.
    !ASSUMES CONSTANT SPACING OF MEASURES
      implicit none
      integer*8, intent(in) :: len,n_params
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, dimension(len-1) ::  dts
      double precision,dimension(n_params) :: params_model
      double precision, intent(out) :: Prob
      character(len=*), intent(in) :: Model
      double precision,dimension(len-1) :: fs,B_2s
      double precision :: dt,sqrt_dt,xi,xip1,D2,mean,var,dummy
      double precision :: mu_OU,K_OU,D_OU,D_OU2
      double precision :: mu_DE,K_DE,D_DE,D_DE2
      double precision :: mu_EN,K_EN,D_EN,D_EN2
      double precision :: mu_CP,K_CP,D_CP,D_CP2
      integer*8 :: ii,meas
      double precision ::  tol
      logical :: is_uniform

      dt = ts(2)-ts(1) !Assuming constant time step
      meas = len-1
      dts = ts(2:len) - ts(1:meas)
      tol = 1.0d-12
      ! Check that all time steps are equal to dt within tolerance
      is_uniform = .true.
      do ii = 1, meas
          if (abs(dts(ii) - dt) > tol) then
              is_uniform = .false.
              exit
          end if
      end do

      if (.not. is_uniform) then
          print *, "Error: time steps are not uniformly spaced."
          do ii = 1, meas
              print *, "dts(", ii, ") =", dts(ii), "  deviation from dt =", abs(dts(ii) - dt)
          end do
          stop
      end if

      sqrt_dt = sqrt(dt)
      xi = xs(1)
      
      if (trim(Model) == "OU") then
          mu_OU = params_model(1)
          K_OU  = params_model(2)
          D_OU  = params_model(3)
          D_OU2 = D_OU * D_OU
          fs    = K_OU * (mu_OU - xs(1:meas))
          B_2s  = D_OU2
      elseif (trim(Model) == "DE") then
          mu_DE = params_model(1)
          K_DE  = params_model(2)
          D_DE  = params_model(3)
          D_DE2 = D_DE * D_DE
          fs    = K_DE * (mu_DE - xs(1:meas))
          B_2s  = D_DE2 * xs(1:meas)
      elseif (trim(Model) == "EN") then
          mu_EN = params_model(1)
          K_EN  = params_model(2)
          D_EN  = params_model(3)
          D_EN2 = D_EN * D_EN
          fs    = K_EN * (mu_EN - xs(1:meas)) * xs(1:meas)
          B_2s  = D_EN2 * xs(1:meas)**2
      elseif (trim(Model) == "CP") then
          mu_CP = params_model(1)
          K_CP  = params_model(2)
          D_CP  = params_model(3)
          D_CP2 = D_CP * D_CP
          fs    = K_CP * (mu_CP - xs(1:meas)) * xs(1:meas)
          B_2s  = D_CP2 * xs(1:meas)
      else
          print *, "Fatal error: model not recognized in Evaluate_log_likelihood_time_series_Gaussian_approx"
          stop
      end if
      prob = 1.0d0
      do ii = 1, meas, 1
        xip1 = xs(ii+1)
        mean = xi+fs(ii)*dt
        var = B_2s(ii)*dt
        dummy = Gauss(mean,var,xip1)
        prob = prob*dummy
        xi = xip1
      end do

    end subroutine Evaluate_likelihood_time_series_Gaussian_approx

    subroutine Evaluate_log_likelihood_time_series_Gaussian_approx(xs,ts,len,lg_prob,Model,params_model,n_params)
    ! Evaluate the normalized log-likelihood of a time series using the Gaussian approximation.
      implicit none
      integer*8, intent(in) :: len, n_params
      double precision, dimension(len), intent(in) :: xs, ts
      double precision, dimension(n_params) :: params_model
      double precision, intent(out) :: lg_prob
      character(len=*), intent(in) :: Model
      double precision, dimension(len-1) :: fs, B_2s, dts
      double precision :: xi, xip1, mean, var, dummy
      double precision :: mu_OU, K_OU, D_OU, D_OU2
      double precision :: mu_DE, K_DE, D_DE, D_DE2
      double precision :: mu_EN, K_EN, D_EN, D_EN2
      double precision :: mu_CP, K_CP, D_CP, D_CP2
      integer*8 :: ii, meas

      meas = len - 1
      dts = ts(2:len) - ts(1:meas)
      xi = xs(1)

      if (trim(Model) == "OU") then
          mu_OU = params_model(1)
          K_OU  = params_model(2)
          D_OU  = params_model(3)
          D_OU2 = D_OU * D_OU
          fs    = K_OU * (mu_OU - xs(1:meas))
          B_2s  = D_OU2
      elseif (trim(Model) == "DE") then
          mu_DE = params_model(1)
          K_DE  = params_model(2)
          D_DE  = params_model(3)
          D_DE2 = D_DE * D_DE
          fs    = K_DE * (mu_DE - xs(1:meas))
          B_2s  = D_DE2 * xs(1:meas)
      elseif (trim(Model) == "EN") then
          mu_EN = params_model(1)
          K_EN  = params_model(2)
          D_EN  = params_model(3)
          D_EN2 = D_EN * D_EN
          fs    = K_EN * (mu_EN - xs(1:meas)) * xs(1:meas)
          B_2s  = D_EN2 * xs(1:meas)**2
      elseif (trim(Model) == "CP") then
          mu_CP = params_model(1)
          K_CP  = params_model(2)
          D_CP  = params_model(3)
          D_CP2 = D_CP * D_CP
          fs    = K_CP * (mu_CP - xs(1:meas)) * xs(1:meas)
          B_2s  = D_CP2 * xs(1:meas)
      else
          print *, "Fatal error: model not recognized in Evaluate_log_likelihood_time_series_Gaussian_approx"
          stop
      end if

      lg_prob = 0.0d0
      do ii = 1, meas
          xip1 = xs(ii+1)
          mean = xi + fs(ii) * dts(ii)
          var  = B_2s(ii) * dts(ii)

          if (var <= 0.0d0) then
              print *, "Fatal error in Subroutine: valuate_log_likelihood_time_series_Gaussian_approx"
              print *, "Error: non-positive variance at step", ii
              print *, "  B_2s(ii) =", B_2s(ii), ", dts(ii) =", dts(ii), ", var =", var
              stop
          end if

          dummy = -((xip1 - mean)**2) / (2.0d0 * var) - 0.5d0 * log(2.0d0 * pi * var)
          lg_prob = lg_prob + dummy
          xi = xip1
      end do

      lg_prob = lg_prob / (ts(len) - ts(1))
    end subroutine Evaluate_log_likelihood_time_series_Gaussian_approx

    subroutine Evaluate_propagator_jump_Gaussian_approx(x0, t0, xf, tf, prob, Model, params_model, n_params)

      implicit none
      !-----------------  inputs / outputs ------------------------
      double precision, intent(in)  :: x0, t0, xf, tf
      integer*8,       intent(in)  :: n_params
      double precision, dimension(n_params), intent(in) :: params_model
      character(len=*), intent(in) :: Model
      double precision, intent(out) :: prob
      !-----------------  locals  --------------------------------
      double precision :: dt, f, B_2, mean, var,lg_prob
      double precision :: mu, K, D, D2
      double precision, parameter :: pi = 3.141592653589793d0
      !------------------------------------------------------------

      dt = tf - t0
      if (dt <= 0.0d0) then
          print *, "Fatal error: tf must be > t0 in Evaluate_propagator_jump_Gaussian_approx"
          stop
      end if

      !-----------  drift  f(x0)  and  noise prefactor  B²(x0)  --
      select case ( trim(Model) )
      case ("OU")
          mu = params_model(1);  K = params_model(2);  D = params_model(3)
          D2 = D*D
          f  = K * (mu - x0)
          B_2 = D2
      case ("DE")
          mu = params_model(1);  K = params_model(2);  D = params_model(3)
          D2 = D*D
          f  = K * (mu - x0)
          B_2 = D2 * x0
      case ("EN")
          mu = params_model(1);  K = params_model(2);  D = params_model(3)
          D2 = D*D
          f  = K * (mu - x0) * x0
          B_2 = D2 * x0**2
      case ("GB")
          K = params_model(1);  D = params_model(2)
          D2 = D*D
          f  = K * x0
          B_2 = D2 * x0**2
      case ("CP")
          mu = params_model(1);  K = params_model(2);  D = params_model(3)
          D2 = D*D
          f  = K * (mu - x0) * x0
          B_2 = D2 * x0
      case default
          print *, "Fatal error: model not recognised in Evaluate_propagator_jump_Gaussian_approx"
          stop
      end select

      !-----------  Gaussian propagator  --------------------------
      mean = x0 + f * dt
      var  = B_2 * dt

      if (var <= 0.0d0) then
          print *, "Fatal error: non-positive variance in Evaluate_propagator_jump_Gaussian_approx"
          stop
      end if

      lg_prob = - ( (xf - mean)**2 ) / (2.0d0*var)            &
                - 0.5d0 * log(2.0d0*pi*var)

      prob = exp(lg_prob) 

    end subroutine Evaluate_propagator_jump_Gaussian_approx

    subroutine Evaluate_likelihood_time_series_exact(xs,ts,len,Prob,Model,params_model,n_params)
    !Evaluate the likelihood of a time series using the exact propagator.
      implicit none
      integer*8, intent(in) :: len,n_params
      double precision,dimension(len), intent(in) :: xs,ts
      double precision,dimension(n_params) :: params_model
      double precision, intent(out) :: Prob
      character(len=*), intent(in) :: Model
      double precision,dimension(len-1) :: fs,B_2s
      double precision :: dt,sqrt_dt,xi,xip1,D2,mean,var,dummy
      double precision :: mu_OU,K_OU,D_OU,D_OU_2
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      integer*8 :: ii,meas

      dt = ts(2)-ts(1) !Assuming constant time step
      sqrt_dt = sqrt(dt)
      xi = xs(1)
      meas = len-1
      if ( trim(Model) == "OU" ) then
        mu_OU = params_model(1)
        K_OU = params_model(2)
        D_OU = params_model(3)
        D_OU_2 = D_OU*D_OU
        prob = 1.0d0
        do ii = 1, meas, 1
          xip1 = xs(ii+1)
          call propagator_OU(mu_OU,K_OU,D_OU_2,xi,xip1,dt,dummy)
          prob = Prob*dummy
          xi = xip1
        end do
      elseif ( trim(Model) == "DE" ) then
        mu_DE = params_model(1)
        K_DE = params_model(2)
        D_DE = params_model(3)
        D_DE_2 = D_DE*D_DE
        prob = 1.0d0
        do ii = 1, meas, 1
          xip1 = xs(ii+1)
          call propagator_dem(mu_DE,K_DE,D_DE_2,xi,xip1,dt,N_series=int8(10000),prop=dummy)
          prob = Prob*dummy
          xi = xip1
        end do
      else
        print *, "Fatal error, model not recognized in Evaluate_likelihood_time_series_Gaussian_approx"
        stop
      end if

    end subroutine Evaluate_likelihood_time_series_exact

    subroutine Evaluate_log_likelihood_time_series_exact(xs,ts,len,lg_prob,Model,params_model,n_params)
    !Evaluate the normalized log-likelihood of a time series using the exact propagator.
      implicit none
      integer*8, intent(in) :: len,n_params
      double precision,dimension(len), intent(in) :: xs,ts
      double precision,dimension(n_params) :: params_model
      double precision, intent(out) :: lg_prob
      character(len=*), intent(in) :: Model
      double precision,dimension(len-1) :: fs,B_2s
      double precision :: dt,sqrt_dt,xi,xip1,D2,mean,var,dummy
      double precision :: mu_OU,K_OU,D_OU,D_OU_2
      double precision :: mu_DE,K_DE,D_DE,D_DE_2
      integer*8 :: ii,meas

      dt = ts(2)-ts(1) !Assuming constant time step
      sqrt_dt = sqrt(dt)
      xi = xs(1)
      meas = len-1
      lg_prob = 0.0d0
      if ( trim(Model) == "OU" ) then
        mu_OU = params_model(1)
        K_OU = params_model(2)
        D_OU = params_model(3)
        D_OU_2 = D_OU*D_OU
        
        do ii = 1, meas, 1
          xip1 = xs(ii+1)
          call propagator_OU(mu_OU,K_OU,D_OU_2,xi,xip1,dt,dummy)
          lg_prob = lg_prob + log(dummy)
          xi = xip1
        end do
        
      elseif ( trim(Model) == "DE" ) then
        mu_DE = params_model(1)
        K_DE = params_model(2)
        D_DE = params_model(3)
        D_DE_2 = D_DE*D_DE
        
        do ii = 1, meas, 1
          xip1 = xs(ii+1)
          call propagator_dem(mu_DE,K_DE,D_DE_2,xi,xip1,dt,N_series=int8(20000),prop=dummy)
          lg_prob = lg_prob + log(dummy)
          xi = xip1
        end do
        
      else
        print *, "Fatal error, model not recognized in Evaluate_likelihood_time_series_Gaussian_approx"
        stop
      end if
      lg_prob = lg_prob/(ts(len)-ts(1))

    end subroutine Evaluate_log_likelihood_time_series_exact

    subroutine Evaluate_log_likelihood_time_series_OU_exact(xs,ts,len,lg_LK,params_model,n_params)
    !Evaluate the likelihood of a time series using the Gaussian approximation.
      implicit none
      integer*8, intent(in) :: len,n_params
      double precision,dimension(len), intent(in) :: xs,ts
      double precision,dimension(n_params) :: params_model
      double precision, intent(out) :: lg_LK
      double precision :: dt,sqrt_dt,xi,xip1,D2,tfmt0,dummy
      double precision :: mu_OU,K_OU,D_OU,D_OU_2
      integer*8 :: ii

      dt = ts(2)-ts(1) !Assuming constant time step
      sqrt_dt = sqrt(dt)
      xi = xs(1)
      tfmt0 = ts(len)-ts(1)
      mu_OU = params_model(1)
      K_OU = params_model(2)
      D_OU = params_model(3)
      D_OU_2 = D_OU*D_OU
      
      lg_LK = 0.0d0
      do ii = 1, len-1, 1
        xip1 = xs(ii+1)
        call propagator_OU(mu_OU,K_OU,D_OU_2,xi,xip1,dt,dummy)
        lg_LK = lg_LK + log(dummy)
        xi = xip1
      end do

      lg_LK = lg_LK/tfmt0

    end subroutine Evaluate_log_likelihood_time_series_OU_exact

    subroutine propagator_dem(mu_dem,K_dem,D2_dem,x0,xt,dt,N_series,prop)
    !Evaluate the propagator of the DEM model using truncated series for the Bessel function and Stirling's approximation for the gamma function.
      implicit none
      double precision, intent(in) :: mu_dem,K_dem,D2_dem,x0,xt,dt
      double precision, intent(out) ::prop
      integer*8, intent(in) :: N_series
      double precision :: a,b,Ebdt,l,w,v,u,z,lg_uv,lg_a_max,mu_old,K_old
      integer*8 :: ii
      double precision, dimension(0:N_series) :: lg_as
      mu_old = mu_dem*K_dem
      K_old = 1.0d0/mu_dem
      a = mu_old
      b = -mu_old*K_old
      Ebdt = exp(b*dt)
      l = 2.0d0*b/(Ebdt-1.0d0)/D2_dem
      w = 2.0d0*a/D2_dem-1.0d0
      v = l*xt
      u = l*x0*Ebdt
      z = 2.0*sqrt(u*v)
      
      lg_as(0) = +log(l) -u-v+w*log(v) -approx_log_gamma(w+1.0d0)
      lg_uv = log(u*v)
      do ii = 1, N_series, 1
          lg_as(ii) = lg_as(ii-1) + lg_uv - log(dble(ii)+w) - log(dble(ii))
      end do
      ! lg_a_max = maxval(lg_as)
      ! lg_prop = lg_a_max+log(sum(exp(lg_as-lg_a_max)))
      prop = sum(exp(lg_as))
      
    end subroutine propagator_dem

    subroutine propagator_OU(mu_OU,K_OU,D_OU_2,x0,xf,dt,prop)
    !Evaluate the propagator of the OU model using exact expression.
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,D_OU_2,x0,xf,dt
      double precision :: var,mean,dumt
      double precision, intent(out) :: prop
      dumt = K_OU*dt
      mean = mu_OU+(x0-mu_OU)*exp(-dumt)
      var = D_OU_2*(1.0d0-exp(-2.0d0*dumt))/K_OU/2.0d0 
      prop = Gauss(mean,var,xf)
      
    end subroutine propagator_OU

    subroutine propagator_GB(k_GB,D_GB_2,x0,xf,dt,prop)
    !Evaluate the propagator of the GB (geometric Brownian) model using exact expression.
    ! GB: dX = k*X*dt + D*X*dW
      implicit none
      double precision, intent(in) :: k_GB,D_GB_2,x0,xf,dt
      double precision :: dummy_2,dummy_1
      double precision, intent(out) :: prop

      if (x0 <= 0.d0 .or. xf <= 0.d0 .or. dt <= 0.d0) then
        print *, "Fatal error in subroutines: propagator_GB"
        print *, 'Error: x0, xf, and dt must be positive.'
        print *, '  x0 =', x0, ', xf =', xf, ', dt =', dt
        stop
      end if

      dummy_1 = (k_GB-0.5d0*D_GB_2)*dt 
      dummy_2 = 2.0d0*D_GB_2*dt 
      prop = exp(-(log(xf/x0)-dummy_1)**2.0d0/(dummy_2))/sqrt(dummy_2*pi)/xf
      
    end subroutine propagator_GB

    subroutine propagator_Unbiased_Brownian(D_B_2,x0,xf,dt,prop)
    !Evaluate the propagator of the OU model using exact expression.
      implicit none
      double precision, intent(in) :: D_B_2,x0,xf,dt
      double precision :: var,mean
      double precision, intent(out) :: prop
      mean = x0
      var = D_B_2*dt 
      prop = Gauss(mean,var,xf)
      
    end subroutine propagator_Unbiased_Brownian


    subroutine estimate_parameters_OU_Girsanov_correct_D(E_D_OU,E_D_OU2,E_mu_OU,E_K_OU,xs,meas,dt,T,m_x,m_x2,xt,x0,n_rows)
    !Estimate parameters of the OU model using Girsanov's theorem.
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_mu_OU,E_K_OU
      double precision, intent(in) ::E_D_OU,E_D_OU2
      double precision, intent(in) :: dt,T,m_x,m_x2,xt,x0
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs
      double precision, dimension(n_rows) :: ys
      double precision :: w_a,w_b,Z,m_y,m_y2,a,b

      w_a = (xt-x0)/E_D_OU/T
      w_b = (xt**2.0-x0**2.0)/2.0d0/E_D_OU2/T-0.5d0

      Z = (m_x2-m_x**2.0)/E_D_OU2
      m_y = m_x/E_D_OU
      m_y2 = m_x2/E_D_OU2

      a = (m_y2*w_a-m_y*w_b)/Z
      b = (m_y*w_a-w_b)/Z

      E_mu_OU = E_D_OU*a/b
      E_K_OU = b
    end subroutine estimate_parameters_OU_Girsanov_correct_D

    subroutine draw_trajectory_log_RN_lamperti_wiener_unconstrained(a,b,len,ts,ys,ls,target_model)
    !Draw a trajectory of the log-RN process of target_model lamperti transformed measure with respect to the Wiener process
      implicit none
      double precision, intent(in) :: a,b
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: ts,ys
      double precision, dimension(len), intent(out) :: ls
      character(len=*), intent(in) :: target_model
      integer*8 :: ii
      
      ls(1) = 0.0d0 !Initial value of the log-RN process is zero
      do ii = 2, len, 1
        call Evaluate_log_RN_OU_W_Girsanov_new(a,b,ii,ts(1:ii),ys(1:ii),ls(ii))
      end do

    end subroutine draw_trajectory_log_RN_lamperti_wiener_unconstrained

    subroutine draw_trajectory_log_RN_lampertis_unconstrained(a_target,b_target,a_ref,b_ref,len,ts,ys,ls,target_model,ref_model)
    !Draw a trajectory of the log-RN process of target_model lamperti transformed measure with respect to another (reference) lamperti transformed process.
    !Rn is computed using chain rule of Radon-Nikodym derivatives.
      implicit none
      double precision, intent(in) :: a_target,b_target,a_ref,b_ref
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: ts,ys
      double precision, dimension(len), intent(out) :: ls
      double precision, dimension(len) :: ls_ref,ls_target
      character(len=*), intent(in) :: target_model,ref_model
      integer*8 :: ii
      
      call draw_trajectory_log_RN_lamperti_wiener_unconstrained(a_target,b_target,len,ts,ys,ls_target,target_model)
      call draw_trajectory_log_RN_lamperti_wiener_unconstrained(a_ref,b_ref,len,ts,ys,ls_ref,ref_model)

      ls = ls_target - ls_ref !Chain rule of Radon-Nikodym derivatives


    end subroutine draw_trajectory_log_RN_lampertis_unconstrained

    subroutine draw_trajectory_log_RN_lamperti_constrained(a,b,len,ts,ys,ls,model)
    !Draw a trajectory of the log-RN process with respect to the lamperti-transformed process
      implicit none
      double precision, intent(in) :: a,b
      integer*8, intent(in) :: len
      double precision, dimension(len), intent(in) :: ts,ys
      double precision, dimension(len), intent(out) :: ls
      character(len=*), intent(in) :: model
      integer*8 :: ii
      double precision :: tfmt0,prop_OU,K_OU,mu_OU,prop_W
      
      ls(1) = 0.0d0 !Initial value of the log-RN process is zero
      if ( trim(model) == "OU" ) then
        mu_OU = a/b
        K_OU = b
        do ii = 2, len, 1
          tfmt0 = ts(ii)-ts(1)
          call propagator_OU(mu_OU,K_OU,1.0d0,ys(1),ys(ii),tfmt0,prop_OU)
          call Evaluate_log_RN_OU_W_Girsanov_new(a,b,ii,ts(1:ii),ys(1:ii),ls(ii)) 
          ls(ii) = ls(ii)- log(prop_OU)/tfmt0
        end do
      else if ( trim(model) == "WI" ) then
        do ii = 2, len, 1
          tfmt0 = ts(ii)-ts(1)
          call propagator_Unbiased_Brownian(1.0d0,ys(1),ys(ii),tfmt0,prop_W)
          ls(ii) = -log(prop_W)/tfmt0
        end do
      endif

    end subroutine draw_trajectory_log_RN_lamperti_constrained

    subroutine Evaluate_log_RN_OU_W_Girsanov_new(a,b,n_rows,ts,ys,l_OU)
    ! Update of Evaluate_log_RN_OU_W_Girsanov_new commented below 28/05/25
    ! Older version computed a, b inside the subroutine, this version receives them as arguments.
    !       ! a = K_OU*mu_OU/D_OU
    !       ! b = K_OU
    !Evaluate normalized log-RN derivative OU/Wiener process 
    !Input data is Lamperti transformed process yt (OU process with D=1)
      implicit none
      double precision, intent(in) :: a,b
      integer*8, intent(in) :: n_rows
      double precision, dimension(n_rows), intent(in) :: ts,ys
      double precision, intent(out) :: l_OU
      double precision :: w_a,w_b,tfmt0,m_y,m_y2,yf,y0
      integer*8 :: meas
      meas = n_rows-1
      m_y = sum(ys(1:meas))/dble(meas)
      m_y2 = sum(ys(1:meas)**2.0d0)/dble(meas)
      y0 = ys(1)
      yf = ys(n_rows)
      tfmt0 = ts(n_rows)-ts(1) 
      
      w_a = (yf-y0)/tfmt0
      w_b = (yf**2.0-y0**2.0)/2.0d0/tfmt0-0.5d0

      l_OU = w_a*a-w_b*b-(a**2.0d0)/2.0d0-(b**2.0d0)*m_y2/2.0d0+a*b*m_y
    end subroutine Evaluate_log_RN_OU_W_Girsanov_new

    subroutine Evaluate_log_RN_DE_W_Girsanov(mu_DE,K_DE,D_DE,n_rows,ts,ys,l_DE)
    !Evaluate normalized log-RN derivative DE/Wiener process 
    !Input data is Lamperti transformed process yt (DE process with D=1)
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: mu_DE,K_DE,D_DE
      integer*8, intent(in) :: n_rows
      double precision, dimension(n_rows), intent(in) :: ts,ys
      double precision, intent(out) :: l_DE
      integer*8 :: meas
      double precision :: a,b,m_inv_y2,m_y2,y0,yf,tfmt0,br_y,br_inv_y

      m_inv_y2 = sum(ys(1:meas)**(-2.0d0))/dble(meas)
      m_y2 = sum(ys(1:meas)**2.0d0)/dble(meas)
      y0 = ys(1)
      yf = ys(n_rows)
      tfmt0 = ts(n_rows)-ts(1) 
      a = 2.0d0*K_DE*mu_DE/D_DE**2.0d0-0.5d0
      b = K_DE/2.0d0
      br_y = (yf**2.0d0 - y0**2.0d0)/2.0d0/tfmt0 -0.5d0
      br_inv_y = log(yf/y0)/tfmt0 + 0.5d0*m_inv_y2


      l_DE = a*br_inv_y - b*br_y -(a**2.0d0)*m_inv_y2/2.0d0-(b**2.0d0)*m_y2/2.0d0+a*b
    end subroutine Evaluate_log_RN_DE_W_Girsanov

    subroutine Evaluate_ell_OU_EM(mu_OU,K_OU,xs,ts,len,ell_OU_EM)
    !Evaluate log-RN derivative in the OU process with D = 1 (Lamperti transformed process yt)
    !No exact evaluation of stochastic integrals, this is the same as evaluating the likelihood ratio process using the Euler-Maruyama approximation of the OU propagator
      implicit none
      double precision, intent(in) :: mu_OU,K_OU
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_EM
      double precision :: a,b,tfmt0,int_y,int_t
      double precision,dimension(len-1) :: dys,dts

      a = K_OU*mu_OU
      b = K_OU
      tfmt0 = ts(len)-ts(1)
  
      dts = ts(2:len)-ts(1:len-1)
      dys = xs(2:len)-xs(1:len-1)

      int_y = sum((a-b*xs(1:len-1))*dys)
      int_t = sum((a-b*xs(1:len-1))**2.0d0*dts)
      ell_OU_EM = (int_y-int_t/2.0d0)/tfmt0
    end subroutine Evaluate_ell_OU_EM

    subroutine Evaluate_ell_OU_OU_no_integrals(mu_OU,K_OU,mu_OU2,K_OU2,xs,ts,len,ell_OU_no_integrals)
    !Evaluate log-RN derivative in the OU/OU process with D = 1 (Lamperti transformed process yt)
    !No exact evaluation of stochastic integrals, this is exactly as evaluating the likelihood ratio process using the Euler-Maruyama approximation of the OU propagator
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,mu_OU2,K_OU2
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_no_integrals
      double precision :: tfmt0,int_y,int_t
      double precision,dimension(len-1) :: dys,dts,drift_1,drift_2

      tfmt0 = ts(len)-ts(1)
      dts = ts(2:len)-ts(1:len-1)
      dys = xs(2:len)-xs(1:len-1)
      drift_1 = -K_OU*(xs(1:len-1)-mu_OU) 
      drift_2 = -K_OU2*(xs(1:len-1)-mu_OU2) 

      int_y = sum((drift_1-drift_2)*dys)
      int_t = sum((drift_1**2.0d0-drift_2**2.0d0)*dts)
      ell_OU_no_integrals = (int_y-int_t/2.0d0)/tfmt0
    end subroutine Evaluate_ell_OU_OU_no_integrals

    subroutine Evaluate_ell_OU_OU_Girsanov_same_K(mu_OU,K_OU,mu_OU2,xs,ts,len,ell_OU_OU_Gir)
    !Evaluate log-RN derivative in the OU process with D = 1 (Lamperti transformed process yt) with respected to another OU process with different mean (same autocorrelation time parametrized here with K_OU).
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,mu_OU2
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_OU_Gir
      double precision :: Dmu,s1,s2,tfmt0,int_t,int_y
      double precision,dimension(len-1) :: dts

      Dmu = mu_OU-mu_OU2
      tfmt0 = ts(len)-ts(1)
  
      dts = ts(2:len)-ts(1:len-1)

      int_y = K_OU*Dmu*(xs(len)-xs(1))
      s1 = (mu_OU**2.0d0-mu_OU2**2.0d0)*tfmt0
      s2 = -2.0*Dmu*sum(xs(1:len-1)*dts)
      int_t = K_OU*K_OU*(s1+s2)
      ! print *, "  int_y=",int_y,"int_t=",int_t
      ell_OU_OU_Gir = (int_y-int_t/2.0d0)/tfmt0
    end subroutine Evaluate_ell_OU_OU_Girsanov_same_K


    subroutine Evaluate_ell_OU_OU_exact(mu_OU,K_OU,mu_OU2,K_OU2,xs,ts,len,ell_OU_exact)
    !Evaluate log-RN derivative in the OU/OU process with D = 1 (Lamperti transformed process yt)
    !Using exact OU propagator
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,mu_OU2,K_OU2
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_exact
      double precision :: RN,prop_OU,prop_OU2,tfmt0,xf,x0,dt
      integer*8 :: ii

      RN = 1.0d0
      dt = ts(2)-ts(1)
      tfmt0 = ts(len)-ts(1)
      do ii = 1, len-1, 1
        x0 = xs(ii)
        xf = xs(ii+1)
        call propagator_OU(mu_OU,K_OU,1.0d0,x0,xf,dt,prop_OU)
        call propagator_OU(mu_OU2,K_OU2,1.0d0,x0,xf,dt,prop_OU2)
        RN = RN*prop_OU/prop_OU2
      end do
      ell_OU_exact = log(RN)/tfmt0
    end subroutine Evaluate_ell_OU_OU_exact

    subroutine Evaluate_ell_OU_OU_EM(mu_OU,K_OU,mu_OU2,K_OU2,xs,ts,len,ell_OU_EM)
    !Evaluate log-RN derivative in the OU/OU process with D = 1 (Lamperti transformed process yt)
    !Using EM approximation of OU propagator
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,mu_OU2,K_OU2
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_EM
      double precision :: RN,prop_OU,prop_OU2,tfmt0,xf,x0,dt,mean,var
      integer*8 :: ii

      RN = 1.0d0
      dt = ts(2)-ts(1)
      tfmt0 = ts(len)-ts(1)
      var = dt
      do ii = 1, len-1, 1
        x0 = xs(ii)
        xf = xs(ii+1)
        mean = x0 -K_OU*(x0-mu_OU)*dt
        prop_OU = Gauss(mean,var,xf)
        mean = x0 -K_OU2*(x0-mu_OU2)*dt
        prop_OU2 = Gauss(mean,var,xf)
        RN = RN*prop_OU/prop_OU2
      end do
      ell_OU_EM = log(RN)/tfmt0
    end subroutine Evaluate_ell_OU_OU_EM

    subroutine Evaluate_ell_Wiener_bridge(D_B,tf,xt,t0,x0,ell_B)
    !Evaluate log-RN derivative between the Wiener process and its bridge. This evaluation doesn't depend on the path.
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: D_B
      double precision, intent(in) :: tf,xt,t0,x0
      double precision, intent(out) :: ell_B
      double precision :: y0,yt

      y0 = x0/D_B
      yt = xt/D_B
      ell_B = (-(yt-y0)**2.0d0/(tf-t0)+log(2*pi*(tf-t0)))/2.0d0/(tf-t0)

    end subroutine Evaluate_ell_Wiener_bridge

    subroutine Estimate_integrals_with_expanded_data(mu_OU,K_OU,D_OU,xs,ts,m_x_refined,m_x2_refined,n_rows,sum_xs,sum_x2s,dt_bridge,dt_integration_bridge,len_bridge)
    !Use OU bridges to refine data and compute integrals appearing in log-likelyhood ratio using this refined paths
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,D_OU,sum_xs,sum_x2s,dt_bridge,dt_integration_bridge
      double precision, intent(out) :: m_x_refined,m_x2_refined
      integer*8, intent(in) :: n_rows, len_bridge
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision :: x0_bridge,xf_bridge,t0_bridge,tf_bridge
      double precision, allocatable :: xs_bridge(:),ts_bridge(:)
      integer*8 :: norm,ii

      m_x_refined = sum_xs
      m_x2_refined = sum_x2s
      norm = n_rows

      allocate(xs_bridge(len_bridge),ts_bridge(len_bridge))
      do ii =1,n_rows-1,1
        t0_bridge = ts(ii)
        tf_bridge = ts(ii +1)
        x0_bridge = xs(ii)
        xf_bridge = xs(ii +1)
        ! call Fill_gaps_with_OU_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,xs_bridge,len_bridge,dt_integration_bridge,mu_OU,K_OU,D_OU)
        call Fill_gaps_with_OU_bridge_multiresolution(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_bridge,D_OU,mu_OU,K_OU,xs_bridge,ts_bridge)
        norm = norm + len_bridge
        m_x_refined = m_x_refined + x0_bridge + xf_bridge + sum(xs_bridge)
        m_x2_refined = m_x2_refined + x0_bridge**2.0 + xf_bridge**2.0 + sum(xs_bridge**2.0)
      enddo
      m_x_refined  = m_x_refined /norm
      m_x_refined  = m_x2_refined /norm
      deallocate(xs_bridge)
    end subroutine Estimate_integrals_with_expanded_data

    subroutine estimate_parameters_OU_gradient_descent_expanded_data(E_D_0,E_mu_0,E_K_0,E_mu_OU,E_K_OU,T,xs,ts,xt,x0,n_rows)
    ! Given an initial guess for mu and K looks for refinement in the estimation using gradient descent. (E_D_0,E_mu_0,E_K_0) --> (E_mu,E_K)
    ! Synthetic data is aggregated using OU bridges
    !Some input variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: E_D_0,E_mu_0,E_K_0
      double precision, intent(out) :: E_mu_OU,E_K_OU
      double precision, intent(in) :: T,xt,x0
      integer*8, intent(in) :: n_rows
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision :: w_a,w_b,Z,m_y,m_y2,a,b,h_a,h_b,E_D_02,D_l_a,D_l_b,errs,an,bn,tolerance,sum_xs,sum_x2s,m_x2,m_x,dt_bridge,dt_integration_bridge
      integer*8 :: iterations, ii,len_bridge 

      sum_xs = sum(xs)
      sum_x2s = sum(xs**2.0)
      
      E_D_02 = E_D_0*E_D_0
      

      h_a = 1.0D-3
      h_b = 1.0D-4
      w_a = (xt-x0)/E_D_0/T
      w_b = (xt**2.0-x0**2.0)/2.0d0/E_D_02/T-0.5d0

      dt_bridge = 0.01d0 ! move out as imput parameter
      dt_integration_bridge = 1.0d-4 ! move out as imput parameter
      len_bridge = int8((ts(2)-ts(1))/dt_bridge)-1
      
      
      errs = 1.0d0
      iterations = 1.0D7
      tolerance = 1.0D-3
      ii = 0
      E_mu_OU = E_mu_0
      E_K_OU = E_K_0
      a = E_mu_0*E_K_0/E_D_0
      b = E_K_0
      do ii =  1, iterations, 1
        ! print *, "  ii/iteration*100=",dble(ii)/iterations*100
        call Estimate_integrals_with_expanded_data(E_mu_0,E_K_0,E_D_0,xs,ts,m_x,m_x2,n_rows,sum_xs,sum_x2s,dt_bridge,dt_integration_bridge,len_bridge)
        m_y = m_x/E_D_0
        m_y2 = m_x2/E_D_02
        Z = (m_x2-m_x**2.0)/E_D_02
        D_l_a = w_a-a+m_y*b
        D_l_b = -w_b+m_y*a-m_y2*b

        a = a + h_a*D_l_a
        b = b + h_b*D_l_b
        
        E_mu_OU = E_D_0*a/b
        E_K_OU = b
      end do
      
    end subroutine estimate_parameters_OU_gradient_descent_expanded_data

    subroutine estimate_parameters_OU_gradient_descent(E_D_0,E_mu_0,E_K_0,E_mu_OU,E_K_OU,T,m_x,m_x2,xt,x0)
    ! Given an initial guess for mu and K looks for refinement in the estimation using gradient descent. (E_D_0,E_mu_0,E_K_0) --> (E_mu,E_K)
    !Some input variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(in) :: E_D_0,E_mu_0,E_K_0
      double precision, intent(out) :: E_mu_OU,E_K_OU
      double precision, intent(in) :: T,m_x,m_x2,xt,x0
      double precision :: w_a,w_b,Z,m_y,m_y2,a,b,h_a,h_b,E_D_02,D_l_a,D_l_b,errs,an,bn,tolerance
      integer*8 :: iterations, ii
      h_a = 1.0D-3
      h_b = 1.0D-4
      E_D_02 = E_D_0*E_D_0
      a = E_mu_0*E_K_0/E_D_0
      b = E_K_0

      w_a = (xt-x0)/E_D_0/T
      w_b = (xt**2.0-x0**2.0)/2.0d0/E_D_02/T-0.5d0

      Z = (m_x2-m_x**2.0)/E_D_02
      m_y = m_x/E_D_0
      m_y2 = m_x2/E_D_02
      
      errs = 1.0d0
      iterations = 1.0D7
      tolerance = 1.0D-3
      ii = 0
      do ii =  1, iterations, 1
      ! do while ((ii<=iterations).and.(errs>=tolerance))
        ! ii = ii + 1 
        D_l_a = w_a-a+m_y*b
        D_l_b = -w_b+m_y*a-m_y2*b
        ! print *
        ! print *, "D_l_a=",D_l_a,"D_l_b=",D_l_b
        an = a + h_a*D_l_a
        bn = b + h_b*D_l_b
        errs = sqrt((an-a)**2.0+(bn-b)**2.0)
        a = an
        b = bn
      end do
      

      E_mu_OU = E_D_0*a/b
      E_K_OU = b
    end subroutine estimate_parameters_OU_gradient_descent

    subroutine Evaluate_ell_OU_OU_Girsanov(mu_OU,K_OU,mu_OU2,K_OU2,xs,ts,len,ell_OU_OU_Gir)
    !Evaluate log-RN derivative in the OU process with D = 1 (Lamperti transformed process yt) with respected to another OU process with different mean and auto-correlation time
      implicit none
      double precision, intent(in) :: mu_OU,K_OU,mu_OU2,K_OU2
      integer*8, intent(in) :: len
      double precision,dimension(len), intent(in) :: xs,ts
      double precision, intent(out) :: ell_OU_OU_Gir
      double precision :: int_t_1,int_t_2,int_y_1,int_y_2,xf,x0,dummy,tfmt0
      double precision,dimension(len-1) :: dts,ys

      tfmt0 = ts(len)-ts(1)  
      dts = ts(2:len)-ts(1:len-1)
      ys = xs(1:len-1)
      xf = xs(len)
      x0 = xs(1)
      dummy = tfmt0 + x0**2.0d0 - xf**2.0d0
      int_y_1 = K_OU* (dummy +2.0d0*mu_OU *(xf-x0))/2.0d0
      int_y_2 = K_OU2*(dummy +2.0d0*mu_OU2*(xf-x0))/2.0d0
      int_t_1 = K_OU*K_OU  *sum(((mu_OU -ys)**2.0d0)*dts)
      int_t_2 = K_OU2*K_OU2*sum(((mu_OU2-ys)**2.0d0)*dts)
      ! print *, "  int_y=",int_y,"int_t=",int_t
      ell_OU_OU_Gir = (int_y_1-int_y_2-int_t_1/2.0d0+int_t_2/2.0d0)/tfmt0
    end subroutine Evaluate_ell_OU_OU_Girsanov

    subroutine Compute_quantities_data(T,m_inv_x,m_x,m_x2,x0,xt,xs,ts,n_rows)
    !Compute quantities used to evaluate log-likelihood from data
      implicit none
      double precision, intent(out) :: T,m_inv_x,m_x,m_x2,x0,xt
      integer*8, intent(in) :: n_rows
      double precision, dimension(n_rows), intent(in) :: xs,ts
      integer*8 :: n_jumps

      T = ts(n_rows)
      n_jumps = n_rows - 1
      m_inv_x = sum(1.0d0/xs)/n_rows
      ! print *, "  m_inv_x= ",m_inv_x,sum(xs**(-1.0d0))/n_jumps
      m_x = sum(xs(1:n_jumps))/n_jumps
      m_x2 = sum((xs(1:n_jumps))**2.0)/n_jumps
      x0 = xs(1)
      xt = xs(n_rows)
    end subroutine Compute_quantities_data

    subroutine estimate_parameters_OU_Girsanov(E_D_OU,E_D_OU2,E_mu_OU,E_K_OU,xs,meas,dt,T,m_x,m_x2,xt,x0,n_rows)
    !Estimate parameters of the OU model using Girsanov's theorem.
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_mu_OU,E_K_OU,E_D_OU,E_D_OU2
      double precision, intent(in) :: dt,T,m_x,m_x2,xt,x0
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs
      double precision, dimension(n_rows) :: ys
      double precision :: w_a,w_b,Z,m_y,m_y2,a,b

      E_D_OU = sqrt(sum((xs(2:) - xs(:meas))**2)/dt/meas)
      E_D_OU2 = E_D_OU*E_D_OU

      w_a = (xt-x0)/E_D_OU/T
      w_b = (xt**2.0-x0**2.0)/2.0d0/E_D_OU2/T-0.5d0

      Z = (m_x2-m_x**2.0)/E_D_OU2
      m_y = m_x/E_D_OU
      m_y2 = m_x2/E_D_OU2

      a = (m_y2*w_a-m_y*w_b)/Z
      b = (m_y*w_a-w_b)/Z

      E_mu_OU = E_D_OU*a/b
      E_K_OU = b
    end subroutine estimate_parameters_OU_Girsanov

    subroutine estimate_parameters_DEM_Girsanov(E_D_dem,E_D_dem2,E_mu_dem,E_K_dem,xs,meas,dt,tf,m_x,m_inv_x,xt,x0,n_rows)
    !Estimate parameters of the DEM model using Girsanov's theorem
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_mu_dem,E_K_dem,E_D_dem,E_D_dem2
      double precision, intent(in) :: dt,tf,m_x,m_inv_x,xt,x0
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs
      double precision, dimension(n_rows) :: ys
      double precision :: m_y2,m_inv_y2,br_y,br_inv_y,Z,dum_mu,dum_K
      double precision :: a,b

      ys = 2.0d0*sqrt(xs)
      E_D_dem = sqrt(sum((ys(2:) - ys(:meas))**2)/dt/meas)
      E_D_dem2 = E_D_dem*E_D_dem
      
      m_y2 = 4.0d0*m_x/E_D_dem2
      m_inv_y2 = E_D_dem2*m_inv_x/4.0d0
      br_y = (4*(xt-x0)/E_D_dem2-tf)/2.0d0/tf
      br_inv_y = 0.5d0*(log(xt/x0)/tf+E_D_dem2*m_inv_x/4.0d0) 

      Z = 1.0d0 - m_y2*m_inv_y2
      a = (br_y - m_y2*br_inv_y)/Z
      b = (br_y*m_inv_y2-br_inv_y)/Z
      dum_mu = E_D_dem2/4.0d0*(a*2.0d0+1.0d0)
      dum_K = 2.0d0*b/dum_mu
      E_K_dem = dum_mu*dum_K
      E_mu_dem = 1.0d0/dum_K
      ! print *, "a= ",a,"b=",b
      ! print *, "dum_mu= ",dum_mu,"dum_k2=",dum_K
    end subroutine estimate_parameters_DEM_Girsanov

    subroutine estimate_parameters_DEM_Girsanov_correct_D(E_D_dem,E_D_dem2,E_mu_dem,E_K_dem,xs,meas,dt,tf,m_x,m_inv_x,xt,x0,n_rows)
    !Estimate parameters of the DEM model using Girsanov's theorem
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_mu_dem,E_K_dem
      double precision, intent(in) :: dt,tf,m_x,m_inv_x,xt,x0,E_D_dem,E_D_dem2
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs
      double precision, dimension(n_rows) :: ys
      double precision :: m_y2,m_inv_y2,br_y,br_inv_y,Z,dum_mu,dum_K
      double precision :: a,b

      ys = 2.0d0*sqrt(xs)
      
      m_y2 = 4.0d0*m_x/E_D_dem2
      m_inv_y2 = E_D_dem2*m_inv_x/4.0d0
      br_y = (4*(xt-x0)/E_D_dem2-tf)/2.0d0/tf
      br_inv_y = 0.5d0*(log(xt/x0)/tf+E_D_dem2*m_inv_x/4.0d0) 

      Z = 1.0d0 - m_y2*m_inv_y2
      a = (br_y - m_y2*br_inv_y)/Z
      b = (br_y*m_inv_y2-br_inv_y)/Z
      dum_mu = E_D_dem2/4.0d0*(a*2.0d0+1.0d0)
      dum_K = 2.0d0*b/dum_mu
      E_K_dem = dum_mu*dum_K
      E_mu_dem = 1.0d0/dum_K
      ! print *, "a= ",a,"b=",b
      ! print *, "dum_mu= ",dum_mu,"dum_k2=",dum_K
    end subroutine estimate_parameters_DEM_Girsanov_correct_D


    subroutine estimate_parameters_ENV_Girsanov(E_D_env,E_D_env2,E_mu_env,E_K_env,xs,ts,n_rows)
    !Estimate parameters of the ENV model using Girsanov's theorem
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_mu_env,E_K_env,E_D_env,E_D_env2
      integer*8, intent(in) :: n_rows
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision, dimension(n_rows) :: ys
      double precision, dimension(n_rows-1) ::dts
      double precision :: a,b,Z,bry,m_x,m_x2,xf,x0,tf,t0,tfmt0
      integer*8 :: meas

      meas = n_rows - 1
      dts = ts(2:n_rows)-ts(1:n_rows-1)
      ys = log(xs)
      E_D_env2 = sum(((ys(2:n_rows) - ys(1:meas))**2) / dts) / meas
      E_D_env = sqrt(E_D_env2)
      m_x = sum(xs)/n_rows
      m_x2 = sum(xs**2.0)/n_rows
      xf = xs(n_rows)
      x0 = xs(1)
      tf = ts(n_rows)
      t0 = ts(1)
      tfmt0 = tf-t0
      Z = m_x2-m_x**2.0d0
      bry = (xf-x0)/tfmt0/E_D_env-m_x*E_D_env/2.0d0

      b = (log(xf/x0)/E_D_env/tfmt0*m_x-bry)/Z
      a = (log(xf/x0)/E_D_env/tfmt0*m_x2-bry*m_x)/Z

      E_mu_env = (a+E_D_env/2.0d0)/b
      E_K_env = E_D_env*b

    end subroutine estimate_parameters_ENV_Girsanov

    subroutine estimate_parameters_CP_Girsanov(E_D_CP,E_D_CP2,E_mu_CP,E_K_CP,xs,ts,n_rows,meas)
    !Estimate parameters of the CP model using Girsanov's theorem
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_D_CP,E_D_CP2,E_mu_CP,E_K_CP
      integer*8, intent(in) :: n_rows,meas
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision, dimension(n_rows) :: ys
      double precision, dimension(meas) :: yss,dys
      double precision :: a,b,Z,m_y_2,m_y_4,m_y_6,yf,y0,tfmt0,dtfmt0,br_y,br_y_3,br_inv_y,c1,c2,dt

      dt = ts(2)-ts(1)
      ys = 2.0d0*sqrt(xs)
      E_D_CP = sqrt(sum((ys(2:) - ys(:meas))**2)/dt/meas)
      E_D_CP2 = E_D_CP*E_D_CP
      ys = ys /E_D_CP
      yss = ys(:meas)
      dys = ys(2:)-ys(:meas)
      tfmt0 = ts(n_rows)-ts(1)
      dtfmt0 = 2.0d0*tfmt0
      yf = ys(n_rows)
      y0 = ys(1)
      m_y_2 = sum(yss**2.0d0)/meas
      m_y_4 = sum(yss**4.0d0)/meas
      m_y_6 = sum(yss**6.0d0)/meas
      br_y = sum(dys*yss)/tfmt0!(yf**2.0d0-y0**2.0d0)/dtfmt0 - 1.0d0/dtfmt0 
      br_y_3 = sum(dys*yss**(3.0d0))/tfmt0!(yf**4.0d0-y0**4.0d0)/4.0d0/tfmt0 - 3.0d0/2.0d0*m_y_2
      br_inv_y = sum(dys/yss)/tfmt0
      Z = 2.0d0*(m_y_4**2.0d0 - m_y_2*m_y_6)
      if (abs(Z) < 1.0d-12) then
          print *, "Warning: Denominator Z is close to zero!"
          return
      end if

      c1 = 2.0d0*br_y_3 - m_y_2
      c2 = 1+2.0d0*br_y
      b = (c1*m_y_4-c2*m_y_6)/Z
      a = (c1*m_y_2-c2*m_y_4)/Z

      E_mu_CP = (a*E_D_CP2/2.0d0)/4.0d0/b
      E_K_CP = 8.0d0*b/E_D_CP2

    end subroutine estimate_parameters_CP_Girsanov

    subroutine estimate_parameters_CP_Girsanov_2(E_D_CP,E_D_CP2,E_mu_CP,E_K_CP,xs,ts,n_rows,meas)
    !Estimate parameters of the CP model using Girsanov's theorem
    !Some inputs variables require calling Compute_quantities_data before.
      implicit none
      double precision, intent(out) :: E_D_CP,E_D_CP2,E_mu_CP,E_K_CP
      integer*8, intent(in) :: n_rows,meas
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision, dimension(n_rows) :: ys
      double precision :: m_x_2,m_x,v,v2,xf,x0,tfmt0,br_x,Z,m_x_3,dt

      dt = ts(2)-ts(1)
      tfmt0 = ts(n_rows)-ts(1)
      ys = 2.0d0*sqrt(xs)
      E_D_CP = sqrt(sum((ys(2:) - ys(:meas))**2)/dt/meas)
      E_D_CP2 = E_D_CP*E_D_CP
      xf = xs(n_rows)
      x0 = xs(1)
      m_x_3 = sum(xs**3.0d0)/n_rows
      m_x_2 = sum(xs**2.0d0)/n_rows
      m_x = sum(xs)/n_rows
      v = (xf-x0)/tfmt0
      v2 = v*v
      br_x = v2/2.0d0 - E_D_CP2*m_x/2.0d0
      Z = v2+m_x_3*m_x- m_x_2**2.0d0

      E_mu_CP = m_x_2/m_x+v
      E_K_CP = -( m_x*br_x - v2 + v*m_x_2)/Z

    end subroutine estimate_parameters_CP_Girsanov_2

    subroutine estimate_parameters_GGM_Girsanov(E_D_GGM,E_D_GGM2,E_mu_GGM,E_K_GGM,xs,meas,dt,tf,n_rows,theta)
    !Estimate parameters of the GGM (general gamma model) model using Girsanov's theorem
    !Theta cannot be equal to one!!
      implicit none
      double precision, intent(out) :: E_mu_GGM,E_K_GGM,E_D_GGM,E_D_GGM2
      double precision, intent(in) :: dt,theta,tf
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs
      double precision, dimension(n_rows) :: ys
      double precision, dimension(meas) :: dys,yss
      double precision :: omt, gamma , br_y_g, m_y_2g, m_inv_y2,m_y_gm1,br_inv_y,dum_mu,dum_k,a,b,Z,W,yf,y0,gp1
      ! print*,"theta=",theta
      omt = (1.0d0-theta)/2.0d0
      ys = xs**(omt)/omt
      E_D_GGM = sqrt(sum((ys(2:) - ys(:meas))**2)/dt/meas)
      
      E_D_GGM2 = E_D_GGM**2.0d0
      ys = ys /E_D_GGM
      yf = ys(n_rows)
      y0 = ys(1)
      dys = ys(2:) - ys(:meas)
      yss = ys(:meas)
      gamma = (1.0d0+theta)/(1.0d0-theta)
      gp1 = gamma + 1

      m_inv_y2 = sum(yss**(-2.0d0))/meas
      m_y_gm1  = sum(yss**(gamma-1.0d0))/meas
      m_y_2g   = sum(yss**(2.0d0*gamma))/meas
      ! br_inv_y = sum(dys/yss)/tf ! Using definition
      br_inv_y = log(yf/y0)/tf+0.5d0*m_inv_y2
      ! br_y_g   = sum(dys * yss**(gamma))/tf ! Using definition
      br_y_g   = (yf**(gp1)-y0**(gp1))/gp1/tf-gamma/2.0d0*m_y_gm1
      Z = m_y_gm1**2.0d0 - m_inv_y2*m_y_2g

      a = (br_y_g*m_y_gm1-br_inv_y*m_y_2g)/Z
      b = (br_y_g*m_inv_y2-br_inv_y**m_y_gm1)/Z
      ! print *, "  Z=",Z,"br_y_g=",br_y_g,"m_y_gm1=",m_y_gm1
      ! print *, "  a=",a,"b=",b
      W = E_D_GGM*(1.0d0-theta)/2.0d0
      dum_mu = E_D_GGM*W*(a+0.5d0*gamma)
      dum_K = E_D_GGM*b/dum_mu*W**(-gamma)

      E_mu_GGM = 1.0d0/dum_K
      E_K_GGM = dum_mu*dum_K

    end subroutine estimate_parameters_GGM_Girsanov

    subroutine estimate_parameters_GGM_Girsanov_2(E_D_GGM,E_D_GGM2,E_mu_GGM,E_K_GGM,xs,ts,meas,n_rows,theta)
    !Estimate parameters of the GGM (general gamma model) model using Girsanov's theorem
    !Theta cannot be equal to one!!
      implicit none
      double precision, intent(out) :: E_mu_GGM,E_K_GGM,E_D_GGM,E_D_GGM2
      double precision, intent(in) :: theta
      integer*8, intent(in) :: meas,n_rows
      double precision, dimension(n_rows), intent(in) :: xs,ts
      double precision, dimension(n_rows) :: ys
      double precision, dimension(meas) :: xss
      double precision :: omt, m_x_t, m_x_tp1, m_x_tm1,v,br_inv_x,dt,tfmt0
      ! print*,"theta=",theta
      dt = ts(2)-ts(1)
      omt = (1.0d0-theta)/2.0d0
      ys = xs**(omt)/omt
      E_D_GGM = sqrt(sum((ys(2:) - ys(:meas))**2)/dt/meas)
      E_D_GGM2 = E_D_GGM**2.0d0
      xss = xs(:meas)
      m_x_t = sum(xss**theta)/meas
      m_x_tp1 = sum(xss**(theta+1.0d0))/meas
      m_x_tm1 = sum(xss**(theta-1.0d0))/meas
      tfmt0 = ts(n_rows)-ts(1)
      v = (xs(n_rows)-xs(1))/tfmt0
      br_inv_x = log(xs(n_rows)/xs(1))/tfmt0 + 0.5d0*E_D_GGM2*m_x_tm1

      E_K_GGM = (v*m_x_tm1-m_x_t*br_inv_x)/(m_x_t**2.0d0 - m_x_tm1*m_x_tp1)
      E_mu_GGM = (v*m_x_t-m_x_tp1*br_inv_x)/(v*m_x_tm1 - m_x_t*br_inv_x)
      

    end subroutine estimate_parameters_GGM_Girsanov_2

    subroutine compute_log_likelihood_OU(l_OU,E_mu_OU,E_K_OU,E_D_OU2,n_rows,T,xs,dt)
        !Compute log-likelihood of the data given the parameters of the OU model
        implicit none
        double precision, intent(out) :: l_OU
        double precision, intent(in) :: E_mu_OU,E_K_OU,E_D_OU2,T,dt
        integer*8, intent(in) :: n_rows
        double precision, dimension(n_rows), intent(in) :: xs
        double precision :: mean,var,x0,xt
        integer*4 :: ii

            l_OU = 0.0d0
            x0 = xs(1)
            var = E_K_OU/E_D_OU2/(1.0d0-exp(-2.0d0*E_K_OU*dt)) !This is actually twice the variance
            do ii = 1, n_rows-1, 1
                xt = xs(ii+1)
                mean = E_mu_OU+(x0-E_mu_OU)*exp(-E_K_OU*dt)
                l_OU = l_OU - (xt-mean)**2.0d0
                x0 = xt
            end do
            l_OU = l_OU*var
            l_OU = l_OU - (n_rows-1)/2.0d0*log(pi/var)
            l_OU = l_OU/T

    end subroutine compute_log_likelihood_OU

    subroutine compute_log_likelihood_DEM(l_dem,E_mu_dem,E_K_dem,E_D_dem2,n_rows,N_Bessel,T,xs,dt)
    !Compute log-likelihood of the data given the parameters of the DEM model
        implicit none
        double precision, intent(out) :: l_dem
        double precision, intent(in) :: E_mu_dem,E_K_dem,E_D_dem2,T,dt
        integer*8, intent(in) :: n_rows,N_Bessel
        double precision, dimension(n_rows), intent(in) :: xs
        double precision :: a,b,Ebdt,l,w,x0,xt,v,u,z,log_prod_Bess
        integer*4 :: ii

        l_dem = 0.0d0
        a = E_mu_dem
        b = -E_mu_dem*E_K_dem
        Ebdt = exp(b*dt)
        l = 2.0d0*b/(Ebdt-1.0d0)/E_D_dem2
        w = 2.0d0*a/E_D_dem2-1.0d0
        
        x0 = xs(1)
        do ii = 1, n_rows-1, 1
            xt = xs(ii+1)
            v = l*xt
            u = l*x0*Ebdt
            z = 2.0*l*sqrt(x0*xt*Ebdt)
            log_prod_Bess = -u-v+Stirlings_log_Bessel_Iv(w,z,N_Bessel)
            ! print *, "  w=",w,"z=",z,"xt=",xt,"x0=",x0,"dt=",dt,"lg_B=",log_Bessel_Iv(w,z,N_Bessel),Stirlings_log_Bessel_Iv(w,z,N_Bessel)
            l_dem = l_dem + log(v/u)*(w/2)+log_prod_Bess
            x0 = xt
        end do
        l_dem = l_dem + (n_rows-1)*log(l)
        l_dem = l_dem/T

    end subroutine compute_log_likelihood_DEM

    subroutine Update_system_OU_bridge(s ,t_initial,x_initial, t_final,x_final,D_OU,mu_OU,K_OU)
        !Sample x(s) of OU Bridge conditioned to X(t_initial) = x_initial and X(t_final) = x_final
        !also advances t
        !t_final < s
        implicit none
        double precision, intent(in) :: s ,t_initial,x_initial, t_final,x_final,D_OU,mu_OU,K_OU
        double precision :: mean,std,u,sinH_tfmt0,sinH_tfms,sinH_smt0,xnew
        double precision :: dran_g

        u=dran_g()
        sinH_tfmt0 = sinh(K_OU*(t_final-t_initial))
        sinH_tfms = sinh(K_OU*(t_final-s))
        sinH_smt0 = sinh(K_OU*(s-t_initial))
        mean = (x_initial-mu_OU)*sinH_tfms/sinH_tfmt0+(x_final-mu_OU)*sinH_smt0/sinH_tfmt0 + mu_OU !mean of Gaussian Bridge propagator
        std = sqrt((D_OU**2.0d0)/K_OU *sinH_smt0*sinH_tfms/sinH_tfmt0) !standard deviation of Gaussian Bridge propagator
        t = s
        xnew = mean + std*u

        if ( xnew<0 ) then !This is not due to numerical errors in this case, the OU bridge can perfectly enter x<0
            print*, "x<0 generated!",t,x 
        end if
        x=xnew

    end subroutine Update_system_OU_bridge

    subroutine Fill_gaps_with_OU_bridge_multiresolution(t0,x0,tf,xf,dt_bridge,len_bridge,D_OU,mu_OU,K_OU,xs,ts)
      !Generate len-1 points as a bridge of the OU process linking the points 
      !X(t0)=x0 and X(tf)=xtf
      !The method uses a sort of multiresolution algorithm
      !Store the points generated and the input extremes in the array xs.
      !Also created the array of times where the process is generated, ts.
      !The number of new points generated needs to be bigger or equal than 1
      !The grid generated by this algorithm is regular, ts(ii+1)-ts(ii) = cte for all ii.
      !The information in dt_bridge and len is redundant and needs to be coherent, len = int((tf-t0)/dt_bridge)+1
      implicit none
      integer*8,intent(in) :: len_bridge
      double precision,intent(in) :: t0,x0,tf,xf,dt_bridge
      double precision,dimension (len_bridge),intent(out) :: xs,ts
      double precision, intent(in) :: D_OU,mu_OU,K_OU
      integer*8 :: ii

      if ( (len_bridge < 1)) then
        print *, 'Error: len must be greater than 1.'
        return
      end if
      do ii = 1, len_bridge, 1
        ts(ii) = t0+(ii-1)*dt_bridge
      end do
      xs(1) = x0
      do ii = 2, len_bridge-1, 1
        call Update_system_OU_bridge(ts(ii) ,ts(ii-1),xs(ii-1), tf,xf,D_OU,mu_OU,K_OU)
        xs(ii) = x
      end do
      xs(len_bridge) = xf

    end subroutine Fill_gaps_with_OU_bridge_multiresolution

    subroutine draw_OU_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,D_OU,mu_OU,K_OU,filename)
    !Save file with trajectory of OU bridge in file named filename.
      implicit none
      double precision, intent(in) :: t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,D_OU,mu_OU,K_OU
      character(len=*), intent(in) :: filename
      integer*8 :: len_bridge
      integer*8 :: ii
      double precision, allocatable :: ts_bridge(:),xs_bridge(:)

      len_bridge = int((tf_bridge-t0_bridge)/dt_bridge)+1
      allocate(ts_bridge(len_bridge),xs_bridge(len_bridge))
      ! call Fill_gaps_with_OU_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,xs_bridge,len_bridge,D_B)
      call Fill_gaps_with_OU_bridge_multiresolution(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_bridge,D_OU,mu_OU,K_OU,xs_bridge,ts_bridge)


      print *, "write data in: ",trim(filename)
      print *
      open(unit=1001, file=trim(filename), status='unknown', action='write')
      do ii = 1, len_bridge, 1
          write(1001,*) ts_bridge(ii),xs_bridge(ii)
      end do
      close(1001)

      deallocate(ts_bridge,xs_bridge)
    end subroutine draw_OU_bridge

    subroutine Update_system_Brownian_bridge(s ,t_initial,x_initial, t_final,x_final,D_B_2)
      !Sample x(s) of Brownian Bridge conditioned to X(t_initial) = x_initial and X(t_final) = x_final
      !also advances t
      !t_final < s
      implicit none
      double precision, intent(in) :: s ,t_initial,x_initial, t_final,x_final,D_B_2
      double precision :: mean,std,u,dumt,xnew
      double precision :: dran_g

      u=dran_g()
      dumt = (t_final-t_initial)
      mean = x_initial+(s-t_initial)*(x_final-x_initial)/dumt !mean of Gaussian Bridge propagator
      std = sqrt(D_B_2*(s-t_initial)*(t_final-s)/dumt) !standard deviation of Gaussian Bridge propagator
      t = s
      xnew = mean + std*u

      if ( xnew<0 ) then !This is not due to numerical errors in this case, the Wiener bridge can perfectly enter x<0
          print*, "x<0 generated!",t,x 
      end if
      x=xnew

    end subroutine Update_system_Brownian_bridge

    subroutine Update_system_Lamperti_bridge(s ,t_initial,x_initial, t_final,x_final,bridge_model,n_params_bridge,params_bridge)
    !Sample x(s) of Bridge conditioned to X(t_initial) = x_initial and X(t_final) = x_final
    !also advances t
    !t_final < s
      implicit none
      integer*8, intent(in) :: n_params_bridge
      double precision, intent(in) :: s ,t_initial,x_initial, t_final,x_final
      character(len=*), intent(in) :: bridge_model
      double precision, dimension(n_params_bridge), intent(in) :: params_bridge
      double precision :: mean,std,u,dumt,xnew
      double precision :: mu_OU,K_OU,sinH_tfmt0,sinH_tfms,sinH_smt0
      double precision :: dran_g

      u=dran_g()
      dumt = (t_final-t_initial)
      if ( trim(bridge_model) == "WI" ) then
        mean = x_initial+(s-t_initial)*(x_final-x_initial)/dumt
        std = sqrt((s-t_initial)*(t_final-s)/dumt) 
        xnew = mean + std*u !Gaussian bridge
      else if ( trim(bridge_model) == "GB" ) then
        mean = log(x_initial)+(s-t_initial)*log(x_final/x_initial)/dumt
        std = sqrt((s-t_initial)*(t_final-s)/dumt) 
        xnew = exp (mean + std*u) !Gaussian bridge
      else if ( trim(bridge_model) == "OU" ) then
        mu_OU = params_bridge(1)
        K_OU = params_bridge(2)
        sinH_tfmt0 = sinh(K_OU*(t_final-t_initial))
        sinH_tfms = sinh(K_OU*(t_final-s))
        sinH_smt0 = sinh(K_OU*(s-t_initial))
        mean = (x_initial-mu_OU)*sinH_tfms/sinH_tfmt0+(x_final-mu_OU)*sinH_smt0/sinH_tfmt0 + mu_OU 
        std = sqrt(1.0d0/K_OU *sinH_smt0*sinH_tfms/sinH_tfmt0) 
        xnew = mean + std*u !Gaussian bridge
      else 
        print *, "Fatal error, bridge model not recognized"
        stop
      end if
      
      t = s
      

      ! if ( xnew<0 ) then !This is not due to numerical errors in this case, the Wiener bridge can perfectly enter x<0
      !     print*, "x<0 generated!",t,x 
      ! end if
      x=xnew

    end subroutine Update_system_Lamperti_bridge

    subroutine Fill_gaps_with_Brownian_bridge_multiresolution(t0,x0,tf,xf,dt_bridge,len_bridge,D_B,xs,ts)
      !Generate len-1 points as a bridge of the brownian process linking the points 
      !X(t0)=x0 and X(tf)=xtf
      !The method uses a sort of multiresolution algorithm
      !Store the points generated and the input extremes in the array xs.
      !Also created the array of times where the process is generated, ts.
      !The number of new points generated needs to be bigger or equal than 1
      !The grid generated by this algorithm is regular, ts(ii+1)-ts(ii) = cte for all ii.
      !The information in dt_bridge and len is redundant and needs to be coherent, len = int((tf-t0)/dt_bridge)+1
      !It is a bit weird to pass both the dt and len of the bridge, but it is easier to define the array "xs" inside the subroutine in this way
      implicit none
      integer*8,intent(in) :: len_bridge
      double precision,intent(in) :: t0,x0,tf,xf,dt_bridge
      double precision,dimension (len_bridge),intent(out) :: xs,ts
      double precision, intent(in) :: D_B
      double precision :: D_B_2
      integer*8 :: ii

      if ( (len_bridge < 1)) then
        print *, 'Error: len must be greater than 1.'
        return
      end if
      do ii = 1, len_bridge, 1
        ts(ii) = t0+(ii-1)*dt_bridge
      end do
      xs(1) = x0
      D_B_2 = D_B*D_B
      do ii = 2, len_bridge-1, 1
        call Update_system_Brownian_bridge(ts(ii) ,ts(ii-1),xs(ii-1), tf,xf,D_B_2)
        xs(ii) = x
      end do
      xs(len_bridge) = xf

    end subroutine Fill_gaps_with_Brownian_bridge_multiresolution

    subroutine draw_Brownian_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,D_B,filename)
    !Save file with trajectory of Brownian bridge in file named filename.
      implicit none
      double precision, intent(in) :: t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,D_B
      character(len=*), intent(in) :: filename
      integer*8 :: len_bridge
      integer*8 :: ii
      double precision, allocatable :: ts_bridge(:),xs_bridge(:)

      len_bridge = int((tf_bridge-t0_bridge)/dt_bridge)+1
      allocate(ts_bridge(len_bridge),xs_bridge(len_bridge))
      ! call Fill_gaps_with_Brownian_bridge(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,xs_bridge,len_bridge,D_B)
      call Fill_gaps_with_Brownian_bridge_multiresolution(t0_bridge,x0_bridge,tf_bridge,xf_bridge,dt_bridge,len_bridge,D_B,xs_bridge,ts_bridge)


      print *, "write data in: ",trim(filename)
      print *
      open(unit=1001, file=trim(filename), status='unknown', action='write')
      do ii = 1, len_bridge, 1
          write(1001,*) ts_bridge(ii),xs_bridge(ii)
      end do
      close(1001)

      deallocate(ts_bridge,xs_bridge)
    end subroutine draw_Brownian_bridge

    subroutine Update_system_OU(tf)
    !Update system from x(t)=x to x(tf), x and t are public variables
        implicit none
        double precision, intent(in) :: tf
        double precision :: mean,std,u,xnew,dt,kdt
        double precision :: dran_g,dran_gbmw

        dt = tf-t
        ! u = dran_g()
        u = dran_gbmw()
        kdt = K*dt
        std = D/sqrt(2.0d0*K)*sqrt(1.0d0-exp(-2.0d0*kdt))
        mean = mu + (x-mu)*exp(-kdt)  
        x = mean + std*u
        t = tf

        if ( x<0 ) then
            print*, "x<0 generated!",t,x
            !x=x reject update
        endif

    end subroutine Update_system_OU

    subroutine Update_system_dem_approx(tf , dt)
        !Update system from x(t)=x to x(tf), x and t are public variables
        !This is the Milstein scheme for the DEM model
        implicit none
        double precision, intent(in) :: tf,dt
        double precision :: drift, noise,u,milstein,xnew,sqdt,D2
        double precision :: dran_g
        integer*8 :: steps,ii

        steps=int8(nint((tf-t)/dt))
        sqdt = sqrt(dt)
        D2 = D*D
        milstein = D2*dt/4.0

        do ii = 1, steps
            t=t+dt
            u=dran_g()
            drift = K*(mu-x)*dt 
            noise = D*sqrt(x)*sqdt
            xnew=x+drift+noise*u+milstein*(u*u-1.0d0)

            if ( xnew<0 ) then
                print*, "Warning in Update_system_dem_approx"
                print*, "x<0 generated!",t,x,xnew
            else
                x=xnew
            end if
        end do

    end subroutine Update_system_dem_approx

    subroutine Update_system_dem_exact(next_t)
        !Update system from x(t)=x to x(next_t), x and t are public variables
        !This is the exact scheme for the DEM model
        !K=0 will not work here, as the limit beta/(eb-1) (see below) is not properly computed (it is possible to use Update_system_dem_approx instead).
        !See Dornic, I., Chaté, H., & Munoz, M. A. (2005). Integration of langevin equations with multiplicative noise and the viability of field theories for absorbing phase transitions. Physical Review Letters, 94(10), 18–21. https://doi.org/10.1103/PhysRevLett.94.100601
        implicit none
        double precision, intent(in) :: next_t
        double precision :: lambda,nu,alpha,beta,D2,dt,eb
        double precision :: dran_gamma
        integer*8 :: steps,ii,n
        integer*8 :: iran_poisson

        dt = next_t-t
        D2 = D*D
        beta = -K
        alpha = mu*K
        eb = exp(beta*dt)
        lambda = 2.0d0*beta/D2/(eb-1.0d0)
        nu = 2.0d0*alpha/D2-1.0d0
        n = iran_poisson(lambda*x*eb)
        x = dran_gamma(n+nu+1.0d0)/lambda
        t = next_t

    end subroutine Update_system_dem_exact
    
    subroutine Update_system_env(tf , dt)
        !Update system from x(t)=x to x(tf), x and t are public variables
        implicit none
        double precision, intent(in) :: tf,dt
        double precision :: drift, noise,u,milstein,xnew,sqdt,D2
        double precision :: dran_g
        integer*8 :: steps,ii

        steps=int8(nint((tf-t)/dt))
        ! print*, " steps=",steps
        sqdt = sqrt(dt)
        D2 = D*D

        do ii = 1, steps
            u=dran_g()
            
            drift    = K * x * (mu - x) * dt
            noise    = D * x * sqdt
            milstein = 0.5d0 * D2 * x * dt

            xnew = x + drift + noise * u + milstein * (u*u - 1.0d0)
            t=t+dt
            if (xnew >= 0.0d0) then
              x  = xnew
            else
                print*, "Warning: x < 0 at t =", t, " with x =", x
                ! reject update, x remains unchanged
            end if
        end do

    end subroutine Update_system_env

    subroutine Update_system_CP(tf , dt)
        !Update system from x(t)=x to x(tf), x and t are public variables
        implicit none
        double precision, intent(in) :: tf,dt
        double precision :: drift, noise,u,milstein,xnew,sqdt,D2
        double precision :: dran_g
        integer*8 :: steps,ii

        steps=int8(nint((tf-t)/dt))
        ! print*, " steps=",steps
        sqdt = sqrt(dt)
        D2 = D*D
        milstein = D2*dt/4.0
        do ii = 1, steps
            t=t+dt
            u=dran_g()
            
            drift = K*x*(mu-x)*dt 
            noise = D*sqrt(x)*sqdt
            xnew=x+drift+noise*u+milstein*(u*u-1.0d0)

            if ( xnew<0 ) then
                print*, "x<0 generated!",t,x
                !x=x reject update
            else
                x=xnew
                ! print *, "x =",x,xnew
            end if
        end do

    end subroutine Update_system_CP

    subroutine Update_system_general_gamma_model(tf , dt, theta)
        !Update system from x(t)=x to x(tf), x and t are public variables
        implicit none
        double precision, intent(in) :: tf,dt,theta
        double precision :: drift, noise,u,milstein,xnew,sqdt,D2,dtheta,xth
        double precision :: dran_g
        integer*8 :: steps,ii

        steps=int8(nint((tf-t)/dt))
        sqdt = sqrt(dt)
        D2 = D*D
        dtheta = (theta+1.0d0)/2.0d0
        
        do ii = 1, steps
            t=t+dt
            u=dran_g()
            xth = x**(theta)
            milstein = D2*dtheta*xth*dt/2.0
            drift = K*(mu-x)*xth*dt 
            noise = D*(x**(dtheta))*sqdt
            xnew=x+drift+noise*u+milstein*(u*u-1.0d0)

            if ( xnew<0 ) then
                print*, "x<0 generated!",t,x
                !x=x reject update
            else
                x=xnew
                ! print *, "x =",x,xnew
            end if
        end do

    end subroutine Update_system_general_gamma_model

    subroutine read_xF_tf()
      !Example of subroutine, read data file
      implicit none
      integer*4 :: ios,i,dum,datapoints
      double precision :: dummy
  
      open(unit=1001, file="data/F_Target_test3.dat", iostat=ios, status="old", action="read")
      if ( ios /= 0 ) stop "Error opening file "
      do i = 1, datapoints, 1
        read(1001,*) dummy,dummy,dummy,dum
      end do
  
  
    end subroutine read_xF_tf
    
    subroutine search_list_binary_algorithm(list,position,p)
      !Rafle event:
      !Given a list of probabilities called "list" such that sum(list)=1 and a probability p.
      !Look for "position" such that C(position)>=p and C(j)<p for all j in [1,position[.
      !Where C is the cumulative of list: C(i)=list(1)+list(2)+...+list(i)-
      !REFERENCE:Brainerd, W. S. (2015). Guide to Fortran 2008 programming (p. 141). Berlin: Springer.
  
      implicit none
      double precision, dimension(:), intent (in) :: list
      double precision, intent (in) :: p
      integer*4, intent(out) :: position
      double precision, dimension(size(list)) :: C
      integer*4 ii,N,first,last,half
  
      N=size(list) !It would be cool to define this as a parameter (constant), I don't know how...
      C(1)=list(1)
      do ii = 2, N, 1 !Compute cumulative of list
        c(ii)=C(ii-1)+list(ii)
      end do
  
      first=1;last=N
      do while ( first.ne.last )
        half=(first+last)/2
        if ( p>C(half) ) then
          first=half+1
        else
          last=half
        end if
      end do
  
      position=first
  
    end subroutine search_list_binary_algorithm
  
    subroutine float_to_string(a, n, result)
      !Convert float number a to string with format int(a)//"d"//dec(a).
      !E.g. 0.00543--> "0d00543"
      !For dec(a) select first n decimals
      !If number of digits in dec(a)>n then add zeros to the left
      implicit none
      double precision, intent(in) :: a
      integer, intent(in) :: n
      character(len=*), intent(out) :: result
      character(len=n) :: dumc
      integer :: digit,zeros
  
    
      integer :: int_part
      double precision :: dec_part
    
      ! Get the integer and decimal parts
      int_part = int(a)
      dec_part = a - int_part
  
      ! Convert the decimal part to a string with n digits
      zeros=0
      digit=int(dec_part*10**(zeros+1))
      do while ((digit==0).and.(zeros<n))
        zeros=zeros+1
        digit=int(dec_part*10**(zeros+1))
      enddo
      dumc=repeat("0",zeros)//trim(str(int(dec_part*10**(n)))) 
  
  
      result=trim(str(int_part))//"d"//trim(dumc)
  
    
    end subroutine float_to_string

    subroutine logspace(start, end, num, base, result)
      ! Inputs
      implicit none
      double precision, intent(in) :: start      ! Start of the range (logarithmic)
      double precision, intent(in) :: end       ! End of the range (logarithmic)
      integer*8, intent(in) :: num        ! Number of points
      double precision, intent(in) :: base       ! Base of the logarithmic scale

      ! Output
      double precision,dimension(num), intent(out) :: result ! Resulting array

      ! Local variables
      double precision :: step, exponent
      integer*8 :: ii

      ! Check for valid number of points
      if (num < 2) then
        print *, "Error: num must be at least 2."
        stop
      end if

      ! Calculate the step size in logarithmic space
      step = (end - start) / real(num - 1)

      ! Generate the values
      do ii = 1, num
        exponent = start + step * real(ii - 1)
        result(ii) = base ** exponent
      end do

    end subroutine logspace

  
end module subroutines
