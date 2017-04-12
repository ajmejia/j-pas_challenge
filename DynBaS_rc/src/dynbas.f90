PROGRAM DynBaS

  USE DYNBAS_PIECES_MODULE

  IMPLICIT NONE

  integer                     :: i, j, k, n, m
  integer                     :: unit, iostatus
  integer                     :: nwef, nwlm, nwlo, nage, nmet, nmod, navs, nrsh, ind(3)
  integer                     :: nhdr, ncol
  integer                     :: it1d(1), it2d(2), it3d(3), ibasis(12)
  integer, allocatable        :: filter_ids(:), iz(:), gen1d(:, :), gen2d(:, :), gen3d(:, :)
  real, allocatable           :: redshifts(:), uni_ages(:)
  real, allocatable           :: t(:), ages(:), wlength(:), mods(:, :, :)
  real, allocatable           :: basis(:, :), fssps(:, :), Zs(:), Zs_(:)
  real, allocatable           :: weff(:), wlo(:), fluxo_s(:), sigma_s(:), dyn_flux(:, :)
  real, allocatable           :: sed(:, :), rflux(:), Av_grid(:)
  real, allocatable           :: models(:, :, :)
  real(8), allocatable        :: fluxo(:), sigma(:), chi_list(:, :)
  real(8), allocatable        :: coeff1d(:, :), coeff2d(:, :), coeff3d(:, :)
  real                        :: t1d(1), t2d(2), t3d(3), m1d(1), m2d(2), m3d(3)
  real                        :: chi(3), z1d(1), z2d(2), z3d(3)
  real                        :: mass_mod(3), t_m_mod(3), t_fr_mod(3), z_m_mod(3), z_fr_mod(3)
  real                        :: Av_mod(3), losvd_mod
  real                        :: z, universe_age(1), basis_ages(12)
  character(:), allocatable   :: models_path
  character(256)              :: record, input_file, output_file
  character(100), allocatable :: long_opts(:), options(:), values(:), arguments(:)
  character(80), allocatable  :: model_ids(:), fname(:)*64
  character(34), allocatable  :: kwords(:), files_ised(:), sel_files_ised(:)
  character(:), allocatable   :: dir_out
  logical, allocatable        :: mask(:), sed_mask(:)
  logical                     :: verbose, given_tu, given_losvd

	CALL INITIALIZE

	CALL PARSE_OPTS_AND_ARGS

	CALL LOAD_SSP_MODELS

  if(verbose) CALL VERBOSE_

  write(6, *)

  allocate(models(nwlo, nmod, navs))
  if(allocated(filter_ids)) then

    write(6, "(A)", advance="no") "Computing model grid..."

    models(:, :, 1) = EFFECTIVE_FLUX(filter_ids, wlength, fssps, z)

    !$OMP PARALLEL DO
    do m=1, navs
      models(:, :, m) = models(:, :, 1) * spread(CCM(wlo/(z+1), Av_grid(m)), 2, nmod)
    end do
    !$OMP END PARALLEL DO

    write(6, *) "         done."
  else
    write(6, "(A)", advance="no") "Computing model grid..."

    models(:, :, 1) = LINEAR_INTER(wlength*(1+z), fssps, wlo)

    if(.not.given_losvd .or. losvd_mod<=40.0) then
      do j = 1, size(ibasis)
        CALL NEAREST_ELEMENT(t(:nage), basis_ages(j), ibasis(j))
      end do

      CALL SUBMATRIX(models(:, :nage, 1), 2, ibasis, basis)

      allocate(fluxo_s(nwlo), sigma_s(nwlo))
      fluxo_s = fluxo
      sigma_s = sigma

      CALL VELOCITY_DISP(wlo, basis, fluxo_s, sigma_s, [5175.-33, 5175.+33], losvd_mod)
    end if

    if(losvd_mod>40.0) CALL APPLY_VD(wlo, models(:, :, 1), losvd_mod)

    !$OMP PARALLEL DO
    do m=1, navs
      models(:, :, m) = models(:, :, 1) * spread(CCM(wlo/(z+1), Av_grid(m)), 2, nmod)
    end do
    !$OMP END PARALLEL DO

    write(6, *) "         done."
  end if

  write(6, "(A)", advance="no") "Computing minimum chi square..."

  !$OMP PARALLEL DO
  do m=1, navs
    CALL MIN_CHI_SQUARE(fluxo/sigma, dble(models(:, :, m))/spread(sigma, 2, nmod), &
                        coeff1d(:, m), coeff2d(:, m), coeff3d(:, m), &
                        gen1d(:, m), gen2d(:, m), gen3d(:, m), chi_list(:, m))
  end do
  !$OMP END PARALLEL DO

  write(6, *) " done."
  write(6, *)

  ind = minloc(chi_list, dim=2)

  Av_mod(1) = Av_grid(ind(1))
  Av_mod(2) = Av_grid(ind(2))
  Av_mod(3) = Av_grid(ind(3))

  m1d  = coeff1d(:, ind(1))
  m2d  = coeff2d(:, ind(2))
  m3d  = coeff3d(:, ind(3))
  it1d = gen1d(:, ind(1))
  it2d = gen2d(:, ind(2))
  it3d = gen3d(:, ind(3))

  chi  = [(chi_list(j, ind(j)), j=1, 3)]

  if(any(it1d<1) .or. any(it1d>nmod)) then
    dyn_flux(:, 1) = 0.0
    mass_mod(1)    = 0.0
    t_m_mod(1)     = 0.0
    t_fr_mod(1)    = 0.0
    z_m_mod(1)     = 0.0
    z_fr_mod(1)    = 0.0
  else
    t1d = t(it1d(1))
    z1d = Zs(it1d(1))

    dyn_flux(:, 1) = m1d(1) * models(:, it1d(1), ind(1))
    mass_mod(1)    = m1d(1)
    t_m_mod(1)     = log10(t1d(1))
    t_fr_mod(1)    = log10(t1d(1))
    z_m_mod(1)     = log10(Zs(it1d(1)) / 0.0199)
    z_fr_mod(1)    = log10(Zs(it1d(1)) / 0.0199)
  end if
  if(any(it2d<1) .or. any(it2d>nmod)) then
    dyn_flux(:, 2) = 0.0
    mass_mod(2)    = 0.0
    t_m_mod(2)     = 0.0
    t_fr_mod(2)    = 0.0
    z_m_mod(2)     = 0.0
    z_fr_mod(2)    = 0.0
  else
    t2d = [(t(it2d(j)),  j=1, 2)]
    z2d = [(Zs(it2d(j)), j=1, 2)]

    CALL MONOCHROMATIC_SEDS(wlength, fssps*spread(CCM(wlength, Av_grid(ind(2))), 2, nmod), 19, 0.0, it2d, rflux)

    dyn_flux(:, 2) = m2d(1) * models(:, it2d(1), ind(2)) + m2d(2) * models(:, it2d(2), ind(2))
    mass_mod(2)    = sum(m2d)
    t_m_mod(2)     = MEAN(log10(t2d), m2d)
    t_fr_mod(2)    = MEAN(log10(t2d), m2d*rflux)
    z_m_mod(2)     = MEAN(log10(z2d/0.0199), m2d)
    z_fr_mod(2)    = MEAN(log10(z2d/0.0199), m2d*rflux)
  end if
  if(any(it3d<1) .or. any(it3d>nmod)) then
    dyn_flux(:, 3) = 0.0
    mass_mod(3)    = 0.0
    t_m_mod(3)     = 0.0
    t_fr_mod(3)    = 0.0
    z_m_mod(3)     = 0.0
    z_fr_mod(3)    = 0.0
  else
    t3d = [(t(it3d(j)),  j=1, 3)]
    z3d = [(Zs(it3d(j)), j=1, 3)]

    CALL MONOCHROMATIC_SEDS(wlength, fssps*spread(CCM(wlength, Av_grid(ind(3))), 2, nmod), 19, 0.0, it3d, rflux)

    dyn_flux(:, 3) = m3d(1) * models(:, it3d(1), ind(3)) + m3d(2) * models(:, it3d(2), ind(3)) + m3d(3) * models(:, it3d(3), ind(3))
    mass_mod(3)    = sum(m3d)
    t_m_mod(3)     = MEAN(log10(t3d), m3d)
    t_fr_mod(3)    = MEAN(log10(t3d), m3d*rflux)
    z_m_mod(3)     = MEAN(log10(z3d/0.0199), m3d)
    z_fr_mod(3)    = MEAN(log10(z3d/0.0199), m3d*rflux)
  end if

  where(fluxo<=0.0) sigma = 0.0

  m = index(input_file, "/", back=.true.)
  n = index(input_file, ".", back=.true.)
  if(n<m) n = len_trim(input_file) + 1

  output_file = dir_out//"dynbasfit_"//input_file(m+1:n-1)//".log"

  unit = AVAILABLE_UNIT()
  open(unit, file=trim(output_file), action="write")

  write(unit,                         "(2A)") "# SED file         = ", trim(input_file)
  write(unit,                          "(A)") "#"
  write(unit, "(A,"//STRING(nage)//"EN13.2)") "# ages             = ", ages
  write(unit,  "(A,"//STRING(nmet)//"G10.3)") "# metallicities    = ", Zs_ / 0.0199
  write(unit,   "(A,"//STRING(navs)//"F5.2)") "# Av sample        = ", Av_grid
  write(unit,                          "(A)") "#"
  write(unit,                      "(A,1I7)") "# generator 1d     = ", it1d
  write(unit,                      "(A,2I7)") "# generator 2d     = ", it2d
  write(unit,                      "(A,3I7)") "# generator 3d     = ", it3d
  write(unit,                  "(A,1EN13.2)") "# ages 1d          = ", t1d
  write(unit,                  "(A,2EN13.2)") "# ages 2d          = ", t2d
  write(unit,                  "(A,3EN13.2)") "# ages 3d          = ", t3d
  write(unit,                   "(A,1G10.3)") "# metallicities 1d = ", z1d / 0.0199
  write(unit,                   "(A,2G10.3)") "# metallicities 2d = ", z2d / 0.0199
  write(unit,                   "(A,3G10.3)") "# metallicities 3d = ", z3d / 0.0199
  write(unit,                   "(A,1E10.2)") "# coefficients 1d  = ", m1d
  write(unit,                   "(A,2E10.2)") "# coefficients 2d  = ", m2d
  write(unit,                   "(A,3E10.2)") "# coefficients 3d  = ", m3d
  write(unit,                          "(A)") "#"
  write(unit,                    "(A,F10.5)") "# redshift         = ", z
  if(allocated(filter_ids)) then
  write(unit,                         "(2A)") "# LOSVD            = ", "None"
  else
  write(unit,                     "(A,F8.2)") "# LOSVD            = ", losvd_mod
  end if
  write(unit,                         "(4A)") "# prop.\generator  = ", "             1d",&
                                                                       "             2d",&
                                                                       "             3d"
  write(unit,                   "(A,3E15.6)") "# M/Mo             = ", mass_mod
  write(unit,                   "(A,3F15.6)") "# <log t>_M        = ", t_m_mod
  write(unit,                   "(A,3F15.6)") "# <log t>_Lr       = ", t_fr_mod
  write(unit,                   "(A,3F15.6)") "# <log Z/Zo>_M     = ", z_m_mod
  write(unit,                   "(A,3F15.6)") "# <log Z/Zo>_Lr    = ", z_fr_mod
  write(unit,                   "(A,3F15.6)") "# Av[mag]          = ", Av_mod
  write(unit,                   "(A,3F15.6)") "# chi square       = ", chi
  write(unit,                          "(A)") "#"
  write(unit,                         "(6A)") "#    wavelength", "       obs_flux",&
                                              "          sigma", "        1d_flux",&
                                              "        2d_flux", "        3d_flux"
  do i = 1, nwlo
  write(unit, "(F15.2,5E15.6)")   wlo(i),          fluxo(i),          sigma(i),    dyn_flux(i, 1), &
                                  dyn_flux(i, 2),    dyn_flux(i, 3)
  end do

	CONTAINS

  SUBROUTINE HELP()
		write(6, *)
		write(6, *) "Usage"
		write(6, *) "-----"
		write(6, *) "  dynbas [options] <input_file>"
		write(6, *)
		write(6, *) "  <input_file>         : a plain text file containing columns:"
    write(6, *) "                           * wlength, flux, sigma; or"
    write(6, *) "                           * wlength, flux, sigma, mask"
		write(6, *) "  [options]            : any of the following ones."
		write(6, *)
		write(6, *) "Options"
		write(6, *) "-------"
		write(6, *)
    write(6, *) "  --output-dir=<path>  : path to the directory where log (output files) shall be"
    write(6, *) "                         placed."
		write(6, *) "  --passbands=<values> : a comma-separated list of filters' IDs as defined in the"
		write(6, *) "                         BC/CB srcs."
    write(6, *) "  --Av-grid=<values>   : dust extinction range (Av) in which the model will be"
    write(6, *) "                         constraint. If a third value is given, is will be taken"
    write(6, *) "                         as the number of values in the grid. All values SHOULD be"
    write(6, *) "                         integers."
		write(6, *) "  --models=<keys>      : a list of comma-separated BC/CB models keys."
		write(6, *) "                         Defaults to m22,..,m72."
    write(6, *) "  -v                   : verbose mode. Use this option to print on screen"
    write(6, *) "                         meaningful information of the current run and SED fitting"
    write(6, *) "                         setups."
		write(6, *) "  -h, --help           : display this message."
		write(6, *)
    stop
  END SUBROUTINE

	SUBROUTINE INITIALIZE()

		dir_out = "./"

		nmet = 20
		allocate(model_ids(nmet))

	  verbose      = .false.
    given_tu     = .false.
    given_losvd  = .false.
		model_ids    = ["Z0.0017", "Z0.0034", "Z0.0051",&
                    "Z0.0067", "Z0.0084", "Z0.0100", "Z0.0117",&
                    "Z0.0133", "Z0.0149", "Z0.0166", "Z0.0183",&
                    "Z0.0199", "Z0.0216", "Z0.0233", "Z0.0249",&
                    "Z0.0265", "Z0.0282", "Z0.0298", "Z0.0315",&
                    "Z0.0332"]
    z            = 0.0
	  universe_age = 13.75e9
    Av_grid      = [0.        , 0.18421053, 0.36842105, 0.55263158, 0.73684211,&
                    0.92105263, 1.10526316, 1.28947368, 1.47368421, 1.65789474,&
                    1.84210526, 2.02631579, 2.21052632, 2.39473684, 2.57894737,&
                    2.76315789, 2.94736842, 3.13157895, 3.31578947, 3.5       ]
    basis_ages   = [1.78e5, 4.17e6, 8.71e6, 1.91e7, 3.60e7, 9.05e7,&
                    3.60e8, 1.28e9, 2.75e9, 4.75e9, 9.25e9, 13.75e9]
	  long_opts    = ["models=    ",&
                    "help       ",&
                    "id=        ",&
                    "output-dir=",&
                    "passbands= ",&
                    "Av-grid=   "]
    kwords       = [character(len(kwords))::"redshift", "age", "LOSVD", "Av"]

	END SUBROUTINE

	SUBROUTINE PARSE_OPTS_AND_ARGS()

		CALL PARSE_ARGS("hivz:", long_opts, options, values, arguments)

    if(.not.(allocated(options) .and. allocated(values) .and. allocated(arguments))) CALL HELP

		do m=1, size(options)
			select case(options(m))
	    case("v")
	      verbose = .true.
			case("models")
				nmet = STRING_COUNT(trim(values(m)), ",") + 1

				deallocate(model_ids)
				allocate(model_ids(nmet))
				read(values(m), *) model_ids
      case("output-dir")
        n       = len_trim(values(m))
        dir_out = values(m)(:n)

        if(dir_out(n:n) /= "/") dir_out = dir_out//"/"
      case("extinction-range")
        n = STRING_COUNT(values(m), ",") + 1
        if(n==2) then
          Av_grid = EVAL(SPLIT(values(m), ","))
          Av_grid = [(Av_grid(1) + (i-1)*(Av_grid(2)-Av_grid(1))/10, i=1, 10+1)]
        else if(n==3) then
          Av_grid = EVAL(SPLIT(values(m), ","))
          Av_grid = [(Av_grid(1) + (i-1)*(Av_grid(2)-Av_grid(1))/int(Av_grid(3)), i=1, int(Av_grid(3))+1)]
        end if
      case("passbands")
        n = index(values(m), ",...,")
        if(n==0) then
          CALL APPEND(filter_ids, EVAL(SPLIT(values(m), ",")))
        else
          do k=1, STRING_COUNT(values(m), ",...,")
            i = index(values(m)(:n-1), ",", back=.true.)
            j = index(values(m)(n+5:), ",")
            if(i+j==0) then
              CALL APPEND(filter_ids, LISTING(values(m)))
            else if(i==0) then
              CALL APPEND(filter_ids, LISTING(values(m)(:j+n+3)))
            else if(j==0) then
              CALL APPEND(filter_ids, EVAL(SPLIT(values(m)(:i-1), ",")))
              CALL APPEND(filter_ids, LISTING(values(m)(i+1:)))
            else
              CALL APPEND(filter_ids, EVAL(SPLIT(values(m)(:i-1), ",")))
              CALL APPEND(filter_ids, LISTING(values(m)(i+1:j+n+3)))
            end if
            values(m) = values(m)(j+n+5:)
            n = index(values(m), ",...,")
          end do
          if(EVAL(values(m))/=filter_ids(size(filter_ids))) CALL APPEND(filter_ids, EVAL(SPLIT(values(m), ",")))
        end if
        nwef = size(filter_ids)
        if(nwef<3) STOP "DynBaS: number of filters must be >= 3 for the system to have one solution."

     		allocate(fname(nwef), weff(nwef))

        do i=1, nwef
          CALL FILTER_INFO(filter_ids(i), mean_wl=weff(i), filter_name=fname(i))
        end do
			case("h", "help")
        CALL HELP
			end select
		end do

		if(size(arguments)/=1) STOP "DynBaS: number of arguments must be 1. Use --help option to see the documentation."

		input_file = arguments(1)
		CALL INQUIRE_FILE(input_file, nhdr, nwlo, ncol, status=iostatus)
		CALL IOERR(iostatus, input_file)

    unit = AVAILABLE_UNIT()
    open(unit, file=input_file, action="read")

    if(nhdr>1) then
      do i=1, nhdr
        read(unit, "(A)") record

        do j=1, size(kwords)
          if(index(record, trim(kwords(j)))==0) cycle

          select case(j)
          case(1)
            read(record(index(record, "=")+1:), *) z
          case(2)
            read(record(index(record, "=")+1:), *) universe_age(1)
            universe_age = universe_age*1e9
            given_tu = .true.
          case(3)
            read(record(index(record, "=")+1:), *) losvd_mod
            given_losvd = .true.
          case(4)
            read(record(index(record, "=")+1:), *) Av_mod(1)
            Av_mod  = Av_mod(1)
            Av_grid = [Av_mod(1)]
          end select
          exit
        end do
      end do
    end if

    select case(ncol)
    case(3)
      allocate(sed(nwlo, ncol))
      do i=1, nwlo
        read(unit, *) sed(i, :)
      end do
    case(4)
      allocate(sed(nwlo, ncol-1), sed_mask(nwlo))
      do i=1, nwlo
        read(unit, *) sed(i, :), sed_mask(i)
      end do
    case default
      STOP "DynBaS: input file must have at most 4 columns. Use --help option to see the documentation."
    end select
    close(unit)

    allocate(wlo(nwlo), fluxo(nwlo), sigma(nwlo))
    wlo   = sed(:, 1)
    fluxo = sed(:, 2)
    sigma = sed(:, 3)
    where(fluxo<=0.0) sigma = huge(1.0)

	END SUBROUTINE

	SUBROUTINE LOAD_SSP_MODELS()

	  CALL GET_PATH("MODELS", path=models_path)
    CALL INQUIRE_FILE(models_path//"challenge.list", nrow=nmod, status=iostatus)
	  CALL IOERR(iostatus, models_path//"challenge.list")
	  unit = AVAILABLE_UNIT()
	  open(unit, file=models_path//"challenge.list", action="read")

	  allocate(files_ised(nmod), sel_files_ised(nmet), iz(nmet), mask(nmod))
          mask = .false.

		do i=1, nmod
			read(unit, *) files_ised(i)
		end do
		close(unit)

    if(.not.given_tu .and. z/=0.0) then
      CALL INQUIRE_FILE("../data/universe_age.dat", nhdr, nrsh, ncol, iostatus)
      CALL IOERR(iostatus, "../data/universe_age.dat")

      allocate(redshifts(nrsh), uni_ages(nrsh))

      open(unit, file="../data/universe_age.dat", action="read", iostat=iostatus)

      do i=1, nhdr
        read(unit, *)
      end do

      do i=1, nrsh
        read(unit, *) redshifts(i), uni_ages(i)
      end do
      close(unit)

      uni_ages = uni_ages*1e9

      universe_age = LINEAR_INTER(redshifts, uni_ages, [z])
    end if

		do i=1, nmet
			CALL STRING_LOC(files_ised, model_ids(i), iz(i))
      if(iz(i)==0) then
			   write(6, "(3A)") "STOP DynBaS: the model labelled '", trim(model_ids(i)), "' was not found in models list 'models.list'."
         stop
      end if
			mask(iz(i)) = .true.
		end do

		sel_files_ised = pack(files_ised, mask)

		CALL READ_SED(models_path//sel_files_ised, wlength, mods, ages, Zs_, .true., .false.)

		n = minloc(abs(ages-universe_age(1)), dim=1)
    if(ages(n)>universe_age(1)) n = n-1

	  CALL SLICE(ages, 1, n)

	  nage = size(ages)
	  nwlm = size(wlength)
	  nmod = nage * nmet
    navs = size(Av_grid)

	  allocate(fssps(nwlm, nmod), t(nmod), Zs(nmod))
	  t     = reshape(spread(ages, 2, nmet), shape=[nmod])
	  Zs    = reshape(spread(Zs_, 1, nage), shape=[nmod])
	  fssps = reshape(mods(:, 1:n, :), shape(fssps))

    deallocate(mods)
    allocate(chi_list(3, navs))
    allocate(coeff1d(1, navs), coeff2d(2, navs), coeff3d(3, navs))
    allocate(gen1d(1, navs), gen2d(2, navs), gen3d(3, navs))
    allocate(dyn_flux(nwlo, 3))

	END SUBROUTINE

	SUBROUTINE VERBOSE_()

		write(6, *)
    write(6, *) "Available ISED models (selected, °):"
    write(6, *) "----------------------------------- "
    do i=1, size(files_ised)
    if(any(trim(files_ised(i))==sel_files_ised)) then
    write(6,             "(2A)") " ° ", trim(files_ised(i))
    else
    write(6,             "(2A)") "   ", trim(files_ised(i))
    end if
    end do
    write(6, *)
    write(6, *) "SED fitting parameters:"
    write(6, *) "----------------------"
    write(6,              "(A)") " * Model grid dimensions:"
    write(6,           "(A,I7)") "     > SSPs               = ", nage*nmet
    write(6,           "(A,I7)") "     > metallicities      = ", nmet
    write(6,           "(A,I7)") "     > ages               = ", nage
    write(6,           "(A,I7)") "     > extinctions        = ", navs
    write(6,           "(A,I7)") "     > wavelengths        = ", nwlo
    write(6,           "(A,I7)") "     > eff. wavelengths   = ", count(fluxo>0.0)
    write(6,              "(A)") " * Physical quantities:"
    write(6,         "(A,F7.3)") "     > z (redshift)       = ", z
    write(6,         "(A,F7.3)") "     > Z_sun              = ", 0.0199
    write(6, "(A,F7.3,A2,F7.3)") "     > age range [log/yr] = ", log10(minval(ages)), "-", log10(maxval(ages))
    write(6, "(A,F7.3,A2,F7.3)") "     > Z   range [Z_sun]  = ", minval(Zs_)/0.0199,    "-", maxval(Zs_)/0.0199
    write(6, "(A,F7.3,A2,F7.3)") "     > Av  range [mag]    = ", minval(Av_grid),     "-", maxval(Av_grid)

    if(allocated(filter_ids)) then
    write(6,              "(A)") " * Passbands:"
    do i=1, nwef
    write(6,             "(2A)") "     > ", trim(fname(i))
    end do
    end if

	END SUBROUTINE

END PROGRAM
