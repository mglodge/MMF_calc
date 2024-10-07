Light scattering code that finds extinction, scattering and absorption efficiencies, as well as asymmetry parameter, for a user-input range of wavelengths/particle radii/refractive indices, adapted from CORAL (see https://github.com/mglodge/CORAL for more information and associated papers).

For this version, I have removed DDA calculations from CORAL. This version calculates Mean-field Theory (MFT), Modified Mean Field Theory (MMF) and Mie Theory only. 

NOTE: This version has a unique parameter input file (i.e. the regular CORAL parameter files will not work with this version, so be sure to use the example given in this directory if switching between codes!).

# Running the Code

1) Download all of the files into a directory of your choice -- these are already set up to run as demonstration. There are 3 key files in this folder:

    - CORAL_MMF_version.c (the main code which we will compile)
    - parameters.txt (sets the parameters for the scattering models)
    - refractive_index_data.txt (describes the wavelength dependence of refractive index)

2) To compile, simply navigate to the directory from within the terminal window, and type: 

         gcc -o CORAL CORAL_MMF_version.c -fopenmp -lm -O3

The flag -fopenmp enables OpenMP and parallel processing, which rapidly speed up the code, and -O3 gives a further increased speed boost. OpenMP is usually necessary 
to enable calculations at a reasonable speed, but users can remove the '#include omp.h' header from within CORAL.c, and remove the '-fopenmp' flag from above if they 
do not want to use it. Further details about OpenMP are [here.](https://www.openmp.org/)

3) To run the code, simply type:

       ./CORAL

# Parameters file

The parameters file details the range of wavelengths to study and the size of the particles etc. There are several options for the user to choose, and these are detailed below.

<b>IMPORTANT NOTE: the formatting has to be kept <i>exactly</i> the same as the original file i.e. if you change "aerosol_radius_min= 0.5" to "aerosol_radius_min=0.5" (if you delete the space after the equals sign), the code will not work. </b>

    aerosol_radius_min: set the minimum radius of aerosols to be examined, in um
    aerosol_radius_max: set the maximum radius of aerosols to be examined, in um
    n_radii: sets the number of radii to be analysed within the range above (max-min)
    radius_increment_type: set equal to 0 for linear-spaced incrememnts, or 1 for log-spaced increments
    requested_N: chooses a rough number of dipoles to create pseudospheres out of, if making our own spheres using shapeswitch=0 (rather than an input file)
    speed_test: Set to 0 to ignore, or set to 1 to record the times that various parts in the code reach. Used for debugging/speed benchmrking.
    
    /* wavelength range */
    wavelength_min: the minimum wavelength of radiation to be examined, in um
    wavelength_max: the maximum wavelength of radiation to be examined, in um
    num_wavelengths: sets the number of wavelengths to be analysed within the range above (max-min)
    wavelength_increment_type= 1  set equal to 0 for linear-spaced incrememnts, or 1 for log-spaced increments
    
    threadnumber: choose the number of threads to distribute the parallelised calculations across
    
    N_monomers: set the number of monomers for the shape that you are analysing (often you can judge this by eye using STAG)
    k_prefactor: set the fractal prefactor k0
    d_f: set the fractal dimension of the aggregate
    MMF_method: select the cut-off function for the structure factor calculation a value of 0 gives a gaussian cutoff and a value of 1 gives a fractal cutoff. We use the fractal cutoff as default to match OPTOOL. For more details, see Eq. (19) of Tazaki, R. and Tanaka, H., 2018. Light scattering by fractal dust aggregates. II. opacity and asymmetry parameter. The Astrophysical Journal, 860(1), p.79.

# Output files

If left at the default settings, CORAL will now calculate scattering properties for IR radiation of wavelength 0.7 um incident on a compact aggregate composed of 23 monomers and a fractal dimension of 2.25. The aggregate is assumed to have a spherical-equivalent-volume radius of 0.5 um. 

MAIN OUTPUT FILES: 

Each of the following files lists calculated values for $Q_{ext}$, $Q_{sca}$, $Q_{abs}$, and $g=<\cos(\theta)>$, as well as the associated cross-sections in $\mu m^2$, for each of the Mie, MMF and MFT models, at each particle radius and wavelength given.

   - mie_results.txt
   - mmf_results.txt
   - mft_results.txt


# Citations/Disclaimer

If CORAL is useful to your research, please cite Lodge et al. (2023). You are free to use and modify this code, but please attribute credit to the author (Matt Lodge) and also the authors below.

I would like to thank the authors of the following two codes, which were really well-written benchmark models. The detailed comments 
throughout were also indredibly helpful for debugging our own models, and we are grateful for their attention to detail:

    OPTOOL: Dominik, C., Min, M. & Tazaki, R. 2021, Optool, 1.9, Astrophysics Source Code Library, ascl:2104.010
