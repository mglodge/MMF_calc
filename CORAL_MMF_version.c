/*

01/12/23 - Intial public release. CORAL (Comparasion Of Radiative AnaLyses) is designed to compare three seperate optical scattering theories:

1) Mie scattering - assumes that particles are spherical
2) MMF (Modified Mean Field theory) - calculates the fractal properties of the input file, and uses these to calcualte scattering properties
3) DDA (Discrete Dipole Approximation) - uses the exact shape file to determine the scattering properties.

A full guide for use is found in the readme file.  All equations used in the code below are explained in Lodge, Wakeford and Leinhardt (2023), MNRAS (attached in git folder). We are incredibly 
grateful to Kouji Adachi, Peter Buseck and Serena Chung for allowing the use of their shape.dat files from:

    Adachi, K., Chung, S.H. and Buseck, P.R., 2010. Shapes of soot aerosol particles and implications for their effects on climate,
    Journal of Geophysical Research: Atmospheres, 115(D15).

If these shape files are used in your own research, please cite their paper. If CORAL is useful, please also cite Lodge et al. (2023).
You are free to use and modify this code, but please attitribute credit to the author (Matt Lodge) and also the authors below.

I would like to thank the authors of the following two codes, which were really well-written benchmark models. The detailed comments 
throughout were also incredibly helpful for debugging our own models, and we are grateful for their attention to detail:

    OPTOOL: Dominik, C., Min, M. & Tazaki, R. 2021, Optool, 1.9, Astrophysics Source Code Library, ascl:2104.010

    DDSCAT: Draine, B.T., & Flatau, P.J., "Discrete dipole approximation for scattering calculations", J. Opt. Soc. Am. A, 11, 1491-1499 (1994)

Although every care has been taken to benchmark the code and test it under a wide range of conditions (see attached paper for details), 
this code is provided without a warantee of any kind. The author cannot be held liable or accept any resposibility for any claim, loss 
or damages that arise in connection with the use of this code.

07/10/24 - Adapted to perform MMF, Mie and MFT calculations only. Note that the parameter input file is different to the regular version of CORAL (DDA inputs removed).

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

void Print_Matrix(double _Complex**A, double _Complex*E, int n){
    int i2, j2;
    printf("\n\n\n");
    for(i2=0;i2<n;i2++){        //Function to print matrix
        printf("\n |");
        for(j2=0;j2<n;j2++){
            if((j2>0)&&(((j2+1)%3)==0)){
                printf(" %5.1f + %5.1f i  | ", creal(A[i2][j2]), cimag(A[i2][j2])); //seperate the 3x3 elements
            }
            else{
                printf(" %5.1f + %5.1f i   ", creal(A[i2][j2]), cimag(A[i2][j2]));
            }
        }
        printf("  | P%d |   =   | %5.1f + %5.1f i |\n", (i2+1), creal(E[i2]), cimag(E[i2]));

        if((i2>0)&&(((i2+1)%3)==0)){
            for(j2=0;j2<(n/3);j2++){
                printf("       -----------------------------------------------    ");
            }
            printf("\n");
        }
    }
}
void Print_back_substitution(double _Complex**A, double _Complex*P, double _Complex*E, int n, int i){
    int i2, j2;
    printf("\n\n");
    for(i2=0;i2<n;i2++){        //Function to print matrix showing the back substitution process at each step
        printf("\n  |");
        for(j2=0;j2<n;j2++){
            printf(" %5.1f + %5.1f i ", creal(A[i2][j2]), cimag(A[i2][j2]));
        }
        if(i2 >= (n -(i+1))){ //only print values of P that we have calculated so far
            printf(" |    | %5.1f + %5.1f i |   =   | %5.1f + %5.1f i |\n", creal(P[i2]), cimag(P[i2]), creal(E[i2]), cimag(E[i2]));
        }
        else{
            printf(" |    | P%d |   =   | %5.1f + %5.1f i |\n", (i2+1), creal(E[i2]), cimag(E[i2]));
        }
    }
}
void Print_Final_Matrix(double _Complex*P, int n){
    int i2;
    for(i2=0;i2<n;i2++){        //Function to print final matrix
        printf("   |  %13g + %13g i  |\n", creal(P[i2]), cimag(P[i2]));
    }
        printf("\n\n");
}
double _Complex get_refractive_index(double wavelength, double _Complex**refractive_index_data, int number_ref_ind_data_points){
    double _Complex get_n=0, gradient_between_points, intercept;
    int i2;

    /* function to get wavelength-dependent refractive indices */

    if(wavelength < creal(refractive_index_data[0][0])){
        printf(" \n\n\n WARNING: The wavelength (%f um) is smaller than the smallest reference wavelength in our table of data (%f um). Values will be estimated. Continue?", wavelength, creal(refractive_index_data[0][0]));
        getchar();
    }

    for(i2=1;i2<number_ref_ind_data_points;i2++){
        if((wavelength < creal(refractive_index_data[i2][0]))||(i2==(number_ref_ind_data_points-1))){ // if the wavelength is smaller than the reference table wavelength, or if we are on the last loop and the wavelength is still larger than our biggest data point
            gradient_between_points = (refractive_index_data[i2][1]-refractive_index_data[i2-1][1])/(refractive_index_data[i2][0]-refractive_index_data[i2-1][0]); //find linear gradient between the points either side of the wavelength

            intercept = refractive_index_data[i2][1] - gradient_between_points*refractive_index_data[i2][0]; // c = y - mx for the second data point
            
            get_n = gradient_between_points*wavelength + intercept; // y = mx + c. if i==13 and we reach this point, wavelength must be larger than the max, and we just assume the same linear gradient as between the last two points.
            break; //exit loop
        }
    }
    
    return get_n;
}
double get_factorial(int i2){
    int i3;
    double factorial_value;

        //SIMPLE FUNCTION THAT RETURNS THE FACTORIAL: i2!

        //evaluate factorial
        factorial_value=1; // because zero factorial is 1
        for(i3=i2; i3>=1; i3--){
            factorial_value = i3*factorial_value;
            //printf("\n i3=%d, factorial = %f", i3, factorial_value);
        }
    return factorial_value;
}
double get_double_factorial(int dbl_fact_max){
    int i3;
    double double_factorial_value;

        // SIMPLE FUNCTION THAT RETURNS DOUBLE FACTORIALS x!! : NOT TO BE CONFUSED WITH EVALUATING THE FACTORIAL OF A FACTORIAL, DENOTED (x!)!
        // where x is variable 'dbl_fact_max'

        //evaluate double factorial (NOTE: THIS ALGORITHM ONLY WORKS FOR ODD INTEGERS, WHICH "2*i + 2*n + 1" ALWAYS IS).
        double_factorial_value=1; //initialise as 1 because it's a multiplicative iteration
        for(i3=2; i3<=dbl_fact_max; i3++){ //no point in evaluating i3=1; it just multiplies by 1x1! Start at i3=2.
            double_factorial_value = double_factorial_value*(2*i3 - 1);
            //printf("\n max= %d, i3=%d, double_factorial = %f", dbl_fact_max, i3, double_factorial_value);
        }
    return double_factorial_value;
}
int calculate_n_max(double x_size_parameter){
    int get_n_max;

    get_n_max = ceil(x_size_parameter + 4*pow(x_size_parameter, (1.0/3.0)) +2);
    
    return get_n_max;
}
double calculate_d_f(int N_monomers, double radius_of_gyration, double r_monomer){
    double best_solution, test_value, test_increment, test_max, min_answer_so_far, function_value;

    // function that finds numerical solution to the equation relating number of monomers to the radius of gyration, monmoer radius and fractal dimension. 
    // runs through values between the initial one set below (1.0) and the max (4.0) with increments of 0.001 and returns the test value that gives the closest solution to zero

    min_answer_so_far=100; //initialise this as a random high value
    test_value=1.0; //set initial test value to try
    test_increment=0.001; //set value to increment test by (we only need 3 sig figs as this function uses an approximation term for fractal prefactor k_0 in it's derivation anyway)
    test_max=4.0; //set max value to test

    //printf(" INPUT TO FUNCTION: \n\t N_monomers: %d \n\t radius of gyration = %f \n\t r_monomer = %f \n", N_monomers, radius_of_gyration, r_monomer);

    while(test_value<=test_max){

        function_value = (0.716 + sqrt(3) - 0.716*test_value)*pow(radius_of_gyration/r_monomer, test_value) - N_monomers;

        //printf("\n test = %f    Function = %f     Best so far = %f    Best x so far = %f ", test_value, function_value, min_answer_so_far, best_solution);

        if(sqrt(pow(function_value,2))<min_answer_so_far){ //if the magnitude of the answer to this equation is the closest to zero so far...
            best_solution = test_value; //...then the current value is the best solution so far - save it and try the next one
            min_answer_so_far = sqrt(pow(function_value,2)); // record the magnitude, because we are only interested in how "close to 0" it is
        }

        test_value += test_increment; //add increment and try the next value
    }

    return best_solution;
}
double calculate_d_f_Wolf_and_Toon_version(int N_monomers, double radius_of_gyration, double r_monomer){
    double best_solution, test_value, test_increment, test_max, min_answer_so_far, function_value;

    // function that finds numerical solution to the equation relating number of monomers to the radius of gyration, monmoer radius and fractal dimension. 
    // runs through values between the initial one set below (1.0) and the max (4.0) with increments of 0.001 and returns the test value that gives the closest solution to zero

    min_answer_so_far=100; //initialise this as a random high value
    test_value=1.0; //set initial test value to try
    test_increment=0.001; //set value to increment test by (we only need 3 sig figs as this function uses an approximation term for fractal prefactor k_0 in it's derivation anyway)
    test_max=4.0; //set max value to test

    //printf(" INPUT TO FUNCTION: \n\t N_monomers: %d \n\t radius of gyration = %f \n\t r_monomer = %f \n", N_monomers, radius_of_gyration, r_monomer);

    while(test_value<=test_max){

        function_value = pow(radius_of_gyration/r_monomer, test_value) - N_monomers;

        //printf("\n test = %f    Function = %f     Best so far = %f    Best x so far = %f ", test_value, function_value, min_answer_so_far, best_solution);

        if(sqrt(pow(function_value,2))<min_answer_so_far){ //if the magnitude of the answer to this equation is the closest to zero so far...
            best_solution = test_value; //...then the current value is the best solution so far - save it and try the next one
            min_answer_so_far = sqrt(pow(function_value,2)); // record the magnitude, because we are only interested in how "close to 0" it is
        }

        test_value += test_increment; //add increment and try the next value
    }

    return best_solution;
}
double _Complex Bessel_1st_kind(int bessel_order, double _Complex bessel_argument){
    int i2, i2_max, downwards_extra, dbl_fact_max; //local iterative variables
    double a1, a2, a3, a4, b1, b2, b3, b4, factorial, double_factorial, ln_x, ln_x_a, ln_x_b;
    double _Complex Bessel_1_value, Bessel_coeff_complex, series_total_complex, complex_sin_x, complex_cos_x;
    double _Complex* local_Bessel_1_x;
    double _Complex* temporary_F_complex;

    // This function calculates Spherical Bessel functions of the 1st kind of type j_n(x) in the form j_order(argument). It can deal with complex arguments.

    //define constant values for assessing which bessel method to use 
    a1 = -1.0*3.6693021*pow(10,-8);
    a2 = -1.0*3.1745158*pow(10,-5);
    a3 = 2.1567720*pow(10,-2);
    a4 = 9.9123380*pow(10,-1);
    b1 = -1.0*9.9921351*pow(10,-8);
    b2 = 8.7822303*pow(10,-6);
    b3 = 1.0238752*pow(10,-2);
    b4 = 3.7588265;
    i2_max=100; // limit set in Tazeki for series expansion method
    downwards_extra=100; // extra terms to calculate in Jablonski downwards recurrence method (set as 100 both in Jablonski + Tazeki)

    // Find the complex versions of sin x and cos x (there is no built in function for this)
    complex_sin_x= sin(creal(bessel_argument))*cosh(cimag(bessel_argument)) + cos(creal(bessel_argument))*sinh(cimag(bessel_argument))*I; // sin is complex because the refractive index is complex (m=a+bi) -- the general rule is sin(a+bi) = sin(a)cosh(b) + cos(a)sinh(b)i
    complex_cos_x= cos(creal(bessel_argument))*cosh(cimag(bessel_argument)) - sin(creal(bessel_argument))*sinh(cimag(bessel_argument))*I; // cos is complex because the refractive index is complex (m=a+bi) -- the general rule is cos(a+bi) = cos(a)cosh(b) + sin(a)sinh(b)i, and so multiplying by x, the rule is: cos(ax+bxi) = cos(ax)cosh(bx) + sin(ax)sinh(bx)i

    //printf("\n Complex sin x = %f + %f i", creal(complex_sin_x), cimag(complex_sin_x));

    /* calculate bessel functions for argument x */
    local_Bessel_1_x= (double _Complex*)malloc((bessel_order+1)*sizeof(double _Complex));  // initialise array to store Bessel 1st kind

    // calculate terms that determine which bessel method to use to find the functions of the 1st kind for x  - these don't need to be in the function but nice to keep them together
    ln_x = log(sqrt(pow(creal(bessel_argument),2) + pow(cimag(bessel_argument),2))); // calculate the log of the magnitude of x (which could be imaginary - this function can deal with complex or real numbers)
    ln_x_a = a1*pow(bessel_order,3) + a2*pow(bessel_order,2) + a3*bessel_order + a4;
    ln_x_b = b1*pow(bessel_order,3) + b2*pow(bessel_order,2) + b3*bessel_order + b4;

    //printf("\n\n ln_x = %g \n ln_x_a = %g \n ln_x_b = %g \n", ln_x, ln_x_a, ln_x_b);

        /* Choose a method based on values above and calculate Bessel function */

    if(ln_x<ln_x_a){
                        /* USE SERIES EXPANSION METHOD */

        //printf(" \n USING SERIES EXPANSION METHOD: ");

        /* set all terms to zero to initialise (important for this one, as we are only calculating the final term) */
        for(i2=0; i2<=bessel_order; i2++){
            local_Bessel_1_x[i2]= 0;
        }  

        //begin series calculation
        series_total_complex=0;
        for(i2=0; i2<=i2_max; i2++){

            //printf("\n\n i2= %d\n", i2);
            factorial = get_factorial(i2); //evaluate factorial
            dbl_fact_max=((2*bessel_order + 2*i2 + 1)+1)/2; //calculate argument for double factorial
            double_factorial=get_double_factorial(dbl_fact_max); //evaluate double factorial (NOTE: THIS ALGORITHM ONLY WORKS FOR ODD INTEGERS, WHICH "2*i + 2*n + 1" ALWAYS IS).
            
            series_total_complex += pow(-1.0,i2)*cpow((cpow(bessel_argument,2)/2),i2)/(factorial*double_factorial);
        }

        local_Bessel_1_x[bessel_order]= cpow(bessel_argument, bessel_order)*series_total_complex; //just evaluate final Bessel term for [i2max], we don't need the others! They are all set to 0 in case we call them/print them during debugging.
        
        //printf(" Spherical_Bessel_j_x(%d, %f + %f i) = %g + %g i", bessel_order, creal(bessel_argument), cimag(bessel_argument), creal(local_Bessel_1_x[bessel_order]), cimag(local_Bessel_1_x[bessel_order]));
        
        //printf("\n bessel_order= %d, Series total complex = %g + %g i", bessel_order, creal(series_total_complex), cimag(series_total_complex));
    }
    else if(ln_x<ln_x_b){
                        /* USE DOWNWARDS RECURRENCE METHOD */

        //printf(" \n USING DOWNWARDS RECURRENCE METHOD: ");

        temporary_F_complex= (double _Complex*)malloc((bessel_order+1+downwards_extra)*sizeof(double _Complex));  // initialise temporary array to store first guess at Bessel values
            
        temporary_F_complex[bessel_order+downwards_extra] = 0; //set final value to 0
        temporary_F_complex[bessel_order+downwards_extra - 1] = 1.0; //set penultimate value to 1

        //printf("\n Temporary_F_complex[%d] = %g + %g i", bessel_order+downwards_extra-1, creal(temporary_F_complex[bessel_order+downwards_extra-1]), cimag(temporary_F_complex[bessel_order+downwards_extra-1]));

        for(i2=(bessel_order+downwards_extra-2); i2>=0; i2--){
            temporary_F_complex[i2] = ((2.0*(i2+2.0-1.0)+1.0)/bessel_argument)*temporary_F_complex[i2+1] - temporary_F_complex[i2+2];
            //printf("\n Temporary_F_complex[%d] = %g + %g i", i2, creal(temporary_F_complex[i2]), cimag(temporary_F_complex[i2]));
            //printf(" coefficient = %g ", ((2.0*(i2+2.0-1.0)+1.0)/bessel_argument));
        }

        //find actual of value of j_0(x)

            /* IMPORTANT - condition to check whether the bessel argument is an integer multiple of pi - if so, element[0] is zero, so use element[1] to find the coefficient instead! Otherwise sin(x)=sin(pi)=0 and the coefficient will be 0 */
        if((creal(complex_sin_x) < 0.000000001) && (creal(complex_sin_x) > -0.000000001)){

            if(bessel_order==0){ //if we only need element zero, return zero (happens to the first term of the bessel function as x -> 0)!
                //we also need this 'if statement to prevent trying to access the second Bessel term (stored as element [1]) in the next 'else' statement when we haven't made an array big enough to store it!
                local_Bessel_1_x[0]= 0; //just set the function to zero
            }
            else{
                //if this is true, the argument is a multiple of pi and we will have to use element [1] to calculate the multiplying coefficient instead
                local_Bessel_1_x[1]= complex_sin_x/cpow(bessel_argument,2) - complex_cos_x/bessel_argument; //set element [1]

                //calculate coefficient of proportionality that links actual values to predicted ones
                Bessel_coeff_complex = local_Bessel_1_x[1]/temporary_F_complex[1]; //use element 1 to find the coefficient
            }
        }
        else{ //use the regular method (element[0])
            local_Bessel_1_x[0] = complex_sin_x / bessel_argument; 

            //calculate coefficient of proportionality that links actual values to predicted ones
            Bessel_coeff_complex = local_Bessel_1_x[0]/temporary_F_complex[0];
        }

        //printf("\n Actual B[0] = %g + %g i, and coeff = %g + %g i \n", creal(local_Bessel_1_x[0]), cimag(local_Bessel_1_x[0]), creal(Bessel_coeff_complex), cimag(Bessel_coeff_complex));

        if( ((creal(complex_sin_x) < 0.000000001) && (creal(complex_sin_x) > -0.000000001)) && (bessel_order==0) ){ //quick check in case we are looking for an order 0 term where the argument is a function of pi, in the region where downwards recurrence will fail (see note above in similar if statement)
            //do nothing - we have accounted for this case above - 'if(bessel_order==0)' - but we just need to make sure we don't try to access the Bessel_coeff in the statement below
        }
        else{
            //REGULAR METHOD - use the coefficient to correct all of the values up to value [bessel_order] (we could correct the others above this, but we don't need them!)
            for(i2=0; i2<=bessel_order; i2++){
                local_Bessel_1_x[i2] = Bessel_coeff_complex*temporary_F_complex[i2];
            }
        }
        //printf("\n Spherical_Bessel_j(%d, %g + %g i)= %g + %g i", bessel_order, creal(bessel_argument), cimag(bessel_argument), creal(local_Bessel_1_x[bessel_order]), cimag(local_Bessel_1_x[bessel_order]));

        free((void*)temporary_F_complex); //free temp array
    }
    else{
                        /* USE UPWARDS RECURRENCE METHOD */

        //printf(" \n USING UPWARDS RECURRENCE METHOD: ");

        /* Set initial value for element [0] */
        local_Bessel_1_x[0] = complex_sin_x / bessel_argument; //these functions are complex too

        /* if bessel_order>0, assign element [1] too */
        if(bessel_order>0){
            local_Bessel_1_x[1]= complex_sin_x/cpow(bessel_argument,2) - complex_cos_x/bessel_argument;
        }

        /* if bessel_order>1, calculate the rest of the Bessel function terms */
        for(i2=2; i2<=bessel_order; i2++){
            local_Bessel_1_x[i2]= ((2*(i2-1)+1)/bessel_argument)*local_Bessel_1_x[i2-1] - local_Bessel_1_x[i2-2];
        }  
        
        //printf("\n Spherical_Bessel_j_x(%d, %f + %f i) = %g + %g i", bessel_order, creal(bessel_argument), cimag(bessel_argument), creal(local_Bessel_1_x[bessel_order]), cimag(local_Bessel_1_x[bessel_order]));
    }

    Bessel_1_value = local_Bessel_1_x[bessel_order]; // assign value to return from function (be sure to do this before freeing the array!)

    //free array
    free((void*)local_Bessel_1_x);

    return Bessel_1_value;
}
double _Complex Bessel_2nd_kind(int bessel_order, double _Complex bessel_argument){
    int i2;
    double _Complex Bessel_2_value, complex_sin_x, complex_cos_x;
    double _Complex* local_Bessel_2_x;

    // This function calculates Spherical Bessel functions of the 2nd kind of type y_n(x) in the form y_order(argument). It can deal with complex arguments.

        complex_sin_x= sin(creal(bessel_argument))*cosh(cimag(bessel_argument)) + cos(creal(bessel_argument))*sinh(cimag(bessel_argument))*I; // sin is complex because the refractive index is complex (m=a+bi) -- the general rule is sin(a+bi) = sin(a)cosh(b) + cos(a)sinh(b)i
        complex_cos_x= cos(creal(bessel_argument))*cosh(cimag(bessel_argument)) - sin(creal(bessel_argument))*sinh(cimag(bessel_argument))*I; // cos is complex because the refractive index is complex (m=a+bi) -- the general rule is cos(a+bi) = cos(a)cosh(b) + sin(a)sinh(b)i, and so multiplying by x, the rule is: cos(ax+bxi) = cos(ax)cosh(bx) + sin(ax)sinh(bx)i

        local_Bessel_2_x= (double _Complex*)malloc((bessel_order+1)*sizeof(double _Complex));  // initialise array to store Bessel 1st kind

        /* Set initial value for Bessel function of 2nd kind */
        local_Bessel_2_x[0]= -complex_cos_x / bessel_argument; 

        /* if bessel_order>0, assign 2nd elements too */
        if(bessel_order>0){
            local_Bessel_2_x[1]= -complex_cos_x/cpow(bessel_argument,2) - complex_sin_x/bessel_argument;
        }

        /* if bessel_order>1, calculate the rest of the Bessel function terms using upwards recurrence (no other methods needed for 2nd kind functions, they are stable!)*/
        for(i2=2; i2<=bessel_order; i2++){ 
            local_Bessel_2_x[i2]=((2*(i2-1)+1)/bessel_argument)*local_Bessel_2_x[i2-1] - local_Bessel_2_x[i2-2]; // Find Bessel_2 terms using upwards recurrence
        }  

        Bessel_2_value = local_Bessel_2_x[bessel_order]; // assign value to return from function (be sure to do this before freeing the array!)

        //free array
        free((void*)local_Bessel_2_x);

    return Bessel_2_value;
}
double _Complex Hankel_x(double _Complex bessel_argument, double _Complex bessel_argument_2){
    double _Complex hankel_value;
    // This function computes the spherical Hankel function, given two spherical Bessel functions (1st and 2nd kind) as inputs (in that order respectively: argument 1, argument 2
    hankel_value = bessel_argument + bessel_argument_2*I;     // Hankel function = j(x) + y(x)*i
    
    return hankel_value;
}
void get_a_n_and_b_n(int i, double x_size_parameter, double _Complex ref_ind_sphere, double magnetic_perm_sphere, double magnetic_perm_medium, double _Complex *a_n, double _Complex *b_n) {
    double diff_xjx;
    double _Complex diff_mxjmx, diff_xhx, Hankel, Hankel_minus_1;

    // This function finds the value of a_n and b_n coefficents. It requires calculating the bessel functions from scratch for each value of 'n', which
    // is inefficient, but it's nicer code-wise to be able to calculate this in a standard way for both the Mie and MMF sections, and it is not a computational
    // bottleneck (even calculating the bessel functions from element [0] each loop, there are only a few (probably <100 loops to do max) -- doesn't take long!)

    /* Calculate coefficients a and b */

    diff_xjx= x_size_parameter*creal(Bessel_1st_kind(i-1, x_size_parameter)) - i*creal(Bessel_1st_kind(i, x_size_parameter)); //just use the real parts as x is always real
    diff_mxjmx= ref_ind_sphere*x_size_parameter*Bessel_1st_kind(i-1, ref_ind_sphere*x_size_parameter) - i*Bessel_1st_kind(i, ref_ind_sphere*x_size_parameter);
    
    Hankel = Hankel_x(Bessel_1st_kind(i, x_size_parameter), Bessel_2nd_kind(i, x_size_parameter)); //combine Bessel 1+2 to create Hankel function
    Hankel_minus_1= Hankel_x(Bessel_1st_kind(i-1, x_size_parameter), Bessel_2nd_kind(i-1, x_size_parameter)); //same as above but for [i-1]th element, which is needed to calculate diff_xhx
    diff_xhx = x_size_parameter*Hankel_minus_1 - i*Hankel;

    //complex_temp_2= x_size_parameter*ref_ind_sphere;

    //printf("\n\n DEBUG: \n i = %d \n x = %f \n mx = %g + %g \n Bessel_1st_kind[%d] = %g + %g i\n Bessel_1st_kind[%d] = %g + %g i\n Bessel_2nd_kind[%d] = %g + %g i\n Hankel[%d] = %g + %g i\n Hankel[%d] = %g + %g i\n\n", i, x_size_parameter, creal(complex_temp_2), cimag(complex_temp_2), i-1, creal(Bessel_1st_kind(i-1, x_size_parameter)), cimag(Bessel_1st_kind(i-1, x_size_parameter)), i, creal(Bessel_1st_kind(i, x_size_parameter)), cimag(Bessel_1st_kind(i, x_size_parameter)), i, creal(Bessel_2nd_kind(i, x_size_parameter)), cimag(Bessel_2nd_kind(i, x_size_parameter)), i, creal(Hankel), cimag(Hankel), i-1, creal(Hankel_minus_1), cimag(Hankel_minus_1));
    //printf("\n\t diff_xjx = %14f \n\t diff_mxjmx = %14f + %14g i \n\t diff_xhx = %14f + %14g i \n", creal(diff_mxjmx), cimag(diff_mxjmx), creal(diff_xhx), cimag(diff_xhx));

    *a_n= (magnetic_perm_medium*cpow(ref_ind_sphere,2)*Bessel_1st_kind(i, ref_ind_sphere*x_size_parameter)*diff_xjx - magnetic_perm_sphere*Bessel_1st_kind(i, x_size_parameter)*diff_mxjmx) / (magnetic_perm_medium*cpow(ref_ind_sphere,2)*Bessel_1st_kind(i, ref_ind_sphere*x_size_parameter)*diff_xhx - magnetic_perm_sphere*Hankel*diff_mxjmx);
    *b_n= (magnetic_perm_sphere*Bessel_1st_kind(i, ref_ind_sphere*x_size_parameter)*diff_xjx - magnetic_perm_medium*Bessel_1st_kind(i, x_size_parameter)*diff_mxjmx) / (magnetic_perm_sphere*Bessel_1st_kind(i, ref_ind_sphere*x_size_parameter)*diff_xhx - magnetic_perm_medium*Hankel*diff_mxjmx);
}
double sigma_integrand(double x_val, double d_f){
    return pow(x_val, -2.0/d_f)*exp(-x_val); //function to calculate the integrand of the sigma term
}
double calculate_geometrical_cross_section(int N_monomers, double r_monomer, double d_f, double k_frac_prefactor, double pi){
    int num_steps;
    double MMF_eta, MMF_eta_th, N_th, get_G, MMF_A, sigma, integrated_sigma, sigma_th, integrated_sigma_th, x_min, x_max, dx; 

    // code to calculate the geometrical cross-section of aggregates, using the method in Tazeki (2021), for use in MMF calculations

    MMF_eta = pow(2.0, d_f-1.0)*k_frac_prefactor/N_monomers;

    N_th= 11*d_f - 8.5;
    if(N_th>8.0){
        N_th= 8.0;   //calculate N_th (should always be <= 8.0)
    }

    MMF_eta_th= pow(2.0, d_f-1.0)*k_frac_prefactor/N_th; //not in the paper, but this is used to find sigma_th

    //printf("\n MMT_eta = %f \n MMT_eta_th = %f \n N_th = %f \n", MMF_eta, MMF_eta_th, N_th);

    // use Simpson's rule to calculate sigma_th (used to find the numerical value that joins the two regions, A)

    x_min = MMF_eta_th;
    x_max = 6; //plotting this function in wolfram, it looks like it is pretty much zero once x > 6 even for the largest fractal dimension value of 3 (so no need to go to infinity!)
    num_steps = 100000; //same value used in optools_fractal.f90
    dx=(x_max-x_min)/num_steps; //calculate step size
    //printf("\n num steps = %d  \n x_min = %f    \n x_max = %f    \n dx = %f   ",  num_steps, x_min, x_max, dx);

    integrated_sigma_th=0; //initialise integral sum

    //begin parallel calculation to integrate sigma_th
    #pragma omp parallel for reduction (+:integrated_sigma_th)  //use "reduction" to ensure that each thread does not access integrated_S_q at the same time
    for(int z=0; z<num_steps; z++){

        double step_a= x_min + z*dx; //set lower integral limit (this will be lowest possible value of integral (0) PLUS "number of steps X step size")
        double step_b= step_a + dx; //set upper integral limit     // NOTE: DEFINE ALL THESE VARIABLES WITHIN PARALLEL SECTION TO KEEP THESE PRIVATE
        double step_c= (step_a+step_b)/2.0; //set middle integral term for Simpson's rule (we integrate between a and b above, but this term is useful)

        //printf("\n\n\n We integrate between %f --> %f  (with a midpoint of %f) and where dx= %f \n", step_a, step_b, step_c, dx);
        //printf("\n f(a) = %f      f(c) = %f        f(b) = %f ", sigma_integrand(step_a, d_f), sigma_integrand(step_c, d_f), sigma_integrand(step_b, d_f));

        //find area underneath curve between points a and b
        double step_area = (dx/6.0) * (sigma_integrand(step_a, d_f) + 4.0*sigma_integrand(step_c, d_f) + sigma_integrand(step_b, d_f)); // Use Simpson's rule to estimate the area under the step between a --> b. And yep: the order is a, c, b!

        integrated_sigma_th += step_area; //add step area to total area under curve

        //printf("\n step area = %f    integrated_sigma_th = %f ", step_area, integrated_sigma_th);
    }
    #pragma omp barrier

    sigma_th = (pow(MMF_eta_th,2.0/d_f)/16.0) * integrated_sigma_th;
    //printf("\n integrated_sigma_th = %f      sigma_th = %f", integrated_sigma_th, sigma_th);

    MMF_A = 12.5*pow(N_th, -0.315)*exp(-2.53/pow(N_th, 0.092))*( 1 + (N_th-1)*sigma_th);

    MMF_A = 1; //OVERRIDE METHOD ABOVE AND ASSUME THAT A =1 (using the above method results in the wrong value for the short wavelength limit; private comm. with Ryo Tazaki)

    //printf("\n A = %f", MMF_A);

        // use Simpson's rule to calculate sigma 

    x_min = MMF_eta;
    x_max = 6; //plotting this function in wolfram, it looks like it is pretty much zero once x > 6 even for the largest fractal dimension value of 3 (so no need to go to infinity!)
    num_steps = 100000; //same value used in optools_fractal.f90
    dx=(x_max-x_min)/num_steps; //calculate step size
    //printf("\n num steps = %d  \n x_min = %f    \n x_max = %f    \n dx = %f   ",  num_steps, x_min, x_max, dx);

    integrated_sigma=0; //initialise integral sum

    //begin parallel calculation to integrate sigma
    #pragma omp parallel for reduction (+:integrated_sigma)  //use "reduction" to ensure that each thread does not access integrated_S_q at the same time
    for(int z=0; z<num_steps; z++){

        double step_a= x_min + z*dx; //set lower integral limit (this will be lowest possible value of integral (0) PLUS "number of steps X step size")
        double step_b= step_a + dx; //set upper integral limit     // NOTE: DEFINE ALL THESE VARIABLES WITHIN PARALLEL SECTION TO KEEP THESE PRIVATE
        double step_c= (step_a+step_b)/2.0; //set middle integral term for Simpson's rule (we integrate between a and b above, but this term is useful)

        //printf("\n\n\n We integrate between %f --> %f  (with a midpoint of %f) and where dx= %f \n", step_a, step_b, step_c, dx);
        //printf("\n f(a) = %f      f(c) = %f        f(b) = %f ", sigma_integrand(step_a, d_f), sigma_integrand(step_c, d_f), sigma_integrand(step_b, d_f));

        //find area underneath curve between points a and b
        double step_area = (dx/6.0) * (sigma_integrand(step_a, d_f) + 4.0*sigma_integrand(step_c, d_f) + sigma_integrand(step_b, d_f)); // Use Simpson's rule to estimate the area under the step between a --> b. And yep: the order is a, c, b!

        integrated_sigma += step_area; //add step area to total area under curve

        //printf("\n step area = %f    integrated_sigma = %f ", step_area, integrated_sigma);
    }
    #pragma omp barrier

    sigma = (pow(MMF_eta,2.0/d_f)/16.0) * integrated_sigma;
    //printf("\n integrated_sigma = %f      sigma = %f", integrated_sigma, sigma);

    if(N_monomers<N_th){
        get_G = N_monomers*pi*pow(r_monomer,2.0) * 12.5*pow(N_monomers, -0.315)*exp(-2.53/pow(N_monomers, 0.092));
        //printf("\n N_monomers = %d and N_th = %f so N_monomers<N_th.       We calculate G_MMF = %f       (normalised value: %f).", N_monomers, N_th, get_G, get_G/(N_monomers*pi*pow(r_monomer,2.0)));
    }
    else{
        get_G = N_monomers*pi*pow(r_monomer,2.0) * (MMF_A / (1 + (N_monomers-1)*sigma));
        //printf("\n N_monomers = %d and N_th = %f so N_monomers>N_th.       We calculate G_MMF = %f.       (normalised value: %f).", N_monomers, N_th, get_G, get_G/(N_monomers*pi*pow(r_monomer,2.0)));
    }

    return get_G;
}
double _Complex S_integrand(double u, double d_f, double x_agg, int p_MMF, int MMF_method){
    double f_c, f_c_constant;
    double _Complex get_S_integrand;

    //printf("\n\t u = %f",u); //print the value of u that the integrand is being found for

    //safety check - if u=0, the first term of S_p, which is u^(d_f-1), becomes zero. Therefore set this and skip the rest of the code (because you also can't calculate the terms of the 2nd Bessel function for u=0)
    if(u<0.000001){
        get_S_integrand=0;
    }
    else{

        /* select appropriate cut-off function f_c (u/x_g) */
        if(MMF_method==0){
            f_c_constant = d_f/4.0; //assign constant
            f_c = ((2.0 * pow(f_c_constant, (d_f/2.0)))/tgamma(d_f/2.0)) * exp(-1.0 * f_c_constant * pow((u/x_agg),2.0));
        }
        else{
            f_c_constant = 0.5; //assign constant
            f_c = f_c_constant*d_f*exp(-1.0 * f_c_constant * pow((u/x_agg),d_f));
        }
        //printf("\n\t\t f_c = %f", f_c);

        //calculate integrand and return it
        get_S_integrand = pow(u, (d_f-1)) * Bessel_1st_kind(p_MMF,u) * Hankel_x(Bessel_1st_kind(p_MMF,u), Bessel_2nd_kind(p_MMF,u)) * f_c; // calculate integrand for this value of 'u' (remember that the Hankel function is called in a slightly different way -- it required Bessel functions as the arguments, not the (order, argument) format that the Bessel functions themselves use)
    }
    return get_S_integrand;
}
double LegendrePoly(int order, double x_LP){
    int i2;
    double Legendre_value;
    double* LegendreP;

    LegendreP = (double*)malloc((order+1)*sizeof(double));  // array to store Legendre Polynomials
    
    // assign initial Legendre Polynomial terms for this value of x
    LegendreP[0] = 1.0;

    /* if p>0, assign 2nd Polynomial too */
    if(order>0){
        LegendreP[1]=x_LP;
    }
    /* if order>1, calculate the rest of the terms */
    for(i2=2; i2<=order; i2++){  //the regular Legendre terms only ever need to calculated for variable p, so only iterate up to the p_MFF term that we need
        LegendreP[i2] = ((2*(i2-1)+1)*x_LP*LegendreP[i2-1] - (i2-1)*LegendreP[i2-2])/i2;
    }
    Legendre_value = LegendreP[order]; //save final value to return from function

    free((void*)LegendreP); //free array
    
    return Legendre_value;
}
double diff_LegendrePoly(int order, double x_LP){
    double diff_Legendre_value;

    if(order==0){
        diff_Legendre_value=0; //because we shouldn't submit 'order-1' as an argument to LegendrePoly function
    }
    else{
        diff_Legendre_value = (order/(pow(x_LP,2)-1)) * (x_LP*LegendrePoly(order, x_LP) - LegendrePoly(order-1, x_LP));
    }
    return diff_Legendre_value;
}
double Associated_LegendrePoly(int order, double x_LP){
    int i2;
    double Associated_Legendre_value;
    double* Associated_LegendreP;

    Associated_LegendreP = (double*)malloc((order+1)*sizeof(double));  // array to store Associated Legendre Polynomials

    // calculate initial Associated Legendre Polynomials for this value of x
    Associated_LegendreP[0]= 0; // this value is actually undefined, but never used anyway (we only calculate the associated Legendre terms for v and n, which are always >= 1 in the summation terms. So we never need to calculate term 0)
    Associated_LegendreP[1]= -1.0*sqrt(1.0-pow(x_LP,2)); // this is always initialised as largest_n_v>=1 because n and v are always >=1
    
    /* if largest_n_v > 1, initialise second term */
    if(order>1){
        Associated_LegendreP[2]= -1.0*sqrt(1.0-pow(x_LP,2))*3.0*x_LP;
    }

    /* if largest_n_v > 2, calculate all other terms too */
    for(i2=3; i2<=order;i2++){
        Associated_LegendreP[i2] = (x_LP*(2.0*i2 - 1)*Associated_LegendreP[i2-1] - i2*Associated_LegendreP[i2-2])/(i2-1);
    }

    Associated_Legendre_value = Associated_LegendreP[order];
        
    free((void*)Associated_LegendreP);  // free arrays
    
    return Associated_Legendre_value;
}
void calculate_GL_roots(int num_points_quadrature, double* GL_roots){
    int i2, k2, midpoint; 
    double x_k, x_k_new;
    double sum_term, precision=0.00000000001;

    /* function that uses Newton-Rhaphson to calculate the roots of Legendre polynomials (to use in Gauss-Legendre quadrature) */
    
    //find midway point
    if((num_points_quadrature%2)==0){
        midpoint=num_points_quadrature/2; //num points is even
    }
    else{
        midpoint=(num_points_quadrature-1)/2; //num points is odd
    }

    for(k2=0; k2<midpoint; k2++){
        x_k=1.0001; //initial guess at root -- this needs to be chosen carefully to make the algorithm numerically stable. x_k = 0 or 1 do NOT work and cause overflows, however x_k = 1.0001 seems a good choice for num_points=400.

        while(1==1){ // loop until we achieve the desired precision
            
            sum_term=0;
            for(i2=0; i2<k2; i2++){
                sum_term = sum_term + (1.0/(x_k-GL_roots[i2]));  //divide function by previously found roots to find the next one
            }

            x_k_new = x_k - ( LegendrePoly(num_points_quadrature, x_k)/(diff_LegendrePoly(num_points_quadrature, x_k) - LegendrePoly(num_points_quadrature, x_k)*sum_term) ); //find next guess
            
            //printf("\n Root %d: Initial guess = %.10f    Updated guess = %.10f", k2, x_k, x_k_new);

            if(sqrt(pow((x_k_new-x_k),2)) < precision){ //compare previous and new guess
                GL_roots[k2] = x_k_new; // if achieved the required precision, store value and exit
                //printf("   ROOT FOUND!");
                break;    
            }

            x_k = x_k_new; // otherwise, set new x_k value as the next x_k and repeat
        }
    }

    //once we find all the positive roots, use symmetry to find the rest of the negative points
    
    i2=0;
    if((num_points_quadrature%2)==0){ //if num points even
        for(k2=(num_points_quadrature-1); k2>=midpoint; k2--){
            GL_roots[k2] = -GL_roots[i2];
            i2++;
        }
    }
    else{ //if num points odd, don't include the midpoint in the for loop (work that out seperately!)
        for(k2=(num_points_quadrature-1); k2>midpoint; k2--){
            GL_roots[k2] = -GL_roots[i2];
            i2++;
        }
        GL_roots[midpoint]=0; //middle root is always 0 for an odd number of points
    }
}
void calculate_GL_weights(int num_points_quadrature, double* GL_roots, double* GL_weights){
    int i2, k2, midpoint;
    
    /* function that calculates the respective weights of each root of the Legendre polynomials (to use in Gauss-Legendre quadrature) */


    //find midway point
    if((num_points_quadrature%2)==0){
        midpoint=num_points_quadrature/2; //num points is even
    }
    else{
        midpoint=(num_points_quadrature-1)/2; //num points is odd
    }

    for(i2=0; i2<midpoint; i2++){
        GL_weights[i2] = 2.0 / ( (1.0-pow(GL_roots[i2],2)) * pow(diff_LegendrePoly(num_points_quadrature,GL_roots[i2]),2) );
    }

    //once we find all the positive roots, use symmetry to find the rest of the negative points
    i2=0;
    if((num_points_quadrature%2)==0){ ///if num points even
        for(k2=(num_points_quadrature-1); k2>=midpoint; k2--){
            GL_weights[k2] = GL_weights[i2];
            i2++;
        }
    }
    else{ //if num points odd, don't include the midpoint in the for loop (work that out seperately!) 
        for(k2=(num_points_quadrature-1); k2>midpoint; k2--){
            GL_weights[k2] = GL_weights[i2];
            i2++;
        }
        GL_weights[midpoint] = 2.0 / pow(diff_LegendrePoly(num_points_quadrature,0),2); // calculate weight of midpoint (at the midpoints for odd Legendre Polynomials, the root x_i = 0 ,so the weight equation simplifies to this)
    }
}
double cmodulus(double _Complex a_n){    //computes the complex modulus (absolute magnitude) of a complex number
    return sqrt( pow(creal(a_n),2) + pow(cimag(a_n),2) );
}
double _Complex cconjugate(double _Complex a_n){    //computes the complex conjugate of a complex number
    return creal(a_n) - cimag(a_n)*I;
}
int is_odd(int p_MMF){
    return p_MMF%2;
}
void calculate_pi_and_tau_function_terms(int n_max, double scatter_angle, double* pi_function, double* tau_function){
    int i2;

    // function to calculate pi and tau functions for a particular scattering angle and store them in an array

    pi_function[0]= 0;
    pi_function[1]= 1;

    tau_function[0]= 0;
    tau_function[1]= cos(scatter_angle);

    for(i2=2;i2<=n_max;i2++){
        // Calculate new function terms (only after the first loop -- we already have terms 0 and 1 from initialisation above
        pi_function[i2]= ((2*i2-1.0)/(i2-1.0))*cos(scatter_angle)*pi_function[i2-1] - (i2/(i2-1.0))*pi_function[i2-2];
        tau_function[i2]= i2*cos(scatter_angle)*pi_function[i2] - (i2+1)*pi_function[i2-1];
    }
}
double _Complex S1(int n_max, double _Complex* d_1_terms, double _Complex* d_2_terms, double* pi_function, double* tau_function){
    int i2;
    double _Complex get_S1;

    // function to calculate S1 term

    get_S1=0;
    for(i2=1; i2<=n_max; i2++){
        get_S1 += ((2.0*i2 + 1.0)/(i2*(i2+1.0))) * (d_1_terms[i2]*pi_function[i2] + d_2_terms[i2]*tau_function[i2]);
        //printf("\n i2= %d   S1 = %f + %f i", i2, creal(get_S1), cimag(get_S1));
        //printf("\n fn= %f   d1 = %f + %f i    d2 = %f + %f i", ((2.0*i2 + 1.0)/(i2*(i2+1.0))),  creal(d_1_terms[i2]), cimag(d_1_terms[i2]),  creal(d_2_terms[i2]), cimag(d_2_terms[i2]));
        //printf("\n pi = %f    tau = %f ", pi_function[i2], tau_function[i2]);
    }

    return get_S1;
}
double _Complex S2(int n_max, double _Complex* d_1_terms, double _Complex* d_2_terms, double* pi_function, double* tau_function){
    int i2;
    double _Complex get_S2;

    // function to calculate S2 term

    get_S2=0;
    for(i2=1; i2<=n_max; i2++){
        get_S2 += ((2.0*i2 + 1.0)/(i2*(i2+1.0))) * (d_1_terms[i2]*tau_function[i2] + d_2_terms[i2]*pi_function[i2]);
    }

    return get_S2;
}
double calculate_S11(int n_max, double _Complex* d_1_terms, double _Complex* d_2_terms, double* pi_function, double* tau_function){
    double find_S11;

    // function to calculate S11 term = 0.5 * ( |S1|^2 + |S2|^2 )

    find_S11 = 0.5*( cpow(cmodulus(S1(n_max, d_1_terms, d_2_terms, pi_function, tau_function)),2) + cpow(cmodulus(S2(n_max, d_1_terms, d_2_terms, pi_function, tau_function)),2) );

    return find_S11;
}
double S_q_integrand(double x_val, double d_f, double q_coeff, double radius_of_gyration){
    
    if(x_val<0.00004){
        return 0; //condition just for the first call of the function -- for some reason, it returns nan...
    }
    else{
        return pow(x_val, d_f-2.0)*sin(q_coeff*radius_of_gyration*x_val)*exp(-0.5*pow(x_val, d_f));
    }
}
double calculate_structure_factor(double d_f, double q_coeff, double radius_of_gyration){
    int num_steps;
    double find_S_q, x_max, integrated_S_q, dx;


    
    /* use simpson's rule to integrate structure factor equation */

    if(q_coeff<0.000001){
        find_S_q=1;     // if q is basically zero, don't bother with the integral and set this value equal to 1 (line 1116 in optools_fractal.f90)
    }
    else{

        // perform integration as normal -- start by finding limits to integral
        
        x_max = pow((25.0/0.5),(1.0/d_f)); // same as in Tazaki + Tanaka
        num_steps = 100000; //same value used in optools_fractal.f90
        dx=x_max/num_steps; //calculate step size
        //printf("\n num steps = %d  \n x_min = 0    \n x_max = %f    \n dx = %f   ",  num_steps, x_max, dx);

        integrated_S_q=0; //initialise integral sum

        //begin parallel calculation to integrate structure factor function
        //#pragma omp parallel for reduction (+:integrated_S_q)  //use "reduction" to ensure that each thread does not access integrated_S_q at the same time
        for(int z=0; z<num_steps; z++){

            double step_a= z*dx; //set lower integral limit (this will be lowest possible value of integral (0) PLUS "number of steps X step size")
            double step_b= step_a + dx; //set upper integral limit     // NOTE: DEFINE ALL THESE VARIABLES WITHIN PARALLEL SECTION TO KEEP THESE PRIVATE
            double step_c= (step_a+step_b)/2.0; //set middle integral term for Simpson's rule (we integrate between a and b above, but this term is useful)

            //printf("\n\n\n We integrate between %f --> %f  (with a midpoint of %f) and where dx= %f \n", step_a, step_b, step_c, dx);
            //printf("\n f(a) = %f      f(c) = %f        f(b) = %f ", S_q_integrand(step_a, d_f, q_coeff, radius_of_gyration), S_q_integrand(step_c, d_f, q_coeff, radius_of_gyration), S_q_integrand(step_b, d_f, q_coeff, radius_of_gyration));

            //find area underneath curve between points a and b
            double step_area = (dx/6.0) * (S_q_integrand(step_a, d_f, q_coeff, radius_of_gyration) + 4.0*S_q_integrand(step_c, d_f, q_coeff, radius_of_gyration) + S_q_integrand(step_b, d_f, q_coeff, radius_of_gyration)); // Use Simpson's rule to estimate the area under the step between a --> b. And yep: the order is a, c, b!

            integrated_S_q += step_area; //add step area to total area under curve

            //printf("\n step area = %f    integrated_S_q = %f ", step_area, integrated_S_q);
        }
        //#pragma omp barrier

        //use the intregal to calculate S(q)
        find_S_q= (0.5*d_f/(q_coeff*radius_of_gyration)) * integrated_S_q;

        //SAFETY CHECK - S(q) should always be >=0
        if(find_S_q<0){
            //printf("\n WARNING: find_S_q < 0. Assuming it should be = 0 instead (this can happen for d_f > 2, and is ok...)");
            find_S_q=0; //if the monomers are far enough apart that they don't affect each other, setting find_S_q = 0 should find a case for two non-interacting spheres. Tazeki made the same note in optools_fractal.f90.
        }
    }

    return find_S_q;
}
double g_integrand(double* phase_function, double pi, int k){
    // function to calculate the integrand of g_asymmetry
    double g_angle;
        g_angle= k*pi/180.0; //covert k from degrees to radians
    return phase_function[k]*sin(g_angle)*cos(g_angle); 
}
double find_g_MMF(double* phase_function, double pi){
    int d_angle;
    double find_g, integrated_g;
    
    // perform integration as normal -- start by finding limits to integral
        
        d_angle= 2; //set step size (d_theta)
        //printf("\n num steps = %d  \n x_min = 0    \n x_max = %f    \n dx = %f   ",  num_steps, x_max, dx);

        integrated_g=0; //initialise integral sum

        //begin parallel calculation to integrate structure factor function
        #pragma omp parallel for reduction (+:integrated_g)  //use "reduction" to ensure that each thread does not access integrated_g at the same time
        for(int z=0; z<90; z++){ //this will loop through array values: phase_function[0->180]

            int int_a= z*d_angle; //set lower integral limit (this will be lowest possible value of integral (0) PLUS "number of steps X step size")
            int int_b= int_a + 2; //set upper integral limit     // NOTE: DEFINE ALL THESE VARIABLES WITHIN PARALLEL SECTION TO KEEP THESE PRIVATE
            int int_c= int_a + 1; //set middle integral term for Simpson's rule (we integrate between a and b above, but this term is useful)

            //printf("\n\n For z = %d, we integrate between %d --> %d  (with a midpoint of %d) and where d_angle= %d \n", z, int_a, int_b, int_c, d_angle);
            //printf("\n f(a) = %f      f(c) = %f        f(b) = %f ", g_integrand(phase_function, pi, int_a), g_integrand(phase_function, pi, int_c), g_integrand(phase_function, pi, int_b));

            //find area underneath curve between points a and b - REMEMBER THAT d_angle NEEDS TO BE CONVERTED FROM DEGREES INTO RADIANS HERE!!
            double step_area = (d_angle*(pi/180.0)/6.0) * (g_integrand(phase_function, pi, int_a) + 4.0*g_integrand(phase_function, pi, int_c) + g_integrand(phase_function, pi, int_b)); // Use Simpson's rule to estimate the area under the step between a --> b. And yep: the order is a, c, b!

            integrated_g += step_area; //add step area to total area under curve

            //printf("\n step area = %f    integrated_g = %f ", step_area, integrated_g);
        }
        #pragma omp barrier

        find_g= 2*pi*integrated_g; //finally, calculate g

    return find_g;
}
double maximum_C_abs(double MFT_C_abs, double G_MMF, double tau_MMF){
    double G_tau_term;

    G_tau_term= G_MMF*(1.0-exp(-1.0*tau_MMF));
    printf("\n G_tau_term = %f", G_tau_term);

    printf("\n Which is bigger?      Modified C_abs (G_tau_term) = %f    or     MFT_C_abs (found using MMF_C_Ext - MFT_C_scatter) = %f", G_tau_term, MFT_C_abs);

    if(MFT_C_abs>G_tau_term){
        return MFT_C_abs;
    }
    else{
        return G_tau_term;
    }
}
FILE* matrixfile;
FILE* parameterfile;
FILE* LDRcriteriafile;
FILE* ddscatfile;
FILE* dipolefile;
FILE* shapedatafile;
FILE* qtablefile;
FILE* results_angles_file_y;
FILE* results_angles_file_z;
FILE* results_averages_file_y;
FILE* results_averages_file_z;
FILE* mieresultsfile;
FILE* mmfresultsfile;
FILE* mftresultsfile;
FILE* wavelengthfile;
FILE* anglesfile;
FILE* timefile;
FILE* refractiveindexfile;
FILE* LDRparameterspacefile;
FILE* polarisationresultsfile;
FILE* refrindfile;
FILE* mie_sphere_mieff_file;
FILE* DDA_mieff_file;
int requested_N, lattice_N, lattice_dim, lattice_radius, x_pos, y_pos, z_pos, sphere_N, N, n, i, j, v, w, k, highest_k=0, b, number_of_spheres, STAG_lattice_dim, temp_int, dipole_number, IX, IY, IZ, ICOMP[3], shapeswitch, marbleswitch, num_wavelengths, N_theta, N_phi, radius_index, wavelength_index, angle_index, N_x, N_y, N_z, e, f, g, orientational_average_resolution, num_vertices, new_vertex, vertex_1, vertex_2, vertex_3, triangular_faces[20][3], edges[30][3], loop_number, speed_test, number_of_threads, wavelength_increment_type, number_ref_ind_data_points, polarisation, num_radii, radius_increment_type;
double _Complex ref_ind_medium, ref_ind_sphere, epsilon_sphere, alpha_CM, alpha_denominator_1, alpha_denominator_2, alpha_denominator_3, alpha_CLDR[3][3], inverse_alpha_CLDR[3][3], complex_temp, conj_inverse_alpha;
double wavelength, wavenumber_k, magnetic_perm_medium, magnetic_perm_sphere, relative_permeability_sphere, accuracy, highest_wavelength, LDR_criteria, LDR_criteria_failpoint, LDR_wavelength_increment, LDR_R_mie, lowest_wavelength_so_far;
double aerosol_radius_min, aerosol_radius_max, aerosol_radius, aerosol_radius_increment, d_to_point, DDA_position_offset_X0, total, previous_total, previous_previous_total, total2;
double d, b_1, b_2, b_3, incident_light[3], initial_polarisation[3], rotated_polarisation[3], DDA_Q_Ext, DDA_Q_Abs, DDA_Q_Scatter, DDA_C_Ext, DDA_C_Abs, DDA_C_Scatter, abs_sum;
double sphere_locations[100][3], min[3], max[3], STAG_odd_even_offset, STAG_offset[3], sphere_odd_even_offset, temp_real, center_of_mass[3], dipole_mass[10], total_mass, radius_to_dipole, radius_of_gyration, radius_of_gyration_unscaled, magnitude_normalise, r_vector_1[3], r_vector_2[3];
double _Complex sum, temp_E, exp_term, rotated_x, rotated_y, rotated_z, ratio1, ratio2, ext_sum, scalar_a, scalar_B, resi_square, scal_alph_denom, old_resi_square, scalar_t_s, scalar_t_t, scalar_w_BiCGSTAB, rho_k_old, rho_k_new;
double wavelength_min, wavelength_max, wavelength_increment, log_increment, theta_3D, phi_3D,  g_asymmetry_DDA, validity_test_DDA, Beta_DDA, angle_aerosol_radius, radius_of_rotation, ensemble_average, hydrodynamic_radius;
double start_time, current_time, end_time;
const double pi = 3.14159265358979323846264;
double* wavelength_list;
double* radius_list;
double _Complex* E;
double _Complex* E_inc;
double _Complex* conj_E;
double _Complex* P;
double _Complex* conj_P;
double _Complex* P_new;
double _Complex* temp_array;
double _Complex* residual;
double _Complex* residual_hat;
double _Complex* search_direction;
double _Complex* A_dk;
double _Complex* h;
double _Complex* S;
double _Complex* t;
double*** results;
double*** Bessel;
double** refractive_index_raw_data;
double** averaged_results;
double** angle_averaged_results_y;
double** angle_averaged_results_z;
double** polarisation_averaged_results;
double** mie_results;
double** MMF_results;
double** MFT_results;
double** lattice;
double** sphere_dipole_positions;
double** original_dipole_positions;
double** dipole_positions;
double** STAG_dipole_positions;
double** vertices;
double _Complex** A;
double _Complex** refractive_index_data;
double _Complex* v_k;
double _Complex* p_k;
double _Complex**** vectors;
char buf[1000];

// Variables used to find g=<cos(theta)> (same names as in DDSCAT)
int ICOSTH, IPHI, NTHETA, NPHI;
double THETA, PHI, THETAL, THETAU, AFAC, SI, TERM, THETA0, DDA_X, DOMEGA, ETASCA, DPHI;
double AK_TF[3], AKK, AK2, E02, DX[3], AKS1[3], AKS2[3], AKS0[3], AKSN[3], EINC, RRR, ES2, COSTH, SINTH, COSPHI, SINPHI, DDA_W, CSCA, CSCAG, FALB, QSCA, QSCAG;
double* SCRRS2;
double** SCRRS1;
double _Complex* CXSCR1;
double _Complex* CXSCR2;
double _Complex CXE01_TF[3], CXES[3];
double _Complex** CXP_TF;

// Mie variables 

int mie_n_max;
double x_size_parameter, mie_Q_scatter, mie_Q_extinction, mie_Q_absorption, mie_C_scatter, mie_C_ext, mie_C_abs, mie_theta, sum_scatter, sum_extinction, sum_absorption, g_asymmetry_mie, sum_cos_Q_1st_term, sum_cos_Q_2nd_term, cos_Q_mie;
double _Complex a_n, b_n;

// MMF variables

int MMF_method, N_monomers, p_MMF, p_max, p_min, v_MMF, v_max, v_min, n_MMF, n_max, n_min, num_steps_integral, bessel_order, num_points_quadrature=400, largest_p_max, largest_n_max;
double r_monomer, k_frac_prefactor, d_f, x_agg, x_mon, du, u_min, u_max, f_c, MMF_Q_extinction, MMF_Q_scatter, MMF_Q_absorption, MMF_C_ext, MMF_C_scatter, MMF_C_abs, RDG_C_abs, MFT_C_abs, largest_size_parameter_x, integrated_C_sca, MFT_C_scatter;
double a_vnp, b_vnp, x_LP, integrated_a_vnp, integrated_b_vnp , tau_MMF, G_MMF, scatter_angle, x_val, q_coeff, S11[181], S_q[181], phase_function[181], dx_scatter;
double MFT_Q_scatter, MFT_Q_absorption, g_asymmetry_MMF, mie_aerosol_area;
double* GL_roots;
double* GL_weights;
double* pi_function;
double* tau_function;
double** LegendreP_at_root;
double** diff_LegendreP_at_root;
double** Associated_LegendreP_at_root;
double _Complex S_p, integrated_S, A_1_v_n, B_1_v_n, A_1_v_n_sum, B_1_v_n_sum, bessel_argument, bessel_argument_2;
double _Complex* structure_factor_terms;
double _Complex* d_terms;
double _Complex* d_1_terms;
double _Complex* d_2_terms;
double _Complex* a_n_terms;
double _Complex* b_n_terms;
double _Complex* final_col;
double _Complex** matrix_MMF;
double _Complex** original_matrix_MMF;

int main()
{
    start_time= omp_get_wtime(); //record start time (so that we can calculate the runtime later)

    printf("\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n\n                             WELCOME TO CORAL (Comparasion Of Radiative AnaLyses) -- MMF and Mie only!          \n\n");
    printf("\n                                                                                                                         ");
    printf("\n                                                           O O         O O                                                   ");
    printf("\n                                       \\                   O           O                                                 ");
    printf("\n                         ~~~~~~~~~~~~   \\                   O O       O        O  O                                     ");
    printf("\n                         ~~~~~~~~~~~~    \\                   O       O          OO       OO                             ");
    printf("\n                         ~~~~~~~~~~~~    /                  OOOO O  OOOOO     O        O                                 ");
    printf("\n                         ~~~~~~~~~~~~   /                  O   OOOOO  OO   OOOO     OOO                                  ");
    printf("\n                                       /                        O        OOOO   OO  O   OO                               ");
    printf("\n                                                               O O        O       OO                                     ");
    printf("\n                                                              O            O        O                                    ");
    printf("\n                                                             O O          O O                                              ");
    printf("\n                                                            O   O            O                                            ");
    printf("\n\n\n ---------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------------------------------------------------------\n\n");

    /* read in parameters from input file */

    if((parameterfile=fopen("parameters.txt","r")) == NULL){
          printf("\n\nError- the parameter file cannot be found!! \n\n\n");
          getchar();
          return 1;
    }
    else{
        parameterfile=fopen("parameters.txt","r");
    }  

    fscanf(parameterfile," %*s %lf %*s %*s %*s %*s ", &aerosol_radius_min);
    fscanf(parameterfile," %*s %lf %*s %*s %*s %*s ", &aerosol_radius_max);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s ", &num_radii);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s %*s %*s ", &radius_increment_type);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s ", &speed_test);

    fscanf(parameterfile," %*s %*s %*s %*s %*s %lf %*s %*s ", &wavelength_min);
    fscanf(parameterfile," %*s %lf %*s %*s ", &wavelength_max);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s ", &num_wavelengths);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s %*s %*s ", &wavelength_increment_type);

    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s %*s %*s ", &number_of_threads);

    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s ", &N_monomers);
    fscanf(parameterfile," %*s %lf %*s %*s %*s ", &k_frac_prefactor);
    fscanf(parameterfile," %*s %lf %*s %*s ", &d_f);
    fscanf(parameterfile," %*s %d %*s %*s %*s %*s %*s %*s %*s %*s ", &MMF_method);


    /* calculate list of wavelengths to analyse */

    wavelength_list=(double*)malloc(num_wavelengths*sizeof(double));
    wavelength_list[0]=wavelength_min; //set initial wavelength

    //printf("\n WAVELENGTH LIST: \n");
    //printf("\n Wavelength[%3d] = %f    ", i, wavelength_list[0]);

    if(num_wavelengths>1){ //only calculate wavelengths if we are looking for more than one wavelength to prevent dividing by 0 (this is particularly important for the batch create program, where we only analyse one at a time)
        if(wavelength_increment_type==0){ //use linearly-spaced wavelengths
            wavelength_increment= (wavelength_max-wavelength_min)/(num_wavelengths-1); //calculate what the increment should be
            for(i=1;i<num_wavelengths;i++){
                wavelength_list[i] =  wavelength_list[i-1] + wavelength_increment;
                //printf("\n Wavelength[%3d] = %f    linear increment = %f ", i, wavelength_list[i], wavelength_list[i]-wavelength_list[i-1]);
            }
        }
        else{ //use log-spaced wavelengths (the logs of any two adjacent wavelengths are seperated by a constant amount)
            log_increment = (log(wavelength_max)-log(wavelength_min))/(num_wavelengths-1);

            for(i=1;i<num_wavelengths;i++){
                wavelength_list[i] =  wavelength_list[i-1] * exp(log_increment);
                //printf("\n Wavelength[%3d] = %f    linear increment = %f    log_increment = %f (should be constant)", i, wavelength_list[i], wavelength_list[i]-wavelength_list[i-1], log_increment);
            }
        }
    }

    /* calculate list of aerosol radii to analyse */
    
    radius_list=(double*)malloc(num_radii*sizeof(double));
    radius_list[0]=aerosol_radius_min; //set initial radius

    if(num_radii>1){ //only calculate radii if we are looking for more than one radius to prevent dividing by 0 (this is particularly important for the batch create program, where we only analyse one at a time)
        if(radius_increment_type==0){ //use linearly-spaced radii
            aerosol_radius_increment= (aerosol_radius_max-aerosol_radius_min)/(num_radii-1); //calculate what the increment should be
            for(i=1;i<num_radii;i++){
                radius_list[i] =  radius_list[i-1] + aerosol_radius_increment;
                //printf("\n radius[%3d] = %f    linear increment = %f ", i, radius_list[i], radius_list[i]-radius_list[i-1]);
            }
        }
        else{
            log_increment = (log(aerosol_radius_max)-log(aerosol_radius_min))/(num_radii-1);

            for(i=1;i<num_radii;i++){
                radius_list[i] =  radius_list[i-1] * exp(log_increment);
                //printf("\n radius[%3d] = %f    linear increment = %f    log_increment = %f (should be constant)", i, radius_list[i], radius_list[i]-radius_list[i-1], log_increment);
            }
        }
    }

    // for the scaling of the visual output files (e.g. STAG), use the smallest aerosol radius in the list.
    aerosol_radius=radius_list[0];

    /* calculate number of vertices required for orientational average */

    if(orientational_average_resolution==0){
        num_vertices= 12;
    }
    else{
        num_vertices= 12 + 20*orientational_average_resolution + 10*pow(orientational_average_resolution,2);
    }

    /* Set other basic parameter values */

    ref_ind_medium = 1.0 + 0.0*I;
    magnetic_perm_medium= 0.0000012566; // 1.26 uH/m
    magnetic_perm_sphere= 0.0000012566; // 1.26 uH/m
    relative_permeability_sphere= magnetic_perm_sphere/0.0000012566; //should come out as "1" usually

    /* Initialise other coefficients */

    b_1= -1.891531;
    b_2= 0.1648469;
    b_3= -1.7700004;

    /* print summary of imported parameters to screen */

    printf(" INITIAL PARAMETERS: \n");
    printf("\n\t medium refractive index = %f + %g i \n\t relative permeability sphere = %g", creal(ref_ind_medium), cimag(ref_ind_medium), relative_permeability_sphere);

    if(radius_increment_type==0){
        printf("\n\n   Radii: %g um -> %g um (linearly-spaced increments of %g um - %d to analyse in total) ", aerosol_radius_min, aerosol_radius_max, aerosol_radius_increment, num_radii);
    }
    else{
        printf("\n\n   Radii: %g um -> %g um (log-spaced increments - %d to analyse in total) ", aerosol_radius_min, aerosol_radius_max, num_radii);
    }

    if(wavelength_increment_type==0){
        printf("\n\n   Wavelength: %g um -> %g um (linearly-spaced increments of %g um - %d to analyse in total) ", wavelength_min, wavelength_max, wavelength_increment, num_wavelengths);
    }
    else{
        printf("\n\n   Wavelength: %g um -> %g um (log-spaced increments - %d to analyse in total) \n", wavelength_min, wavelength_max, num_wavelengths);
    }

    if(number_of_threads>=100000){
        printf("\n   Working from any cores that are available. \n");
    }
    else{
         omp_set_num_threads(number_of_threads);
         printf("\n   Working from %d cores. \n", number_of_threads);
    }


    printf("\n MMF SECTION: \n\t Number of monomers: %d \n\t Fractal prefactor k0: %f \n\t Fractal dimension d_f: %f \n\t MMF_method = %d (Gaussian cutoff = 0, Fractal cutoff = 1)\n", N_monomers, k_frac_prefactor, d_f, MMF_method);

    if(speed_test==1){
        timefile=fopen("times.txt","w"); //open file for saving times of key processes throughout program
    }

    /* import refactive index data */

    if((refractiveindexfile=fopen("refractive_index_data.txt", "r")) == NULL){
        printf("\n\nError- the refractive index data cannot be found!! Please check that the file is titled 'refractive_index_data.txt' and that it is located in the current working directory...\n\n\n");
        getchar();
        return 1;
    }
    else{

        // open file and go through it once to count number of data points

        refractiveindexfile=fopen("refractive_index_data.txt", "r");
        number_ref_ind_data_points= -1; // start at -1 so we skip the header line in the data file
        while (fgets(buf,1000, refractiveindexfile)!=NULL){ //if something is read
            number_ref_ind_data_points++;
        }
        fclose(refractiveindexfile);

        printf("\n Refractive index data loaded: %d data points detected.", number_ref_ind_data_points);

        if(number_ref_ind_data_points<2){
            printf("\n\nError- not enough data points in refractive index file. Add at least 2 to be able to extrapolate, and try again! \n\n\n");
            getchar();
            return 1;
        }

        /* now that we know number of data points, initialise a "raw data" matrix to store data from input file */
        refractive_index_raw_data=(double**)malloc(number_ref_ind_data_points*sizeof(double*));
        for (i=0; i<number_ref_ind_data_points; i++) {
            refractive_index_raw_data[i]=(double*)malloc(4*sizeof(double));  //e.g. refractive_index_raw_data[n][0,1,2,3] where [n] is line number and 0,1,2,3 = column
        }

        // reopen original file and load the data into "raw_data" matrix

        refractiveindexfile=fopen("refractive_index_data.txt", "r");
        number_ref_ind_data_points=0;
        fscanf(refractiveindexfile, " %*s %*s %*s %*s "); //remove header line
        while(fscanf(refractiveindexfile, " %lf %lf %lf %lf ", &refractive_index_raw_data[number_ref_ind_data_points][0], &refractive_index_raw_data[number_ref_ind_data_points][1], &refractive_index_raw_data[number_ref_ind_data_points][2], &refractive_index_raw_data[number_ref_ind_data_points][3])>0){
            number_ref_ind_data_points++; // increment number of data points in the file
        }
        fclose(refractiveindexfile);


            /* merge data into a single complex number array */
        refractive_index_data=(double _Complex**)malloc(number_ref_ind_data_points*sizeof(double _Complex*));
        for (i=0; i<number_ref_ind_data_points; i++) {
            refractive_index_data[i]=(double _Complex*)malloc(2*sizeof(double _Complex));  //e.g. '[n][0] = wavelength' and '[n][1] = complex n + ki value'
        }


        //combine refractive index data into a complex matrix that stores both wavelength and complex refractive index
        for(i=0; i<number_ref_ind_data_points; i++){
            refractive_index_data[i][0] = refractive_index_raw_data[i][1]; //wavelength
            refractive_index_data[i][1] = refractive_index_raw_data[i][2] + refractive_index_raw_data[i][3]*I; //n + ki
            //printf("\n %10f %10f %10f ", creal(refractive_index_data[i][0]),creal(refractive_index_data[i][1]),cimag(refractive_index_data[i][1]));
        }
    }






        /* CALCULATE MIE SCATTERING FOR A SPHERE OF EQUIVALENT VOLUME FOR ALL WAVELENGTHS AND RADII ANALYSED */

    printf("\n\n\n\n -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n -------------------------------------------------------------------------------     MIE SCATTERING    -------------------------------------------------------------------------------");
    printf("\n -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

    mie_results=(double**)malloc(num_wavelengths*sizeof(double*));
    for (i=0; i<num_wavelengths; i++) {
        mie_results[i]=(double*)malloc(7*sizeof(double));  //makes a 2D matrix with index [a][b] where [a]= wavelength index and [b]= 6 elements [0,1,2] = Q_ext, Q_scatter, Q_abs and elements  [3,4,5] = C_ext, C_scatter, C_abs and [6] is g=<cos(theta)>
    }

    /* open file for saving results */
    mieresultsfile=fopen("mie_results.txt","w"); //open file for saving results
    fprintf(mieresultsfile,"    Radius(um)       Wavelength(um)         Q_Ext          Q_Scatter        Q_Abs             C_Ext         C_Scatter         C_Abs      g=<cos(theta)>\n");

    
    /* BEGIN RADIUS LOOPS HERE */

    for(radius_index=0; radius_index<num_radii; radius_index++){

        aerosol_radius=radius_list[radius_index]; //set new aerosol radius

        /* BEGIN WAVELENGTH LOOPS HERE */
        for(wavelength_index=0; wavelength_index<num_wavelengths; wavelength_index++){

            wavelength = wavelength_list[wavelength_index]; //look up wavelength values from list made at beginning of program (depends on whether they are chosen to be linearly-spaced or log-spaced)

            wavenumber_k= 2.0*pi/wavelength;
            ref_ind_sphere=get_refractive_index(wavelength, refractive_index_data, number_ref_ind_data_points); //remember to do this in the Mie section too - refractive index is wavelength-dependent!

            x_size_parameter= 2*pi*aerosol_radius/wavelength; // wavelength and sphere radius were both in m (rather than um) in the original mie program, but this is the only calculation that uses them (and it's a ratio), so we can leave them both in um.
            mie_n_max = calculate_n_max(x_size_parameter);

            printf("\n\n INITIAL PARAMETERS: \n\n Wavelength = %g um \n Mie sphere radius = %g um \n Sphere refractive index = %g + %g i \n Medium refractive index = %g + %g i \n Size parameter is: %f \n mie_n_max= %d \n", wavelength, aerosol_radius, creal(ref_ind_sphere), cimag(ref_ind_sphere), creal(ref_ind_medium) ,cimag(ref_ind_medium), x_size_parameter, mie_n_max);
            
            // arrays to store a_n and b_n (standard monomer mie coefficients)
            a_n_terms = (double _Complex*)malloc((mie_n_max+1)*sizeof(double _Complex));  // initialise array to store a_n terms
            b_n_terms = (double _Complex*)malloc((mie_n_max+1)*sizeof(double _Complex));  // initialise array to store b_n terms

            a_n_terms[0]=0;
            b_n_terms[0]=0; //initialise these array elements (we will never use them, but just to be safe)

            //initialise summation terms
            sum_scatter=0;
            sum_extinction=0;
            for(i=1;i<=mie_n_max;i++){
                
                /* Calculate coefficients a and b */
                get_a_n_and_b_n(i, x_size_parameter, ref_ind_sphere, magnetic_perm_sphere, magnetic_perm_medium, &a_n, &b_n); //i represents term n (e.g. i=1 finds a_1, and i=2 finds a_2...)
                
                a_n_terms[i]= a_n; //save both of these terms in an array for use in the calculation of g_asymmetry_mie
                b_n_terms[i]= b_n;

                //printf("\n a_%d is: %f + %f i         b_%d is: %f + %f i\n", i, creal(a_n), cimag(a_n), i, creal(b_n), cimag(b_n));

                /* calculate summation terms for Q_sca and Q_ext */
                sum_scatter= sum_scatter + (2*i+1)*(pow(creal(a_n),2) + pow(cimag(a_n),2) + pow(creal(b_n),2) + pow(cimag(b_n),2));
                sum_extinction= sum_extinction + (2*i+1)*(creal(a_n+b_n));
            }

            mie_Q_scatter= (2/pow(x_size_parameter,2))*sum_scatter;
            mie_Q_extinction = (2/pow(x_size_parameter,2))*sum_extinction;
            mie_Q_absorption = mie_Q_extinction - mie_Q_scatter;

            mie_C_scatter= mie_Q_scatter*pi*pow(aerosol_radius,2);
            mie_C_ext= mie_Q_extinction*pi*pow(aerosol_radius,2);
            mie_C_abs= mie_Q_absorption*pi*pow(aerosol_radius,2);

            /* calculate summation terms for g_asymmetry_mie */
            
            sum_cos_Q_1st_term=0;
            sum_cos_Q_2nd_term=0;
            for(i=1;i<mie_n_max;i++){ //only go up to mie_n_max-1, because we require term i+1
                sum_cos_Q_1st_term += (i*(i+2.0)/(i+1.0)) * creal(a_n_terms[i]*cconjugate(a_n_terms[i+1]) + b_n_terms[i]*cconjugate(b_n_terms[i+1]));
                sum_cos_Q_2nd_term += ((2.0*i+1.0)/(i*(i+1.0))) * creal(a_n_terms[i]*cconjugate(b_n_terms[i]));
            }

            // calculate g_asymmetry_mie
            cos_Q_mie = (4.0/pow(x_size_parameter,2.0)) * (sum_cos_Q_1st_term + sum_cos_Q_2nd_term);
            g_asymmetry_mie = cos_Q_mie/mie_Q_scatter;
            printf("\n cos_Q_Mie = %f    g_asymmetry_Mie = %f \n", cos_Q_mie, g_asymmetry_mie);

            /* save results from mie analysis for this particular wavelength */

            mie_results[wavelength_index][0]=mie_Q_extinction;
            mie_results[wavelength_index][1]=mie_Q_scatter;
            mie_results[wavelength_index][2]=mie_Q_absorption;
            mie_results[wavelength_index][3]=mie_C_ext;
            mie_results[wavelength_index][4]=mie_C_scatter;
            mie_results[wavelength_index][5]=mie_C_abs;
            mie_results[wavelength_index][6]=g_asymmetry_mie;

            printf("\n ------------------------------- FOR WAVELENGTH %5.3f um ----------------------------------------------------------", wavelength);
            printf("\n\n Q EFFICIENCIES: \n\n                Q (Extinction)          Q (Scatter)          Q (Absorption) \n");
            printf("\n     Mie \t  %f \t         %f \t       %f \n",  mie_Q_extinction, mie_Q_scatter, mie_Q_absorption);
            printf("\n -------------------------------------------------------------------------------------------------------------------");
            printf("\n\n CROSS-SECTIONS (um^2): \n\n                C (Extinction)          C (Scatter)          C (Absorption) \n");
            printf("\n     Mie \t  %f \t         %f \t       %f \n",  mie_C_ext, mie_C_scatter, mie_C_abs);
            printf("\n ------------------------------------------------------------------------------------------------------------------- \n\n");

            //free arrays on each loop
            free((void*)a_n_terms); //NOTE: these have to be created/freed WITHIN the loops because they depend on aize parameter which is a function of wavelength
            free((void*)b_n_terms);
        }

        /* save mie results */
        for(i=0;i<num_wavelengths;i++){
            fprintf(mieresultsfile, "    %10e      %10e      %10e     %10e    %10e      %10e    %10e    %10e    %10e\n", aerosol_radius, wavelength_list[i], mie_results[i][0], mie_results[i][1], mie_results[i][2], mie_results[i][3], mie_results[i][4], mie_results[i][5], mie_results[i][6]);
        }
    

    } //return to top of radius loop and repeat for next value

    fclose(mieresultsfile);








        /* CALCULATE MMF SCATTERING FOR ALL RADII AND WAVELENGTHS */

    printf("\n\n\n\n -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n ---------------------------------------------------------------------    MODIFIED MEAN FIELD THEORY (MMF)   -------------------------------------------------------------------------");
    printf("\n -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");

    if(speed_test==1){
        current_time= omp_get_wtime();
        fprintf(timefile, " %12f   Beginning MMF analysis.\n", current_time-start_time);
    }

    //set up matrix to store results of MMF
    MMF_results=(double**)malloc(num_wavelengths*sizeof(double*));
    for (i=0; i<num_wavelengths; i++) {
        MMF_results[i]=(double*)malloc(7*sizeof(double));  //makes a 2D matrix with index [a][b] where [a]= wavelength index and [b]= 6 elements [0,1,2] = Q_ext, Q_scatter, Q_abs and elements  [3,4,5] = C_ext, C_scatter, C_abs and [7] is g=<cos(theta)>
    }

    //set up matrix to store results of MFT
    MFT_results=(double**)malloc(num_wavelengths*sizeof(double*));
    for (i=0; i<num_wavelengths; i++) {
        MFT_results[i]=(double*)malloc(7*sizeof(double));  //makes a 2D matrix with index [a][b] where [a]= wavelength index and [b]= 6 elements [0,1,2] = Q_ext, Q_scatter, Q_abs and elements  [3,4,5] = C_ext, C_scatter, C_abs and [7] is g=<cos(theta)>
    }



        /* PREPARE MATRICES AND VALUES FOR GAUSS-LEGENDRE INTEGRATION */

    num_points_quadrature=400;

    if(speed_test==1){
        current_time= omp_get_wtime();
        fprintf(timefile, " %12f   Beginning Gauss-Legendre quadrature.\n", current_time-start_time);
    }

    // find the positions of the Gauss-Legendre points (roots of the poynomials) and weights. Place this outside the main MMF wavelength loop, because we only need to calculate the values once */
            
    GL_roots = (double*)malloc(num_points_quadrature*sizeof(double));  //initialise array to hold the roots of the Legendre Polynomial
    GL_weights = (double*)malloc(num_points_quadrature*sizeof(double));  //initialise array to hold the roots of the Legendre Polynomial

    printf("\n Initiating Gauss-Legendre quadrature. Calculating Legendre Polynomial values at %d roots...", num_points_quadrature);
    
    calculate_GL_roots(num_points_quadrature, GL_roots); //find the roots of the Legendre polynomials
    calculate_GL_weights(num_points_quadrature, GL_roots, GL_weights); //find the weights of the Legendre polynomials

    if(speed_test==1){
        current_time= omp_get_wtime();
        fprintf(timefile, " %12f   Gauss-Legendre roots and weights calculated (%d of them).\n", current_time-start_time, num_points_quadrature);
    }

    //for(i=0; i<num_points_quadrature; i++){
        //printf("\n %d th root: %12f    Weight: %12f ", i, GL_roots[i], GL_weights[i]);
    //}

    /* open results files */
    mmfresultsfile=fopen("mmf_results.txt","w"); //open file
    fprintf(mmfresultsfile,"   Radius(um)        Wavelength(um)         Q_Ext          Q_Scatter        Q_Abs             C_Ext         C_Scatter         C_Abs     g=<cos(theta)>\n");
        
    mftresultsfile=fopen("mft_results.txt","w"); //open file
    fprintf(mftresultsfile,"   Radius(um)        Wavelength(um)         Q_Ext          Q_Scatter        Q_Abs             C_Ext         C_Scatter         C_Abs     g=<cos(theta)>\n");
        


    /* BEGIN RADIUS LOOPS */
    for(radius_index=0; radius_index<num_radii; radius_index++){

        aerosol_radius=radius_list[radius_index]; //set new aerosol radius

        
        /* calculate other required parameters from the number of monomers */

        mie_aerosol_area= pi*aerosol_radius*aerosol_radius; //calculate area of spherical particle equivalent (used to find Q_Ext, Q_Sca and Q_abs)
        
        /* calculate radius of monomers (volume is conserved between the spherical and non-spherical MMF particle equivalents */
        r_monomer = aerosol_radius/pow(N_monomers,(1.0/3.0));

        /* calculate radius of gyration */
        radius_of_gyration = r_monomer*pow(N_monomers/k_frac_prefactor, (1.0/d_f));

        /* calculate the geometrical cross-section of the fractal aggregate using Tazeki's "Analytic expressions for geometric cross-sections of fractal dust aggregates". 
        Have benchmarked against the code in optools_geofractal.f90 */

        G_MMF = calculate_geometrical_cross_section(N_monomers, r_monomer, d_f, k_frac_prefactor, pi);        
        
        
        printf("\n For MMF, we consider an aerosol shape composed of %d individual monomers each with a radius of %f um, and a geometrical cross-section of %f um^2.\n", N_monomers, r_monomer, G_MMF);        
        printf(" The monomer radius is derived by ensuring that volume is conserved between the spherical and non-spherical models.\n");        
        printf("\n\t SPHERICAL MODEL (Mie): Radius is %f um, and volume is given by (4/3)*pi*R^3 = %f um^3 \n\t NON-SPHERICAL MODEL (MMF): volume given by N_mon*(4/3)*pi*r_mon^3 = %f um^3\n", aerosol_radius, (4.0/3.0)*pi*pow(aerosol_radius,3.0), N_monomers*(4.0/3.0)*pi*pow(r_monomer,3.0));        
        printf("\n The shape has a fractal dimension of %f and fractal prefactor of %f, read from the input file. \n", d_f, k_frac_prefactor);
        if(MMF_method==0){
            printf("\n A gaussian cutoff function is used to calculate the structure factor term (MMF_method = %d).\n", MMF_method);
        }
        else{
            printf("\n A fractal dimension cutoff function is used to calculate the structure factor term (MMF_method = %d).\n", MMF_method);
        }

        // find absolute largest p_max and n_max from largest size_parameter that will be studied (smallest wavelength, largest particle size)
        largest_size_parameter_x = 2.0*pi*r_monomer/wavelength_min; // the MMF size parameter is based on the monomer radius
        largest_p_max = 2*calculate_n_max(largest_size_parameter_x); // p_max = n_max+v+max, and v_max = n_maxn, so just multiply by 2
        largest_n_max = calculate_n_max(largest_size_parameter_x); //demo version for now! add a function to calculate using note above, using n_max for the largest size parameter
        
        // initialise matrices to store Legendre Polynomial values at each Gauss-Legendre root P(i, x_i) up to P(p_max, x_i)
        LegendreP_at_root=(double**)malloc((largest_p_max+1)*sizeof(double*));
        for (i=0; i<(largest_p_max+1); i++) {
            LegendreP_at_root[i]=(double*)malloc(num_points_quadrature*sizeof(double)); 
        }  

        // initialise matrices to store differentiated Legendre Polynomial values at each Gauss-Legendre root P'(i, x_i) up to P'(p_max, x_i)
        diff_LegendreP_at_root=(double**)malloc((largest_p_max+1)*sizeof(double*));
        for (i=0; i<(largest_p_max+1); i++) {
            diff_LegendreP_at_root[i]=(double*)malloc(num_points_quadrature*sizeof(double));  
        }  

        // initialise matrices to store 'degree 1' Associated Legendre Polynomial values at each Gauss-Legendre root P(i, 1, x_i) up to P(n_max, 1, x_i) - note smaller matrix size needed for these times
        Associated_LegendreP_at_root=(double**)malloc((largest_n_max+1)*sizeof(double*));
        for (i=0; i<(largest_n_max+1); i++) {
            Associated_LegendreP_at_root[i]=(double*)malloc(num_points_quadrature*sizeof(double)); 
        }  
        
            /* store the values of each of the three Legendre types (regular, differentiated and 1st degree associated) at each of the Gauss-Legendre root values found above */
        
        // Store values for regular and differentiated Legendre Polynomials first (their matrix is larger)
        for(i=0; i<(largest_p_max+1); i++){
            for(j=0; j<num_points_quadrature; j++){
                LegendreP_at_root[i][j]= LegendrePoly(i,GL_roots[j]); // Calculate P(i, x_i) for each Gauss-Legendre root up to polynomial p_max (the maximum one that will be needed in calculations)
                diff_LegendreP_at_root[i][j]= diff_LegendrePoly(i,GL_roots[j]); // Calculate P'(i, x_i) for each Gauss-Legendre root up to polynomial p_max (the maximum one that will be needed in calculations)
                //printf("\n Legendre P(%d, %f) = %g       diff_Legendre P(%d, %f) = %g", i, GL_roots[j], LegendreP_at_root[i][j], i, GL_roots[j], diff_LegendreP_at_root[i][j]);
            }
        }

        // Store values for associated Legendre Polynomials (smaller matrix, we only ever need to go up to n_max for this type)
        for(i=0; i<(largest_n_max+1); i++){
            for(j=0; j<num_points_quadrature; j++){
                Associated_LegendreP_at_root[i][j]= Associated_LegendrePoly(i,GL_roots[j]); // Calculate P(i, 1, x_j) for each Gauss-Legendre root up to polynomial n_max (the maximum one that will be needed in calculations)
                //printf("\n Associated P(%d, 1, %f) = %g", i, GL_roots[j], Associated_LegendreP_at_root[i][j]);
            }
        }




        if(speed_test==1){
            current_time= omp_get_wtime();
            fprintf(timefile, " %12f   Legendre Polynomial values calculated and stored for up to %d orders for each of the %d roots.\n", current_time-start_time, largest_p_max+1, num_points_quadrature);
        }

            /* BEGIN WAVELENGTH LOOPS FOR MMF */

        for(wavelength_index=0; wavelength_index<num_wavelengths; wavelength_index++){

            wavelength = wavelength_list[wavelength_index]; //look up wavelength values from list made at beginning of program (depends on whether they are chosen to be linearly-spaced or log-spaced)
            
            // derived parameters from inputs
            ref_ind_sphere=get_refractive_index(wavelength, refractive_index_data, number_ref_ind_data_points); //remember to do this in the Mie/MMF section too - refractive index uses the same wavelength-dependent laws throughout each analysis!
            wavenumber_k= 2.0*pi/wavelength;

            x_agg= 2.0*pi*radius_of_gyration/wavelength;
            x_mon= 2.0*pi*r_monomer/wavelength; // different to Mie sphere model - size parameter of a SINGLE MONOMER (uses the monomer radius, not mie sphere of equivalent volume to bulk fractal)   
            n_max= calculate_n_max(x_mon); //x represents the monomer size parameter in Tazeki, not x_c (the aggregate size parameter)
            n_min= 1;
            n_MMF= n_min; //initialise first term

            printf(" \n MMF inputs: \n\n \t number of monomers = %d \n\t monomer radius =  = %f \n\t wavelength = %f \n\t k_prefactor = %f \n\t refractive index = %f + %f i \n\t wavenumber_k = %f \n\t d_f = %f \n\t radius of gyration = %f \n\t x_agg = %f \n\t x_monomer = %f \n\t ", N_monomers, r_monomer, wavelength, k_frac_prefactor, creal(ref_ind_sphere), cimag(ref_ind_sphere), wavenumber_k, d_f, radius_of_gyration, x_agg, x_mon);

            printf("\n n_min= %d\n n_max= %d\n n_MMF= %d\n", n_min, n_max, n_MMF);

            //find min and max limits for v
            v_min= 1;
            v_max= calculate_n_max(x_mon); //this is the same equation as for the Mie terms
            v_MMF= v_min; //initialise first term

            //printf("\n v_min= %d\n v_max= %d\n v_MMF= %d\n", v_min, v_max, v_MMF);

            // initiliase arrays to store matrix and calculate d_terms -- do this once per wavelength/particle radius (as size parameter changes, n_max changes, and therefore the matrix size changes)
            matrix_MMF=(double _Complex**)malloc(2*n_max*sizeof(double _Complex*));
            for (i=0; i<(2*n_max); i++) {
                matrix_MMF[i]=(double _Complex*)malloc(2*n_max*sizeof(double _Complex));  //makes an n_max x n_max 2D matrix to store the terms in the MMF matrix
            }    

            // arrays to store a_n and b_n (standard monomer mie coefficients)
            a_n_terms = (double _Complex*)malloc((n_max+1)*sizeof(double _Complex));  // initialise array to store a_n terms
            b_n_terms = (double _Complex*)malloc((n_max+1)*sizeof(double _Complex));  // initialise array to store b_n terms

            a_n_terms[0]=0;
            b_n_terms[0]=0; //initialise these array elements (we will never use them, but just to be safe)

                /* CALCULATE STRUCTURE FACTOR TERMS */

            // find limits to integral
            u_min= x_agg*exp(-40.0/d_f);

            if(MMF_method==0){
                u_max = 2.0* x_agg * sqrt(25.0/d_f); // Tazaki + Tanaka 2016 (GAUSS)
            }
            else{
                u_max = x_agg*pow((2.0*25.0),(1.0/d_f)); // Tazaki + Tanaka 2016 (FLDIM)
            }
            num_steps_integral = 10000; //same value used in optools_fractal.f90
            du=(u_max-u_min)/num_steps_integral; //calculate step size
            //printf("\n num steps = %d  \n u_min = %f    \n u_max = %f",  num_steps_integral, u_min, u_max);

            if(speed_test==1){
                current_time= omp_get_wtime();
                fprintf(timefile, " %12f   MMF variables set-up and ready for calculations.\n", current_time-start_time);
            }

            printf("\n Calculating structure factor terms...");

            structure_factor_terms= (double _Complex*)malloc((n_max+v_max+1)*sizeof(double _Complex));  // we will need all terms between 0 -> (n_max + v_max)

            for(p_MMF=0; p_MMF <= (n_max+v_max); p_MMF++){

                //printf(" \n p = %3d/%d ", p_MMF, n_max+v_max);
                
                /* calculate integral in S_p(k*Rg) term using Simpson's rule */

                integrated_S=0; //initialise integral sum

                //begin parallel calculation to integrate Legendre functions and find a(v,n,p)
                #pragma omp parallel for reduction (+:integrated_S)  //use "reduction" to ensure that each thread does not access integrated_a_vnp at the same time
                for(int z=0; z<num_steps_integral; z++){

                    double step_a= u_min + z*du; //set lower integral limit (this will be lowest possible value of integral (u_min) PLUS "number of steps X step size")
                    double step_b= step_a + du; //set upper integral limit     // NOTE: DEFINE ALL THESE VARIABLES WITHIN PARALLEL SECTION TO KEEP THESE PRIVATE
                    double step_c= (step_a+step_b)/2.0; //set middle integral term for Simpson's rule (we integrate between a and b above, but this term is useful)

                    //printf("\n\n\n We integrate between %f --> %f  (with a midpoint of %f) and where du= %f \n", step_a, step_b, step_c, du);

                    //find area underneath curve between points a and b
                    double _Complex step_area = (du/6.0) * (S_integrand(step_a, d_f, x_agg, p_MMF, MMF_method) + 4.0*S_integrand(step_c, d_f, x_agg, p_MMF, MMF_method) + S_integrand(step_b, d_f, x_agg, p_MMF, MMF_method)); // Use Simpson's rule to estimate the area under the step between a --> b. And yep: the order is a, c, b!

                    integrated_S += step_area; //add step area to total area under curve

                    //printf("\n step area = %f + %f i     integrated_S = %f + %f i ", creal(step_area), cimag(step_area), creal(integrated_S), cimag(integrated_S));
                }
                #pragma omp barrier

                //printf("\n p_MMF = %d \n d_f = %f \n x_agg = %f \n u_min= %f \n u_max = %f \n integrated_S = %f + %f i \n\n", p_MMF, d_f, x_agg, u_min, u_max, creal(integrated_S), cimag(integrated_S));

                structure_factor_terms[p_MMF] = (1.0/(2.0*pow(x_agg,d_f)))*integrated_S; // calculate and store structure factor term for each value of p
            
                //SAFETY CHECK - S(p) should always be >=0
                if(creal(structure_factor_terms[p_MMF])<0){
                    printf("\n WARNING: S[%d] < 0. Assuming it should be = 0 instead (this can happen for d_f > 2, and is ok...)", p_MMF);
                    structure_factor_terms[p_MMF]=0; //if the monomers are far enough apart that they don't affect each other, setting S_p = 0 should find a case for two non-interacting spheres. I made this decision independently, but Tazeki made the same note in optools_fractal.f90!
                }
            }

            /* Debugging - print structure factor terms 

            for(p_MMF=0; p_MMF <= (n_max+v_max); p_MMF++){
                printf("\n For p = %d    Sp: %14g + %14g i ", p_MMF, creal(structure_factor_terms[p_MMF]), cimag(structure_factor_terms[p_MMF]));
            } */

            /* store pi and tau function terms in arrays outside main loop */

            pi_function = (double*)malloc((n_max+1)*sizeof(double));
            tau_function = (double*)malloc((n_max+1)*sizeof(double));

            if(speed_test==1){
                current_time= omp_get_wtime();
                fprintf(timefile, " %12f   Structure factor array found.\n", current_time-start_time);
            }

            /* Initialise arrays for d_terms (two seperate arrays for easy reading) */

            d_1_terms = (double _Complex*)malloc((n_max+1)*sizeof(double _Complex));  // initialise array to store d_1 terms
            d_2_terms = (double _Complex*)malloc((n_max+1)*sizeof(double _Complex));  // initialise array to store d_2 terms

            d_1_terms[0]=0;
            d_2_terms[0]=0; //initialise these array elements (we will never use them, but just to be safe)

                
            if(N_monomers==1){
                /* simple case - the d_terms are just the a_terms */
                for(n_MMF=n_min; n_MMF <= n_max ; n_MMF++){

                    get_a_n_and_b_n(n_MMF, x_mon, ref_ind_sphere, magnetic_perm_sphere, magnetic_perm_medium, &a_n, &b_n); //n_MMF represents term n (e.g. n_MMF=1 finds a_1, and n_MMF=2 finds a_2...)

                    a_n_terms[n_MMF]= a_n; //save both of these terms in an array for use in the calculation of C_abs
                    b_n_terms[n_MMF]= b_n;

                    d_1_terms[n_MMF] = a_n;  //save all d_1 terms in their own array
                    d_2_terms[n_MMF] = b_n;   //save all d_2 terms in their own array
                }
            }
            else{
                /* begin looping through n and v and assigning A and B elements to the MMF_matrix */
                for(n_MMF=n_min; n_MMF <= n_max ; n_MMF++){

                    /* Calculate Mie coefficients a_n and b_n for the MONOMERS */
                    get_a_n_and_b_n(n_MMF, x_mon, ref_ind_sphere, magnetic_perm_sphere, magnetic_perm_medium, &a_n, &b_n); //n_MMF represents term n (e.g. n_MMF=1 finds a_1, and n_MMF=2 finds a_2...)

                    a_n_terms[n_MMF]= a_n; //save both of these terms in an array for use in the calculation of C_abs
                    b_n_terms[n_MMF]= b_n;

                    if(speed_test==1){
                        current_time= omp_get_wtime();
                        fprintf(timefile, " %12f   Mie coefficients a(%d) and b(%d) calculated.\n", current_time-start_time, n_MMF, n_MMF);
                    }
                    
                    //printf("\n a(%d) is: %f + %f i\n b(%d) is: %f + %f i\n", n_MMF, creal(a_n), cimag(a_n), n_MMF, creal(b_n), cimag(b_n));


                    for(v_MMF=v_min; v_MMF <= v_max; v_MMF++){   //this loop is only <= v_max (NOT 2*v_max) because we assign 4 terms (both A and B, in both the top and bottom "halves" of the matrix) each loop

                        //printf("\n Calculating A and B terms for [n][v]=[%d][%d]...", n_MMF,v_MMF);

                                /* CALCULATE A_1_v_n TERMS */

                        //find min and max limits for p
                        p_min= sqrt(pow((n_MMF - v_MMF),2));
                        p_max= n_MMF + v_MMF;
                        p_MMF= p_min; //initialise first term
                        //printf("\n p_min= %d\n p_max= %d\n p_MMF= %d\n", p_min, p_max, p_MMF);
                        
                        /* CALCULATE A_1_v_n TERMS */

                        A_1_v_n_sum=0; //initialise summation terms
                        B_1_v_n_sum=0;
                        for(p_MMF = p_min;  p_MMF <= p_max;  p_MMF++){

                                    /* STRUCTURE FACTOR */

                            S_p = structure_factor_terms[p_MMF]; //recover structure factor terms (calculated outside main loop)
                            
                                /* CALCULATE a(v,n,p) and b(v,n,p) */

                            if( ((is_odd(p_MMF)==0)&&(is_odd(n_MMF+v_MMF)==0)) || ((is_odd(p_MMF)==1)&&(is_odd(n_MMF+v_MMF)==1)) ){ //if p = even and (n+v) = even, or p= odd and (n+v) = odd, then B = 0 for these values, and we need to calculate A below
                                
                                    /* CALCULATE a(v,n,p) (and the B term is zero for this value of p, so no need to add to B_1_v_n_sum this time) */

                                /* find a_vnp term by adding terms of Gauss-Legendre series, the terms of which are calculated outside of the wavelength loop */
                                integrated_a_vnp = 0;
                                for(i=0; i<num_points_quadrature; i++){
                                    //printf(" \n i = %3d    weight = %9f     P1 = %14f   P2 = %14f   P3 = %14f   Adding %14f to %14f", i, GL_weights[i], Associated_LegendreP_at_root[v_MMF][i], Associated_LegendreP_at_root[n_MMF][i], LegendreP_at_root[p_MMF][i], GL_weights[i] * Associated_LegendreP_at_root[v_MMF][i]*Associated_LegendreP_at_root[n_MMF][i]*LegendreP_at_root[p_MMF][i], integrated_a_vnp);
                                    integrated_a_vnp +=  GL_weights[i] * Associated_LegendreP_at_root[v_MMF][i]*Associated_LegendreP_at_root[n_MMF][i]*LegendreP_at_root[p_MMF][i];  // we sum the weights x the integrand: w_i * P(v,1,x_i) * P(n,1,x_i) * P(p,x_i)
                                }
                                a_vnp = ((2.0*p_MMF + 1.0)/2.0)*integrated_a_vnp;
                                //printf("\n integrated_a_vnp = %f     a_vnp = %f \n", integrated_a_vnp, a_vnp);

                                A_1_v_n_sum += (n_MMF*(n_MMF+1) + v_MMF*(v_MMF+1) - p_MMF*(p_MMF+1))*a_vnp*S_p; 
                            }
                            else{ // otherwise, p and (n+v) have different parities (one is odd an one is even) -- so A = 0 for these values, and we need to calculate B below
                                
                                    /* CALCULATE b(v,n,p) (and the A term is zero for this value of p, so no need to add to A_1_v_n_sum this time) */
                                
                                /* find b_vnp term by adding terms of Gauss-Legendre series, the terms of which are calculated outside of the wavelength loop */
                                integrated_b_vnp = 0;
                                for(i=0; i<num_points_quadrature; i++){
                                    //printf(" \n i = %3d    weight = %9f     P1 = %14f   P2 = %14f   P3 = %14f   Adding %14f to %14f", i, GL_weights[i], Associated_LegendreP_at_root[v_MMF][i], Associated_LegendreP_at_root[n_MMF][i], diff_LegendreP_at_root[p_MMF][i], GL_weights[i] * Associated_LegendreP_at_root[v_MMF][i]*Associated_LegendreP_at_root[n_MMF][i]*diff_LegendreP_at_root[p_MMF][i], integrated_b_vnp);
                                    integrated_b_vnp +=  GL_weights[i] * Associated_LegendreP_at_root[v_MMF][i]*Associated_LegendreP_at_root[n_MMF][i]*diff_LegendreP_at_root[p_MMF][i];  // we sum the weights x the integrand: w_i * P(v,1,x_i) * P(n,1,x_i) * P(p,x_i)
                                }
                                b_vnp = ((2.0*p_MMF + 1.0)/2.0)*integrated_b_vnp;
                                //printf("\n integrated_b_vnp = %f     b_vnp = %f \n", integrated_b_vnp, b_vnp);
                                
                                B_1_v_n_sum += b_vnp*S_p;
                            }

                            //printf("\n A_1_v_n_sum = %f + %f i \n B_1_v_n_sum = %f + %f i", creal(A_1_v_n_sum), cimag(A_1_v_n_sum), creal(B_1_v_n_sum), cimag(B_1_v_n_sum));
                        }

                        if(speed_test==1){
                            current_time= omp_get_wtime();
                            fprintf(timefile, " %12f   Found A and B terms for [n][v]=[%d][%d]...\n", current_time-start_time, n_MMF,v_MMF);
                        }

                        //printf("\n\n\n\n A_1_v_n_sum = %f + %f i \n B_1_v_n_sum = %f + %f i", creal(A_1_v_n_sum), cimag(A_1_v_n_sum), creal(B_1_v_n_sum), cimag(B_1_v_n_sum));

                        A_1_v_n = ((2.0*v_MMF+1.0)/(n_MMF*(n_MMF+1.0)*v_MMF*(v_MMF+1.0)))*A_1_v_n_sum;
                        B_1_v_n = (2.0*(2.0*v_MMF+1.0)/(n_MMF*(n_MMF+1.0)*v_MMF*(v_MMF+1.0)))*B_1_v_n_sum;

                        //printf(" \n\n A[v][n] = A[%d][%d] = %f + %f i  \n\n B[v][n] = B[%d][%d] = %f + %f i", v_MMF, n_MMF, creal(A_1_v_n), cimag(A_1_v_n), v_MMF, n_MMF,  creal(B_1_v_n), cimag(B_1_v_n));

                        //assign terms to the matrix here! Remember to -1 from each element so that n=1 is stored as element [0] in both rows and columns

                        if(v_MMF==n_MMF){ // add extra term along the main diagonal
                            matrix_MMF[n_MMF-1][v_MMF-1] = (1.0/((N_monomers-1.0)*a_n)) + A_1_v_n; // top-left term (which is along main diagonal) - (Remember to -1 from ALL ELEMENTS so that n=1 v=1 is stored as element [0][0] etc
                            matrix_MMF[n_MMF-1][v_MMF+v_max-1] = B_1_v_n; //top-right term)
                            matrix_MMF[n_MMF+v_max-1][v_MMF-1] = B_1_v_n; // bottom-left term
                            matrix_MMF[n_MMF+v_max-1][v_MMF+v_max-1] = (1.0/((N_monomers-1.0)*b_n)) + A_1_v_n; //bottom-right term (which is along main diagonal)
                        }
                        else{ //all other terms follow same pattern. Add them 4 at a time (the n/v loops over 'v_max x v_max, which only a quarter of the full matrix '2v_max x 2v_max').
                            matrix_MMF[n_MMF-1][v_MMF-1] = A_1_v_n; // top-left term
                            matrix_MMF[n_MMF-1][v_MMF+v_max-1] = B_1_v_n; //top-right term
                            matrix_MMF[n_MMF+v_max-1][v_MMF-1] = B_1_v_n; // bottom-left term
                            matrix_MMF[n_MMF+v_max-1][v_MMF+v_max-1] = A_1_v_n; //bottom-right term
                        }
                    }
                }

                /* Optional print statements to print entire matrix to screen for debugging
                printf(" \n\n ");
                for(n_MMF=0; n_MMF < (2*n_max) ; n_MMF++){
                    for(v_MMF=0; v_MMF < (2*v_max); v_MMF++){  

                        printf(" %f ", creal(matrix_MMF[n_MMF][v_MMF]));

                    }
                    printf(" \n\n ");
                } */

                        
                        
                        /* SOLVE LINEAR EQUATIONS USING CONJUGATE GRADIENT METHOD */

                // Here we aim to find unknown values 'd' from the linear matrix equation:       Ad = F     ( A * d_terms = final_col)

                n = 2*n_max; // set the numnber of equations to solve

                if(speed_test==1){
                    current_time= omp_get_wtime();
                    fprintf(timefile, " %12f   Solving Linear equations for %d terms...\n", current_time-start_time, n);
                }

                /* initialise column vectors "d_terms" (unknowns) and "final_col" (simply filled with constants) */
                d_terms = (double _Complex*)malloc(2*n_max*sizeof(double _Complex));  // initialise array to store d_terms
                final_col = (double _Complex*)malloc(2*n_max*sizeof(double _Complex));  // initialise array to store final column (initial solution to the matrix equation, where all terms are '1/(N-1)'

                // set initial values for final column and initialise d_terms
                for(i=0; i<n; i++){
                    d_terms[i] = 0;
                    final_col[i]= 1.0/(N_monomers-1.0);
                }

                //save values of original MMF matrix for debugging/checking solution has worked

                original_matrix_MMF=(double _Complex**)malloc(2*n_max*sizeof(double _Complex*));
                for (i=0; i<(2*n_max); i++) {
                    original_matrix_MMF[i]=(double _Complex*)malloc(2*n_max*sizeof(double _Complex));  //makes an n_max x n_max 2D matrix to store the terms in the MMF matrix
                }    

                for(n_MMF=0; n_MMF < (2*n_max) ; n_MMF++){
                    for(v_MMF=0; v_MMF < (2*v_max); v_MMF++){  
                        original_matrix_MMF[n_MMF][v_MMF] = matrix_MMF[n_MMF][v_MMF];
                    }
                }

                /* Solve equations using Gaussian elimination */

                //printf("\n\n\n--------------------------------------------------------------------CONVERTING TO UPPER TRIANGULAR MATRIX------------------------------------------------------------------------------- \n\n\n");

                /* Create upper triangular matrix */

                for(k=0;k<(n-1);k++){ // k represents the top row of the reduced matrix

                    #pragma omp parallel for
                    for(int p=(k+1);p<n;p++){ // p represents the row that we are comparing to the top row of the reduced matrix (declare within parallel section to keep private)

                        double _Complex coeff= matrix_MMF[p][k] / matrix_MMF[k][k]; //find coefficient (declare within parallel section to keep private)
                        //printf(" \n THREAD %d: top row k= %d, manipulated row p = %d, coeff = %f", omp_get_thread_num(), k, p, coeff);

                        for(int t=0;t<n;t++){ //t represents the column number (declare within parallel section to keep private)
                            matrix_MMF[p][t]= matrix_MMF[p][t] - matrix_MMF[k][t]*coeff; //subtract "top row x coefficient" from each of the elements of the current row
                        }

                        /*adjust E matrix */
                        final_col[p]= final_col[p] - final_col[k]*coeff;
                    }
                    #pragma omp barrier
                    //printf("\n Row %d complete...", k+2);
                }

                /* Use back substitution to find elements of P */

                //printf("\n\n\n\n -------------------------------------------------------------------------------- BACK SUBSTITUTION --------------------------------------------------------------------------------\n\n\n");
                for(i=0;i<n;i++){
                    double _Complex sum=0;
                    for(j=0;j<i;j++){
                        sum = sum + matrix_MMF[n-(i+1)][n-(j+1)]*d_terms[n-(j+1)];
                        //printf("\n i= %d, j= %d, sum= %f, n-(i+1)= %d, n-(j+1)= %d", i, j, sum, n-(i+1), n-(j+1));
                    }
                    d_terms[n-(i+1)]= (final_col[n-(i+1)] - sum) / matrix_MMF[n-(i+1)][n-(i+1)];
                    //printf("\n Back substitution for row %d complete...", n-i);
                }

                //printf("\n\n -----------------------------------------------------------------------------------    FINAL MATRIX    -----------------------------------------------------------------------------------\n\n\n");
                
                if(speed_test==1){
                    current_time= omp_get_wtime();
                    fprintf(timefile, " %12f   Linear equations solved. Calculating cross-sections...\n", current_time-start_time);
                }
                
                /* check MARBLES solution has worked by multiplying out Ad = E (Important: be sure to use the ORIGINAL MMF matrix!) and seeing what it equals
                for(i=0;i<n;i++){ // moves DOWN one row in the column
                    complex_temp=0;
                    for(j=0;j<n;j++){ // go ACROSS each element one at a time 
                        
                        complex_temp += original_matrix_MMF[i][j]*d_terms[j];

                    }
                    printf(" \n Multiplication of A*d for row [%3d] = %14g + %14g i          (it should equal %14g + 0 i)", i, creal(complex_temp), cimag(complex_temp), 1.0/(N_monomers-1.0));
                } */

                /* print d terms 
                for(i=0; i<n; i++){
                    printf("\n d_term[%d] = %g + %g i ", i, creal(d_terms[i]), cimag(d_terms[i]));
                } */
                
                /* store d_terms into two seperate arrays for easy reading */

                for(i=0; i<n_max; i++){
                    d_1_terms[i+1] = d_terms[i];  //save all d_1 terms in their own array
                    d_2_terms[i+1] = d_terms[i+n_max];   //save all d_2 terms in their own array
                }
            
                /* free matrices used in this conditional (whenever N>1 monomers) */

                for (i=0; i<(2*n_max); i++) {
                    free((void*)original_matrix_MMF[i]);
                }
                free((void*)original_matrix_MMF);

                free((void*)d_terms);
                free((void*)final_col);
            }

            /* print an and dn terms in their new arrays
            printf("\n\n a terms:  \n");
            for(i=1; i<=n_max; i++){
                printf("\n\t a_n_terms[%d] = %12g + %12g i         b_n_terms[%d] = %12g + %12g i ", i, creal(a_n_terms[i]), cimag(a_n_terms[i]), i, creal(b_n_terms[i]), cimag(b_n_terms[i]));
            }

            printf("\n\n d terms:  \n");
            for(i=1; i<=n_max; i++){
                printf("\n\t d_1_terms[%d] = %12g + %12g i         d_2_terms[%d] = %12g + %12g i ", i, creal(d_1_terms[i]), cimag(d_1_terms[i]), i, creal(d_2_terms[i]), cimag(d_2_terms[i]));
            } */



            /* Calculate MMF Extinction cross-section using:
            
        !                    n_max
        !            2*pi*N  ---
        !   C_ext =  ------  \    (2n+1) Re ( d_1(n) + d_2(n) )
        !             k^2    /--
        !                    n=1
            
            */

            sum_extinction=0;
            for(n_MMF=1; n_MMF<=n_max; n_MMF++){
                sum_extinction += (2.0*n_MMF+1.0) * creal(d_1_terms[n_MMF] + d_2_terms[n_MMF]); 
            } 
            
            MMF_C_ext = (2*pi/pow(wavenumber_k,2))*N_monomers*sum_extinction; // extinction cross-section
            
            /* Calculate Absorption cross-section using RDG theory (simply absorption by individual monomers X number of monomers) to compare to C_abs_MFT later (this is the MODIFIED section of MODIFIED mean field theory) :
            

        !                    n_max
        !            2*pi*N  ---                      1                    1
        !   C_abs =  ------  \    (2n+1) Re{|an|^2 ( --- - 1 ) + |bn|^2 ( --- - 1 )}
        !             k^2    /--                     an*                  bn*
        !                    n=1
        !   */

            sum_absorption=0;
            for(n_MMF=1; n_MMF<=n_max; n_MMF++){
                sum_absorption += (2.0*n_MMF+1.0) * creal( pow(cmodulus(a_n_terms[n_MMF]),2)*( (1.0/cconjugate(a_n_terms[n_MMF])) - 1.0 ) +  pow(cmodulus(b_n_terms[n_MMF]),2)*( (1.0/cconjugate(b_n_terms[n_MMF])) - 1.0 ) );
            }

            RDG_C_abs = (2*pi/pow(wavenumber_k,2))*N_monomers*sum_absorption; // calculation of absorption cross-section from RDG Theory

                /* to check this value, calculate MFT scattering cross-section */

            // calculate S11 and S(q) for each angle between 0 and 180 degrees

            for(i=0; i<=180; i++){ //z repesents angle in degrees (explore 0 -> 180 degrees)

                scatter_angle = i*pi/180.0; //covert i from degrees to radians

                //printf ("\n i = %d    angle = %f", i, scatter_angle);
                
                /* find and store the pi and tau functions for this scattering angle (stored in the arrays "pi_function" and "tau_function")*/
                calculate_pi_and_tau_function_terms(n_max, scatter_angle, pi_function, tau_function);

                //calculate S11 scattering term for this angle and store in an array with 181 elements (one for each angle between 0-> 180 degrees)
                S11[i] = calculate_S11(n_max, d_1_terms, d_2_terms, pi_function, tau_function);
                
                //calculate q and integrate structure factor for this angle
                q_coeff = 2.0*wavenumber_k*sin(scatter_angle/2.0);
                x_val=0; // initialise this value to zero (it will be changed within the functions)
                S_q[i] = calculate_structure_factor(d_f, q_coeff, radius_of_gyration);
            }

            /* key values to print if debugging 
            for(i=0; i<=180; i++){ //z repesents angle in degrees (explore 0 -> 180 degrees)

                scatter_angle = i*pi/180.0; //covert i from degrees to radians
                q_coeff = 2.0*wavenumber_k*sin(scatter_angle/2.0);

                //printf("\n S11[%d] = %f    S(%f)[%d] = %f", i, S11[i], q_coeff, i, S_q[i]);
            } */



                /* use simpson's rule to integrate C_sca equation */

            dx_scatter = 2.0*pi/180.0; //set dx as 2 degrees in radians
            integrated_C_sca=0; //initialise integral sum

            //begin calculation to integrate structure factor function

            for(i=0; i<179; i+=2){ //i repesents angle in degrees. IMPORTANT - the step size is two degrees, so increase the angle by two degrees each time we iterate, using i+=2. The only reason we calculate the middle angles are because they are useful for step c. Also, we only iterate up to 178, because there is a [i+2] element below (and our max is 180 degrees)

                scatter_angle = i*pi/180.0; //covert i from degrees to radians

                //printf("\n i = %d, Smat_1 (optools version) = %f", i, N_monomers * S11[i] * (1.0 + (N_monomers - 1.0) * S_q[i]));

                //find area underneath curve between points a and b, which are every 2 degrees. This rule is in the same format as the other integrals, step_a = 0 point, step_c = 1 degree higher (where 1 degree in radians = dx_scatter/2) and step_b = 2 degrees higher (2 degrees = dx_scatter)
                double step_area = (dx_scatter/6.0) * (N_monomers*S11[i]*(1.0 + (N_monomers - 1.0)*S_q[i])*sin(scatter_angle) + 4.0*N_monomers*S11[i+1]*(1.0 + (N_monomers - 1.0)*S_q[i+1])*sin(scatter_angle+(dx_scatter/2.0)) + N_monomers*S11[i+2]*(1.0 + (N_monomers - 1.0)*S_q[i+2])*sin(scatter_angle+dx_scatter)); // Use Simpson's rule to estimate the area of S11*S(q)*sin(theta) under the step between a --> b. And yep: the order is a, c, b!

                integrated_C_sca += step_area; //add step area to total area under curve
            }

            //complete the calculation by using the integrated part to find C_sca below
            MFT_C_scatter = (2.0*pi/pow(wavenumber_k,2.0))*integrated_C_sca;

            printf("\n\n For %f um, integrated C_sca = %f and RDG C_scatter = %f \n", wavelength, integrated_C_sca, MFT_C_scatter);
            
            printf("\n Initial values:     C_ext = %f     MFT_C_sca = %f     RDG_C_abs = %f", MMF_C_ext, MFT_C_scatter, RDG_C_abs);


                /* Calculate Asymmetry parameter */
            
            for(i=0; i<=180; i++){ //z repesents angle in degrees (explore 0 -> 180 degrees)
                phase_function[i] = N_monomers*S11[i]*(1.0 + (N_monomers - 1.0)*S_q[i]) / (MFT_C_scatter*pow(wavenumber_k,2)); //calculate the phase function - remember to use S11,agg here, given by:   N_monomers*S11[i]*(1.0 + (N_monomers - 1.0)*S_q[i])
                //printf("\n i= %d        S11_agg = %f         phase_function = %12f ", i, N_monomers*S11[i]*(1.0 + (N_monomers - 1.0)*S_q[i]), phase_function[i]);
            }

            g_asymmetry_MMF= find_g_MMF(phase_function, pi); //calculate asymmetry parameter g

            printf(" \n\n\n g_asymmetry= %f", g_asymmetry_MMF);



                    /* MODIFIED MEAN FIELD CALCULATIONS - adjust C_Abs */

            // Calculate tau value
            tau_MMF = RDG_C_abs/G_MMF;

            printf("\n tau_MMF = %f", tau_MMF);

            MFT_C_abs = MMF_C_ext-MFT_C_scatter; //calculate Mean Field Theory version of C_abs

            MMF_C_abs = maximum_C_abs(MFT_C_abs, G_MMF, tau_MMF); //compare MFT C_abs with modified term -- use whichever is larger

            printf("\n MMF_C_abs = %f", MMF_C_abs);

            /* Finally, calculate MMF scatter from our original C_Ext and our corrected C_Abs values */

            MMF_C_scatter = MMF_C_ext - MMF_C_abs;


            // Calculate Q efficiencies
            MMF_Q_extinction = MMF_C_ext/mie_aerosol_area;
            MMF_Q_scatter = MMF_C_scatter/mie_aerosol_area;
            MMF_Q_absorption = MMF_C_abs/mie_aerosol_area;

            printf("\n ------------------------------- FOR WAVELENGTH %5.3f um ----------------------------------------------------------", wavelength);
            printf("\n\n Q EFFICIENCIES: \n\n                Q (Extinction)          Q (Scatter)          Q (Absorption) \n");
            printf("\n     Mie \t  %f \t         %f \t       %f \n",  mie_results[wavelength_index][0], mie_results[wavelength_index][1], mie_results[wavelength_index][2]);
            printf("\n     MMF \t  %f \t         %f \t       %f \n",  MMF_Q_extinction, MMF_Q_scatter, MMF_Q_absorption);
            printf("\n -------------------------------------------------------------------------------------------------------------------");
            printf("\n\n CROSS-SECTIONS (um^2): \n\n                C (Extinction)          C (Scatter)          C (Absorption) \n");
            printf("\n     Mie \t  %f \t         %f \t       %f \n",  mie_results[wavelength_index][3], mie_results[wavelength_index][4], mie_results[wavelength_index][5]);
            printf("\n     MMF \t  %f \t         %f \t       %f \n",  MMF_C_ext, MMF_C_scatter, MMF_C_abs);
            printf("\n ------------------------------------------------------------------------------------------------------------------- \n\n");

            MMF_results[wavelength_index][0]=MMF_Q_extinction;
            MMF_results[wavelength_index][1]=MMF_Q_scatter;
            MMF_results[wavelength_index][2]=MMF_Q_absorption;
            MMF_results[wavelength_index][3]=MMF_C_ext;
            MMF_results[wavelength_index][4]=MMF_C_scatter;
            MMF_results[wavelength_index][5]=MMF_C_abs;
            MMF_results[wavelength_index][6]=g_asymmetry_MMF;


            /* CODE TO PRODUCE MFT INSTEAD OF MMF - NOT NORMALLY USED, BUT IF SO, IT WILL BE SAVED IN THE MFT_results.txt FILE! */


            MMF_Q_extinction = MMF_C_ext/mie_aerosol_area; //note - MMF_ext is the same as MFT_ext!
            MFT_Q_scatter = MFT_C_scatter/mie_aerosol_area;
            MFT_Q_absorption = MFT_C_abs/mie_aerosol_area;

            MFT_results[wavelength_index][0]=MMF_Q_extinction; //note - MMF_ext is the same as MFT_ext!
            MFT_results[wavelength_index][1]=MFT_Q_scatter;
            MFT_results[wavelength_index][2]=MFT_Q_absorption;
            MFT_results[wavelength_index][3]=MMF_C_ext; //note - MMF_ext is the same as MFT_ext!
            MFT_results[wavelength_index][4]=MFT_C_scatter;
            MFT_results[wavelength_index][5]=MFT_C_abs;
            MFT_results[wavelength_index][6]=g_asymmetry_MMF;

            /* free memory for MMF matrix; do this every time the wavelength or size parameter changes */

            for (i=0; i<(2*n_max); i++) {
                free((void*)matrix_MMF[i]);
            }
            free((void*)matrix_MMF);

            free((void*)a_n_terms);
            free((void*)b_n_terms);

            free((void*)pi_function);
            free((void*)tau_function);

            free((void*)d_1_terms);
            free((void*)d_2_terms);

            free((void*)structure_factor_terms);
        } //increase wavelength and repeat

        /* save MMF results for all wavelengths analysed at this radius */
        for(i=0;i<num_wavelengths;i++){
            fprintf(mmfresultsfile, "    %10e      %10e      %10e     %10e    %10e      %10e    %10e    %10e    %10e\n", aerosol_radius, wavelength_list[i], MMF_results[i][0], MMF_results[i][1], MMF_results[i][2], MMF_results[i][3], MMF_results[i][4], MMF_results[i][5], MMF_results[i][6]);
        }

        /* save MFT results for all wavelengths analysed at this radius */
        for(i=0;i<num_wavelengths;i++){
            fprintf(mftresultsfile, "    %10e      %10e      %10e     %10e    %10e      %10e    %10e    %10e    %10e\n", aerosol_radius, wavelength_list[i], MFT_results[i][0], MFT_results[i][1], MFT_results[i][2], MFT_results[i][3], MFT_results[i][4], MFT_results[i][5], MFT_results[i][6]);
        }


        //free all matrices at the end of the MMF program (outside of the wavelength loop)

        for (i=0; i<(largest_p_max+1); i++) {
            free((void*)LegendreP_at_root[i]);
        }
        free((void*)LegendreP_at_root);

        for (i=0; i<(largest_p_max+1); i++) {
            free((void*)diff_LegendreP_at_root[i]);
        }
        free((void*)diff_LegendreP_at_root);

        for (i=0; i<(largest_n_max+1); i++) {
            free((void*)Associated_LegendreP_at_root[i]);
        }
        free((void*)Associated_LegendreP_at_root);
    
    } // return to top of loop and analyse next radius 

    /* close results files */
    fclose(mmfresultsfile);
    //fclose(mftresultsfile);

    
    //free Gauss-Legendre arrays (calculated outside of the radius loop above)
    free((void*)GL_roots);
    free((void*)GL_weights);

    /* find total time that program took to run */
    end_time= omp_get_wtime();
    if(speed_test==1){
        fprintf(timefile, " %12f   Program finished!\n", end_time-start_time);
    }

    printf("\n\n PROCESS COMPLETE - SUCCESS!\n\n");

    printf(" Program took %f secs to complete.\n\n", end_time-start_time);


    for (i=0; i<num_wavelengths; i++) {
        free((void*)mie_results[i]);
    }
    free((void*)mie_results);

    for (i=0; i<num_wavelengths; i++) {
        free((void*)MMF_results[i]);
    }
    free((void*)MMF_results);

    for (i=0; i<num_wavelengths; i++) {
        free((void*)MFT_results[i]);
    }
    free((void*)MFT_results);

    /* free matrice that store refractive index data */

    for (i=0; i<number_ref_ind_data_points; i++) {
        free((void*)refractive_index_raw_data[i]);
    }
    free((void*)refractive_index_raw_data);

    for (i=0; i<number_ref_ind_data_points; i++) {
        free((void*)refractive_index_data[i]);
    }
    free((void*)refractive_index_data);

    /* free other matrices */
    
    free((void*)wavelength_list);
    free((void*)radius_list);

    if(speed_test==1){
        fclose(timefile);
    }

    return 0;
}

