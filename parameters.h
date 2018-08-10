#pragma once

static const double PI= 3.14159265359;

//###########################################################################################################//
//                                Choose the lattice volume with T, X, Y and Z:                              //
//###########################################################################################################//

#define T 8
#define X 6
#define Y 6
#define Z 6

/*
#define T 48
#define X 14 //24
#define Y 14 //24
#define Z 14 //24
*/
 
#define volume T*X*Y*Z

// Parameters of theory in the continuum:
#define m2_0 -4.9  //(A)
#define lambda_c 10.0 //(A)

// Parameters of the action S:
extern double LAMBDA;
extern double KAPPA;
 
// Parameters needed in the Metropolis algorithm. In each step one makes one update in the magnitude
// and phase of the field
static const double deltaphi = 3.14159/16; //3.14159265359/8;
static const double deltarho = 1;//1; //0.45;

#define nprint_field 1000
#define n_update_phi 1


//###########################################################################################################//
//      In main.cpp "start_random_conf" is needed in (II.J) in order to print a specific configuration:      //
//     "fprint_field(phi, (long long) (i+1)*n_term_save + (long long) start_random_conf*(1-start_random))"   //
//###########################################################################################################//

#define start_random_conf 0


//###########################################################################################################//
//       "n_save" is the number of saved phi field configurations called "scalar_X_Y_Z_T_n_(conf).txt":      //
//###########################################################################################################//

#define n_save 20


//###########################################################################################################//
//          If start_random 0: The the file "start_conf" is chosen to be the start configuration.            //
//          If start_random 1: Hot start, where a disordered, random configuration is utilized.              //
//###########################################################################################################//

#define start_random 0


//###########################################################################################################//
//      In case of "start_random 0" one has to choose a start configuration file which must be included      //
//   in "start_config". If "start_conf" is available, the program continues with the Metropolis algorithm:   //
//###########################################################################################################//

static const char *start_conf = "./start_config/scalar_6_6_6_8_6000.txt";


//###########################################################################################################//
//    HOT START (start_random 1):                                                                            //
//    If we start with a disordered, random configuration (start_random 1), n_term_field is the number of    //
//    steps until thermalization:                                                                            //
//###########################################################################################################//

#define n_term_field 1000


//###########################################################################################################//
//    COLD START:                                                                                            //
//    Maybe later...                                                                                         //
//###########################################################################################################//


//###########################################################################################################//
//               Number of local updates in a single Metropolis step (see metropolis.cpp):                   //
//###########################################################################################################//

#define n_metropolis  250*10*4//650 //for L=18


//###########################################################################################################//
//                   After n_term_save Metropolis steps a field configuration is saved, so                   //
//    in "main.c" the acceptance = metropolis() is calculated for n_field = n_term_save and the endings of   //
//               the field configurations computed in "main.c" are multiplied with n_term_save:              //
//###########################################################################################################//

#define n_term_save 1*1000 //10000


//###########################################################################################################//
//  In "main.c" this is the path where the folder including the scalar field configurations phi is created,  //
//  in "calculate_corr.c" this path is used to read in the fields for the creation of correlation functions: //
//###########################################################################################################//

static const char *path_read = "./field_configs/";


//###########################################################################################################//
//      In "main.c" and "calculate_corr.c" (II.A) this is the path where the folder "analysis" (in which     //
//        the correlation functions will be printed) is created (the path must end with "/name-of-folder/"): //
//###########################################################################################################//

static const char *path_corr = "./corr_analysis/";


//###########################################################################################################//
//          This parameter determines the particle number n of the n particle correlation function           //
//                      calculated in "operators.c" and "calculate_corr.c", respectively:                    //
//###########################################################################################################//

#define n_fields 5


//###########################################################################################################//
//   Number of field configurations phi that are read in in "calculate_corr.c" to calculate the correlation  // 
//functions (number must be chosen n_analyse = n_save-1 where n_save is the number of created field configs)://
//###########################################################################################################//

#define n_analyse 20//100000 //36000


//###########################################################################################################//
//   In "calculate_corr.c" (II.G) only every (n_restrict)th correlation function is printed into the ".tsv"  // 
//  file in folder "analysis" for the 1, 2, 3, 4, 5 particle case, respectively. If n_restrict is set to 1,  //
//                                   all correlation functions are printed:                                  //
//###########################################################################################################//

#define n_restrict 1


//###########################################################################################################//
// If correlator == 0: n particle correlation function is utilized in "calculate_corr.cpp"                   //
// If correlator == 1: derivative of n particle correlation function is utilized in "calculate_corr.cpp"     //
//###########################################################################################################//

#define correlator 1

