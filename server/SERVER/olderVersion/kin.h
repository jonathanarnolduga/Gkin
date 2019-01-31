// first try of fortran -> c++ translation (by Aleman)
// this is the .def file 
// converted to a .h file

#ifndef __kin_h__
#define __kin_h__


// #define OPTIMIZE 


//
//c=================================================
//c Reaction Kinetics Network Modeling Code "kin.for"
//c Header File "kin.def"
//c=================================================
//c H.-B. Schuttler
//c Department of Physics and Astronomy
//c University of Georgia
//c Athens, GA 30602
//c E-Mail hbs@physast.uga.edu
//c
//c Version 1.1, 00-05-17
//c

      const int MAX_NUMBER_SPECIES   = 100;
      const int MAX_NUMBER_REACTIONS = 100;
      const int MAX_NUMBER_PARTICIPANT_REACTIONS = 10;
      const int MAX_NUMBER_TIME_STEPS = 50000; //20000


      char tmpLine[ 72 ] = { 0 };

      char *bline[ 10 + MAX_NUMBER_REACTIONS ] = { 0 };
      char *nameOfSpecies[ MAX_NUMBER_SPECIES + 1 ] = { 0 };               
      char *nameOfInputSpecies[ MAX_NUMBER_PARTICIPANT_REACTIONS + 1 ][ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
      char *nameOfOutputSpecies[ MAX_NUMBER_PARTICIPANT_REACTIONS + 1 ][ MAX_NUMBER_REACTIONS + 1 ] = { 0 };

      int numberOfDataSets;
      int numberOfSpecies, numberOfReactions;
      int numberOfTimeSteps;
      int integrationOption;
      int ntskip;

      int iispec[ MAX_NUMBER_PARTICIPANT_REACTIONS + 1 ][ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
      int iospec[ MAX_NUMBER_PARTICIPANT_REACTIONS + 1 ][ MAX_NUMBER_REACTIONS + 1 ] = { 0 };

      int *jfix = new int[ MAX_NUMBER_SPECIES + 1 ];
      int *numberInputParticipants = new int[ MAX_NUMBER_REACTIONS + 1 ];
      int *numberOutputParticipants = new int[ MAX_NUMBER_REACTIONS + 1 ];
      int *jkin = new int[ MAX_NUMBER_REACTIONS + 1 ];

//      int jfix[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
//      int numberInputParticipants[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
//      int numberOutputParticipants[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
//      int jkin[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };

      double initialTime, finalTime, dtime;

      double *initialConcentration = new double[ MAX_NUMBER_SPECIES + 1 ];
//      double initialConcentration[ MAX_NUMBER_SPECIES + 1 ] = { 0 };

//      double xspec[ MAX_NUMBER_SPECIES + 1 ][ MAX_NUMBER_TIME_STEPS + 1 ] = { 0 };
//      double xreac[ MAX_NUMBER_REACTIONS + 1 ][ MAX_NUMBER_TIME_STEPS + 1 ] = { 0 };
      double **xspec;
#ifndef OPTIMIZE
      double **xreac;
#endif

      double *xtime;
      int xtime_index;

//      double forwardReactionRates[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
//      double backwardReactionRates[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
//      double forwardReactionRates2[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
//      double backwardReactionRates2[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
      double *forwardReactionRates = new double[ MAX_NUMBER_REACTIONS + 1 ];
      double *backwardReactionRates = new double[ MAX_NUMBER_REACTIONS + 1 ];
      double *forwardReactionRates2 = new double[ MAX_NUMBER_REACTIONS + 1 ];
      double *backwardReactionRates2 = new double[ MAX_NUMBER_REACTIONS + 1 ];

   const char *kin_i01 = "kin.i01";
   const char *kin_o01 = "kin.o01";
   const char *kin_o02 = "kin.o02";
   const char *kin_o03 = "kin.o03";

   const double c20 = 0.25, c21 = 0.25;
   const double c30 = 0.375, c31 = 0.09375, c32 = 0.28125;
   const double c40 = 12.0 / 13.0, c41 = 1932.0 / 2197.0;
   const double c42 = -7200.0 / 2197.0, c43 = 7296.0 / 2197.0;
   const double c51 = 439.0 / 216.0, c52 = -8.0;
   const double c53 = 3680.0 / 513.0, c54 = -845.0 / 4104.0;
   const double c60 = 0.5, c61 = -8.0 / 27.0, c62 = 2.0;
   const double c63 = -3544.0 / 2565.0, c64 = 1859.0 / 4104.0;
   const double c65 = -0.275;
   const double a1 = 25.0 / 216.0, a2 = 0.0, a3 = 1408.0 / 2565.0;
   const double a4 = 2197.0 / 4104.0, a5 = -0.2;
   const double b1 = 16.0 / 135.0, b2 = 0.0, b3 = 6656.0 / 12825.0;
   const double b4 = 28561.0 / 56430.0, b5 = -0.18, b6 = 2.0 / 55.0;

   struct Evaluate {
      int       participants;
      double    multiplier;
      int      *indexes;
      int       appears;
      Evaluate *next;
   };


// functions:

   int  readInputData();
   void setfix();
   void setnet();
   void runkin();
   int outputDataFile1();
   int outputDataFile2();
   int outputDataFile3();
   char* trim( const char *str );
   int ispec4name( const char *name2look );
   void xrate( double vspec[], double vfor[], double vbak[], int t );

   double myabs( double number );
   void freeMemory();
   void gauss( double **a, int d[], int n );
   void gauss_solve( double **a, int d[], double b[], double x[], int n );
   void prepareJacobian( Evaluate *** jac );
   void eval_prep_jacobian( int time, double ** jac, Evaluate *** prepared_jac);



#endif

