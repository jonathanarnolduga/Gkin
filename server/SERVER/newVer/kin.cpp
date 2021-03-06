/*
c========================================
c Reaction Kinetics Network Modeling Code
c========================================
c 
c Yihai Yu
c
c 04-11-10
c   Initialze the big matrix takes about 10 sec, which
c   is the most of the time, and for water, integration
c   is just 0.4 sec including multiple integrations for
c   external control. But we only need to declare the
c   the biggest matrix once for MINE code, and reuse it.
c   But, we need to delete the malloc matrix for multiple
c   integrations for each one complete ens.o02
c   Still sometthing wrong with RK4 method for multiple
c   integrations
c
c   Note: need to optimize number of steps for pulses for
c         external control
c
c 04-11-09
c   Finished integrating multiple integrations, 
c   Note, only outputDataFile2() are updated
c
c 04-11-09
c   Now, we can slice the time into multiple integrations, 
c   Note, this version hasn't incorporated all integrations.
c
c 04-09-30
c   Add external control for MINE
c
c 04-01-10: 
c   Add lsodes as option 6
c
*/

/*
c
c Aleman
c
c 00-08-01:
c   Added time evol. option "jtime=2": Modified Euler alg.
c   (orig.: "jtime=1": Standard Euler alg.)
c 00-07-29:
c   Added input parameter "ntskip"
c   to skip time steps on output
c 00-05-17:
c   Changed output file protocol
c   to overwrite pre-existing "kin.o0*"
c 00-05-13: 
c   Added warn./stop for excessive array size
c 00-03-26:
c   Added 2nd + 3rd output file for graphics
c 00-03-24: 
c   Added MM kinetics option
c
*/

#include "kin.h"
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <string.h>
#include <time.h>

#ifndef DEBUG_TIMESTEP
#define DEBUG_TIMESTEP 0
#endif

#define DEBUG 0
#define DEBUG_XRATE 0
#define DEBUG_EXECUTION_TIME 1
#define LINE 72
#define DEC8 8
#define FRAC 16


/*
*
**
* Initializes some of the 'big' matrices used
*/
int initializeBigMatrices()
{
xspec = new double* [ MAX_NUMBER_TIME_STEPS + 10 ];
if( xspec == 0 ) {
cerr << "ERROR: Unable to allocate memory for xspec!" << endl;
return 1;
}
for( int i = 0; i < ( MAX_NUMBER_TIME_STEPS + 10 ); i++ ) {
xspec[ i ] = new double[ MAX_NUMBER_SPECIES + 1 ];
if( xspec[ i ] == 0 ) {
cerr << "ERROR: Unable to allocate memory for xspec[ " << i;
cerr << " ]!" << endl;
return 1;
}
}

#ifndef OPTIMIZE
//xreac = new (double*) [ MAX_NUMBER_REACTIONS + 1 ];
xreac = new (double*) [ MAX_NUMBER_TIME_STEPS + 10 ];
if( xreac == 0 ) {
cerr << "ERROR: Unable to allocate memory for xreac!" << endl;
return 1;
}
//for( int i = 0; i < ( MAX_NUMBER_REACTIONS + 1 ); i++ ) {
for( int i = 0; i < ( MAX_NUMBER_TIME_STEPS + 10 ); i++ ) {
//xreac[ i ] = new double[ MAX_NUMBER_TIME_STEPS + 10 ];
xreac[ i ] = new double[ MAX_NUMBER_REACTIONS + 1 ];
if( xreac[ i ] == 0 ) {
cerr << "ERROR: Unable to allocate memory for xreac[ " << i;
cerr << " ]!" << endl;
return 1;
}
}
#endif

xtime = new double[ MAX_NUMBER_TIME_STEPS + 10 ];
if( xtime == 0 ) {
cerr << "ERROR: Unable to allocate memory for xtime !" << endl;
return 1;
}

return 0;
}

void setnumofintegrations(){
double remaintime;
if(!extFlag){
return;
}//endof if
else{
numofcycle = (int)((true_finalTime - true_initialTime - extOfSpecies[whichExt][1])/oneCycle[whichExt]);
//cout << numofcycle <<endl;
numOfIntegrations = numofcycle*npulse[whichExt] + 1;
//cout << numOfIntegrations <<endl;
remaintime = true_finalTime - true_initialTime - oneCycle[whichExt]*numofcycle;
for(int i=3; i<=2*npulse[whichExt]+2; i=i+2){
numOfIntegrations++;
remaintime = remaintime - extOfSpecies[whichExt][i];
if(remaintime<0.0) break;
}//endof while
}//endof else
//cout << numOfIntegrations <<endl;
}

void setintg_(){
intg_initialTime = new double[numOfIntegrations+1];
intg_finalTime = new double[numOfIntegrations+1];
intg_numberOfTimeSteps = new int[numOfIntegrations+1];
intg_xtime_index = new int[numOfIntegrations+1];
intg_dtime = new double[numOfIntegrations+1];

intg_xspec = new double**[numOfIntegrations+1];
intg_xtime = new double*[numOfIntegrations+1];

if(!extFlag){
intg_initialTime[1] = true_initialTime;
intg_finalTime[1] = true_finalTime;
//todo for various time steps
intg_numberOfTimeSteps[1] = true_numberOfTimeSteps;
intg_dtime[1] = (intg_finalTime[1]-intg_initialTime[1])/intg_numberOfTimeSteps[1];
return;
}
else{

intg_initialTime[1] = true_initialTime;
intg_finalTime[1] = true_initialTime + extOfSpecies[whichExt][1];
//todo for various time steps
intg_numberOfTimeSteps[1] = true_numberOfTimeSteps;
intg_dtime[1] = (intg_finalTime[1]-intg_initialTime[1])/intg_numberOfTimeSteps[1];

for(int i=2; i<=numOfIntegrations; i++){
intg_initialTime[i] = intg_finalTime[i-1];
intg_finalTime[i] = intg_initialTime[i] + extOfSpecies[whichExt][2*((i-1)-((i-2)/npulse[whichExt])*npulse[whichExt])+1];
if(intg_finalTime[i]>true_finalTime){
intg_finalTime[i] = true_finalTime;
}//endof if
//cout << "final time " << true_finalTime <<endl;
intg_numberOfTimeSteps[i] = true_numberOfTimeSteps;         
//cout << intg_numberOfTimeSteps[i] <<endl;
intg_dtime[i] = (intg_finalTime[i] - intg_initialTime[i])/intg_numberOfTimeSteps[i];
}//endof for
}//endof else

/*
for(int i=1; i<=numOfIntegrations; i++){
cout << intg_initialTime[i] << " " << intg_finalTime[i] << " " << intg_dtime[i] <<endl;
}
*/

}

void setinitial(int i){
initialTime = intg_initialTime[i];
}

void setfinal(int i){
finalTime = intg_finalTime[i];
}

void setfix(int i)
{
//dtime = ( finalTime - initialTime ) / ( 1.0 * numberOfTimeSteps );
numberOfTimeSteps = intg_numberOfTimeSteps[i];
dtime = intg_dtime[i];   
}

void setinitialdata(int i){

if(i==1){
for( int i = 1; i <= numberOfSpecies; i++ ) {
xspec[ 0 ][ i ] = initialConcentration[ i ];
}//endof for
}//endof if
else{
int top = numberOfTimeSteps;
if( integrationOption == 4 ) { // RK Adaptive
top = xtime_index;
}

for( int i = 1; i <= numberOfSpecies; i++ ) {
xspec[ 0 ][ i ] = xspec[ top ][ i ];
}//endof for
}//endof else

}

void storedata(int i){

int top = numberOfTimeSteps;
if( integrationOption == 4 ) { // RK Adaptive
top = xtime_index;
intg_xtime_index[i] = top;
intg_xtime[i] = new double[ top+1 ];
}//endof if

intg_xspec[i] = new double* [ top + 1 ];
if( intg_xspec[i] == 0 ) {
cerr << "ERROR: Unable to allocate memory for intg_xspec!" << endl;  
//return 1;
}//endof if
for( int j = 0; j < ( top + 1 ); j++ ) {
intg_xspec[ i ][ j ] = new double[ numberOfSpecies + 1 ];
if( intg_xspec[ i ][ j ] == 0 ) {
cerr << "ERROR: Unable to allocate memory for intg_xspec[ " << i;
cerr << " ]! " << endl;
//return 1;
}
for(int k=1; k<=numberOfSpecies; k++){
intg_xspec[ i ][ j ][ k ] = xspec[ j ][ k ];
}//endof for
if(integrationOption == 4){
intg_xtime[ i ][ j ] = xtime[ j ];
}//endof if
}//endof for

}

// solves kinetic equations for arbitrary reaction networks
int main(int argc, char **argv)
{

lsodesArgv = argv[1];
outputArgv = argv[2];
includeArgv = argv[3]; 

clock_t start, end;
double elapsed = 0.0;

int retValue = 0;

numberOfDataSets = 0;

start = clock();
if( initializeBigMatrices() != 0 ) {
return 1;
}
end = clock();
elapsed = ( (double) ( end - start ) ) / CLOCKS_PER_SEC;
if( DEBUG_EXECUTION_TIME ) {
cerr << " initializeBigMatrices time: " << elapsed <<endl;
}

do {
start = clock();
retValue = readInputData();

if( retValue != 0 ) {
break;
}

end = clock();
elapsed = ( (double) ( end - start ) ) / CLOCKS_PER_SEC;
if( DEBUG_EXECUTION_TIME ) {
cerr << " readInputData time: " << elapsed <<endl;
}

setnet();
setnumofintegrations();
setintg_();

start = clock();

for(int i=1; i<=numOfIntegrations; i++){
cout << i << endl;
setinitial(i);
setfinal(i);
setfix(i);
setinitialdata(i);
runkin();
storedata(i);
}

end = clock();
elapsed = ( (double) ( end - start ) ) / CLOCKS_PER_SEC;
if( DEBUG_EXECUTION_TIME ) {
cerr << " total integration time:  " << elapsed <<endl;
}

//for lsodesMethod
if( integrationOption == 6 ){
//just stop
}//endof if
//for all other methods

else{
retValue = 1;//stop it

retValue = outputDataFile1();

if( retValue != 0 ) {
cerr << "error at outputDataFile1()" << endl;
break;
}

retValue = outputDataFile2();

if( retValue != 0 ) {
cerr << "error at outputDataFile2()" << endl;
break;
}

retValue = outputDataFile3();

if( retValue != 0 ) {
cerr << "error at outputDataFile3()" << endl;
break;
}

}
} while( retValue == 0 );

freeMemory();
}

/**
* Reads input data from "kin.i01" and writes them into "kin.o01"
* @return  zero if everything was ok
*/
int readInputData()
{

int iset = 0;
char buf[ LINE ];

// open files kin_i01 and kin_o01 
ifstream inputFile( kin_i01, ios::in );
if( ! inputFile ) {
cerr << "Couldn't open input file " << kin_i01 << endl;
return 1;
}

ofstream outputFile1( kin_o01, ( numberOfDataSets == 0 ? ios::out : ios::app ) );
if( ! outputFile1 ) {
cerr << "Couldn't open input file " << kin_i01 << endl;
return 1;
}

// find first/next input data set:
int dataSet = 0;
while( dataSet <= numberOfDataSets ) {
inputFile.getline( buf, LINE );
if( strstr( buf, "data set" ) != NULL ) {
dataSet++;
}
if( inputFile.eof() ) {
break;
}
}

// verify that a data set was actually read
if( dataSet <= numberOfDataSets ) {
return 1;
}
if( DEBUG ) {
cerr << "dataSet = " << dataSet << endl;
}
numberOfDataSets = dataSet;

// read in input parameters:
int copiedLines = 0;
bline[ copiedLines ] = new char[ strlen( buf ) + 1 ];
strcpy( bline[ copiedLines ], buf );

// number of species and reactions:
copiedLines++;
inputFile.getline( buf, LINE );
bline[ copiedLines ] = new char[ strlen( buf ) + 1 ];
strcpy( bline[ copiedLines ], buf );
inputFile >> numberOfSpecies;
inputFile >> numberOfReactions;
inputFile.getline( buf, LINE );  // skip the new line
if( DEBUG ) {
cerr << "numberOfSpecies=" << numberOfSpecies;
cerr << "  numberOfReactions=" << numberOfReactions << endl;
}

if( numberOfSpecies > MAX_NUMBER_SPECIES ) { 
cerr << endl << endl << " ERROR" << endl;
cerr << "   Increase array size MAX_NUMBER_SPECIES" << endl;
cerr << "  EXECUTION TERMINATED" << endl;
return 1;
}

if( numberOfReactions > MAX_NUMBER_REACTIONS ) { 
cerr << endl << endl << " ERROR" << endl;
cerr << "   Increase array size MAX_NUMBER_SPECIES" << endl;
cerr << "  EXECUTION TERMINATED" << endl;
return 1;
}

// initial + final time, number of time steps, time integrtn. option
copiedLines++;
inputFile.getline( buf, LINE );
bline[ copiedLines ] = new char[ strlen( buf ) + 1 ];
strcpy( bline[ copiedLines ], buf );
inputFile >> initialTime;
true_initialTime = initialTime;
inputFile >> finalTime;
true_finalTime = finalTime;
inputFile >> numberOfTimeSteps;
true_numberOfTimeSteps = numberOfTimeSteps;
inputFile >> ntskip;
inputFile >> integrationOption;
inputFile.getline( buf, LINE );  // skip the new line
if( DEBUG ) {
cerr << "initialTime=" << initialTime << "  finalTime=";
cerr << finalTime << "  numberOfTimeSteps=" << numberOfTimeSteps;
cerr << "  integrationOption=" << integrationOption << endl;
}

double h = 0.0;
double hMin = 0.0;
double errorBound = 0.0;
h = ( finalTime - initialTime ) / numberOfTimeSteps;
if( numberOfTimeSteps > MAX_NUMBER_TIME_STEPS ) { 
cerr << endl << endl << " ERROR" << endl;
cerr << "   Increase array size MAX_NUMBER_TIME_STEPS" << endl;
cerr << "  EXECUTION TERMINATED" << endl;
return 1;
}

// names and initial concentration of reacting species:
copiedLines++;
inputFile.getline( buf, LINE );
bline[ copiedLines ] = new char[ strlen( buf ) + 1 ];
strcpy( bline[ copiedLines ], buf );

char *abc = new char[75];

for( int i = 1; i <= numberOfSpecies; i++ ) {
inputFile.getline( tmpLine, LINE );
//cerr << "OK " << tmpLine;
nameOfSpecies[ i ] = new char[ strlen( tmpLine ) + 1 ];
strcpy( nameOfSpecies[ i ], tmpLine );
inputFile >> initialConcentration[ i ];
inputFile >> jfix[ i ];
if((jfix[i]==10)||(jfix[i]==20)){
extFlag = true;
whichExt = i;
inputFile >> npulse[ i ];  
// process the external control information
inputFile.getline( buf, LINE ); //skip the new line
inputFile.getline( buf, LINE ); //skip one line

for(int j = 1; j <= 2*npulse[i]+2; j++){
inputFile >> abc; //actually ipm
inputFile >> extOfSpecies[i][j]; //pmspec
if((j!=1)&&(j%2==1)){
oneCycle[i] = oneCycle[i] + extOfSpecies[i][j];
slopeTrap[i][j] = (extOfSpecies[i][j+1]-extOfSpecies[i][j-1])/extOfSpecies[i][j];
}
}//endof for j

// get slopes
for(int j = 1; j <= 2*npulse[i]+2; j++){
if((j!=1)&&(j%2==1)){
slopeTrap[i][j] = (extOfSpecies[i][j+1]-extOfSpecies[i][j-1])/extOfSpecies[i][j];
//cout << extOfSpecies[i][j+1] << " " << extOfSpecies[i][j-1] << " " << extOfSpecies[i][j] <<endl;
//cout << slopeTrap[i][j] << endl;
}
}//endof for j

//cerr << " ok " << oneCycle[i];
}   
inputFile.getline( buf, LINE );  // skip the new line
if( DEBUG ) {
cerr << "nameOfSpecies[" << i << "]=" << nameOfSpecies[ i ];
cerr << "  initialConcentration=" << initialConcentration[ i ];
cerr << " jfix=" << jfix[ i ] << endl;
}
}

// reaction network by raction number:
// reaction kinetics and rates, number of reactant and product species:
//  jkin= 1: standard single-step kinetics with
//    rkfor, rkbak   = forward, backward rate constants
//  jkin=11: Michaelis-Menten (MM) two-step kinetics with
//    nameispec(1,.)=nameospec(1,.) = enzyme species
//    rkfor,  rkbak  = step-1 rate constants
//    rkfor2, rkbak2 = step-2 rate constants

for( int r = 1; r <= numberOfReactions; r++ ) {
copiedLines++;
inputFile.getline( buf, LINE );
bline[ copiedLines ] = new char[ strlen( buf ) + 1 ];
strcpy( bline[ copiedLines ], buf );
inputFile >> forwardReactionRates[ r ];
inputFile >> backwardReactionRates[ r ];
inputFile >> numberInputParticipants[ r ];
inputFile >> numberOutputParticipants[ r ];
inputFile >> jkin[ r ];
inputFile.getline( buf, LINE );  // skip the new line
if( DEBUG ) {
cerr << bline[ copiedLines ] << endl;
cerr << "forwardReactionRates=" << forwardReactionRates[ r ];
cerr << "  backwardReactionRates=" << backwardReactionRates[ r ];
cerr << "  numberInputParticipants=" << numberInputParticipants[ r ];
cerr << "  numberOutputParticipants=" << numberOutputParticipants[ r ];
cerr << "  jkin=" << jkin[ r ] << endl;
}
if( jkin[ r ] == 11 ) {
inputFile >> forwardReactionRates2[ r ];
inputFile >> backwardReactionRates2[ r ];
inputFile.getline( buf, LINE );  // skip the new line
if( DEBUG ) {
cerr << "forwardReactionRates2=" << forwardReactionRates2[ r ];
cerr << "  backwardReactionRates2=" << backwardReactionRates2[ r ] << endl;
}
}
if( numberInputParticipants[ r ] > MAX_NUMBER_PARTICIPANT_REACTIONS 
||numberOutputParticipants[ r ] > MAX_NUMBER_PARTICIPANT_REACTIONS) {
cerr << endl << endl << " ERROR" << endl;
cerr << "   At reaction : " << r << " : " << bline[ copiedLines ];
cerr << endl << "   nipart=" << numberInputParticipants[ r ];
cerr << ", nopart=" << numberOutputParticipants[ r ] << endl;
cerr << "   Increase array size MAX_NUMBER_PARTICIPANT_REACTIONS";
cerr << endl << "  EXECUTION TERMINATED" << endl;
return 1;
}
// reactants:
for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
inputFile.getline( tmpLine, LINE );
nameOfInputSpecies[ q ][ r ] = new char[ strlen( tmpLine ) + 1 ];
strcpy( nameOfInputSpecies[ q ][ r ], tmpLine );
if( DEBUG ) {
cerr << "nameOfInputSpecies[" << q << "]=";
cerr << nameOfInputSpecies[ q ][ r ] << endl;
}
}
// products:
for( int p = 1; p <= numberOutputParticipants[ r ]; p++ ) {
inputFile.getline( tmpLine, LINE );
nameOfOutputSpecies[ p ][ r ] = new char[ strlen( tmpLine ) + 1 ];
strcpy( nameOfOutputSpecies[ p ][ r ], tmpLine );
if( DEBUG ) {
cerr << "nameOfOutputSpecies[" << p << "]=";
cerr << nameOfOutputSpecies[ p ][ r ] << endl;
}
}
}
inputFile.close();

// write out input parameters ( to file kin_o01 )
copiedLines = 0;
outputFile1 << " " << endl;
outputFile1 << " " << endl;
outputFile1 << bline[ copiedLines ] << endl;

copiedLines++;
outputFile1 << bline[ copiedLines ] << endl;
outputFile1 << setw( DEC8 ) << numberOfSpecies;
outputFile1 << setw( DEC8 ) << numberOfReactions << endl;

copiedLines++;
outputFile1 << bline[ copiedLines ] << endl;
outputFile1 << setiosflags( ios::scientific | ios::uppercase );
outputFile1 << setw( FRAC ) << initialTime;
outputFile1 << setw( FRAC ) << finalTime;
outputFile1 << setw( DEC8 ) << numberOfTimeSteps;
outputFile1 << setw( DEC8 ) << ntskip;
outputFile1 << setw( DEC8 ) << integrationOption << endl;

copiedLines++;
outputFile1 << bline[ copiedLines ] << endl;

for( int i = 1; i <= numberOfSpecies; i++ ) {
outputFile1 << nameOfSpecies[ i ] << endl;
outputFile1 << setw( FRAC ) << initialConcentration[ i ];
outputFile1 << setw( DEC8 ) << jfix[ i ] << endl;
}

for( int r = 1; r <= numberOfReactions; r++ ) {
copiedLines++;
outputFile1 << bline[ copiedLines ] << endl;
outputFile1 << setw( FRAC ) << forwardReactionRates[ r ];
outputFile1 << setw( FRAC ) << backwardReactionRates[ r ];
outputFile1 << setw( DEC8 ) << numberInputParticipants[ r ];
outputFile1 << setw( DEC8 ) << numberOutputParticipants[ r ];
outputFile1 << setw( DEC8 ) << jkin[ r ] << endl;
if( jkin[ r ] == 11 ) {
outputFile1 << forwardReactionRates2[ r ] << " ";
outputFile1 << backwardReactionRates2[ r ] << endl;
}
for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
outputFile1 << nameOfInputSpecies[ q ][ r ] << endl;
}
for( int p = 1; p <= numberOutputParticipants[ r ]; p++ ) {
outputFile1 << nameOfOutputSpecies[ p ][ r ] << endl;
}
}

outputFile1.close();
return 0;
}


// tabulate species number "ispec" for each named
// reactant and product species "nispec", "nospec" 
// of each reaction "ireac"
void setnet()
{
char * name2look = 0;
int index = 0;

for( int r = 1; r <= numberOfReactions; r++ ) {

for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
iispec[ q ][ r ] = ispec4name( trim( nameOfInputSpecies[ q ][ r ] ) );
}

for( int q = 1; q <= numberOutputParticipants[ r ]; q++ ) {
iospec[ q ][ r ] = ispec4name( trim( nameOfOutputSpecies[ q ][ r ] ) );
}
}

if( DEBUG ) {
cerr << "iispec[][] :" << endl;
for( int r = 1; r <= numberOfReactions; r++ ) {
for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
cerr << trim( nameOfInputSpecies[ q ][ r ] ) << "=";
cerr << iispec[ q ][ r ] << "\t ";
}
cerr << " ." << endl;
}
cerr << "iospec[][] :" << endl;
for( int r = 1; r <= numberOfReactions; r++ ) {
for( int q = 1; q <= numberOutputParticipants[ r ]; q++ ) {
cerr << trim( nameOfOutputSpecies[ q ][ r ] ) << "=";
cerr << iospec[ q ][ r ] << "\t ";
}
cerr << " ." << endl;
}
}

}

/*
* e v a l _ f x 
*
* Evaluates the r.h.s.
* version by Aleman
* March 19, 2001
*/
void eval_fx( double x_vector[], double result[], 
double leftProduct[], double rightProduct[] )
{
for( int r = 1; r <= numberOfReactions; r++ ) {
leftProduct[ r ] = forwardReactionRates[ r ];
for( int input = 1; input <= numberInputParticipants[ r ]; input++ ) {
leftProduct[ r ] *= x_vector[ iispec[ input ][ r ] ];
}

rightProduct[ r ] = backwardReactionRates[ r ];
for( int out = 1; out <= numberOutputParticipants[ r ]; out++ ) {
rightProduct[ r ] *= x_vector[ iospec[ out ][ r ] ];
}
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
result[ i ] = 0.0;
}
for( int r = 1; r <= numberOfReactions; r++ ) {
for( int input = 1; input <= numberInputParticipants[ r ]; input++ ) {
result[ iispec[ input ][ r ] ] = result[ iispec[ input ][ r ] ] 
- leftProduct[ r ] + rightProduct[ r ];
}
for( int out = 1; out <= numberOutputParticipants[ r ]; out++ ) {
result[ iospec[ out ][ r ] ] = result[ iospec[ out ][ r ] ] 
+ leftProduct[ r ] - rightProduct[ r ];
}
}
}

/*************************************************************
*   t e s t M e t h o d _ a l e m a n ()
*************************************************************
* test method by Aleman - N O T - U S E D
* March 19, 2001
*/
void testMethod_aleman()
{

Evaluate ***prepared_jac = 0;
prepared_jac = new Evaluate**[ numberOfSpecies + 1 ];
if( prepared_jac == 0 ) {
cerr << "ERROR: allocation for prepared_jac failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
prepared_jac[ i ] = new Evaluate*[ numberOfSpecies + 1 ];
if( prepared_jac[ i ] == 0 ) {
cerr << "ERROR: allocation for prepared_jac[" << i << "] failed!";
cerr << endl;
}
}
prepareJacobian( prepared_jac );

double **jac = 0;
jac = new double*[ numberOfSpecies + 1 ];
if( jac == 0 ) {
cerr << "ERROR: allocation for jac failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
jac[ i ] = new double[ numberOfSpecies + 1 ];
if( jac[ i ] == 0 ) {
cerr << "ERROR: allocation for jac[" << i << "] failed!";
cerr << endl;
}
}

eval_prep_jacobian( 1, jac, prepared_jac );

}

/*************************************************************
*lsodesMethod, first generate tmp.f file
************************************************************
* generate tmp.f file for lsodes Fortran code
* Added 02/12/04
*by Yihai Yu
*/
int lsodesMethod()
{
ofstream fortranFile( tmp_f);
if( ! fortranFile ) {
cerr << "Couldn't open input file " << fortranFile << endl;
return 1;
}

fortranFile << "      include \'" << includeArgv << "\'" << endl;
fortranFile << "*DECK PROGTMP" << endl;
fortranFile << "      PROGRAM PROGTMP" << endl;
fortranFile << endl;
fortranFile << "      external fex, jex" << endl;
fortranFile << "      double precision atol, rtol, rwork, t, tout, y" << endl;
//   fortranFile << "      dimension y("<< numberOfSpecies <<"), rwork(500), iwork(30)" << endl;
//   fortranFile << "      data lrw/500/, liw/30/" << endl;
int rwork = 20 + 3*numberOfSpecies*numberOfSpecies + 20*numberOfSpecies;
int iwork = 31 + numberOfSpecies + 2;
fortranFile << "      dimension y("<< numberOfSpecies <<"), rwork(" << rwork << "), iwork(" << iwork << ")" << endl;
fortranFile << "      data lrw/" << rwork << "/, liw/" << iwork << "/" << endl;


fortranFile << "      neq = "<< numberOfSpecies << endl;
for(int i=1; i<=numberOfSpecies; i++){
fortranFile << "      y("<< i <<") = " << initialConcentration[i] << endl;
}//endof for
fortranFile << endl;
fortranFile << "      t = " << initialTime << endl;
double tout = (finalTime - initialTime)/(double)(numberOfTimeSteps/ntskip);
fortranFile << "      tout = " << initialTime << endl;
fortranFile << "      itol = 1" << endl;
fortranFile << "      rtol = " << rtolLSODES <<endl;
fortranFile << "      atol = " << atolLSODES <<endl;
fortranFile << "      itask = 1" <<endl;
fortranFile << "      istate = 1" <<endl;
fortranFile << "      iopt = 0" <<endl;
fortranFile << "      mf = "<< mfLSODES <<endl;
fortranFile << "      OPEN (UNIT=20, FILE=" << "\n     1" << "'" << outputArgv << "')" << endl;

//   fortranFile << "      OPEN (UNIT=20, FILE='columns.data')" << endl;

fortranFile << "      do 40 iout = 0, "<< numberOfTimeSteps/ntskip <<endl;
fortranFile << "        call lsodes (fex, neq, y, t, tout, itol, rtol, atol, " <<endl;
fortranFile << "     1     itask, istate, iopt, rwork, lrw, iwork, liw, jex, mf)" << endl;
fortranFile << "      write(20, *)iout, t, y(" << lsodesArgv << ")"<< endl;
fortranFile << endl;
fortranFile << "        if (istate .lt. 0) go to 80" << endl;
fortranFile << "        tout = tout + " << tout <<endl;
fortranFile << "  40    continue" << endl;
fortranFile << "      lenrw = iwork(17)" << endl;
fortranFile << "      leniw = iwork(18)" << endl;
fortranFile << "      nst = iwork(11)" << endl;
fortranFile << "      nfe = iwork(12)" << endl;
fortranFile << "      nje = iwork(13)" << endl;
fortranFile << "      nlu = iwork(21)" << endl;
fortranFile << "      nnz = iwork(19)" << endl;
fortranFile << "      nnzlu = iwork(25) + iwork(26) + neq" << endl;
fortranFile << "      stop" << endl;
fortranFile << "  80  write(6,90)istate" << endl;
fortranFile << "  90  format(///22h error halt.. istate =,i3)" << endl;
fortranFile << "      stop" << endl;
fortranFile << "      end" << endl;

//fex
fortranFile << "      subroutine fex (neq, t, y, ydot)" << endl;
fortranFile << "      double precision t, y, ydot" << endl;
fortranFile << "      dimension y(" << numberOfSpecies << "), ydot(" << numberOfSpecies << ")" << endl;


//   cout << "OK1" << endl;

//   cout << "1 of 8: " << nameOfOutputSpecies[1][8] << endl;

//Pre-process all the reactions to get the index of all the participants
for(int i=1; i<=numberOfReactions; i++){
//cout << "i = " << i << endl;
for(int j=1; j<=numberInputParticipants[i]; j++){
for(int k=1; k<=numberOfSpecies; k++){
if(mystrcmp(nameOfInputSpecies[j][i], nameOfSpecies[k])==0){
//cout << strlen(nameOfInputSpecies[j][i]) << "  " << strlen(nameOfSpecies[k]) <<endl;
//cout << nameOfInputSpecies[j][i] << "  " << nameOfSpecies[k] <<endl;
//if(nameOfInputSpecies[j][i][17]=='2')   cout<<"this is space" <<endl;
indexOfInputSpecies[j][i] = k;
}//endof if
}//endof k
}//endof j
for(int j=1; j<=numberOutputParticipants[i]; j++){
for(int k=1; k<=numberOfSpecies; k++){
if(mystrcmp(nameOfOutputSpecies[j][i], nameOfSpecies[k])==0){
indexOfOutputSpecies[j][i] = k;
}//endof if
}//endof k
}//endof j
}//endof for i

//   cout << "OK2" << endl;

for(int i=1; i<=numberOfSpecies; i++){

int notHere = 0;//1 here, 0 not here, for this specie

//for each specie, calculate the dy/dt
//jfix

if(jfix[i]==1){
fortranFile << "      ydot(" << i << ")  = 0" << endl;
}//endof if

else{

//LHS
fortranFile << "      ydot(" << i << ")  = ";

//RHS
//all the reactions
for(int j=1; j<=numberOfReactions; j++){
//all the input participansts
for(int k=1; k<=numberInputParticipants[j]; k++){
if(mystrcmp(nameOfSpecies[i], nameOfInputSpecies[k][j])==0){
notHere = 1;
//forward reaction
fortranFile << "\n" << "     1   -" << forwardReactionRates[j];
for(int l=1; l<=numberInputParticipants[j]; l++){
fortranFile << "*" << "y(" << indexOfInputSpecies[l][j] << ")";
}//endof l
//backward reaction
fortranFile << "\n" << "     1   +" << backwardReactionRates[j];
for(int l=1; l<=numberOutputParticipants[j]; l++){
fortranFile << "*" << "y(" << indexOfOutputSpecies[l][j] << ")";
}//endof l
}//endof if
}//endof for k
//all the output participants
for(int k=1; k<=numberOutputParticipants[j]; k++){
if(mystrcmp(nameOfSpecies[i], nameOfOutputSpecies[k][j])==0){
notHere = 1;
//forward reaction
fortranFile << "\n" << "     1   +" << forwardReactionRates[j];
for(int l=1; l<=numberInputParticipants[j]; l++){
fortranFile << "*" << "y(" << indexOfInputSpecies[l][j] << ")";
}//endof l
//backward reaction
fortranFile << "\n" << "     1   -" << backwardReactionRates[j];
for(int l=1; l<=numberOutputParticipants[j]; l++){
fortranFile << "*" << "y(" << indexOfOutputSpecies[l][j] << ")";
}//endof l
}//endof if
}//endof for k
}//endof for j
if(notHere==0){
fortranFile << "0";
}
}//endof else      
fortranFile << endl;
}//endof for i
fortranFile << "      return" << endl;
fortranFile << "      end" << endl;

//jex
fortranFile << "      subroutine jex (neq, t, y, j, ia, ja, pdj)" << endl;
fortranFile << "      double precision t, y, pdj" << endl;
fortranFile << "      dimension y(1), ia(1), ja(1), pdj(1)" << endl;
fortranFile << "      end" << endl;


fortranFile.close();
return 0;      

}



/*************************************************************
*   e v a l _ p r e p _ j a c o b i a n
*************************************************************
* Evals the jacobian matrix with xspec[ time ][]
* Added Mar 20, 2001
* this version is supposed to work faster because it uses
* a prepared jacobian structure
* by Aleman-Meza, Boanerges
*/
void eval_prep_jacobian( int time, double ** jac, Evaluate *** prepared_jac )
{
double sum   = 0.0;
double multi = 0.0;

for(  int i = 1; i <= numberOfSpecies; i++ ) {

for( int j = 1; j <= numberOfSpecies; j++ ) {

sum = 0.0;
for( Evaluate *eval = prepared_jac[ i ][ j ]; eval != 0; 
eval = eval->next ) {
multi = eval->multiplier * eval->appears;
for( int m = 0; m < eval->participants; m++ ) {
multi *= xspec[ time ][ eval->indexes[ m ] ];
}
for( int q = 1; q < eval->appears; q++ ) {
multi *= xspec[ time ][ j ];
}
sum += multi;
}
jac[ i ][ j ] = sum;
}
}
}

double con_jfix_10(int i, double currentTime){

double startTime;
double startCon, endCon, inCycleStart, inCycleEnd;
double inCycleTime;

startTime = true_initialTime+extOfSpecies[i][1];
if(integrationOption!=4){
currentTime = initialTime + currentTime;
}
if(currentTime <= startTime){
return initialConcentration[i];
}//endof if
else{
inCycleTime = (currentTime-startTime) - oneCycle[i]*((int)((currentTime-startTime)/oneCycle[i]));
//cerr << inCycleTime;
inCycleStart = 0.0;
inCycleEnd = 0.0;
for(int j=2; j<2*npulse[i]+2; j=j+2){
inCycleEnd = inCycleEnd + extOfSpecies[i][j+1];
if(inCycleTime<=inCycleEnd){
return slopeTrap[i][j+1]*(inCycleTime-inCycleStart)+extOfSpecies[i][j];
break;
}
inCycleStart = inCycleEnd;
}//endof for j
}//endof else
}

double con_jfix_20(int i, double currentTime){

double startTime;
double startCon, endCon, inCycleStart, inCycleEnd;
double inCycleTime;

startTime = true_initialTime+extOfSpecies[i][1];
if(integrationOption!=4){
currentTime = initialTime + currentTime;
}
if(currentTime <= startTime){
return extOfSpecies[i][2];
}//endof if
else{
inCycleTime = (currentTime-startTime) - oneCycle[i]*((int)((currentTime-startTime)/oneCycle[i]));
//cerr << inCycleTime;
inCycleStart = 0.0;
inCycleEnd = 0.0;
for(int j=2; j<2*npulse[i]+2; j=j+2){
inCycleEnd = inCycleEnd + extOfSpecies[i][j+1];
if(inCycleTime<=inCycleEnd){
return extOfSpecies[i][j+2];
break;
}
inCycleStart = inCycleEnd;
}//endof for j
}//endof else
}


/*************************************************************
*   e u l e r M e t h o d ()
*************************************************************
* jtime=1: Original Euler method (integrationOption)
*/
void eulerMethod()
{

double vspec[ MAX_NUMBER_SPECIES ]  = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };

double currentTime;

for( int t = 1; t <= numberOfTimeSteps; t++ ) {

currentTime = dtime*t;

xrate( vspec, vfor, vbak, t );

for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + vspec[ i ] * dtime;
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);

}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);

}//endof if jfix = 20
} // for i
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ]
xreac[ t ][ r ] = xreac[ t - 1 ][ r ]
+ ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif
}

}

/*************************************************************
*   m o d i f i e d E u l e r M e t h o d ()
*************************************************************
* T. Taha
* Department of Computer Science
* University of Georgia
* Athens, GA 30602
* E-Mail thiab@cs.uga.edu
* added 00-08-01
*/
void modifiedEulerMethod()
{

double vspec[ MAX_NUMBER_SPECIES ]  = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };
double k1[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
#ifndef OPTIMIZE
double k1_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
#endif

double currentTime;

for( int t = 1; t <= numberOfTimeSteps; t++ ) {

currentTime = dtime*t;

// evaluate f( t-1, xpec[ t-1 ] )
xrate( vspec, vfor, vbak, t );
// put h * f( t-1, xpec[ t-1 ] ) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
k1[ i ] = vspec[ i ] * dtime;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k1_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

// put x[t-1] + k1 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + k1[ i ];
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);
}//endof if jfix = 20

}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ] + k1_reac[ r ];
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + k1_reac[ r ];
}
#endif

// evaluate f( t, xspec[ t ] + k1 )
xrate( vspec, vfor, vbak, t + 1 );
// put xspec[ t-1 ] + (1/2) ( k0 + h * f(...) into xspec[ t ]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + 0.5 * ( k1[ i ] 
+ vspec[ i ] * dtime );
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + 0.5 * ( k1_reac[ r ] 
+ ( vfor[ r ] - vbak[ r ] ) * dtime );
}
#endif


}
}

void modifiedEulerMethod_original()
{ 

double vspec[ MAX_NUMBER_SPECIES ]  = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };

for( int t = 1; t <= numberOfTimeSteps; t++ ) {

xrate( vspec, vfor, vbak, t );

for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + vspec[ i ] * dtime;
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}

} // for i

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ]
xreac[ t ][ r ] = xreac[ t - 1 ][ r ]
+ ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

xrate( vspec, vfor, vbak, t + 1 );

for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = 0.5 * ( xspec[ t ][ i ] + xspec[ t - 1 ][ i ]
+ vspec[ i ] * dtime );
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
} // for i

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = 0.5 * ( xreac[ r ][ t ] + xreac[ r ][ t - 1 ]
xreac[ t ][ r ] = 0.5 * ( xreac[ t ][ r ] + xreac[ t - 1 ][ r ]
+ ( vfor[ r ] - vbak[ r ] ) * dtime );
}
#endif

}
}

/*************************************************************
*   r u n g e K u t t a O r d e r 4
*************************************************************
* Runge-Kutta Method order h^4
* Added Jan 30, 2001
* by Aleman-Meza, Boanerges
* x[t+1] = x[t] + (1/6) * ( k1 + 2 * k2 + 2 * k3 + k4 )
* where:
* k1 = h * f( t,       x[t] )
* k2 = h * f( t + h/2, x[t] + (1/2) * k1 )
* k3 = h * f( t + h/2, x[t] + (1/2) * k2 )
* k4 = h * f( t + 1,   x[t] + k3 )
*/
void rungeKuttaOrder4()
{
double vspec[ MAX_NUMBER_SPECIES ]  = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };
double k1[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k2[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k3[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k4[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
#ifndef OPTIMIZE
double k1_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k2_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k3_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k4_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
#endif

double currentTime;

for( int t = 1; t <= numberOfTimeSteps; t++ ) {

currentTime = dtime*t;

// evaluate f( t-1, xpec[ t-1 ] )
xrate( vspec, vfor, vbak, t );
// put h * f( t-1, xpec[ t-1 ] ) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
k1[ i ] = vspec[ i ] * dtime;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k1_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

// put x[t-1] + k1 / 2.0 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + k1[ i ] / 2.0;
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);

}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);

}//endof if jfix = 20

}

}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ] + k1_reac[ r ] / 2.0;
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + k1_reac[ r ] / 2.0;
}
#endif

// evaluate f( t, xspec[ t ] + (1/2) * k1 )
xrate( vspec, vfor, vbak, t + 1 );
// put h * f( t, xspec[ t ] + (1/2) * k1 ) into k2
for( int i = 1; i <= numberOfSpecies; i++ ) {
k2[ i ] = vspec[ i ] * dtime;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k2_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

// put x[t-1] + k2 / 2.0 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + k2[ i ] / 2.0;
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);

}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);

}//endof if jfix = 20

}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ] + k2_reac[ r ] / 2.0;
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + k2_reac[ r ] / 2.0;
}
#endif

// evaluate f( t, xspec[ t ] + (1/2) * k2 )
xrate( vspec, vfor, vbak, t + 1 );
// put h * f( t, xspec[ t ] + (1/2) * k2 ) into k3
for( int i = 1; i <= numberOfSpecies; i++ ) {
k3[ i ] = vspec[ i ] * dtime;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k3_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

// put x[t-1] + k3  into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + k3[ i ];
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);

}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);

}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ] + k3_reac[ r ];
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + k3_reac[ r ];
}
#endif

// evaluate f( t, xspec[ t ] + k3 )
xrate( vspec, vfor, vbak, t + 1 );
// put h * f( t, xspec[ t ] + k3 ) into k4
for( int i = 1; i <= numberOfSpecies; i++ ) {
k4[ i ] = vspec[ i ] * dtime;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k4_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * dtime;
}
#endif

// put x[t-1] + (1/6)*(k1 + 2*k2 + 2*k3 + k4) into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ] + ( k1[ i ] + 2.0 * k2[ i ]
+ 2.0 * k3[ i ] + k4[ i ] ) / 6.0;
if( xspec[ t ][ i ] < 0.0 ) {
xspec[ t ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t][i] = con_jfix_10(i, currentTime);

}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t][i] = con_jfix_20(i, currentTime);

}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t ] = xreac[ r ][ t - 1 ] + ( k1_reac[ r ]
xreac[ t ][ r ] = xreac[ t - 1 ][ r ] + ( k1_reac[ r ]
+ 2.0 * k2_reac[ r ] + 2.0 * k3_reac[ r ]
+ k4_reac[ r ] ) / 6.0;
}
#endif

} 

} // end method: rungeKuttaOrder4

/*************************************************************
*   r u n g e K u t t a 4 5
*************************************************************
* Runge-Kutta 45 Method
* Added Feb 07, 2001
* by Aleman-Meza, Boanerges
* x[t+h] = x[t] + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4
*               - (9/50) * k5 + (2/55)
* where:
* k1 = h * f( t, x[t] )
* k2 = h * f( t, x[t] + (1/4) * k1 )
* k3 = h * f( t, x[t] + (3/32) * k1 + (9/32) * k2 )
* k4 = h * f( t, x[t] + (1932/2197) * k1 - (7200/2197) * k2 
*                     + (7296/2197) * k3 )
* k5 = h * f( t, x[t] + (439/216) * k1 - 8 * k2
*                     + (3680/513) * k3 - (845/4104) * k4 )
* k6 = h * f( t, x[t] - (8/27) * k1 + 2 * k2
*                     - (3544/2565) * k3 + (1859/4104) * k4 - (11/40) * k5 )
*/
void rungeKutta45()
{

double vspec[ MAX_NUMBER_SPECIES ]  = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };

double k1[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k2[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k3[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k4[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k5[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k6[ MAX_NUMBER_SPECIES + 1 ] = { 0 };

#ifndef OPTIMIZE
double k1_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k2_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k3_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k4_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k5_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k6_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
#endif

double h = dtime;

double currentTime = initialTime;

for( int t_index = 1; t_index <= numberOfTimeSteps; t_index++ ) {

currentTime = currentTime + h;

// k1
// evaluate f( t-1, xpec[ t-1 ] )
xrate( vspec, vfor, vbak, t_index );
// put h * f(...) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
k1[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k1_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k2
// put x[t-1] + (1/4) * k1 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c21 * k1[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c21 * k1_reac[ r ];
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c21 * k1_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (1/4) * k1 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k2
for( int i = 1; i <= numberOfSpecies; i++ ) {
k2[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k2_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k3
// put x[t-1] + (3/32) * k1 + (9/32) * k2 ) into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c31 * k1[ i ]
+ c32 * k2[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c31 * k1_reac[ r ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c31 * k1_reac[ r ]
+ c32 * k2_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (3/32) * k1 + (9/32) * k2 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k3
for( int i = 1; i <= numberOfSpecies; i++ ) {
k3[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k3_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k4
// put x[t-1] + (1932/2197) * k1 - (7200/2197) * k2
//            + (7296/2197) * k3 ) into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c41 * k1[ i ]
+ c42 * k2[ i ] + c43 * k3[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c41 * k1_reac[ r ]
+ c42 * k2_reac[ r ] + c43 * k3_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (1932/2197) * k1 - (7200/2197) * k2
//            + (7296/2197) * k3 ) )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k4
for( int i = 1; i <= numberOfSpecies; i++ ) {
k4[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k4_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k5
// put x[t-1] + (439/216) * k1 - 8 * k2
//            + (3680/513) * k3 - (845/4104) * k4 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c51 * k1[ i ]
+ c52 * k2[ i ] + c53 * k3[ i ]
+ c54 * k4[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c51 * k1_reac[ r ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c51 * k1_reac[ r ]
+ c52 * k2_reac[ r ] + c53 * k3_reac[ r ]
+ c54 * k4_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (439/216) * k1 - 8 * k2
//            + (3680/513) * k3 - (845/4104) * k4 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k5
for( int i = 1; i <= numberOfSpecies; i++ ) {
k5[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k5_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k6
// put x[t-1] - (8/27) * k1 + 2 * k2 - (3544/2565) * k3
//            + (1859/4104) * k4 - (11/40) * k5 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c61 * k1[ i ]
+ c62 * k2[ i ] + c63 * k3[ i ]
+ c64 * k4[ i ] + c65 * k5[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c61 * k1_reac[ r ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c61 * k1_reac[ r ]
+ c62 * k2_reac[ r ] + c63 * k3_reac[ r ]
+ c64 * k4_reac[ r ] + c65 * k5_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] - (8/27) * k1 + 2 * k2
//                              - (3544/2565) * k3 + (1859/4104) * k4
//                              - (11/40) * k5 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k6
for( int i = 1; i <= numberOfSpecies; i++ ) {
k6[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k6_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// put x[t-1] + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4
//            - (9/50) * k5 + (2/55)  into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + b1 * k1[ i ]
+ b3 * k3[ i ] + b4 * k4[ i ]
+ b5 * k5[ i ] + b6 * k6[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ]
+ b1 * k1_reac[ r ]
+ b3 * k3_reac[ r ] + b4 * k4_reac[ r ]
+ b5 * k5_reac[ r ] + b6 * k6_reac[ r ];
}
#endif
}

} // end method: rungeKutta45



/*************************************************************
*   r u n g e K u t t a A d a p t i v e 
*************************************************************
* Adaptive Runge-Kutta Fehlberg Method
* Added Feb 07, 2001
* by Aleman-Meza, Boanerges
* x[t+h] = x[t] + (25/216) * k1 + (1408/2565) * k3 + (2197/4104) * k4 
*               - (1/5) * k5
* where:
* k1 = h * f( t,               x[t] )
* k2 = h * f( t + (1/2) * h,   x[t] + (1/4) * k1 )
* k3 = h * f( t + (3/8) * h,   x[t] + (3/32) * k1 + (9/32) * k2 )
* k4 = h * f( t + (12/13) * h, x[t] + (1932/2197) * k1 - (7200/2197) * k2
*                                   + (7296/2197) * k3 )
* k5 = h * f( t + h,           x[t] + (439/216) * k1 - 8 * k2 
*                                   + (3680/513) * k3 - (845/4104) * k4 )
* k6 = h * f( t,               x[t] - (8/27) * k1 + 2 * k2 
*                                   - (3544/2565) * k3 + (1859/4104) * k4 
*                                   - (11/40) * k5 )
*/
void rungeKuttaAdaptive()
{
double vspec[ MAX_NUMBER_SPECIES ] = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };

double k1[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k2[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k3[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k4[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k5[ MAX_NUMBER_SPECIES + 1 ] = { 0 };
double k6[ MAX_NUMBER_SPECIES + 1 ] = { 0 };

#ifndef OPTIMIZE
double k1_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k2_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k3_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k4_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k5_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
double k6_reac[ MAX_NUMBER_REACTIONS + 1 ] = { 0 };
#endif

int iterations = 0;

double epsilon_max = 1e-6;  //<- these 2 work with accuracy 1e-6
double epsilon_min = 1e-13;


double h_max = 1.0 / 4;   //dtime * 10.0;
double h_min = 0.001;

//change h 
//double h = h_max; //( finalTime - initialTime ) * 1.0 / 1000.0;
double h = h_min;

double diff = 0.0;
double x4 = 0.0;
double x4_reac = 0.0;
double epsilon = 0.0;

// initialy xpec[][ 0 ] holds the i.c.

int exit_loop = 0;
int t_index = 0;
double time = initialTime;
double time_save = 0.0;

xtime_index = 0;
xtime[ t_index ] = time;

while( ( t_index <= numberOfTimeSteps ) && ( exit_loop < 2 ) ) {

/*        
if( h < h_min ) {
h = h_min;
}
if( h > h_max ) {
h = h_max;
}
*/  

diff = myabs( finalTime - time );
if( diff <= myabs( h ) ) {
if( diff < ( 0.5e-5 * finalTime ) ) {
break;
}
h = diff;
exit_loop++;
}      

t_index++;
iterations++;
time_save = time;
time = time + h;

// k1 
// evaluate f( t-1, xpec[ t-1 ] )
xrate( vspec, vfor, vbak, t_index );
// put h * f(...) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
k1[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k1_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k2 
// put x[t-1] + (1/4) * k1 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c21 * k1[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] +  c21 * k1_reac[ r ];
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] +  c21 * k1_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (1/4) * k1 ) 
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k2
for( int i = 1; i <= numberOfSpecies; i++ ) {
k2[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k2_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k3 
// put x[t-1] + (3/32) * k1 + (9/32) * k2 ) into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c31 * k1[ i ]
+ c32 * k2[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20

}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] 
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] 
+ c31 * k1_reac[ r ] + c32 * k2_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (3/32) * k1 + (9/32) * k2 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k3
for( int i = 1; i <= numberOfSpecies; i++ ) {
k3[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k3_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k4
// put x[t-1] + (1932/2197) * k1 - (7200/2197) * k2
//            + (7296/2197) * k3 ) into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c41 * k1[ i ]
+ c42 * k2[ i ] + c43 * k3[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c41 * k1_reac[ r ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c41 * k1_reac[ r ]
+ c42 * k2_reac[ r ] + c43 * k3_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (1932/2197) * k1 - (7200/2197) * k2
//            + (7296/2197) * k3 ) )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k4
for( int i = 1; i <= numberOfSpecies; i++ ) {
k4[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k4_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k5
// put x[t-1] + (439/216) * k1 - 8 * k2
//            + (3680/513) * k3 - (845/4104) * k4 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c51 * k1[ i ]
+ c52 * k2[ i ] + c53 * k3[ i ] + c54 * k4[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] 
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] 
+ c51 * k1_reac[ r ] + c52 * k2_reac[ r ] 
+ c53 * k3_reac[ r ] + c54 * k4_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] + (439/216) * k1 - 8 * k2
//            + (3680/513) * k3 - (845/4104) * k4 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k5
for( int i = 1; i <= numberOfSpecies; i++ ) {
k5[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k5_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// k6
// put x[t-1] - (8/27) * k1 + 2 * k2 - (3544/2565) * k3
//            + (1859/4104) * k4 - (11/40) * k5 into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + c61 * k1[ i ]
+ c62 * k2[ i ] + c63 * k3[ i ]
+ c64 * k4[ i ] + c65 * k5[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + c61 * k1_reac[ r ]
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + c61 * k1_reac[ r ]
+ c62 * k2_reac[ r ] + c63 * k3_reac[ r ]
+ c64 * k4_reac[ r ] + c65 * k5_reac[ r ];
}
#endif

// evaluate f( t-1, xpec[ t-1 ] - (8/27) * k1 + 2 * k2 
//                              - (3544/2565) * k3 + (1859/4104) * k4 
//                              - (11/40) * k5 )
xrate( vspec, vfor, vbak, t_index + 1 );
// put h * f( ... ) into k6
for( int i = 1; i <= numberOfSpecies; i++ ) {
k6[ i ] = vspec[ i ] * h;
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
k6_reac[ r ] = ( vfor[ r ] - vbak[ r ] ) * h;
}
#endif

// put x[t-1] + (16/135) * k1 + (6656/12825) * k3 + (28561/56430) * k4
//            - (9/50) * k5 + (2/55)  into x[t]
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
xspec[ t_index ][ i ] = xspec[ t_index - 1 ][ i ] + b1 * k1[ i ] 
+ b3 * k3[ i ] + b4 * k4[ i ] 
+ b5 * k5[ i ] + b6 * k6[ i ];
x4 = xspec[ t_index - 1 ][ i ] + a1 * k1[ i ]
+ a3 * k3[ i ] + a4 * k4[ i ] + a5 * k5[ i ];
if( xspec[ t_index ][ i ] < 0.0 ) {
xspec[ t_index ][ i ] = 0.0;
}
xtime[ t_index ] = time;
if( x4 < 0.0 ) {
x4 = 0.0;
}
diff = myabs( xspec[ t_index ][ i ] - x4 );
if( diff > epsilon_min ) {
epsilon = diff;
if( epsilon > epsilon_max ) {
break;
}
}
else {
epsilon = diff;
}
}
else if( jfix[ i ] == 10){
// trapezoidal polygon pulse
xspec[t_index][i] = con_jfix_10(i, time); 
}//endof if jfix=10
else if( jfix[ i ] == 20){
// rectangle pulse
xspec[t_index][i] = con_jfix_20(i, time);
}//endof if jfix = 20

}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ t_index ] = xreac[ r ][ t_index - 1 ] + b1 * k1_reac[ r ] 
xreac[ t_index ][ r ] = xreac[ t_index - 1 ][ r ] + b1 * k1_reac[ r ] 
+ b3 * k3_reac[ r ] + b4 * k4_reac[ r ] 
+ b5 * k5_reac[ r ] + b6 * k6_reac[ r ];
//x4_reac = xreac[ r ][ t_index - 1 ] + a1 * k1_reac[ r ]
x4_reac = xreac[ t_index - 1 ][ r ] + a1 * k1_reac[ r ]
+ a3 * k3_reac[ r ] + a4 * k4_reac[ r ] + a5 * k5_reac[ r ];
//diff = myabs( xreac[ r ][ t_index ] - x4_reac );
diff = myabs( xreac[ t_index ][ r ] - x4_reac );
if( diff > epsilon_min && diff > epsilon_max ) {
epsilon = diff;
break;
}
}
#endif

if( epsilon < epsilon_min ) {
h = 2.0 * h;
}
if( epsilon > epsilon_max ) {
h = h / 2.0;
time = time_save;
t_index = t_index - 1;
}

// avoid to stay forever
if( iterations > ( numberOfTimeSteps * 2 ) ) {
break;
}

} 

// cerr << "iterations = " << iterations << endl;
// cerr << "timesteps = " << t_index << endl;
xtime_index = t_index;

for( int i = 1; i <= numberOfSpecies; i++ ) {
xspec[ numberOfTimeSteps ][ i ] = xspec[ t_index ][ i ];
}
#ifndef OPTIMIZE
for( int i = 1; i <= numberOfReactions; i++ ) {
//xreac[ i ][ numberOfTimeSteps ] = xreac[ i ][ t_index ];
xreac[ numberOfTimeSteps ][ i ] = xreac[ t_index ][ i ];
}
#endif

} // end method: rungeKuttaAdaptive


/*************************************************************
*   e v a l _ j a c o b i a n
*************************************************************
* Evals the jacobian matrix with xspec[ time ][]
* Added Mar 01, 2001
* by Aleman-Meza, Boanerges
*/
void eval_jacobian( int time, double ** jac )
{

double multi   = 1.0;
double sum     = 0.0;
double fwdPart = 0.0;
double bckPart = 0.0;
int existsInLeft  = 0;
int existsInRight = 0;
int respecttoAppears = 0;

for( int equation = 1; equation <= numberOfSpecies; equation++ ) {

for( int respectto = 1; respectto <= numberOfSpecies; respectto++ ) {

sum = 0.0;

for( int r = 1; r <= numberOfReactions; r++ ) {

multi = 1.0;
respecttoAppears = 0;
existsInLeft  = 0;
existsInRight = 0;

for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
if( equation == iispec[ q ][ r ] ) {
existsInLeft++;
}
if( respectto == iispec[ q ][ r ] ) {
respecttoAppears++;
}
else {
multi *= xspec[ time ][ iispec[ q ][ r ] ];
}
}
// for the case x^3 ( which derivative is 3*x^2), compute x^2
for( int k = 1; k < respecttoAppears; k++ ) {
multi *= xspec[ time ][ respectto ];
}
fwdPart = multi * forwardReactionRates[ r ] * ( (double) respecttoAppears );

multi = 1.0;
respecttoAppears = 0;

for( int q = 1; q <= numberOutputParticipants[ r ]; q++ ) {
if( equation == iospec[ q ][ r ] ) {
existsInRight++;
}
if( respectto == iospec[ q ][ r ] ) {
respecttoAppears++;
}
else {
multi *= xspec[ time ][ iospec[ q ][ r ] ];
}
}
// for the case x^3 ( which derivative is 3*x^2), compute x^2
for( int k = 1; k < respecttoAppears; k++ ) {
multi *= xspec[ time ][ respectto ];
}
bckPart = multi * backwardReactionRates[ r ] * ( (double) respecttoAppears );

if( existsInLeft ) {
sum = sum - fwdPart + bckPart;
}
if( existsInRight ) {
sum = sum + fwdPart - bckPart;
}

}

jac[ equation ][ respectto ] = sum;

}

}
}

/*************************************************************
*   p r e p a r e J a c o b i a n
*************************************************************
* prepares the jacobian matrix ( will be evaluated later many times )
* Added Mar 19, 2001
* by Aleman-Meza, Boanerges
*/
void prepareJacobian( Evaluate *** jac )
{

double   *tMulti = new double[ MAX_NUMBER_PARTICIPANT_REACTIONS * 2 ];
int      *tIndex = new    int[ MAX_NUMBER_PARTICIPANT_REACTIONS * 2 ];
int      tPartic = 0;
Evaluate *left   = 0;
Evaluate *right  = 0;
Evaluate *eval   = 0;
Evaluate *list   = 0;
int existsInLeft = 0;
int existsInRight = 0;
int respecttoAppears = 0;

for( int equation = 1; equation <= numberOfSpecies; equation++ ) {

for( int respectto = 1; respectto <= numberOfSpecies; respectto++ ) {

list = 0;

for( int r = 1; r <= numberOfReactions; r++ ) {

tPartic = 0;
respecttoAppears = 0;
existsInLeft  = 0;
left = 0;

for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
if( equation == iispec[ q ][ r ] ) {
existsInLeft++;
}
if( respectto == iispec[ q ][ r ] ) {
respecttoAppears++;
}
else {
tIndex[ tPartic++ ] = iispec[ q ][ r ];
}
}
if( respecttoAppears ) {
left = new Evaluate;
left->appears = respecttoAppears;
left->multiplier = forwardReactionRates[ r ];
left->participants = tPartic;
left->indexes = new int[ tPartic ];
for( int i = 0; i < tPartic; i++ ) {
left->indexes[ i ] = tIndex[ i ];
}
left->next = 0;
}

tPartic = 0;
respecttoAppears = 0;
existsInRight = 0;
right = 0;

for( int q = 1; q <= numberOutputParticipants[ r ]; q++ ) {
if( equation == iospec[ q ][ r ] ) {
existsInRight++;
}
if( respectto == iospec[ q ][ r ] ) {
respecttoAppears++;
}
else {
tIndex[ tPartic++ ] = iospec[ q ][ r ];
}
}
if( respecttoAppears ) {
right = new Evaluate;
right->appears = respecttoAppears;
right->multiplier = backwardReactionRates[ r ];
right->participants = tPartic;
right->indexes = new int[ tPartic ];
for( int i = 0; i < tPartic; i++ ) {
right->indexes[ i ] = tIndex[ i ];
}
right->next = 0;
}
if( ( ! existsInLeft && ! existsInRight )
|| ( existsInLeft && existsInRight ) ) {
// they cancel each other, do not store them
if( left != 0 ) {
delete [] left->indexes;
delete left;
}
if( right != 0 ) {
delete [] right->indexes;
delete right;
}
continue;
}
if( existsInLeft && left != 0 ) {
left->multiplier = -left->multiplier;
// add it to the list
if( list == 0 ) {
list = left;
}
else {
for( Evaluate *tmp = list; tmp != 0; tmp = tmp->next ) {
if( tmp->next == 0 ) {
tmp->next = left;
break;
}
}
}
} 
if( existsInRight && right != 0 ) {
right->multiplier = -right->multiplier;
// add it to the list
if( list == 0 ) {
list = right;
}
else {
for( Evaluate *tmp = list; tmp != 0; tmp = tmp->next ) {
if( tmp->next == 0 ) {
tmp->next = right;
break;
}
}
}
} 

}
jac[ equation ][ respectto ] = list;

}

}
delete [] tMulti;
delete [] tIndex;
}


/*************************************************************
*   s t i f f S o l v e r
*************************************************************
*
*/
void stiffSolver()
{
Evaluate ***prepared_jac = 0;
prepared_jac = new Evaluate**[ numberOfSpecies + 1 ];
if( prepared_jac == 0 ) {
cerr << "ERROR: allocation for prepared_jac failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
prepared_jac[ i ] = new Evaluate*[ numberOfSpecies + 1 ];
if( prepared_jac[ i ] == 0 ) {
cerr << "ERROR: allocation for prepared_jac[" << i << "] failed!";
cerr << endl;
}
}
prepareJacobian( prepared_jac );

/*
cerr << "h = " << dtime << ", initialTime = " << initialTime;
cerr << ", finalTime = " << finalTime << endl;
*/

clock_t start, end;
double elapsed_jac = 0.0;
double elapsed_gauss = 0.0;
double elapsed_xrate = 0.0;

double vspec[ MAX_NUMBER_SPECIES ] = { 0 };
double vfor[ MAX_NUMBER_REACTIONS ] = { 0 };
double vbak[ MAX_NUMBER_REACTIONS ] = { 0 };

const int maxIts     = 100;
//const double epsilon = 1.0E-6;
const double epsilon = 1.0E-3;
int Y0 = 0;
int Y1 = 1;
int Y2 = 2;

double time                = 0.0;
double maxerror            = 0.0;
double newvalue            = 0.0;
double *initialConditions = new double[ numberOfSpecies + 1 ];
double *bigF      = new double[ numberOfSpecies + 1 ];
double *m         = new double[ numberOfSpecies + 1 ];
double *temp      = new double[ numberOfSpecies + 1 ];
int    *tempInt   = new int[ numberOfSpecies + 1 ];
int    its        = 0;
double totalIts   = 0;

double **bigFprime = 0;
bigFprime = new double*[ numberOfSpecies + 1 ];
if( bigFprime == 0 ) {
cerr << "ERROR: allocation for bigFprime failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
bigFprime[ i ] = new double[ numberOfSpecies + 1 ];
if( bigFprime[ i ] == 0 ) {
cerr << "ERROR: allocation for bigFprime[" << i << "] failed!";
cerr << endl;
}
}

double **jac = 0;
jac = new double*[ numberOfSpecies + 1 ];
if( jac == 0 ) {
cerr << "ERROR: allocation for jac failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
jac[ i ] = new double[ numberOfSpecies + 1 ];
if( jac[ i ] == 0 ) {
cerr << "ERROR: allocation for jac[" << i << "] failed!";
cerr << endl;
}
}

/*
double **jac2 = 0;
jac2 = new (double*)[ numberOfSpecies + 1 ];
if( jac2 == 0 ) {
cerr << "ERROR: allocation for jac2 failed!" << endl;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
jac2[ i ] = new double[ numberOfSpecies + 1 ];
if( jac2[ i ] == 0 ) {
cerr << "ERROR: allocation for jac2[" << i << "] failed!";
cerr << endl;
}
}
*/

// put y0 into y1 as a first guess
for( int i = 1; i <= numberOfSpecies; i++ ) {
xspec[ Y1 ][ i ] = xspec[ Y0 ][ i ];
}
#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ Y1 ] = xreac[ r ][ Y0 ];
xreac[ Y1 ][ r ] = xreac[ Y0 ][ r ];
}
#endif


// for time = 1 :
// approximate y1 with BDF1
// F( y1 ) = y1 - y0 - h * f( t1, y1 ) = 0
// xspec[ 0 ] actually holds y0, x[ 1 ] will hold y1
time = initialTime + dtime;

its = 0;
int t_index = 1;

do {
// eval fx( time, y1, temp, numberOfSpecies );
// evaluate f( t, xpec[ t ] )
start = clock();
xrate( vspec, vfor, vbak, t_index + 1 );
end = clock();
elapsed_xrate += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

// put h * f(...) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
bigF[ i ] = xspec[ Y1 ][ i ] - xspec[ Y0 ][ i ] 
- dtime * vspec[ i ];
}

start = clock();
//eval_jacobian( t_index, jac2 ); 
eval_prep_jacobian( t_index, jac , prepared_jac ); 
end = clock();
elapsed_jac += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

for( int i = 1; i <= numberOfSpecies; i++ ) {
for( int j = 1; j <= numberOfSpecies; j++ ) {
bigFprime[ i ][ j ] = ( i == j ? 1.0 : 0.0 ) 
- dtime * jac[ i ][ j ];
}
}

start = clock();
gauss( bigFprime, tempInt, numberOfSpecies );
gauss_solve( bigFprime, tempInt, bigF, m, numberOfSpecies );
end = clock();
elapsed_gauss += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

maxerror = 0.0;
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
newvalue = xspec[ Y1 ][ i ] - m[ i ];
maxerror = maxerror > myabs( xspec[ Y1 ][ i ] - newvalue ) 
? maxerror : myabs( xspec[ Y1 ][ i ] - newvalue );
xspec[ Y1 ][ i ] = newvalue;
}
}

its++;
if( its > maxIts ) {
cerr << "ERROR 1:  Max iterations reached!" << endl;
return;
}
} while( maxerror > epsilon );

//cerr << "time= " << time << ", h= " << dtime << "\tnumberOfTimeSteps = ";

totalIts += its;
//(end) approximate y1 with BDF1

// now go on with BDF2 for the remaining timeSteps
// F( y2 ) = y2 - (4/3) y1 + (1/3) y0 - (2/3) * f( t2, y2) = 0
const double a1 =  4.0 / 3.0;
const double a2 =  1.0 / 3.0;
const double b0 =  2.0 / 3.0;

for( int i = 1; i <= numberOfSpecies; i++ ) {
xspec[ Y2 ][ i ] = xspec[ Y1 ][ i ];
}

for( t_index = 2; t_index <= numberOfTimeSteps; t_index++ ) {

time = initialTime + dtime * (double) t_index;
its = 0;
do {
// eval_fx( time, y2, temp, numberOfSpecies );
// evaluate f( t, xpec[ t ] )
start = clock();
xrate( vspec, vfor, vbak, t_index + 1 );
end = clock();
elapsed_xrate += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

// put h * f(...) into k1
for( int i = 1; i <= numberOfSpecies; i++ ) {
bigF[ i ] = xspec[ Y2 ][ i ] - a1 * xspec[ Y1 ][ i ] 
+ a2 * xspec[ Y0 ][ i ] - b0 * dtime * vspec[ i ];
}

start = clock();
//eval_jacobian( t_index, jac );
eval_prep_jacobian( t_index, jac, prepared_jac ); 
end = clock();
elapsed_jac += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

for( int i = 1; i <= numberOfSpecies; i++ ) {
for( int j = 1; j <= numberOfSpecies; j++ ) {
bigFprime[ i ][ j ] = ( i == j ? 1.0 : 0.0 ) 
- b0 * dtime * jac[ i ][ j ];
}
}

start = clock();
gauss( bigFprime, tempInt, numberOfSpecies );
gauss_solve( bigFprime, tempInt, bigF, m, numberOfSpecies );
end = clock();
elapsed_gauss += ( (double) ( end - start ) ) / CLOCKS_PER_SEC;

maxerror = 0.0;
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 0 ) {
newvalue = xspec[ Y2 ][ i ] - m[ i ];
maxerror = maxerror > myabs( xspec[ Y2 ][ i ] - newvalue ) 
? maxerror : myabs( xspec[ Y2 ][ i ] - newvalue );
xspec[ Y2 ][ i ] = newvalue;
}
}

its++;
if( its > maxIts ) {
cerr << "ERROR: Max iterations reached!" << endl;
return;
}
} while( maxerror > epsilon );

totalIts += its;

Y0++;
Y1++;
Y2++;

}

/*
cerr << "time = " << time << ", h = " << dtime << "\ttimesteps = ";
cerr << numberOfTimeSteps << endl;
//      cerr << "Y0 = " << Y0 << ", Y1 = " << Y1 << ", Y2 = " << Y2 << endl;
//      for( int i = 1; i <= numberOfSpecies; i++ ) {
//         cerr << trim( nameOfSpecies[ i ] ) << "\t";
//         cerr << xspec[ Y0 ][ i ] << " \t" << xspec[ Y1 ][ i ] << endl;
//      }
*/
cerr << " its = " << totalIts;
cerr << ", jac_time = " << elapsed_jac;
cerr << ", xrate_time = " << elapsed_xrate;
cerr << ", gauss_time = " << elapsed_gauss << endl;

// free up allocated memory
for( int i = 1; i <= numberOfSpecies; i++ ) {
delete [] bigFprime[ i ];
}
delete [] bigFprime;
for( int i = 1; i <= numberOfSpecies; i++ ) {
delete [] jac[ i ];
}
delete [] jac;
delete [] bigF;
delete [] m;
delete [] temp;
delete [] tempInt;

} // end method: stiffSolver

/*************************************************************
*   r u n k i n
*************************************************************
* Runs the integration depening on integrationOption(jtime)
*/
void runkin()
{


#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
//xreac[ r ][ 0 ] = 0.0;
xreac[ 0 ][ r ] = 0.0;
}
#endif

for( int t = 1; t <= numberOfTimeSteps; t++ ) {
for( int i = 1; i <= numberOfSpecies; i++ ) {
if( jfix[ i ] == 1 ) {
xspec[ t ][ i ] = xspec[ t - 1 ][ i ];
}
}
}

if( integrationOption == 1 ) {
eulerMethod();
} 

if( integrationOption == 101 ) {
testMethod_aleman();
} 


if( integrationOption == 2 ) {
modifiedEulerMethod();
} 


if( integrationOption == 3 ) {
rungeKuttaOrder4();
} 


if( integrationOption == 45 ) {
rungeKutta45();
} 


if( integrationOption == 4 ) {
rungeKuttaAdaptive();
}


if( integrationOption == 5 ) {
stiffSolver();
}

if( integrationOption == 6 ) {
lsodesMethod();
}

}

/* Writes out the second part of kin_o01
* returns 0 when ok
*/
int outputDataFile1()
{
ofstream outputFile1( kin_o01, ios::app );
if( ! outputFile1 ) {
cerr << "Couldn't open input file " << kin_i01 << endl;
return 1;
}

outputFile1 << endl << endl;
outputFile1 << setiosflags( ios::scientific | ios::uppercase );
outputFile1 << " initial and final species concentrations" << endl << endl;
outputFile1 << "   ispec         timei         xspec" << endl << endl;

for( int i = 1; i <= numberOfSpecies; i++ ) {

outputFile1 << nameOfSpecies[ i ] << endl;

outputFile1 << setw( DEC8 ) << i;
outputFile1 << setw( FRAC ) << initialTime;
outputFile1 << setw( FRAC ) << xspec[ 0 ][ i ] << endl;

outputFile1 << setw( DEC8 ) << i;
outputFile1 << setw( FRAC ) << finalTime;
outputFile1 << setw( FRAC ) << xspec[ numberOfTimeSteps ][ i ] << endl;
}

outputFile1 << endl << endl;
outputFile1 << " final reaction concentrations" << endl << endl;
outputFile1 << "   ireac         timei         xreac" << endl;

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
outputFile1 << setw( DEC8 ) << r;
outputFile1 << setw( FRAC ) << finalTime;
//outputFile1 << setw( FRAC ) << xreac[ r ][ numberOfTimeSteps ] << endl;
outputFile1 << setw( FRAC ) << xreac[ numberOfTimeSteps ][ r ] << endl;
}
#endif

outputFile1 << endl << endl;
outputFile1 << " time-dep. species concentrations" << endl;

int mtime = -1;
int top = numberOfTimeSteps;
if( integrationOption == 4 ) { // RK Adaptive
top = xtime_index;
}
for( int i = 1; i <= numberOfSpecies; i++ ) {
mtime = -1;
outputFile1 << endl << "   ispec" << endl;
outputFile1 << setw( DEC8 ) << i << endl;

outputFile1 << " namespec:" << endl;
outputFile1 << nameOfSpecies[ i ] << endl;
outputFile1 << "   itime         timei         xspec" << endl;
for( int t = 0; t <= top; t++ ) {
mtime = mtime + 1;
if( mtime == ntskip ) {
mtime = 0;
}
if( mtime == 0 || t == top ) {
outputFile1 << setw( DEC8 ) << t;
if( integrationOption == 4 ) {
outputFile1 << setw( FRAC ) << xtime[ t ];
}
else {
outputFile1 << setw( FRAC ) << ( initialTime + t * dtime );
}
outputFile1 << setw( FRAC ) << xspec[ t ][ i ] << endl;
}
}
}

outputFile1 << endl << endl << " time-dep. reaction concentrations" << endl;

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {
outputFile1 << endl << "   ireac" << endl;
outputFile1 << setw( DEC8 ) << r << endl;
outputFile1 << "   itime         timei         xreac" << endl;
mtime = -1;
for( int t = 0; t <= top; t++ ) {
mtime = mtime + 1;
if( mtime == ntskip ) {
mtime = 0;
}
if( mtime == 0 || t == top ) {
outputFile1 << setw( DEC8 ) << t;
if( integrationOption == 4 ) {
outputFile1 << setw( FRAC ) << xtime[ t ];
}
else {
outputFile1 << setw( FRAC ) << ( initialTime + t * dtime );
}
//outputFile1 << setw( FRAC ) << xreac[ r ][ t ] << endl;
outputFile1 << setw( FRAC ) << xreac[ t ][ r ] << endl;
}
}
}
#endif

outputFile1.close();

return 0;
}

/**
* Writes out file kin_o02 
* returns 0 when ok
*/
int outputDataFile2()
{
ofstream outputFile2( kin_o02, ( numberOfDataSets == 1 ? ios::out : ios::app ) );
int mtime = -1;
int top = numberOfTimeSteps;
if( integrationOption == 4 ) { // RK Adaptive
top = xtime_index;
}

if( ! outputFile2 ) {
cerr << "Couldn't open output file " << kin_o02 << endl;
return 1;
}

outputFile2 << setiosflags( ios::scientific | ios::uppercase );

outputFile2 << "#" << endl;
outputFile2 << "#" << endl;
outputFile2 << "#" << bline[ 0 ] << endl;
outputFile2 << "#  time-dep. species concentrations" << endl;
outputFile2 << "#" << endl;

for( int i = 1; i <= numberOfSpecies; i++ ) {
//cout <<i<<endl;      
outputFile2 << "#" << endl;
outputFile2 << "#   ispec" << endl;
outputFile2 << "#" << setw( DEC8 ) << i << endl;
outputFile2 << "# namespec:" << endl;
outputFile2 << "#" << nameOfSpecies[ i ] << endl;
outputFile2 << "#" << "  itime         timei         xspec" << endl;

mtime = -1;

//cout << i << endl;

for( int i_intg = 1; i_intg <= numOfIntegrations; i_intg++){
//cout << " i_intg " << i_intg <<endl;
if(integrationOption == 4){
top = intg_xtime_index[i_intg];
//cout << " top " << top <<endl;
}
else{
top = intg_numberOfTimeSteps[i_intg];
//cout << " top " << top <<endl;
}
for( int t = 0; t <= top; t++ ) {
mtime = mtime + 1;
if( mtime == ntskip ) {
mtime = 0;
}
if( mtime == 0 || t == top ) {
outputFile2 << setw( DEC8 ) << t;
if( integrationOption == 4 ) {
outputFile2 << setw( FRAC ) << intg_xtime[ i_intg ][ t ];
} 
else {
outputFile2 << setw( FRAC ) << ( intg_initialTime[ i_intg ] + t * intg_dtime[ i_intg ] );
}
outputFile2 << setw( FRAC ) << intg_xspec[ i_intg ][ t ][ i ] << endl;
}
}//endof for t
}//endof for i_intg

outputFile2 << " " << endl;
outputFile2 << " " << endl;


}

outputFile2.close();

return 0;
}

/**
* Writes out file kin_o03
* returns 0 when ok
*/
int outputDataFile3()
{
ofstream outputFile3( kin_o03, ( numberOfDataSets == 1 ? ios::out : ios::app ) );
int mtime = -1;
int top = numberOfTimeSteps;
if( integrationOption == 4 ) { // RK Adaptive
top = xtime_index;
}

if( ! outputFile3 ) {
cerr << "Couldn't open output file " << kin_o03 << endl;
return 1;
}

outputFile3 << setiosflags( ios::scientific | ios::uppercase );

outputFile3 << "#" << endl;
outputFile3 << "#" << endl;
outputFile3 << "#  time-dep. reaction concentrations" << endl;
outputFile3 << "#" << endl;

#ifndef OPTIMIZE
for( int r = 1; r <= numberOfReactions; r++ ) {

outputFile3 << "#" << endl;
outputFile3 << "#   ireac" << endl;
outputFile3 << "#" << setw( DEC8 ) << r << endl;
outputFile3 << "#" << "  itime         timei         xreac" << endl;

mtime = -1;
for( int t = 0; t <= top; t++ ) {

mtime = mtime + 1;
if( mtime == ntskip ) {
mtime = 0;
}
if( mtime == 0 || t == top ) {
outputFile3 << setw( DEC8 ) << t;
if( integrationOption == 4 ) {
outputFile3 << setw( FRAC ) << xtime[ t ];
}
else {
outputFile3 << setw( FRAC ) << ( initialTime + t * dtime );
}
//outputFile3 << setw( FRAC ) << xreac[ r ][ t ] << endl;
outputFile3 << setw( FRAC ) << xreac[ t ][ r ] << endl;
}
}

outputFile3 << " " << endl;
outputFile3 << " " << endl;

}
#endif

outputFile3.close();

return 0;
}

// Output: 
// --reaction forward + backward reaction 
//   rates "vfor", "vbak" for each reaction "ireac"
// --net production rate "vspec" for each species "ispec",
void xrate( double vspec[], double vfor[], double vbak[], int t )
{
int nzyme = 0;
double fmm = 0.0;

// calculate each reaction's forward and backward reaction rates
for( int r = 1; r <= numberOfReactions; r++ ) {

nzyme = jkin[ r ] == 11 ? 1 :0;

vfor[ r ] = forwardReactionRates[ r ];
for( int q = ( 1 + nzyme ); q <= numberInputParticipants[ r ]; q++ ) {
vfor[ r ] = vfor[ r ] * xspec[ t - 1 ][ iispec[ q ][ r ] ];
}

vbak[ r ] = backwardReactionRates[ r ];
if( jkin[ r ] == 11 ) {
vbak[ r ] = backwardReactionRates2[ r ];
}

for( int p = ( 1 + nzyme ); p <= numberOutputParticipants[ r ]; p++ ) {
vbak[ r ] = vbak[ r ] * xspec[ t - 1 ][ iospec[ p ][ r ] ];
}

if( jkin[ r ] == 11 ) {
fmm = xspec[ t - 1 ][ iispec[ 1 ][ r ] ]
/ ( vfor[ r ] + vbak[ r ] + backwardReactionRates[ r ]
+ forwardReactionRates2[ r ] );
vfor[ r ] = fmm * forwardReactionRates2[ r ] * vfor[ r ];
vbak[ r ] = fmm * backwardReactionRates[ r ] * vbak[ r ];
}
}

// calculate each species' net production rate
// (= total creation rate - total annihilation rate)
for( int i = 1; i <= numberOfSpecies; i++ ) {
vspec[ i ] = 0.0;
}

for( int r = 1; r <= numberOfReactions; r++ ) {

for( int q = 1; q <= numberInputParticipants[ r ]; q++ ) {
vspec[ iispec[ q ][ r ] ] = vspec[ iispec[ q ][ r ] ] - vfor[ r ]
+ vbak[ r ];
}

for( int p = ( 1 + nzyme ); p <= numberOutputParticipants[ r ]; p++ ) {
vspec[ iospec[ p ][ r ] ] = vspec[ iospec[ p ][ r ] ] + vfor[ r ]
- vbak[ r ];
}

}

}

//removes the spaces of a string of characters
char* trim( const char *str )
{
static char ret[ LINE ] = { 0 };
int q = 0;

for( int i = 0; i < strlen( str ); i++ ) {
if( str[ i ] != ' ' ) {
ret[ q++ ] = str[ i ];
}
}
ret[ q ] = 0;
char * value = new char[ q ];
strcpy( value, ret );
return value;
}

// returns the index where the specie name is 
int ispec4name( const char *name2look )
{
int index = 0;
if( name2look != 0 ) {
for( int i = 1; i <= numberOfSpecies && index == 0; i++ ) {
if( strstr( nameOfSpecies[ i ], name2look ) != NULL 
&& strstr( name2look, trim( nameOfSpecies[ i ] ) ) != NULL ) {
index = i;
}
}
}
if( index == 0 ) {
cerr << " ERROR" << endl;
cerr << " Could not analyze following species name " << endl;
cerr << " name2look=" << name2look << endl;
}
return index;
}

/*************************************************************
* m y a b s
*************************************************************
* returns the absolute value of a double number
*/
double myabs( double number )
{
if( number < 0.0 )
return -1.0 * number;
return number;
}

/*************************************************************
* f r e e M e m o r y
*************************************************************
* Clears up memory allocated with 'new'
*/
void freeMemory()
{
for( int i = 0; i < ( 10 + MAX_NUMBER_REACTIONS ); i++ ) {
if( bline[ i ] != 0 ) { 
delete [] bline[ i ];
}
}
for( int i = 0; i < ( 1 + MAX_NUMBER_SPECIES ); i++ ) {
if( nameOfSpecies[ i ] != 0 ) {
delete [] nameOfSpecies[ i ];
}
}
for( int i = 0; i < ( MAX_NUMBER_PARTICIPANT_REACTIONS + 1 ); i++ ) {
for( int j = 0; j < ( MAX_NUMBER_REACTIONS + 1 ); j++ ) {
if( nameOfOutputSpecies[ i ][ j ] != 0 ) {
delete [] nameOfOutputSpecies[ i ][ j ];
}
if( nameOfInputSpecies[ i ][ j ] != 0 ) {
delete [] nameOfInputSpecies[ i ][ j ];
}
}
}
for( int i = 0; i < ( MAX_NUMBER_TIME_STEPS + 10 ); i++ ) {
delete [] xspec[ i ];
}
delete [] xspec;
#ifndef OPTIMIZE
for( int i = 0; i < ( MAX_NUMBER_TIME_STEPS + 10 ); i++ ) {
delete [] xreac[ i ];
}
delete [] xreac;
#endif
delete [] xtime;

delete [] jfix;
delete [] numberInputParticipants;
delete [] numberOutputParticipants;
delete [] jkin;
delete [] initialConcentration;
delete [] forwardReactionRates;
delete [] backwardReactionRates;
delete [] forwardReactionRates2;
delete [] backwardReactionRates2;

}

/*************************************************************
*   g a u s s
*************************************************************
* gauss elimination
* d will hold the scaled index array
* a values are modified
* page 260. Ward Cheney, David Kincaid
* the arrays and matrix a are used with indices [ 1..numEquations+1 ]
*/
void gauss( double **a, int d[], int n )
{
int    temp  = 0;
double smax  = 0.0;
double rmax  = 0.0;
double xmult = 0.0;
double r     = 0.0;
double *s = new double[ n + 1 ];   // scale array

for( int i = 1; i <= n; i++ ) {
d[ i ] = i;
smax = 0.0;
for( int j = 1; j <= n; j++ ) {
smax = smax > myabs( a[ i ][ j ] ) ? smax : myabs( a[ i ][ j ] );
}
s[ i ] = smax;
}

int j = 1;

for( int k = 1; k <= ( n - 1 ); k++ ) {

rmax = 0.0;
j = 1;

for( int i = k; i <= n; i++ ) {
r = myabs( a[ d[ i ] ][ k ] / s[ d[ i ] ] );
if( r > rmax ) {
rmax = r;
j = i;
}
}
temp = d[ k ];
d[ k ] = d[ j ];
d[ j ] = temp;

for( int i = k + 1; i <= n; i++ ) {
xmult = a[ d[ i ] ][ k ] / a[ d[ k ] ][ k ];
a[ d[ i ] ][ k ] = xmult;
for( j = k + 1; j <= n; j++ ) {
a[ d[ i ] ][ j ] = a[ d[ i ] ][ j ] - xmult * a[ d[ k ] ][ j ];
}
}
}
delete [] s;
} // end method: gauss


/*************************************************************
*   g a u s s
*************************************************************
* gauss substitution Ax=b
* d has the scaled index array
* page 262. Ward Cheney, David Kincaid
* the arrays and matrix a are used with indices [ 1..numEquations+1 ]
*/
void gauss_solve( double **a, int d[], double b[], double x[], int n )
{
double sum = 0.0;
for( int k = 1; k <= ( n - 1 ); k++ ) {
for( int i = k + 1; i <= n; i++ ) {
b[ d[ i ] ] = b[ d[ i ] ] - a[ d[ i ] ][ k ] * b[ d[ k ] ];
}
}
x[ n ] = b[ d[ n ] ] / a[ d[ n ] ][ n ];
for( int i = n - 1; i >= 1; i-- ) {
sum = b[ d[ i ] ];
for( int j = i + 1; j <= n; j++ ) {
sum -= a[ d[ i ] ][ j ] * x[ j ];
}
x[ i ] = sum / a[ d[ i ] ][ i ];
}
}

int mystrcmp(const char *str1, const char *str2)
{

//   cout << str1 << "  " << str2 <<endl;
int cmpresult = 1;   
int i = 0;
int j = 0;

while(str1[i]==' ') i++;
while(str2[j]==' ') j++;

for(;(str1[i]==str2[j])&&(str1[i]!='\0')&&(str2[j]!='\0')&&(str1[i]!=' ')&&(str2[j]!=' ');){
i++;j++;
}

if((str1[i]=='\0')&&(str2[j]=='\0')){   
cmpresult = 0;
}
if((str1[i]=='\0')&&(str2[j]==' ')){
cmpresult = 0;
}
if((str1[i]==' ')&&(str2[j]=='\0')){
cmpresult = 0;
}
if((str1[i]==' ')&&(str2[j]==' ')){
cmpresult = 0;
}
return cmpresult;

}





