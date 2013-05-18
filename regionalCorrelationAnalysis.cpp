/* regionalCorrelationAnalysis.cpp
 * 
 * 
 *
 * This program implements 2D regional correlation analysis, an algorithm proposed in the paper:
 *
 * Xuefei Wang, H. Peter Lu (2008) "2D Regional Correlation Analysis of Single-Molecule Trajectories," J. Phys. Chem. B. 112, 1492-14926
 *
 *
 * Edited:
 * 	Kristofer Gryte - 2012-06-25 - Created.
 * 
 * Compilation notes:
 * 	- Need to compile using -lm option so sqrt() can take variable input.
 * 
 * 	- LINUX:
 * 	- ...to create the executable...
 * 	 	- g++ -o regionalCorrelationAnalysis regionalCorrelationAnalysis.cpp -lm
 * 	- ...to run...
 * 		- ./regionalCorrelationAnalysis <inputDataFilename> <outputFilename>
 * 			e.g., ./regionalCorrelationAnalysis intputData.dat output.dat
 * 
 * 
 * Contact:
 * 	Kristofer Gryte
 * 	University of Oxford
 * 	Oxford, UK
 * 	k.gryte1(at)physics.ox.ac.uk
 * 
 * 
 * 	Copyright (c) 2012
 * 
 * 
 ************************************************************************************************** */

/* HEADER FILES */
#include <iostream>

using std::cout;
using std::ios;
using std::endl;
using std::cerr;

#include <iomanip>

using std::setw;

#include <cstdlib> 	// Need for exit()

using std::exit;

#include <cmath>  	// Need for math stuff; sqrt()

using std::sqrt;

#include <fstream> 	// Need for handling files

using std::fstream;
using std::ofstream;

//#include "mex.h"	// MATLAB Mex header file






/* PROTOTYPES */

void checkInputs(int x, char* y, char* z); // prototype for checkInputs()

int getData(fstream& s); // prototype for getData()

double mean(double x[], int y); // prototype for mean()

void fluctuationAmplitude(double x[], int y, int z); // prototype for fluctuationAmplitude()

double pixelAverage(int x, int y, int z); //prototype for pixelAverage()

double crossCorrelation(int x, int y, int z); // prototype for crossCorrelation()

void writeToFile(ofstream& s, double** x, int y, int z); // prototype for writeToFile()





/* Definitions */
#define MAXNUMDATAPTS 50000

typedef double* doubleArrayPtr; // pointer to a double array



/* GLOBAL VARIABLES */
double data[MAXNUMDATAPTS][2]; // DataX, DataY ---> Fluctation Amp X, Fluctuation Amp Y


// Constants:
const int pixelSize = 3;
  
// Create object for the files to read and write:
fstream inDataFile;
ofstream outDataFile;
  







/* MAIN PROGRAM */
int main(int argc, char *argv[])
{ 
  
  cout << "Beginning analysis.\n";
  
  /////////////////////////////////////////////////////////////////////////
  /* CHECKS */
    
  checkInputs(argc, argv[1], argv[2]);
  
  
  //////////////////////////////////////////////////////////////////////////
  /* READ IN DATA */

  int numData = getData(inDataFile); // the number of data pts is equal to the number of lines in the input data file
  
  cout << "Number of Data Points: " << numData << endl;
  
  
  //////////////////////////////////////////////////////////////////////////
  /* VARIABLE DECLARATIONS */
    
  int numRows = numData - pixelSize - 1;
    
  double meanVals[2]; 
  
  int startTime, stopTime;  
  
   
  // Build our cross-correlation array dynamically:
  doubleArrayPtr *crossCorr;
  crossCorr = new doubleArrayPtr[numRows];
  int i;
  for (i = 0; i < numRows; i++) {
    crossCorr[i] = new double[numData]; 
  } // end FOR crossCorr --> crossCorr is a numRows x numData array. This will be our Cross-correlation array (our image)
  
  cout << "Finished dynamically creating image array.\n";
  
  
  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////
  /* ANALYSIS */
  
  /* Cross Correlation Calculation */ 
  cout << "Starting Cross-Correlation Calculation.\n";
  
  // SAMPLE MEANS
  
  for (i=0; i<2; i++) {    
      
    int n;
    double tempData[numData];
    
    for (n=0; n < numData; n++) {
     tempData[n] = data[n][i];           
    } // end FOR n
    
    meanVals[i] = mean(tempData, numData); // 
    
    cout << "Mean of Data Set " << i << ": " << meanVals[i] << endl;    
    
  } // end FOR i
    
    
    
    
  // FLUCTUATION AMPLITUDE:
  
  // Get the fluctuation amplitude (subtract the mean):
  fluctuationAmplitude(meanVals, numData, 2); // 2 data columns 
  
  
   
  
 
  /* Outer loop for our StartTimes */
  for (startTime = 0; startTime < numRows; startTime++) { 
   
    /* Inner loop for our StopTimes */
    for (stopTime = startTime + pixelSize; stopTime < startTime + 11; stopTime++) { // ;stopTime < numData; 

      if (stopTime == (numData-1)) {
	break; // break the loop as stopTime will exceed the allocated memory for crossCorr
      }
      
      /* Calculate the cross-correlation for the region of interest (pixel average), building our cross-correlation image */
      crossCorr[startTime][stopTime] = pixelAverage(startTime, stopTime, pixelSize);    
            
    } // end FOR stopTime
   
  } // end FOR startTime
  
  cout << "End of Cross-Correlation Calculation. \n";
  
  
  
  
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  /* WRITE TO FILE: */
  
  
  writeToFile(outDataFile, crossCorr, numRows, numData); 
  
  
  
  
  /////////////////////////////////////////////////////////////////////////////////////
  /* END ANALYSIS */
  
  cout << "Closing files. \n";
  
  // Close the Data Files:
  inDataFile.close( );
  outDataFile.close( );
  
  cout << "End of Regional Correlation Analysis. \n";
  
  return 0;
  
  
} // end function MAIN










/* ************************************************************************************* */
/* SUB-FUNCTIONS:*/



// FUNCTION: mean()
double mean(double dataVector[], int numData)
{
  
 /* Calculate the mean of data */
 int i;
 double sum = 0;
   
 for(i=0; i<numData; i++) {   
   // Compuete sum:
   sum +=  dataVector[i];
   
 } // end FOR
  
 /* Sample mean */
 return sum / numData;
 
 
} // end function MEAN



// FUNCTION: fluctuationAmplitude()
void fluctuationAmplitude(double meanVals[], int numRows, int numCols)
{
  int i,j;
  
  // Subtract the mean...
  for (i=0; i < numRows; i++) {
    
    for (j=0; j < numCols; j++) {      
      
      data[i][j] = (data[i][j] - meanVals[j]); // subtract the mean, getting the fluctuation amplitude      
      
    } // end FOR j
    
  } // end FOR i  
  
} // end function FLUCTUATIONAMPLITUDE




// FUNCTION: pixelAverage()
double pixelAverage(int startTime, int stopTime, int pixelSize)
{
  int deltaPixels;
  double localAvg[pixelSize+1]; 	
  
  /* Loop for average (pixel) correlation */
  for (deltaPixels = 0; deltaPixels < pixelSize+1; deltaPixels++) {
    
    /* Compute the pixel (correlation) average: */
    localAvg[deltaPixels] = crossCorrelation(startTime, stopTime, deltaPixels);  
	    
  } // end FOR deltaPixels
  
  return mean(localAvg, pixelSize);

} // end function PIXELAVERAGE




// FUNCTION: crossCorrelation()
double crossCorrelation(int startTime, int stopTime, int deltaPixels)
{
  /* Loop for actual correlation calculations */
  double autoCorrX=0, autoCorrY=0, corrXY=0;
  int i;
  for (i = startTime; i < stopTime-deltaPixels+1; i++) { 
    
    int j = i + deltaPixels;
    
    autoCorrX += data[i][0] * data[i][0]; // square the fluctuation amplitude in X
    autoCorrY += data[j][1] * data[j][1]; // square the fluctuation amplitude in Y
    
    corrXY += data[i][0] * data[j][1]; // calculate the cross fluctuation amplitude XY
    
  } // end FOR i
  
  //cout << autoCorrX << " " << autoCorrY << "\n";
  
  //cout << corrXY << endl;
  
  
  return corrXY / sqrt(autoCorrX * autoCorrY);
  
} // end function CROSSCORRELATION






// FUNCTION: checkInputs()
void checkInputs(int numArgs, char* inputDataFilename, char* outputDataFilename)
{
 
  /* Parse the input arguments */
  if (numArgs<2) {
    cerr << "ERROR: insufficient input arguments. \n";
    exit(1); /* Kill the program */
  }
  
  
  inDataFile.open(inputDataFilename, ios::in); // Open the file to be read
  
  if ( inDataFile.fail( ) ) {
    cerr << "Input data file could not be opened \n";
    exit( 1 );
  }
  
  outDataFile.open(outputDataFilename, ios::out); // Open the file to written to
  
  if ( outDataFile.fail( ) ) {
    cerr << "Output data file could not be opened \n";
    exit( 1 );
  }
  
  cout << "Opened files successfully.\n";
  
} // end function CHECKINPUTS





// FUNCTION: getData()
int getData(fstream& inDataFile) 
{
  int line_counter; // keep track of the number of lines of data in the input file

  line_counter = 0;
  while (!inDataFile.eof( ) && line_counter < MAXNUMDATAPTS) {
    
    // Get the data from the file:
    inDataFile >> data[line_counter][0]; // first number in the row    
    inDataFile >> data[line_counter][1]; // second number in the row
        
    // Increment the line counter:
    line_counter++;
    
  }
  
  // Decrement the line counter by one, as we have over counted:
  line_counter--;
  
  cout << "Finished getting data.\n";
  
  return line_counter;
  
} // end function GETDATA




// FUNCTION: writeToFile()
void writeToFile(ofstream& filename, double** dataArray, int numRows, int numCols){
  
  cout << "Writing to output file. \n";
  
  int i,j;
  
  for (i=0; i < numRows; i++) {
      
   for (j=0; j < numCols; j++) {
      filename << setw(12) << dataArray[i][j];
   } // end FOR j
   
   filename << endl;
  } // end FOR i
   
} // end function WRITETOFILE
