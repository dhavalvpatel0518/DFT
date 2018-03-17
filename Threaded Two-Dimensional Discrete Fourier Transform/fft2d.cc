// Distributed two-dimensional Discrete FFT transform
// EDITIED by Dhaval Patel

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <cmath>

#include "Complex.h"
#include "InputImage.h"

constexpr unsigned int NUMTHREADS = 4;

using namespace std;

//undergrad students can assume NUMTHREADS will evenly divide the number of rows in tested images
//graduate students should assume NUMTHREADS will not always evenly divide the number of rows in tested images.
// I will test with a different image than the one given
void Transform1D(Complex* h, int w, Complex* H)
{
    // Implement a simple 1-d DFT using the double summation equation
    // given in the assignment handout.  h is the time-domain input
    // data, w is the width (N), and H is the output array.
    Complex _W;
    Complex s(0,0);

    for (int n = 0; n <= w-1; n++)
    {
        for (int k = 0; k <= w-1; k++)
        {
            _W = Complex(cos(2*M_PI*n*k/w),-sin(2*M_PI*n*k/w));
            s=s+(_W*h[k]);
        }
        H[n].imag=s.imag;
        H[n].real=s.real;
        s= Complex(0,0);
    }

}

void Transform1D_inv(Complex* h, int w, Complex* H)
{
    //run under fist Transform then this function
    //Implementing th column transform
    Complex _W;
    Complex s(0,0);

    for (int n = 0; n <= w-1; n++)
    {
        for (int k = 0; k <= w-1; k++)
        {
            _W = Complex(cos(2*M_PI*n*k/w),sin(2*M_PI*n*k/w));
            s=s+_W*h[k];
        }
        H[n]=s;
        H[n].imag=s.imag/w;
        H[n].real=s.real/w;
        s= Complex(0,0);
    }
}

void transpose(Complex* _i, Complex* o, int w, int h)
{
    int k = 0;
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++)
        {
         o[k]= _i[i+j*w];
         k++;
        }

    }
}

void Transform2D(const char* inputFN) 
{ // Do the 2D transform here.
  // 1) Use the InputImage object to read in the Tower.txt file and
  //    find the width/height of the input image.
  // 2) Create a vector of complex objects of size width * height to hold
  //    values calculated
  // 3) Do the individual 1D transforms on the rows assigned to each thread
  // 4) Force each thread to wait until all threads have completed their row calculations
  //    prior to starting column calculations
  // 5) Perform column calculations
  // 6) Wait for all column calculations to complete
  // 7) Use SaveImageData() to output the final results

  // 1)
  InputImage image(inputFN);  // Create the helper object for reading the image
  // Step (1) in the comments is the line above.
  // Your code here, steps 2-7
  // 2)
  int w; int h; int t;
  w = image.GetWidth();
  h = image.GetHeight();
  t=w*h;
  Complex* _1Dw= new Complex[t];
  Complex* _1Dwt= new Complex[t];
  Complex* _1Dh= new Complex[t];
  Complex* _1Dht= new Complex[t];

  //DFT

  // 3) & 4) have each thread join after the previous one is done.
    int rows= (h/NUMTHREADS);
    int last=rows+rows+rows-1;
    Complex* _in = image.GetImageData(); //Get pointer of the first value in array


    //calculate rows
    for (int i = 0; i <= rows; i++) {
        thread _1(&Transform1D,(_in+(i)*w), w, (_1Dw+(i)*w));
        thread _2(&Transform1D,(_in+(i+rows)*w), w, (_1Dw+(i+rows)*w));
        thread _3(&Transform1D,(_in+(i+rows+rows)*w), w, (_1Dw+(i+rows+rows)*w));
        thread _4(&Transform1D,(_in+(i+last)*w), w, (_1Dw+(i+last)*w));
        _1.join();
        _2.join();
        _3.join();
        _4.join();
    }

    //image.SaveImageData("../MyAfter1D.txt",_1Dw,w,h);
    transpose(_1Dw,_1Dwt,w,h);
    //image.SaveImageData("../MyAfter1Dt.txt",_1Dwt,w,h);
    // 5)
    //calculate column
    thread _1(&Transform1D,(_1Dwt), w, (_1Dh));
    _1.join();
    for (int i = 0; i <= rows-1; i++) {
        thread _1(&Transform1D,(_1Dwt+(i*w)+w)          , w,   (_1Dh+(i*w)+w));
        thread _2(&Transform1D,(_1Dwt+(i+rows)*w+w)     , w, (_1Dh+(i+rows)*w+w));
        thread _3(&Transform1D,(_1Dwt+(i+rows+rows)*w+w), w, (_1Dh+(i+rows+rows)*w+w));
        thread _4(&Transform1D,(_1Dwt+(i+last)*w+w)     , w, (_1Dh+(i+last)*w+w));
        _1.join();
        _2.join();
        _3.join();
        _4.join();
    }

    //The DFT is completed

    transpose(_1Dh,_1Dht,w,h);
    image.SaveImageData("../MyAfter2D.txt",_1Dht,w,h);
    // Inverse DFT
    Complex* _1Dwi= new Complex[t];
    Complex* _1Dwti= new Complex[t];
    Complex* _1Dhi= new Complex[t];



    //do the top of inverse
    thread _2(&Transform1D_inv,_1Dh, w, (_1Dwi));
    _2.join();
    //take care of imaginary precision error
    for (int j = 0; j <255 ; j++) {
        _1Dwi[j].imag=0;
    }
    //Get the 1D iDFT.
    //This part of the code should result in same values of a transposed 1d DFT output
    for (int i = 0; i <= rows-1; i++)
    {
        thread _1(&Transform1D_inv,_1Dh+(i*w)+w         , w, (_1Dwi+(i*w)+w));
        thread _2(&Transform1D_inv,(_1Dh+(i+rows)*w+w)     , w, (_1Dwi+(i+rows)*w+w));
        thread _3(&Transform1D_inv,(_1Dh+(i+rows+rows)*w+w), w, (_1Dwi+(i+rows+rows)*w+w));
        thread _4(&Transform1D_inv,(_1Dh+(i+last)*w+w)     , w, (_1Dwi+(i+last)*w+w));
        _1.join();
        _2.join();
        _3.join();
        _4.join();
    }

    //image.SaveImageData("../MyAfter2Di.txt",_1Dwi,w,h);

    //take the transpose of column inverse output
    transpose(_1Dwi,_1Dwti,w,h);

    //image.SaveImageData("../MyAfter2Dti.txt",_1Dwti,w,h);

    //find the original row values
    for (int i = 0; i <= rows; i++)
    {
        thread _1(&Transform1D_inv,_1Dwti+(i)*w             , w, _1Dhi+(i)*w);
        thread _2(&Transform1D_inv,(_1Dwti+(i+rows)*w)      , w, (_1Dhi+(i+rows)*w));
        thread _3(&Transform1D_inv,(_1Dwti+(i+rows+rows)*w) , w, (_1Dhi+(i+rows+rows)*w));
        thread _4(&Transform1D_inv,(_1Dwti+(i+last)*w)      , w, (_1Dhi+(i+last)*w));
        _1.join();
        _2.join();
        _3.join();
        _4.join();
    }

    image.SaveImageData("../MyAfterInverse.txt",_1Dhi,w,h);


}


int main(int argc, char** argv)
{
  string fn("../Tower.txt"); // default file name
  if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
  Transform2D(fn.c_str()); // Perform the transform.
}  
  

  
