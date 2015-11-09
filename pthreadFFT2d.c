// Threaded two-dimensional Discrete FFT transform
// Xiaofei Qiu
// ECE4122 Project 2


#include <iostream>
#include <string>
#include <math.h>
#include <pthread.h>
#include <cstdlib>
#include <algorithm>

#include "Complex.h"
#include "InputImage.h"


using namespace std;


/////////      Routines     //////////

unsigned ReverseBits(unsigned v);
void Transform1D(Complex* h, int N);
void* Transform2DTHread(void* v);
void Transform2D(const char* inputFN);
void Reorder(Complex*);                         // reverse binary indexing
Complex* preCalW();                             // cals W, and return a array
Complex  getW(const int&, const int&);          // get W value depens on n,N
void transpose(Complex* source,const int& rows, const int& cols);
void Print(const Complex*,const int&);          // for debug print
/////////  Global Variables //////////

int Npt;        // N point FFT
int width;      // image width
int height;     // image height
int row_per_cpu; // row per cpu
Complex *data;  // image data

Complex* W;     // contain for pre calculated W

static const int num_threads = 16;      // number of threads
int stopCount = 0;                      // tracking num of stopped thread
int transposeReady = 0;                 // signal for transpose ready

pthread_mutex_t exitMutex;              // exit mutex
pthread_cond_t  exitCond;               // exit condition
pthread_mutex_t stopCountMutex;         // mutex for stop counter
pthread_mutex_t transposeMutex;         // ued for lock transpose
pthread_barrier_t barrier;              // used for waiting all threads to finish 1d;

Complex *garbage[2];

int main(int argc, char** argv)
{
        string fn("Tower.txt"); // default file name
        if (argc > 1) fn = string(argv[1]);  // if name specified on cmd line
        Transform2D(fn.c_str()); // Perform the transform.
}

// Function to reverse bits in an unsigned integer
// This assumes there is a global variable N that is the
// number of points in the 1D transform.
unsigned ReverseBits(unsigned v)
{ //  Provided to students
  unsigned n = Npt; // Size of array (which is even 2 power k value)
  unsigned r = 0; // Return value

  for (--n; n > 0; n >>= 1)
    {
      r <<= 1;        // Shift return value
      r |= (v & 0x1); // Merge in next bit
      v >>= 1;        // Shift reversal value
    }
  return r;
}

void Transform1D(Complex* h, int N)
{
        Complex result[width];          // temperart result container;
        Complex w;                      // W

        Reorder(h);                     // reverse binary indexing;

        //Print(h,8);
        int x=2;

        while(x<=N)                     // star from 2 point FFT, loop until finish N point FFT
        {
                for(int i = 0; i < N/x; i++ )                   // for each x point FFT
                {
                        for(int j = 0; j < x/2; j++)            // for each element in x point  FFT
                        {
                                w = getW(j,x);                  // get W of jth element, x point FFT
                                result[i*x+j]   = h[i*x+j] + w * h[i*x+j+x/2];  // FFT formula
                                result[i*x+j+x/2] = h[i*x+j] - w * h[i*x+j+x/2];

                        }
                }
                copy(result,result+width,h);
                x*=2;
        }

}

void* Transform2DTHread(void* v)
{       
        long id =(long) v;

        int start_row = id * row_per_cpu;

        for(int i=0;i<row_per_cpu;i++)          // first 1D fft
        {
                Transform1D(&data[(start_row+i)*width],Npt);
        }

        pthread_barrier_wait(&barrier);         // wait all threads finish 1D;  

        pthread_mutex_lock(&transposeMutex);    // only single cpu does transpose
        if(!transposeReady)
        {
                transpose(data,height,width);           // transpose back
        }
        transposeReady = 1;

        pthread_mutex_unlock(&transposeMutex);

        for(int i=0;i<row_per_cpu;i++)          // second 1D fft
        {
                Transform1D(&data[(start_row+i)*width],Npt);
        }

        pthread_mutex_lock(&stopCountMutex);    // Lock here to update sopt count
        stopCount++;

        if(stopCount==num_threads)              // if all finished
        {

                transpose(data,height,width);           // transpose back
                pthread_mutex_unlock(&stopCountMutex);
                pthread_mutex_lock(&exitMutex);
                pthread_cond_signal(&exitCond);
                pthread_mutex_unlock(&exitMutex);
        }
        else
        {
                pthread_mutex_unlock(&stopCountMutex);
        }

  return 0;
}

void Transform2D(const char* inputFN)
{
        InputImage image(inputFN);              // Create the helper object for reading the image
        data   = image.GetImageData();          // Put data in global buffer
        height = image.GetHeight();             // Get h
        width  = image.GetWidth();              // Get w

        garbage[0] = data;                      // we will delete data once program finished
        garbage[1] = W;                         // we will delete data once program finished

        Npt = width;                            // init value N
        W = preCalW();                          // pre calculate W

        row_per_cpu = height / num_threads;     // init row_per_cpu

        pthread_mutex_init(&exitMutex,0);       // init exitMutex
        pthread_mutex_init(&stopCountMutex,0);  // init stop counter
        pthread_mutex_init(&transposeMutex,0);  // init transpose
        pthread_cond_init(&exitCond,0);         // init exit condition
        pthread_barrier_init(&barrier,0,16);    // init barrier

        pthread_t threads[num_threads];         // threads' ids
        int rc;                                 // return code
        for(long i=0;i<num_threads;i++)
        {
                rc = pthread_create(&threads[i],0,Transform2DTHread,(void*)i);  // create thread
                if(rc)
                {
                        cout<<"Fail to create thread. rc: "<<rc<<endl;
                        exit(1);
                }
        }

        pthread_cond_wait(&exitCond,&exitMutex);                // wait untill all thread finished
        image.SaveImageData("Tower-DFT2D.txt",data,width,height);       // write to file
        delete [] garbage[0];                                   // delete garbage
        delete [] garbage[1];
}

void Reorder(Complex* input)
{
        Complex temp[Npt];
        copy(input,input+Npt,temp);
        for(int i = 1; i < Npt; i++)
        {
                int j = ReverseBits(i);
                swap(temp[j],input[i]);
        }
}

Complex* preCalW()
{
        Complex* w = new Complex[Npt/2];                // allc N/2 for W

        for(int n = 0; n<Npt/2; n++)
        {
                w[n].real = cos(2*M_PI*n/Npt);
                w[n].imag = -sin(2*M_PI*n/Npt);
        }
        return w;
}


Complex  getW(const int& n,const int& x)        // n = nth sample, x = x point fft transform
{
        int index = n*Npt/x;                    // cal index
       if(index<Npt/2) return W[index];        // in index is within the range, return correct W
        else                                    // else return negative W
        {
                Complex temp;
                temp.real = -W[index - Npt/2].real;
                temp.imag = -W[index - Npt/2].imag;
                return temp;
        }
}

void transpose(Complex* source,const int& rows, const int& cols)
{
    Complex *temp = new Complex[cols*rows];
    int r = 0; int c = 0; int index = 0;

    for(int i =0;i<cols*rows;i++)
    {
        temp[i] = source[i];
    }

    for(int j = 0; j<cols*rows;j++)
    {
        r = j/cols;
        c = j%cols;
        index = c*rows + r;
        source[index] = temp[j];
    }
    delete [] temp;
}

void Print(const Complex* c,const int& l)
{
        for(int i=0;i<l;i++)
        {
                c[i].Print();
                cout<<endl;
        }
}

