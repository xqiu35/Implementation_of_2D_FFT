CXX      = /usr/bin/g++
CXXFLAGS = -Wall -g

pthreadFFT2d:	pthreadFFT2d.o Complex.o InputImage.o
	$(CXX) -g -o pthreadFFT2d pthreadFFT2d.o Complex.o InputImage.o -lpthread

clean:
	@rm *.o pthreadFFT2d
