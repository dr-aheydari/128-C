//
//
// This program calculates the error bound for problem 7
//

#include <iostream>
#include <math.h>
using namespace std;

int main()

{
// setting all the values from the problem

double h = 0.001; // 10^-3
int M = 562000; // calculated for a = 5
int L = 1000; // a = 5
int a = 5;
int b = -5;
// finding N
int N = (a - b) / h;

double err_bd = 0;

int t;
int counter = 0;
int modulo = 100;

for (int i=0 ; i <= N;i++)

{
// finding t_i

t = a + (i * h);

err_bd = ((h*M)/(2*L)) * (exp(L*(t - a)) - 1);

// to show every 100th value
if ( i % modulo == 0)
{

cout << " error is " << err_bd << endl;
counter++;

}

}

cout << "the program showed only " << counter << " values" << endl ;
cout << " for every " << modulo << " t_i 's" << endl ;

}

