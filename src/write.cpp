#include <iostream>
#include <stdlib.h>
#include <thread>
#include "plplot/plplot.h"
#include "plplot/plstream.h"
//#include "plc++demos.h"


class x00 {
public:
    x00( int, const char ** );

private:
    // Class data
    plstream         *pls;

    static const int NSIZE;
};

const int x00::NSIZE = 101;

x00::x00( int argc, const char **argv )
{
    PLFLT x[NSIZE], y[NSIZE];
    PLFLT xmin = 0., xmax = 1., ymin = 0., ymax = 100.;
    int   i;

    // Prepare data to be plotted.
    for ( i = 0; i < NSIZE; i++ )
    {
        x[i] = (PLFLT) ( i ) / (PLFLT) ( NSIZE - 1 );
        y[i] = ymax * x[i] * x[i];
    }

    pls = new plstream();


    // Parse and process command line arguments
    pls->parseopts( &argc, argv, PL_PARSE_FULL );

    // Initialize plplot
    //	pls->sdev("qtwidget");
    pls->init();

    // Create a labelled box to hold the plot.
    pls->env( xmin, xmax, ymin, ymax, 0, 0 );
    pls->lab( "x", "y=100 x#u2#d", "Simple PLplot demo of a 2D line plot" );

    // Plot the data that was prepared above.
    pls->line( NSIZE, x, y );

	x[20] = 0;
	x[21] = 0;
	y[20] = 0;
	y[21] = 0;

	int j = 0;
    std::chrono::milliseconds dura( 200 );
	while(true){
		j++;


		std::this_thread::sleep_for(dura);
		std::cout << "test" << std::endl;

		for (uint i = 0; i < NSIZE; i++ )
		{
			x[i] = (PLFLT) ( i ) / (PLFLT) ( NSIZE - 1 );
			y[i] = j * ymax * x[i] * x[i];
		}

		pls->line( NSIZE, x, y );
		if (j==5){ break;}

		pls->ptex(0.5, 0.5, 0, 0, 0, "test");

		pls->flush();
	}


    // In C++ we don't call plend() to close PLplot library
    // this is handled by the destructor
    delete pls;
}

int main( int argc, const char ** argv )
{
    x00 *x = new x00( argc, argv );
    delete x;
}


//--------------------------------------------------------------------------
//                              End of x00.cc
//--------------------------------------------------------------------------
