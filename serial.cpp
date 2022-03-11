#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );
    double max_velocity = getCutoff();
    double size2 = getSize();
    int nrElem = (int) (size2/max_velocity) + 1;
    struct matrixElement matrix[nrElem][nrElem];

    struct matrixElement **arg;
    struct matrixElement *arg2[nrElem];
    arg = arg2;
    for(int i = 0; i < nrElem; i++){
        arg2[i] = matrix[i];
        
        //printf("ma i . %p, nr: %d \n", matrix[i], i);
    }
/*
    void *argd1;
    argd1 = malloc(sizeof(struct matrixElement*) * nrElem*nrElem);
    struct matrixElement * arg3 = (struct matrixElement *) argd1;

 */
   



    printf("coreDumped here?");
    fflush(stdout);

    initMatrix(arg, particles, n , nrElem, max_velocity); //Not core dumped in initmatrix


     //Check if pointers point right, it's initialized right
    int nrRight = 0;
    for(int i = 0; i < n;i++){
        particle_t *ps = particles + i;
            int y = (int) ((*(particles + i)).y/max_velocity);
            int xz = (int) ((*(particles + i)).x/max_velocity);

            struct matrixElement *el = (*(arg + y) + xz);

        int rightPlace = 0;
        for(int x = 0; x < el -> size; x++){
            if(el -> list[x] == ps){
                rightPlace = 1;
            }
        }
        if(rightPlace){
            nrRight++;
        }
    }
    printf("nrRight: %d\n", nrRight);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ ) //Tror den core dumpar eftersom antalet partiklar på en plats byggs upp tills det är över gränsen, kolla.
    {
        printf("New step \n");
        //
        //  compute forces
        //
        for(int y = 0; y < nrElem; y++){
            for(int x = 0; x < nrElem; x++){//Goes through every matrix element
                for(int z = 0; z < matrix[x][y].size; z++){ //Goes through each particle in each matrix
                    particle_t *part = matrix[x][y].list[z];
                    for(int g = -1; g <= 1; g++){ //Applies forces from every particle on adjacet and the same matrix element
                        for(int h = -1; h <= 1; h++){
                            if(x+g >= 0 && x+g < nrElem && y+h >= 0 && y+h < nrElem){
                                for(int gh = 0; gh < matrix[x+g][y+h].size; gh++){ 
                                if(g != 0 || h != 0 || gh != z){// an if-statement so it doesn't apply force to itself
                                    apply_force(*part, *(matrix[x+g][y+h].list[gh]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        /*
        for( int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( particles[i], particles[j] );
        } */
        
        //
        //  move particles
        //
        for( int i = 0; i < n; i++ ) 
            move( particles[i] );
        

    reposition(arg, particles, n, nrElem, max_velocity);
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
