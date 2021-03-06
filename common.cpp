#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

//My functions

double getCutoff(){
    return cutoff;
}
double getSize(){
    return size;
}

void initMatrix(struct matrixElement **matrixP, particle_t *particles, int n, int nrElem, double max_velocity){
    for(int x = 0; x < n; x++){
        //gets the matrix element to which particle[x] belongs
            int y = (int) ((*(particles + x)).y/max_velocity);
            int xz = (int) ((*(particles + x)).x/max_velocity);

            struct matrixElement *element = (*(matrixP + y) + xz);

            element -> list[element -> size++] = particles + x;
    }
    int nrOver = 0;
    int nrUnd = 0;
    for(int y = 0; y < nrElem; y++){
        for(int x = 0; x < nrElem; x++){
            struct matrixElement *element = (*(matrixP + y) +x);
            if(element -> size){
                nrOver++;
            } else{
                nrUnd++;
            }
        }
    }
    printf("nrOver: %d, nrUnder: %d", nrOver, nrUnd);
    printf("coreDumped in initMatrix");
    fflush(stdout);
}

void cleanList(struct matrixElement *element){//Cleans the list so that there aren??t nul pointer laying around in the middle of the list
    for(int i = 0; i < element -> size; i++){
        if(element -> list[i] == 0){
            element -> list[i] = element -> list[--(element -> size)];
        }
    }
}

void reposition(struct matrixElement **matrixP, particle_t *particles, int n, int nrElem, double max_velocity){
    for(int y = 0; y < nrElem; y++){
        for(int x = 0; x < nrElem; x++){

            //Gets the matrix element to which particle[x] belongs
            struct matrixElement *element = *(matrixP + y) + x;
            int cleanNeeded = 0;
            for(int h = 0; h < element -> size; h++){
                if((int) ((element -> list[h] -> x/max_velocity)) != x || (int) ((element -> list[h] -> y/max_velocity)) != y){
                    cleanNeeded = 1;
                    particle_t *moveParticle = element -> list[h];
                    //Getting the new element
                    int my = (int) (moveParticle -> y/max_velocity);
                    int mx = (int) (moveParticle -> x/max_velocity);

                    struct matrixElement *element2 = *(matrixP + my) + mx;
                    element -> list[h] == 0; //Removes it from previous matrix element
                    
                    element2 -> list[(element2 -> size)++] = moveParticle; //Puts it in the new matrix elem
                    if(cleanNeeded){
                        cleanList(element);
                        cleanNeeded = 0;
                    }
                }
            }
        }
    }
    printf("coreDumped in reposition\n");
    fflush(stdout);
}