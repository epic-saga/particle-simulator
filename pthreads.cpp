#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include <vector>

//
//  global variables
//
int n, n_threads;
particle_t *particles;
FILE *fsave;
pthread_barrier_t barrier;

double max_velocity;
double size2;
int nrElem;
std::vector<std::vector<std::vector<particle_t *>>> mat;

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (n + n_threads - 1) / n_threads;
    int first = min(  thread_id    * particles_per_thread, n );
    int last  = min( (thread_id+1) * particles_per_thread, n );
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        for(int y = 0; y < nrElem; y++){        
            for(int x = 0; x < nrElem; x++){//Goes through every matrix element
            std::vector<particle_t *> element = mat.at(y).at(x);
                for(int z = 0; z < element.size(); z++){ //Goes through each particle in each matrix
                    particle_t *part = element.at(z);
                    for(int g = -1; g <= 1; g++){ //Applies forces from every particle on adjacet and the same matrix element
                        for(int h = -1; h <= 1; h++){
                            if(x+g >= 0 && x+g < nrElem && y+h >= 0 && y+h < nrElem){
                                std::vector<particle_t *> checkAgainst = mat.at(y+h).at(x+g);
                                for(int gh = 0; gh < checkAgainst.size(); gh++){ 
                                if(g != 0 || h != 0 || gh != z){// an if-statement so it doesn't apply force to itself
                                    apply_force(*part, *checkAgainst.at(gh));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        pthread_barrier_wait( &barrier );
        
        //
        //  move particles
        //
        for( int i = first; i < last; i++ ) {
            move( particles[i] );
            
        }

        reposition(mat, particles, n, nrElem, max_velocity); 
        
        pthread_barrier_wait( &barrier );
        
        //
        //  save if necessary
        //
        if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
    }
    
    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    
    //
    //  allocate resources
    //
    fsave = savename ? fopen( savename, "w" ) : NULL;

    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );


    max_velocity = getCutoff();
    size2 = getSize();
    nrElem = (int) (size2/max_velocity) + 1;

    //std::vector<std::vector<std::vector<particle_t *>>> mat;
    //mat.reserve(nrElem);

    for(int i = 0; i < nrElem; i++){
        std::vector<std::vector<particle_t *>> another = {};
        for(int x = 0; x < nrElem; x++){
            std::vector<particle_t *> an;
            another.push_back(an);
        }
        mat.push_back(another);
        
    }
    fflush(stdout);

    initMatrix(mat, particles, n , nrElem, max_velocity); //Not core dumped in initmatrix


     //Check if pointers point right, it's initialized right
    int nrRight = 0;
    for(int i = 0; i < n;i++){
        particle_t *ps = particles + i;
            int y = (int) ((*(particles + i)).y/max_velocity);
            int xz = (int) ((*(particles + i)).x/max_velocity);

            std::vector<particle_t *> el = mat.at(i/nrElem).at(i%nrElem);

        int rightPlace = 0;
        for(int x = 0; x < el.size(); x++){
            if(el.at(x) == ps){
                rightPlace = 1;
            }
        }
        if(rightPlace){
            nrRight++;
        }
    }

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );
    
    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, n_threads = %d, simulation time = %g seconds\n", n, n_threads, simulation_time );
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
