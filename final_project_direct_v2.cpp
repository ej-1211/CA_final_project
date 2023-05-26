#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>


// define the constants and user specified parameters
const int  N = 1024;                    // number of particles (**must be a multiple of BLOCK_SIZE**)
const float BOX_SIZE = 1.0;             // simulation box size
const float MAX_R = 0.5;                // the max radius of the stars
const float V_MAX = 0.1;                // maximum initial velocity
const int END_STEP = 500;                   // total number of evolution steps
const float Soften_Length = 2e-1;     // soften length for calculating the gravitational acceleration
const double G = 1;               // gravitational constant
const float dt = 0.0005;                  // time step



// Function to measure execution time
double getExecutionTime(const std::chrono::steady_clock::time_point& start, const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
}

// Define the struct for a star
struct Star {
    double mass;
    double x;
    double y;
    double vx;
    double vy;
    double ax;
    double ay;
    double energy;
};


// Function to calculate the gravitaional force for all stars
void calculateGravityForce(std::vector<Star>& stars){
    for (auto& star : stars) {
        // reset the acceleration and energy of the star
        star.ax = 0;
        star.ay = 0;
        star.energy = 0;
        // Calculate the gravitational force between stars
        for (const auto& otherStar : stars) {
            if (&star != &otherStar) {
                double dx = otherStar.x - star.x;
                double dy = otherStar.y - star.y;
                double distance = std::sqrt(dx * dx + dy * dy + Soften_Length * Soften_Length);
                double force = G*(star.mass * otherStar.mass) / (distance * distance);
                double fx = force * (dx / distance);
                double fy = force * (dy / distance);
                star.ax += fx/star.mass;
                star.ay += fy/star.mass;
                // Calculate the energy of the star
                star.energy += -G * star.mass * otherStar.mass / distance;
                star.energy += 0.5 * star.mass * (star.vx * star.vx + star.vy * star.vy);
            }
        }
    }
}


double updateStars(std::vector<Star>& stars, const double dt) {
    double Energy=0;
    for (auto& star : stars) {
        // Update positions using first-order Euler integration
        star.x += star.vx * dt;
        star.y += star.vy * dt;
        // Update velocities based on the gravitational force using first-order Euler integration
        star.vx += star.ax * dt;
        star.vy += star.ay * dt;
        Energy += star.energy;
    }
    return Energy;
}

// -----------------------------------------------------------
// output data
// -----------------------------------------------------------
void DumpData( const int Step, std::vector<Star>& stars)
{

   char FileName[100];
   sprintf( FileName, "Data_%d%d%d%d", Step/1000, (Step%1000)/100, (Step%100)/10, Step%10 );

   FILE *File = fopen( FileName, "w" );

   fprintf( File, "#%12s  %13s  %13s  %13s  %13s  %13s\n",
            "x", "y", "vx", "vy","ax", "ay");

   for (auto& star : stars){
      fprintf( File, "%13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e \n",
               star.x, star.y, star.vx,star.vy,star.ax,star.ay);
   }
   fclose( File );
} // FUNCTION : DumpData

int main() {
    float t = 0;
    float r;

    // Create random distribution for star positions and velocities
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> posDist(0, BOX_SIZE);
    std::uniform_real_distribution<double> velDist(-V_MAX, V_MAX);

    // Create stars
    std::vector<Star> stars(N);
    for (auto& star : stars) {
        star.mass = 1;  // Choose an appropriate mass for the stars
        r = MAX_R + 1;
        while ( r > MAX_R)
        {
            star.x = posDist(gen)-BOX_SIZE*0.5;
            star.y = posDist(gen)-BOX_SIZE*0.5;
            r = std::sqrt(star.x * star.x + star.y * star.y);
        }
        star.vx = velDist(gen);
        star.vy = velDist(gen);
    }



    // std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    // Evolve the stars
    for (int step=0;step<END_STEP;step++){
        double return_energy = updateStars(stars, dt);
        calculateGravityForce(stars);
        fprintf( stdout, "Step %4d -> %4d: t = %13.7e -> %13.7e (dt = %13.7e) Energy = %13.7e\n", step, step+1, step*dt, (step+1)*dt, dt ,return_energy);
        DumpData( step, stars );
    }
    // // Clean up and exit
    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    // double elapsedSeconds = getExecutionTime(start, end);
    // std::cout << "Simulation time: " << elapsedSeconds << " seconds" << std::endl;

    return 0;
}
