#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>
#include <SDL2/SDL.h>

const int SCREEN_WIDTH = 1000;
const int SCREEN_HEIGHT = 1000;
const double G = 6.67e-8;
// Function to measure execution time
double getExecutionTime(const std::chrono::steady_clock::time_point& start, const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
}

struct Star {
    double mass;
    double x;
    double y;
    double vx;
    double vy;
    double energy;
};

SDL_Window* window = nullptr;
SDL_Renderer* renderer = nullptr;

void updateStars(std::vector<Star>& stars, double dt, const Star& sun, double end_time, bool applyBoundaryConditions) {
    for (auto& star : stars) {
        // Update positions
        star.x += star.vx * dt;
        star.y += star.vy * dt;

        if (applyBoundaryConditions) {
            // Apply periodic boundary conditions
            star.x = std::fmod(star.x, SCREEN_WIDTH);
            if (star.x < 0.0)
                star.x += SCREEN_WIDTH;

            star.y = std::fmod(star.y, SCREEN_HEIGHT);
            if (star.y < 0.0)
                star.y += SCREEN_HEIGHT;
        }

        // Calculate the gravitational force between stars
        for (const auto& otherStar : stars) {
            if (&star != &otherStar) {
                double dx = otherStar.x - star.x;
                double dy = otherStar.y - star.y;
                double distance = std::sqrt(dx * dx + dy * dy);
                
                // Add a softening factor to prevent extremely large forces at close distances
                double softeningFactor = 1;  // Adjust the value to control the softening effect
                double softeningDistance = distance + softeningFactor;

                double force = G * (star.mass * otherStar.mass) / (softeningDistance * softeningDistance);
                double fx = force * (dx / distance);
                double fy = force * (dy / distance);

                // Update velocities based on the gravitational force
                star.vx += fx / star.mass * dt;
                star.vy += fy / star.mass * dt;

                //Update the potential energy for every star
                star.energy += -force*distance;
            }
        }

        // Calculate the gravitational force from the sun
        double dx_sun = sun.x - star.x;
        double dy_sun = sun.y - star.y;
        double distance_sun = std::sqrt(dx_sun * dx_sun + dy_sun * dy_sun);
        double softeningFactor = 100;
        double softeningDistance_sun = distance_sun + softeningFactor;
        double force_sun = G*(star.mass * sun.mass) / (softeningDistance_sun * softeningDistance_sun);
        double fx_sun = force_sun * (dx_sun / distance_sun);
        double fy_sun = force_sun * (dy_sun / distance_sun);

        // Update velocities based on the gravitational force from the sun
        star.vx += fx_sun / star.mass * dt;
        star.vy += fy_sun / star.mass * dt;

        // add the updated kinetic energy
        star.energy += (0.5*star.mass*star.vx*star.vx+0.5*star.mass*star.vy*star.vy);
        // add the potential energy with the sun
        star.energy -= force_sun * distance_sun;
        

    }

    // Check if end time is reached
    if (SDL_GetTicks() / 1000.0 >= end_time) {
        SDL_Event quitEvent;
        quitEvent.type = SDL_QUIT;
        SDL_PushEvent(&quitEvent);
    }
}

void renderStars(const std::vector<Star>& stars, double cellSize, const Star& sun) {
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    // Render the sun as a filled rectangle
    SDL_SetRenderDrawColor(renderer, 255, 255, 0, 255);
    SDL_Rect sunRect = { static_cast<int>(sun.x - cellSize / 2.0), static_cast<int>(sun.y - cellSize / 2.0), static_cast<int>(cellSize), static_cast<int>(cellSize) };
    SDL_RenderFillRect(renderer, &sunRect);

    // Render stars as filled rectangles
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    for (const auto& star : stars) {
        SDL_Rect starRect = { static_cast<int>(star.x - cellSize / 2.0), static_cast<int>(star.y - cellSize / 2.0), static_cast<int>(cellSize), static_cast<int>(cellSize) };
        SDL_RenderFillRect(renderer, &starRect);
    }

    SDL_RenderPresent(renderer);
}

SDL_Point calculateGravityForce(const Star& a, const Star& b) {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double distance = std::sqrt(dx * dx + dy * dy);
    double force = (a.mass * b.mass) / (distance * distance);
    double fx = force * (dx / distance);
    double fy = force * (dy / distance);

    SDL_Point gravityForce;
    gravityForce.x = static_cast<int>(fx);
    gravityForce.y = static_cast<int>(fy);

    return gravityForce;
}

int main() {
    // Set simulation parameters
    int N;  // Number of stars
    double dt;  // Time step
    double cellSize=1.5;  // Cell size
    double end_time;  // End time in seconds
    int step = 0;
    // Check if SDL should be used
    bool useSDL = true;
    bool applyBoundaryConditions = true;
    float t = 0;

    // Get user input
    std::cout << "Enter the number of stars: ";
    std::cin >> N;

    std::cout << "Enter the time step (dt): ";
    std::cin >> dt;

    std::cout << "Enter the end time: ";
    std::cin >> end_time;

    // std::cout << "Enter the particle size (only for simulation): ";
    // std::cin >> cellSize;

    std::cout << "Apply boundary conditions? (1 for yes, 0 for no): ";
    std::cin >> applyBoundaryConditions;

    std::cout << "Use SDL for visualization? (true or false): ";
    std::cin >> useSDL;


    // Initialize SDL
    if (useSDL) {
        if (SDL_Init(SDL_INIT_VIDEO) < 0) {
            std::cerr << "Failed to initialize SDL: " << SDL_GetError() << std::endl;
            return 1;
        }

        // Create SDL window and renderer
        window = SDL_CreateWindow("2D Gravity Simulation", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
        if (!window) {
            std::cerr << "Failed to create SDL window: " << SDL_GetError() << std::endl;
            SDL_Quit();
            return 1;
        }
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
        if (!renderer) {
            std::cerr << "Failed to create SDL renderer: " << SDL_GetError() << std::endl;
            SDL_DestroyWindow(window);
            SDL_Quit();
            return 1;
        }
    }

    // Create random distribution for star positions and velocities
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> posDist(1*SCREEN_WIDTH/4, 3*SCREEN_WIDTH/4);
    std::uniform_real_distribution<double> velDist(-1, 1);

    // Create stars
    std::vector<Star> stars(N);
    for (auto& star : stars) {
        star.mass = 10000.0;  // Choose an appropriate mass for the stars
        star.x = posDist(gen);
        star.y = posDist(gen);
        star.vx = velDist(gen);
        star.vy = velDist(gen);
        // star.vx = 0;
        // star.vy = 0;
    }

    // Create the sun at the center
    Star sun;
    sun.mass = 100000000000000;  // Choose an appropriate mass for the sun
    sun.x = SCREEN_WIDTH / 2.0;
    sun.y = SCREEN_HEIGHT / 2.0;
    sun.vx = 0.0;
    sun.vy = 0.0;

    // Simulation loop
    bool quit = false;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    double startTime = SDL_GetTicks() / 1000.0;

    //initial energy
    double E0;
    double E;
    double E_error_percent;
    for (auto& star : stars) {
                E0 += star.energy;
                }
    // long startTime_point = std::chrono::steady_clock::now();
    while (!quit) {
        if (useSDL) {
            SDL_Event event;
            while (SDL_PollEvent(&event)) {
                if (event.type == SDL_QUIT) {
                    quit = true;
                }
            }

            updateStars(stars, dt, sun, startTime + end_time ,applyBoundaryConditions);
            
            renderStars(stars, cellSize, sun);
            SDL_Delay(10);

        } else {
            if (t >= end_time) {
                quit = true;
            } else {

                
                updateStars(stars, dt, sun, startTime + end_time,applyBoundaryConditions);
                
                
               
                
            }
        }
        double t_buf = t;
        t += dt;
        // Calculate the energy error
        E = 0;
        for (auto& star : stars) {
            if (step==0){
                E0 += star.energy;
                E += star.energy;
            }
        else{
            E += star.energy;
            }
        }
        E_error_percent = (E-E0)/E0*100;
        std::printf("t:%.3f -> %.3f\t Energy error = %.3f %% \n", t_buf,t,E_error_percent);
        step += 1;



        if (t >= end_time) {
                quit = true;
        }
    }

    // Clean up and exit
    if (useSDL) {
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
    } else {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        double elapsedSeconds = getExecutionTime(start, end);
        std::cout << "Simulation time: " << elapsedSeconds << " seconds" << std::endl;
    }

    return 0;
}
