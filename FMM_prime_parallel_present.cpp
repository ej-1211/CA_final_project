//
//  main.cpp
//  FMM_prime
//
//  Created by 劉柏賢 on 2023/5/28.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <random>
#include <chrono>
#include <omp.h> //openMP


// include the final_project_direct_v2.cpp file


using namespace std;


// define the constants and user specified parameters
const int  N = 1024;                    // number of particles (**must be a multiple of BLOCK_SIZE**)
const float BOX_SIZE = 1.0;             // simulation box size
const float MAX_R = 0.5;                // the max radius of the stars
const float V_MAX = 0.1;                // maximum initial velocity
// const int END_STEP = 50;                   // total number of evolution steps
const float Soften_Length = 2e-1;     // soften length for calculating the gravitational acceleration
const double G = 1;               // gravitational constant
const float dt = 0.005;                  // time step

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

double updateStars(std::vector<Star>& stars, const double dt) {
    double Energy=0;
    for (auto& star : stars) {
        
        // Update velocities based on the gravitational force using first-order Euler integration
        star.vx += star.ax * dt;
        star.vy += star.ay * dt;
        // Update positions using first-order Euler integration
        star.x += star.vx * dt;
        star.y += star.vy * dt;
        Energy += star.energy;
    }
    return Energy;
}

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



// Domain of the x and y grid
const double x_min = 0;                                // Minimum x grid point.
const double y_min = x_min;                            // Minimum y grid point.
const double x_max = 1;                                // Maximum x grid point.
const double y_max = x_max;                            // Maximum y grid point.

// Function to measure execution time
double getExecutionTime(const std::chrono::steady_clock::time_point& start, const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
}
bool compareRows(const vector<double>& row1, const vector<double>& row2) {
    if (row1[3] < row2[3])
        return true;
    else if (row1[3] > row2[3])
        return false;
    else
        return row1[4] < row2[4];
}


void Initialize( int N,vector<vector<double> >& m,vector<double>& x, vector<double>& y) {
    // Seed the random number generator
    srand(static_cast<unsigned int>(time(0)));
    // Generate random values for m x and y
    // Loop through every particle
    for (int i = 0; i < N; i++) {
        // Set the mas and position to be between 0 and 1
        for (int j = 0; j < 3; j++) {
            m[i][j] = static_cast<double>(rand()) / RAND_MAX;
        }
        // Extract and store x and y from m
        x[i] = m[i][0];
        y[i] = m[i][1];
    }
}

vector<vector<vector<double> > > GenerateBoxes(vector<double> &x_edg,vector<double>& y_edg, int N_BOX){
    float L = sqrt(N_BOX);
    double step_size = (x_max - x_min) / round(L);
    vector<double> x_cen1(round(L),0);
    vector<double> y_cen1(round(L),0);
    for (int i = 0; i <= round(L); i++) {
        x_edg[i] = x_min + i * step_size;
        y_edg[i] = y_min + i * step_size;
    }

    // Indexing box centers
    for (int aa = 0; aa < x_edg.size() - 1; aa++) {
        x_cen1[aa] = 0.5 * (x_edg[aa] + x_edg[aa + 1]);
        y_cen1[aa] = 0.5 * (y_edg[aa] + y_edg[aa + 1]);
    }

    vector<vector<double> > x_cen2(y_cen1.size(), vector<double>(x_cen1.size()));
    vector<vector<double> > y_cen2(y_cen1.size(), vector<double>(x_cen1.size()));

    // Generating 2D meshgrid with dimensions y_cen1.size() * x_cen1.size()
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            x_cen2[i][j] = x_cen1[j];
            y_cen2[i][j] = y_cen1[i];
        }
    }
    vector<vector<double> > xy_cen(x_cen1.size() * y_cen1.size(), vector<double>(2));
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            xy_cen[i * x_cen1.size() + j][0] = x_cen2[i][j];
            xy_cen[i * x_cen1.size() + j][1] = y_cen2[i][j];
        }
    }
    vector<vector<vector<double> > > centers(L, vector<vector<double> >(L));
    int count = 0;
    for (int cc = 0; cc < L; cc++) {
        for (int dd = 0; dd < L; dd++) {
            // !! Something wrong with the following line !!
            // centers[cc][dd] ={ xy_cen[count][0],xy_cen[count][1]};
            centers[cc][dd].assign(xy_cen[count].begin(),xy_cen[count].end()) ;
            count++;
        }
    }
    
    /*
    for (int cc = 0; cc < L; cc++) {
        for (int dd = 0; dd < L; dd++) {
            cout << "(" << centers[cc][dd][0] << ", " << centers[cc][dd][1] << ") ";
        }
        cout << endl;
    }
    */
    return centers;
}


vector<vector<vector<vector<double> > > > SortParticlestoBox(int N,vector<double> &x_edg,vector<double>& y_edg,vector<vector<double> >& m,vector<double>& x,vector<double>& y, int N_BOX){
    float L = sqrt(N_BOX);
    vector<int> n_x(N), M_x(N);
    vector<int> n_y(N), M_y(N);

    for (size_t i = 0; i < N; i++) {
        auto it = upper_bound(x_edg.begin(), x_edg.end(), x[i]);
        n_x[i] = static_cast<int>(distance(x_edg.begin(), it));
        M_x[i] = n_x[i] - 1;
    }
    for (size_t i = 0; i < N; i++) {
        auto it = upper_bound(y_edg.begin(), y_edg.end(), y[i]);
        n_y[i] = static_cast<int>(distance(y_edg.begin(), it));
        M_y[i] = n_y[i] - 1;
    }
    vector<vector<double> > M;
    for (size_t i = 0; i < m.size(); i++) {
        vector<double> row = m[i];
        row.push_back(static_cast<double>(M_y[i]));
        row.push_back(static_cast<double>(M_x[i]));
        M.push_back(row);
    }
    
    sort(M.begin(), M.end(), compareRows);
    vector<vector<vector<vector<double> > > > boxes(L, vector<vector<vector<double> > >(L));
    for (int i = 0; i < M.size(); i++) {
        int ee = M[i][3];
        int ff = M[i][4];
        boxes[ee][ff].push_back(M[i]);
    }
    return boxes;
}

int main() {

    // Set up for the constant and variables
    int NThread = 64;
    omp_set_num_threads( NThread );

    // Set up for the constant and variables

    int N; // Number of sources, N.
    cout << "Number of sources, N, to be allocated in the unit box [def = 10,000]: ";
    cin >> N;
    const double G = 1; // Gravitational constant, G.
    int N_BOX = round(sqrt(N)); // Number of boxes, M, based on N = O(M^2)
    int L = round(sqrt(N_BOX)); // Rounding the box size to obtain the number of boxes in each direction
    N_BOX = L*L;                // Update the Number of generated boxes after confirming # boxes in each direction
    int N_step = 1;
    float t = 0;
    float r;

    std::vector<Star> stars(N);
    vector<vector<double> > m(N, vector<double>(3));  // Sources matrix, m = [x, y, mass].
    vector<double> x(N);                             // Sources x-position vector.
    vector<double> y(N);                             // Sources y-position vector.

    vector<double> xx(N);
    vector<double> yy(N);
    vector<double> mass(N);
    vector<double> ax(N);
    vector<double> ay(N);
    vector<double> vx(N);
    vector<double> vy(N);

    Initialize(N,m,x,y);

    // Check the initialization
    /*
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            cout << m[i][j] << " ";
        }
        cout << endl;
    }
    */
    
    // Mesh(boxes) generation
    cout << endl;
    cout << "Generating the mesh..." << endl;
    cout << "# of generated boxes, based on N = O(M^2): " << N_BOX << endl;
    cout << endl;


    vector<double> x_edg(round(L) + 1);
    vector<double> y_edg(round(L) + 1);
    vector<vector<vector<double> > > centers(L, vector<vector<double> >(L));
    centers = GenerateBoxes(x_edg,y_edg,N_BOX);
    



    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    for (int step=0;step<N_step;step++)
    {
    vector<vector<vector<vector<double> > > > boxes(L, vector<vector<vector<double> > >(L));
    
    if (step!=0){
        for (int i=0; i<N;i++){
            x[i] = xx[i];
            y[i] = yy[i];
            m[i][0] = x[i];
            m[i][1] = y[i];
            m[i][2] = mass[i];
        }
    }

    // cout << "step="<< step << " x[0]=" << x[0] << endl;
    // cout << "step="<< step << " m[0][0]=" << m[0][0] << endl;


    boxes = SortParticlestoBox(N,x_edg,y_edg,m,x,y,N_BOX);

    
    
    double eps = 1e-1;

    
    if (eps <= 0) {
        cout << endl;
        cout << "   This is not an acceptable value!" << endl;
        return 0;
    }
    
    // cout << "   Epsilon = " << eps << endl;
    int p = ceil(-log2(eps));
    vector<int> p_vec(p + 1);
    for (int i = 0; i <= p; ++i) {
        p_vec[i] = i;
    }
    
    
    // cout << endl;
    // cout << " - Constructing the Laurent series at each box center..." << endl;

    vector<double> dist;
    vector<vector<vector<vector<double> > > >  dist_pow(L, vector<vector<vector<double> > >(L));
    vector<vector<vector<double> > > Q_k(L, vector<vector<double> > (L));
# pragma omp parallel for private(dist)
    for (int gg = 0; gg < L; gg++) {
        for (int hh = 0; hh < L; hh++) {
            auto numPoints = boxes[gg][hh].size();
            auto numPowers = p_vec.size();
            dist.clear();
            for (int i = 0; i < numPoints; i++) {
                dist.push_back(sqrt(pow(boxes[gg][hh][i][0] - centers[gg][hh][0], 2) + pow(boxes[gg][hh][i][1] - centers[gg][hh][1], 2)));
            }
           
            vector<double> temp(numPowers);
            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPowers; j++){
                    temp[j] = pow(dist[i], p_vec[j])*boxes[gg][hh][i][2];
                }
                dist_pow[gg][hh].push_back(temp);
            }
           
            double Q;
            for (int i = 0; i < numPowers; i++) {
                Q = 0;
                for (int j = 0; j < numPoints; j++){
                    Q += dist_pow[gg][hh][j][i];
                }
                Q_k[gg][hh].push_back(Q);
            }
            
        }
    }
    
    // cout << " - Evaluating the potential..." << endl;
    
    vector<vector<int> >  pr(L, vector<int>(L));
    vector<vector<int> >  qr(L, vector<int>(L));
# pragma omp parallel for //collapse(2) ((will slow down
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            pr[i][j] = j;
            qr[i][j] = i;
        }
    }
    
    vector<vector<vector<int> > > pairs(L, vector<vector<int> > (L));
# pragma omp parallel for //collapse(2) ((will slow down
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            pairs[i][j].push_back(qr[i][j]);
            pairs[i][j].push_back(pr[i][j]);
        }
    }
    
    
    
    vector<vector<vector<vector<vector<int> > > > > ind(L, vector<vector<vector<vector<int> > > > (L, vector<vector<vector<int> > >(pairs.size(), vector<vector<int> > (pairs.size()))));
    vector<vector<vector<vector<int> > > >  flag(L, vector<vector<vector<int> > >(L, vector<vector<int> > (pairs.size(), vector<int>(pairs.size()))));
    vector<vector<vector<vector<int> > > >  int_list(L, vector<vector<vector<int> > >(L));
    vector<vector<vector<vector<int> > > >  n_neigh(L, vector<vector<vector<int> > >(L));
# pragma omp parallel for //collapse(2)
    for (int ii = 0; ii < L; ii++) {
        for (int jj = 0; jj < L; jj++) {
            for (int k = 0; k < pairs.size(); k++) {
                for (int l = 0; l < pairs.size(); l++) {
                    ind[ii][jj][k][l].push_back(abs(ii - pairs[k][l][0]));
                    ind[ii][jj][k][l].push_back(abs(jj - pairs[k][l][1]));
                    
                    if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] >= 2) {
                        flag[ii][jj][k][l] = 2;
                    }
                    else if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] < 2) {
                        flag[ii][jj][k][l] = 1;
                    }
                    else if (ind[ii][jj][k][l][0] < 2 && ind[ii][jj][k][l][1] >= 2) {
                        flag[ii][jj][k][l] = 1;
                    }
                    else {
                        flag[ii][jj][k][l] = 0;
                    }
                }
            }
            
            for (int k = 0; k < pairs.size(); k++) {
                for (int l = 0; l < pairs.size(); l++) {
                    if (flag[ii][jj][k][l] != 0) {
                        vector<int> list = {k,l};
                        int_list[ii][jj].push_back(list);
                    }
                }
            }
            
            for (int k = 0; k < pairs.size(); k++) {
                for (int l = 0; l < pairs.size(); l++) {
                    if (flag[ii][jj][k][l] == 0) {
                        vector<int> neigh = {k,l};
                        n_neigh[ii][jj].push_back(neigh);
                    }
                }
            }
            
        }
    }
    vector<vector<vector<double> > > phi_INT(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_INT_x(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_INT_y(L, vector<vector<double> > (L));
    // vector<double> dist_pow2(p_vec.size());
    
# pragma omp parallel for //private(dist_pow2)

    for (int kk = 0; kk < L; kk++) {
        for (int ll = 0; ll < L; ll++) {
            int n_par = static_cast<int>(boxes[kk][ll].size());
            int temp = static_cast<int>(int_list[kk][ll].size());
            phi_INT[kk][ll] = vector<double>(n_par, 0.0);
            // field_INT_x[kk][ll] = vector<double>(n_par, 0.0);
            // field_INT_y[kk][ll] = vector<double>(n_par, 0.0);
            
            //cout << "box {" << kk << "," << ll << "}: " << endl;
            
            for (int mm = 0; mm < temp; mm++) {
                
                //cout << "int_list" << mm << endl;
                
                for (int mmm = 0; mmm < n_par; mmm++) {
                    double dist2 = sqrt(pow(boxes[kk][ll][mmm][0] - centers [int_list[kk][ll][mm][0]] [int_list[kk][ll][mm][1]] [0], 2) + pow(boxes[kk][ll][mmm][1] - centers [int_list[kk][ll][mm][0]] [int_list[kk][ll][mm][1]] [1], 2));
                    // double dist2_x = centers [int_list[kk][ll][mm][0]] [int_list[kk][ll][mm][1]] [0] - boxes[kk][ll][mmm][0];
                    // double dist2_y = centers [int_list[kk][ll][mm][0]] [int_list[kk][ll][mm][1]] [1] - boxes[kk][ll][mmm][1];
                    
                    //cout << "particle" << mmm << ":";
                    //cout << dist2 << endl;
                    
                    for (int i = 0; i < p_vec.size(); i++) {
                        vector<double> dist_pow2(p_vec.size());
                        dist_pow2[i] = pow(dist2, p_vec[i]+1);
                        
                        //cout << dist_pow2[i] << " ";
                        double ans = Q_k [int_list[kk][ll][mm][0]] [int_list[kk][ll][mm][1]] [i] / dist_pow2[i];
                        // double field = ans*(p_vec[i]+1)/dist_pow2[i];
                        
                        phi_INT[kk][ll][mmm] += ans;
                        // field_INT_x[kk][ll][mmm] += field*dist2_x/dist2;
                        // field_INT_y[kk][ll][mmm] += field*dist2_y/dist2;
                        
                        //cout << "i=" << i << " mm=" << mm << " mmm=" << mmm << " phi_INT="<<phi_INT[kk][ll][mmm] << " dist=" << dist_pow2[i] << endl;
                    }
                }
            }
            /*
            cout << "phi_INT {" << kk << "," << ll << "}:" << endl;
            for (int i = 0; i < phi_INT[kk][ll].size(); i++) {
                cout << phi_INT[kk][ll][i] << " ";
            }
            cout << endl;
            cout << "field_INT_x {" << kk << "," << ll << "}:" << endl;
            for (int i = 0; i < field_INT_x[kk][ll].size(); i++) {
                cout << field_INT_x[kk][ll][i] << " ";
            }
            cout << endl;
            cout << "field_INT_y {" << kk << "," << ll << "}:" << endl;
            for (int i = 0; i < field_INT_y[kk][ll].size(); i++) {
                cout << field_INT_y[kk][ll][i] << " ";
            }
            cout << endl;
            */
        }
    }
    
    vector<vector<vector<double> > > phi_NN(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_NN_x(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_NN_y(L, vector<vector<double> > (L));
# pragma omp parallel for 
    for (int nn = 0; nn < L; nn++) {
        for (int oo = 0; oo < L; oo++) {
            int n_part = static_cast<int>(boxes[nn][oo].size());
            int temp3 = static_cast<int>(n_neigh[nn][oo].size());
            phi_NN[nn][oo] = vector<double>(n_part, 0.0);
            // field_NN_x[nn][oo] = vector<double>(n_part, 0.0);
            // field_NN_y[nn][oo] = vector<double>(n_part, 0.0);
            
            for (int pp = 0; pp < temp3; pp++) {
                int n = static_cast<int>(boxes[n_neigh[nn][oo][pp][0]][n_neigh[nn][oo][pp][1]].size());
                vector<double> pair_dist(n, 0.0);
                // vector<double> pair_dist_x(n, 0.0);
                // vector<double> pair_dist_y(n, 0.0);
                vector<double> pair_used(n, 0.0);
                // vector<double> field(n, 0.0);
                
                for (int qq = 0; qq < n_part; qq++) {
                    for (int i = 0; i < n; i++) {
                        pair_dist[i] = sqrt( pow((boxes[nn][oo][qq][0] - boxes [n_neigh[nn][oo][pp][0]] [n_neigh[nn][oo][pp][1]] [i] [0]),2) + pow((boxes[nn][oo][qq][1] - boxes [n_neigh[nn][oo][pp][0]] [n_neigh[nn][oo][pp][1]] [i] [1]),2) );
                        // pair_dist_x[i] = boxes [n_neigh[nn][oo][pp][0]] [n_neigh[nn][oo][pp][1]] [i] [0] - boxes[nn][oo][qq][0];
                        // pair_dist_y[i] = boxes [n_neigh[nn][oo][pp][0]] [n_neigh[nn][oo][pp][1]] [i] [1] - boxes[nn][oo][qq][1];
                        pair_used[i] = boxes[n_neigh[nn][oo][pp][0]] [n_neigh[nn][oo][pp][1]] [i] [2] / pair_dist[i];
                        
                        // if (0 < pair_dist[i] < 0.01) {
                        //     pair_dist[i] = 0.01;
                        // }
                        // if (0 < pair_dist_x[i] && pair_dist_x[i] < 0.01) {
                        //     pair_dist_x[i] = 0.01;
                        // }
                        // if (0 > pair_dist_x[i] && pair_dist_x[i] > -0.01) {
                        //     pair_dist_x[i] = -0.01;
                        // }
                        // if (0 < pair_dist_y[i] && pair_dist_y[i] < 0.01) {
                        //     pair_dist_y[i] = 0.01;
                        // }
                        // if (0 > pair_dist_y[i] && pair_dist_y[i] > -0.01) {
                        //     pair_dist_y[i] = -0.01;
                        // }
                        
                        if (isinf(pair_used[i])) {
                            pair_used[i] = 0;
                        }
                        
                        // field[i] = pair_used[i]/pair_dist[i];
                        
                        // if (isnan(field[i])) {
                        //     field[i] = 0;
                        // }
                        
                        //cout << "{" << nn << "," << oo << "}:" << " neigh " << pp << " particle " << qq << " num " << i << " field is " << field[i] << endl;
                        
                        // double ans_x = field[i]*pair_dist_x[i]/pair_dist[i];
                        // double ans_y = field[i]*pair_dist_y[i]/pair_dist[i];
                        
                        // if (isnan(ans_x)) {
                        //     ans_x = 0;
                        // }
                        // if (isnan(ans_y)) {
                        //     ans_y = 0;
                        // }
                        
                        //cout << "{" << nn << "," << oo << "}:" << " neigh " << pp << " particle " << qq << " num " << i << " is " << ans_x << " " << ans_y << endl;
                        
                        phi_NN[nn][oo][qq] += pair_used[i];
                        // field_NN_x[nn][oo][qq] += ans_x;
                        // field_NN_y[nn][oo][qq] += ans_y;
                    }
                }
            }
            /*
            cout << "phi_NN {" << nn << "," << oo << "}:" << endl;
            for (int i = 0; i < phi_NN[nn][oo].size(); i++) {
                cout << phi_NN[nn][oo][i] << " ";
            }
            cout << endl;
            cout << "field_NN_x {" << nn << "," << oo << "}:" << endl;
            for (int i = 0; i < field_NN_x[nn][oo].size(); i++) {
                cout << field_NN_x[nn][oo][i] << " ";
            }
            cout << endl;
            cout << "field_NN_y {" << nn << "," << oo << "}:" << endl;
            for (int i = 0; i < field_NN_y[nn][oo].size(); i++) {
                cout << field_NN_y[nn][oo][i] << " ";
            }
            cout << endl;
            */
        }
    }
    
    vector<vector<vector<double> > > phi_FMM(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_FMM_x(L, vector<vector<double> > (L));
    // vector<vector<vector<double> > > field_FMM_y(L, vector<vector<double> > (L));
# pragma omp parallel for
    for (int rr = 0; rr < L; rr++) {
        for (int ss = 0; ss < L; ss++) {
            int num = static_cast<int>(boxes[rr][ss].size());
            phi_FMM[rr][ss] = vector<double>(num, 0.0);
            // field_FMM_x[rr][ss] = vector<double>(num, 0.0);
            // field_FMM_y[rr][ss] = vector<double>(num, 0.0);
            
            for (int i = 0; i < num; i++) {
                phi_FMM[rr][ss][i] = phi_INT[rr][ss][i] + phi_NN[rr][ss][i];
                // field_FMM_x[rr][ss][i] = field_INT_x[rr][ss][i] + field_NN_x[rr][ss][i];
                // field_FMM_y[rr][ss][i] = field_INT_y[rr][ss][i] + field_NN_y[rr][ss][i];
            }
            /*
            cout << "phi_FMM {" << rr << "," << ss << "}:" << endl;
            for (int i = 0; i < phi_FMM[rr][ss].size(); i++) {
                cout << phi_FMM[rr][ss][i] << " ";
            }
            cout << endl;
            cout << "field_FMM_x {" << rr << "," << ss << "}:" << endl;
            for (int i = 0; i < field_FMM_x[rr][ss].size(); i++) {
                cout << field_FMM_x[rr][ss][i] << " ";
            }
            cout << endl;
            cout << "field_FMM_y {" << rr << "," << ss << "}:" << endl;
            for (int i = 0; i < field_FMM_y[rr][ss].size(); i++) {
                cout << field_FMM_y[rr][ss][i] << " ";
            }
            cout << endl;
            */
        }
    }

    
//     int count = 0;
// # pragma omp parallel for //(slow down?)
//     for (int q = 0;q<L;q++){
//         for (int w = 0;w<L;w++){
//             for (int n_par = 0;n_par<boxes[q][w].size();n_par++){
//                 xx[count] = boxes[q][w][n_par][0];
//                 yy[count] = boxes[q][w][n_par][1];
//                 mass[count] = boxes[q][w][n_par][2];
//                 ax[count] = field_FMM_x[q][w][n_par];
//                 ay[count] = field_FMM_y[q][w][n_par];
                
//                 // cout << "count" << count << "x: " << xx[count] << " y: " << yy[count] << " mass: " << mass[count] << " ax: " << ax[count] << " ay: " << ay[count] << endl;
//                 # pragma omp critical 
//                 count++;
//             }
//         }
//     }
    // cout << "step="<< step << " xx[0]=" << xx[0] << endl;

    // int count1 = 0;
    // for (auto& star : stars) {
    //     star.mass = mass[count1];  // Choose an appropriate mass for the stars
    //     star.x = xx[count1];
    //     star.y = yy[count1];
    //     star.ax = ax[count1];
    //     star.ay = ay[count1];
    //     if (count1==0){
    //         // cout << "step="<< step << " ax[0]=" << ax[count1] << endl;
    //     }
    //     count1++;
    // }


    // DumpData( step, stars );
    // double return_energy = updateStars(stars, dt);
    // fprintf( stdout, "Step %4d -> %4d: t = %13.7e -> %13.7e (dt = %13.7e) Energy = %13.7e\n", step, step+1, step*dt, (step+1)*dt, dt ,return_energy);
    // t += dt;
    }



    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsedSeconds = getExecutionTime(start, end);
    std::cout << "Simulation time: " << elapsedSeconds << " seconds" << std::endl;


    return 0;
}
