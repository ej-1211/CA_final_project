//
//  main.cpp
//  FMM
//
//  Created by 劉柏賢 on 2023/5/25.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;


// Domain of the x and y grid
const double x_min = 0;                                // Minimum x grid point.
const double y_min = x_min;                            // Minimum y grid point.
const double x_max = 1;                                // Maximum x grid point.
const double y_max = x_max;                            // Maximum y grid point.


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
        row.push_back(static_cast<double>(M_x[i]));
        row.push_back(static_cast<double>(M_y[i]));
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

    int N; // Number of sources, N.
    cout << "Number of sources, N, to be allocated in the unit box [def = 10,000]: ";
    cin >> N;
    const double G = 1; // Gravitational constant, G.
    int N_BOX = round(sqrt(N)); // Number of boxes, M, based on N = O(M^2)
    int L = round(sqrt(N_BOX)); // Rounding the box size to obtain the number of boxes in each direction
    N_BOX = L*L;                // Update the Number of generated boxes after confirming # boxes in each direction

    vector<vector<double> > m(N, vector<double>(3));  // Sources matrix, m = [x, y, mass].
    vector<double> x(N);                             // Sources x-position vector.
    vector<double> y(N);                             // Sources y-position vector.

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
    
    
    double step_size = (x_max - x_min) / round(L);

    vector<double> x_edg(round(L) + 1);
    vector<double> y_edg(round(L) + 1);
    vector<vector<vector<double> > > centers(L, vector<vector<double> >(L));
    centers = GenerateBoxes(x_edg,y_edg,N_BOX);
    vector<vector<vector<vector<double> > > > boxes(L, vector<vector<vector<double> > >(L));
    boxes = SortParticlestoBox(N,x_edg,y_edg,m,x,y,N_BOX);

    // Check for the center coordinates
    /*
    for (int cc = 0; cc < L; cc++) {
        for (int dd = 0; dd < L; dd++) {
            cout << "(" << centers[cc][dd][0] << ", " << centers[cc][dd][1] << ") ";
        }
        cout << endl;
    }
    */

    // Final checking for the particle and the boxes -> x,y,mass,#box_x,#box_y    
    
    for (int ee = 0; ee < L; ee++) {
        for (int ff = 0; ff < L; ff++) {
            cout << "boxes[" << ee << "][" << ff << "](" << centers[ee][ff][0] << ","<< centers[ee][ff][1] << "):" << endl;
            for (int i = 0; i < boxes[ee][ff].size(); i++) {
                for (int j = 0; j < boxes[ee][ff][i].size(); j++) {
                    cout << boxes[ee][ff][i][j] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }
    
    
    cout << endl;
    cout << "Using FMM to compute all pairwise interactions..." << endl;
    cout << " - Input the desired multipole expansion accuracy [def=1e-12]: ";
    
    double eps;
    cin >> eps; 
    cout << "   Epsilon = " << eps << endl;
    int p = ceil(-log2(eps));
    vector<int> p_vec(p + 1);
    for (int i = 0; i <= p; ++i) {
        p_vec[i] = i;
    }
    
    cout << endl;
    cout << " - Constructing the Laurent series at each box center..." << endl;
    
    vector<double> dist;
    vector<vector<vector<vector<double> > > > dist_pow(L, vector<vector<vector<double> > >(L));
    vector<vector<vector<double> > > Q_k(L, vector<vector<double> >(L));

    for (int gg = 0; gg < L; gg++) {
        for (int hh = 0; hh < L; hh++) {
            auto numPoints = boxes[gg][hh].size();
            auto numPowers = p_vec.size();
            dist.clear();
            for (int i = 0; i < numPoints; i++) {
                dist.push_back(sqrt(pow(boxes[gg][hh][i][0] - centers[gg][hh][0], 2) + pow(boxes[gg][hh][i][1] - centers[gg][hh][1], 2)));
            }
            /*
            for (int i = 0; i < dist.size(); i++) {
                cout << dist[i] << " ";
            }
            cout << endl;
            */
            vector<double> temp(numPowers);
            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPowers; j++){
                    temp[j] = pow(dist[i], p_vec[j])*boxes[gg][hh][i][2];
                }
                dist_pow[gg][hh].push_back(temp);
            }
            /*
            cout << "dist_pow[" << gg << "][" << hh << "]: ";
            for (int i = 0; i < numPoints; i++) {
                for (int j = 0; j < numPowers; j++) {
                    cout << dist_pow[gg][hh][i][j] << " ";
                }
                cout << endl;
            }
            cout << endl;
            */
            double Q = 0;
            for (int i = 0; i < numPowers; i++) {
                for (int j = 0; j < numPoints; j++){
                    Q += dist_pow[gg][hh][j][i];
                }
                Q_k[gg][hh].push_back(Q);
            }
            /*
            for (int i = 0; i < Q_k[gg][hh].size(); i++) {
                cout << Q_k[gg][hh][i] << " ";
            }
            cout << endl;
            */
        }
    }
    
    cout << " - Evaluating the potential..." << endl;


    
    return 0;
}

