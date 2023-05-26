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

bool compareRows(const vector<double>& row1, const vector<double>& row2) {
    if (row1[3] < row2[3])
        return true;
    else if (row1[3] > row2[3])
        return false;
    else
        return row1[4] < row2[4];
}

int main() {
    int N;
    cout << "Number of sources, N, to be allocated in the unit box [def = 10,000]: ";
    cin >> N;

    // Number of sources, N.
    if (N < 8) {
        cout << endl;
        cout << "Try a bigger N to generate at least four boxes!" << endl;
        cout << endl;
        return 0;
    }

    int N_b = round(sqrt(N));
    // Number of boxes, M, based on N = O(M^2).
    int L = round(sqrt(N_b));
    // Rounding the box size to fit the unit box.
    int M_disp = L*L;
    // Number of generated boxes after rounding the box size.

    vector<vector<double>> m(N, vector<double>(3));
    // Sources matrix, m = [x, y, q].
    vector<double> x(N);
    // Sources x-position vector.
    vector<double> y(N);
    // Sources y-position vector.

    double x_min = 0;
    // Minimum x grid point.
    double y_min = x_min;
    // Minimum y grid point.
    double x_max = 1;
    // Maximum x grid point.
    double y_max = x_max;
    // Maximum y grid point.

    // Generate random values for m
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            m[i][j] = static_cast<double>(rand()) / RAND_MAX;
        }
    }
    
    /*
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            cout << m[i][j] << " ";
        }
        cout << endl;
    }
    */
    
    // Extract x and y from m
    for (int i = 0; i < N; i++) {
        x[i] = m[i][0];
        y[i] = m[i][1];
    }
    
    /*
    cout << "x: ";
    for (int i = 0; i < N; i++) {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "y: ";
    for (int i = 0; i < N; i++) {
        cout << y[i] << " ";
    }
    cout << endl;
    */
    
    // Mesh generation
    cout << endl;
    cout << "Generating the mesh..." << endl;
    cout << "# of generated boxes M, based on N = O(M^2): " << M_disp << endl;
    cout << endl;
    
    double step_size = (x_max - x_min) / round(L);

    vector<double> x_edg(round(L) + 1);
    vector<double> y_edg(round(L) + 1);
    vector<double> x_cen1(round(L));
    vector<double> y_cen1(round(L));

    for (int i = 0; i <= round(L); i++) {
        x_edg[i] = x_min + i * step_size;
        y_edg[i] = y_min + i * step_size;
    }
    
    /*
    cout << "x_edg: ";
    for (int i = 0; i <= L; i++) {
        cout << x_edg[i] << " ";
    }
    cout << endl;
     
    cout << "y_edg: ";
    for (int i = 0; i <= L; i++) {
        cout << y_edg[i] << " ";
    }
    cout << endl;
    */

    x_cen1.assign(round(L), 0);
    y_cen1.assign(round(L), 0);
    
    // Indexing box centers
    for (int aa = 0; aa < x_edg.size() - 1; aa++) {
        x_cen1[aa] = 0.5 * (x_edg[aa] + x_edg[aa + 1]);
    }
    for (int bb = 0; bb < y_edg.size() - 1; bb++) {
        y_cen1[bb] = 0.5 * (y_edg[bb] + y_edg[bb + 1]);
    }
    
    /*
    cout << "x_cen1: ";
    for (double val : x_cen1) {
        cout << val << " ";
    }
    cout << endl;

    cout << "y_cen1: ";
    for (double val : y_cen1) {
        cout << val << " ";
    }
    cout << endl;
    */
    
    vector<vector<double>> x_cen2(y_cen1.size(), vector<double>(x_cen1.size()));
    vector<vector<double>> y_cen2(y_cen1.size(), vector<double>(x_cen1.size()));

    // Indexing box centers using meshgrid
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            x_cen2[i][j] = x_cen1[j];
            y_cen2[i][j] = y_cen1[i];
        }
    }
    
    /*
    cout << "x_cen2:" << endl;
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            cout << x_cen2[i][j] << " ";
        }
        cout << endl;
    }
    
    cout << "y_cen2:" << endl;
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            cout << y_cen2[i][j] << " ";
        }
        cout << endl;
    }
    */
    
    // Combine x_cen2 and y_cen2 to obtain xy_cen
    vector<vector<double>> xy_cen(x_cen1.size() * y_cen1.size(), vector<double>(2));
    for (int i = 0; i < y_cen1.size(); i++) {
        for (int j = 0; j < x_cen1.size(); j++) {
            xy_cen[i * x_cen1.size() + j][0] = x_cen2[i][j];
            xy_cen[i * x_cen1.size() + j][1] = y_cen2[i][j];
        }
    }
    
    /*
    for (int i = 0; i < xy_cen.size(); i++) {
        cout << xy_cen[i][0] << " " << xy_cen[i][1] << endl;
    }
    */
    
    vector<vector<vector<double>>> centers(L, vector<vector<double>>(L));

    int count = 0;

    for (int cc = 0; cc < L; cc++) {
        for (int dd = 0; dd < L; dd++) {
            centers[cc][dd] = {xy_cen[count][0], xy_cen[count][1]};
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
    
    // Indexing boxes and source positions
    vector<int> n_x(x.size()), M_x(x.size());
    vector<int> n_y(y.size()), M_y(y.size());

    for (size_t i = 0; i < x.size(); i++) {
        auto it = upper_bound(x_edg.begin(), x_edg.end(), x[i]);
        n_x[i] = static_cast<int>(distance(x_edg.begin(), it));
        M_x[i] = n_x[i] - 1;
    }
    for (size_t i = 0; i < y.size(); i++) {
        auto it = upper_bound(y_edg.begin(), y_edg.end(), y[i]);
        n_y[i] = static_cast<int>(distance(y_edg.begin(), it));
        M_y[i] = n_y[i] - 1;
    }
    
    vector<vector<double>> M;
    
    for (size_t i = 0; i < m.size(); i++) {
        vector<double> row = m[i];
        row.push_back(static_cast<double>(M_x[i]));
        row.push_back(static_cast<double>(M_y[i]));
        M.push_back(row);
    }
    
    sort(M.begin(), M.end(), compareRows);
    
    /*
    for (int i = 0; i < M.size(); i++) {
        cout << "Data " << i+1 << ": ";
        for (int j = 0; j < M[i].size(); j++) {
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
    */
    
    vector<vector<vector<vector<double>>>> boxes(L, vector<vector<vector<double>>>(L));

    for (int i = 0; i < M.size(); i++) {
        int ee = M[i][3];
        int ff = M[i][4];
        boxes[ee][ff].push_back(M[i]);
    }
    
    /*
    for (int ee = 0; ee < L; ee++) {
        for (int ff = 0; ff < L; ff++) {
            cout << "boxes[" << ee << "][" << ff << "]:" << endl;
            for (int i = 0; i < boxes[ee][ff].size(); i++) {
                for (int j = 0; j < boxes[ee][ff][i].size(); j++) {
                    cout << boxes[ee][ff][i][j] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }
    */
    
    cout << endl;
    cout << "Using FMM to compute all pairwise interactions..." << endl;
    cout << " - Input the desired multipole expansion accuracy [def=1e-12]: ";
    
    double eps;
    cin >> eps;
    
    if (eps <= 0) {
        cout << endl;
        cout << "   This is not an acceptable value!" << endl;
        return 0;
    }
    
    cout << "   Epsilon = " << eps << endl;
    int p = ceil(-log2(eps));
    vector<int> p_vec(p + 1);
    for (int i = 0; i <= p; ++i) {
        p_vec[i] = i;
    }
    
    cout << endl;
    cout << " - Constructing the Laurent series at each box center..." << endl;
    
    vector<double> dist;
    vector<vector<vector<vector<double>>>> dist_pow(L, vector<vector<vector<double>>>(L));
    vector<vector<vector<double>>> Q_k(L, vector<vector<double>>(L));

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

