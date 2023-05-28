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
            centers[cc][dd] ={ xy_cen[count][0],xy_cen[count][1]};
            // centers[cc][dd].assign(xy_cen[count].begin(),xy_cen[count].end()) ;
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



vector<vector<vector<double> > >  CalculateQk(vector<vector<vector<vector<double> > > > &boxes,vector<vector<vector<double> > > &centers, vector<int> p_vec,int L){
    vector<double> dist;
    vector<vector<vector<vector<double> > > >  dist_pow(L, vector<vector<vector<double> > >(L));
    vector<vector<vector<double> > > Q_k(L, vector<vector<double> > (L));
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
            // Check for Qk
            cout << "Q_k {" << gg << "," << hh << "}: ";
            for (int i = 0; i < Q_k[gg][hh].size(); i++) {
                cout << Q_k[gg][hh][i] << " ";
            }
            cout << endl;
            
        }
    }
    return Q_k;
}


void CalaulateFarorNear(vector<vector<vector<vector<int> > > > & far_list,vector<vector<vector<vector<int> > > > & near_list,int L){

    vector<vector<int> >  pr(L, vector<int>(L));
    vector<vector<int> >  qr(L, vector<int>(L));
    vector<vector<vector<int> > > pairs(L, vector<vector<int> > (L));
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            pr[i][j] = j;
            pairs[i][j].push_back(pr[i][j]);
            qr[i][j] = i;
            pairs[i][j].push_back(qr[i][j]);
        }
    }
    vector<vector<vector<vector<vector<int> > > > > ind(L, vector<vector<vector<vector<int> > > > (L, vector<vector<vector<int> > >(pairs.size(), vector<vector<int> > (pairs.size()))));
    vector<vector<vector<vector<int> > > >  flag(L, vector<vector<vector<int> > >(L, vector<vector<int> > (pairs.size(), vector<int>(pairs.size()))));
    
    for (int ii = 0; ii < L; ii++) {
        for (int jj = 0; jj < L; jj++) {
            for (int k = 0; k < pairs.size(); k++) {
                for (int l = 0; l < pairs.size(); l++) {
                    vector<int> list = {k,l};
                    ind[ii][jj][k][l].push_back(abs(ii - pairs[k][l][0]));
                    ind[ii][jj][k][l].push_back(abs(jj - pairs[k][l][1]));
                    
                    if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] >= 2) {
                        flag[ii][jj][k][l] = 2;
                        far_list[ii][jj].push_back(list);
                    }
                    else if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] < 2) {
                        flag[ii][jj][k][l] = 1;
                        far_list[ii][jj].push_back(list);
                    }
                    else if (ind[ii][jj][k][l][0] < 2 && ind[ii][jj][k][l][1] >= 2) {
                        flag[ii][jj][k][l] = 1;
                        far_list[ii][jj].push_back(list);
                    }
                    else {
                        flag[ii][jj][k][l] = 0;
                        near_list[ii][jj].push_back(list);
                    }
                }
            }
            
            
        }
    }



}


void CalculatePhiFar(vector<vector<vector<double> > > &phi_far,vector<vector<vector<vector<int> > > > & far_list,vector<vector<vector<vector<double> > > > &boxes,vector<vector<vector<double> > > &centers,vector<vector<vector<double> > > &Q_k, vector<int>& p_vec,int L){
    vector<double> dist_pow2(p_vec.size());
    for (int kk = 0; kk < L; kk++) {
        for (int ll = 0; ll < L; ll++) {
            int n_par = static_cast<int>(boxes[kk][ll].size());
            int temp = static_cast<int>(far_list[kk][ll].size());
            phi_far[kk][ll] = vector<double>(n_par, 0.0);
            
            // cout << "box {" << kk << "," << ll << "}: " << endl;
            
            for (int mm = 0; mm < temp; mm++) {
                
                cout << "far_list" << mm << endl;
                
                for (int mmm = 0; mmm < n_par; mmm++) {
                    double dist2 = sqrt(pow(boxes[kk][ll][mmm][0] - centers [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [0], 2) + pow(boxes[kk][ll][mmm][1] - centers [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [1], 2));
                    
                    //cout << "particle" << mmm << ":";
                    //cout << dist2 << endl;
                    
                    for (int i = 0; i < p_vec.size(); i++) {
                        dist_pow2[i] = pow(dist2, p_vec[i]+1);
                        
                        // cout << dist_pow2[i] << " ";
                        // if (kk==0 && ll==2){
                        // cout << " phi_far="<<phi_far[kk][ll][mmm]  << endl;
                        // }
                        phi_far[kk][ll][mmm] += Q_k [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [i] / dist_pow2[i];
                        if (kk==0 && ll==2){
                        cout << "i=" << i << " mm=" << mm << " mmm=" << mmm << " phi_far="<<phi_far[kk][ll][mmm] << " Qk="<< Q_k [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [i]  <<" dist=" << dist_pow2[i] << endl;
                        }
                    }
                }
            }
            
            
            cout << "phi_far {" << kk << "," << ll << "}:" << endl;
            for (int i = 0; i < phi_far[kk][ll].size(); i++) {
                cout << phi_far[kk][ll][i] << " ";
            }
            cout << endl;
            
            
        }
    }
}

void CalculatePhiNear(vector<vector<vector<double> > > &phi_near,vector<vector<vector<vector<int> > > > & near_list,vector<vector<vector<vector<double> > > > &boxes,int L){
    for (int nn = 0; nn < L; nn++) {
        for (int oo = 0; oo < L; oo++) {
            int n_part = static_cast<int>(boxes[nn][oo].size());
            int temp3 = static_cast<int>(near_list[nn][oo].size());
            phi_near[nn][oo] = vector<double>(n_part, 0.0);
            
            for (int pp = 0; pp < temp3; pp++) {
                int n = static_cast<int>(boxes[near_list[nn][oo][pp][0]][near_list[nn][oo][pp][1]].size());
                vector<double> pair_dist(n, 0.0);
                vector<double> pair_used(n, 0.0);
                for (int qq = 0; qq < n_part; qq++) {
                    for (int i = 0; i < n; i++) {
                        pair_dist[i] = sqrt( pow((boxes[nn][oo][qq][0] - boxes [near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [0]),2) + pow((boxes[nn][oo][qq][1] - boxes [near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [1]),2) );
                        pair_used[i] = boxes[near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [2] / pair_dist[i];
                        if (isinf(pair_used[i])) {
                            pair_used[i] = 0;
                        }
                        phi_near[nn][oo][qq] += pair_used[i];
                    }
                }
            }
            
            // cout << "phi_near {" << nn << "," << oo << "}:" << endl;
            // for (int i = 0; i < phi_near[nn][oo].size(); i++) {
            //     cout << phi_near[nn][oo][i] << " ";
            // }
            // cout << endl;
            
        }
    }
}

void Final_Phi_FMM(vector<vector<vector<double> > > &phi_FMM,vector<vector<vector<double> > > &phi_far,vector<vector<vector<double> > > &phi_near,vector<vector<vector<vector<double> > > > &boxes,int L){
    for (int rr = 0; rr < L; rr++) {
        for (int ss = 0; ss < L; ss++) {
            int num = static_cast<int>(boxes[rr][ss].size());
            phi_FMM[rr][ss] = vector<double>(num, 0.0);
            
            for (int i = 0; i < num; i++) {
                phi_FMM[rr][ss][i] = phi_far[rr][ss][i] + phi_near[rr][ss][i];
            }
            /*
            cout << "phi_FMM {" << rr << "," << ss << "}:" << endl;
            for (int i = 0; i < phi_FMM[rr][ss].size(); i++) {
                cout << phi_FMM[rr][ss][i] << " ";
            }
            cout << endl;
            */
        }
    }
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
    vector<double> x(N);                              // Sources x-position vector.
    vector<double> y(N);                              // Sources y-position vector.
    double step_size = (x_max - x_min) / round(L);


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
    
    
    // double step_size = (x_max - x_min) / round(L);

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
    
    /*
    for (int i = 0; i <= p; ++i) {
        cout << p_vec[i] << " ";
    }
    cout << endl;
    */
    
    cout << endl;
    cout << " - Constructing the Laurent series at each box center..." << endl;
    
    vector<vector<vector<double> > > Q_k(L, vector<vector<double> > (L));
    Q_k = CalculateQk(boxes,centers,p_vec,L);
    
    cout << " - Evaluating the potential..." << endl;
    
    // vector<vector<int> >  pr(L, vector<int>(L));
    // vector<vector<int> >  qr(L, vector<int>(L));
    // for (int i = 0; i < L; i++) {
    //     for (int j = 0; j < L; j++) {
    //         pr[i][j] = j;
    //         qr[i][j] = i;
    //     }
    // }
    
    // vector<vector<vector<int> > > pairs(L, vector<vector<int> > (L));
    // for (int i = 0; i < L; i++) {
    //     for (int j = 0; j < L; j++) {
    //         pairs[i][j].push_back(qr[i][j]);
    //         pairs[i][j].push_back(pr[i][j]);
    //     }
    // }
    
    /*
    for (int i = 0; i < pairs.size(); i++) {
        for (int j = 0; j < pairs[i].size(); j++) {
            cout << "(" << i << "," << j << "): ";
            for (int k = 0; k < pairs[i][j].size(); k++) {
                cout << pairs[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
    */
    
    vector<vector<vector<vector<int> > > >  far_list(L, vector<vector<vector<int> > >(L));
    vector<vector<vector<vector<int> > > >  near_list(L, vector<vector<vector<int> > >(L));
    CalaulateFarorNear(far_list,near_list,L);

    // vector<vector<vector<vector<vector<int> > > > > ind(L, vector<vector<vector<vector<int> > > > (L, vector<vector<vector<int> > >(pairs.size(), vector<vector<int> > (pairs.size()))));
    // vector<vector<vector<vector<int> > > >  flag(L, vector<vector<vector<int> > >(L, vector<vector<int> > (pairs.size(), vector<int>(pairs.size()))));
    // vector<vector<vector<vector<int> > > >  far_list(L, vector<vector<vector<int> > >(L));
    // vector<vector<vector<vector<int> > > >  near_list(L, vector<vector<vector<int> > >(L));

    // for (int ii = 0; ii < L; ii++) {
    //     for (int jj = 0; jj < L; jj++) {
    //         for (int k = 0; k < pairs.size(); k++) {
    //             for (int l = 0; l < pairs.size(); l++) {
    //                 ind[ii][jj][k][l].push_back(abs(ii - pairs[k][l][0]));
    //                 ind[ii][jj][k][l].push_back(abs(jj - pairs[k][l][1]));
                    
    //                 if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] >= 2) {
    //                     flag[ii][jj][k][l] = 2;
    //                 }
    //                 else if (ind[ii][jj][k][l][0] >= 2 && ind[ii][jj][k][l][1] < 2) {
    //                     flag[ii][jj][k][l] = 1;
    //                 }
    //                 else if (ind[ii][jj][k][l][0] < 2 && ind[ii][jj][k][l][1] >= 2) {
    //                     flag[ii][jj][k][l] = 1;
    //                 }
    //                 else {
    //                     flag[ii][jj][k][l] = 0;
    //                 }
    //             }
    //         }
    //         /*
    //         cout << "flag {" << ii << "," << jj << "}: ";
    //         for (int i = 0; i < pairs.size(); i++) {
    //             for (int j = 0; j < pairs.size(); j++) {
    //                 cout << flag[ii][jj][i][j] << " ";
    //             }
    //             cout << endl;
    //         }
    //         cout << endl;
    //         */
    //         for (int k = 0; k < pairs.size(); k++) {
    //             for (int l = 0; l < pairs.size(); l++) {
    //                 if (flag[ii][jj][k][l] != 0) {
    //                     vector<int> list = {k,l};
    //                     far_list[ii][jj].push_back(list);
    //                 }
    //             }
    //         }
    //         /*
    //         cout << "far_list{" << ii << "," << jj << "}: ";
    //         for (int i = 0; i < far_list[ii][jj].size(); i++) {
    //             vector<int> sublist = far_list[ii][jj][i];
    //             for (int j = 0; j < sublist.size(); j++) {
    //                 cout << sublist[j] << " ";
    //             }
    //             cout << endl;
    //         }
    //         */
    //         for (int k = 0; k < pairs.size(); k++) {
    //             for (int l = 0; l < pairs.size(); l++) {
    //                 if (flag[ii][jj][k][l] == 0) {
    //                     vector<int> neigh = {k,l};
    //                     near_list[ii][jj].push_back(neigh);
    //                 }
    //             }
    //         }
    //         /*
    //         cout << "near_list{" << ii << "," << jj << "}: ";
    //         for (int i = 0; i < near_list[ii][jj].size(); i++) {
    //             vector<int> sublist = near_list[ii][jj][i];
    //             for (int j = 0; j < sublist.size(); j++) {
    //                 cout << sublist[j] << " ";
    //             }
    //             cout << endl;
    //         }
    //         */
    //     }
    // }
  
    vector<vector<vector<double> > > phi_far(L, vector<vector<double> > (L));
    CalculatePhiFar(phi_far,far_list,boxes,centers,Q_k,p_vec,L);

    // vector<double> dist_pow2(p_vec.size());
    
    // for (int kk = 0; kk < L; kk++) {
    //     for (int ll = 0; ll < L; ll++) {
    //         int n_par = static_cast<int>(boxes[kk][ll].size());
    //         int temp = static_cast<int>(far_list[kk][ll].size());
    //         phi_far[kk][ll] = vector<double>(n_par, 0.0);
            
    //         cout << "box {" << kk << "," << ll << "}: " << endl;
            
    //         for (int mm = 0; mm < temp; mm++) {
                
    //             //cout << "far_list" << mm << endl;
                
    //             for (int mmm = 0; mmm < n_par; mmm++) {
    //                 double dist2 = sqrt(pow(boxes[kk][ll][mmm][0] - centers [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [0], 2) + pow(boxes[kk][ll][mmm][1] - centers [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [1], 2));
                    
    //                 //cout << "particle" << mmm << ":";
    //                 //cout << dist2 << endl;
                    
    //                 for (int i = 0; i < p_vec.size(); i++) {
    //                     dist_pow2[i] = pow(dist2, p_vec[i]+1);
                        
    //                     //cout << dist_pow2[i] << " ";
                        
    //                     phi_far[kk][ll][mmm] += Q_k [far_list[kk][ll][mm][0]] [far_list[kk][ll][mm][1]] [i] / dist_pow2[i];
                        
    //                     //cout << "i=" << i << " mm=" << mm << " mmm=" << mmm << " phi_far="<<phi_far[kk][ll][mmm] << " dist=" << dist_pow2[i] << endl;
    //                 }
    //             }
    //         }
            
    //         cout << "phi_far {" << kk << "," << ll << "}:" << endl;
    //         for (int i = 0; i < phi_far[kk][ll].size(); i++) {
    //             cout << phi_far[kk][ll][i] << " ";
    //         }
    //         cout << endl;
            
    //     }
    // }
    
    

    vector<vector<vector<double> > > phi_near(L, vector<vector<double> > (L));
    CalculatePhiNear(phi_near,near_list,boxes,L);

    
    // for (int nn = 0; nn < L; nn++) {
    //     for (int oo = 0; oo < L; oo++) {
    //         int n_part = static_cast<int>(boxes[nn][oo].size());
    //         int temp3 = static_cast<int>(near_list[nn][oo].size());
    //         phi_near[nn][oo] = vector<double>(n_part, 0.0);
            
    //         for (int pp = 0; pp < temp3; pp++) {
    //             int n = static_cast<int>(boxes[near_list[nn][oo][pp][0]][near_list[nn][oo][pp][1]].size());
    //             vector<double> pair_dist(n, 0.0);
    //             vector<double> pair_used(n, 0.0);
    //             for (int qq = 0; qq < n_part; qq++) {
    //                 for (int i = 0; i < n; i++) {
    //                     pair_dist[i] = sqrt( pow((boxes[nn][oo][qq][0] - boxes [near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [0]),2) + pow((boxes[nn][oo][qq][1] - boxes [near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [1]),2) );
    //                     pair_used[i] = boxes[near_list[nn][oo][pp][0]] [near_list[nn][oo][pp][1]] [i] [2] / pair_dist[i];
    //                     if (isinf(pair_used[i])) {
    //                         pair_used[i] = 0;
    //                     }
    //                     phi_near[nn][oo][qq] += pair_used[i];
    //                 }
    //             }
    //         }
    //         /*
    //         cout << "phi_near {" << nn << "," << oo << "}:" << endl;
    //         for (int i = 0; i < phi_near[nn][oo].size(); i++) {
    //             cout << phi_near[nn][oo][i] << " ";
    //         }
    //         cout << endl;
    //         */
    //     }
    // }
    
    vector<vector<vector<double> > > phi_FMM(L, vector<vector<double> > (L));
    Final_Phi_FMM(phi_FMM,phi_far,phi_near,boxes,L);
    // for (int rr = 0; rr < L; rr++) {
    //     for (int ss = 0; ss < L; ss++) {
    //         int num = static_cast<int>(boxes[rr][ss].size());
    //         phi_FMM[rr][ss] = vector<double>(num, 0.0);
            
    //         for (int i = 0; i < num; i++) {
    //             phi_FMM[rr][ss][i] = phi_far[rr][ss][i] + phi_near[rr][ss][i];
    //         }
    //         /*
    //         cout << "phi_FMM {" << rr << "," << ss << "}:" << endl;
    //         for (int i = 0; i < phi_FMM[rr][ss].size(); i++) {
    //             cout << phi_FMM[rr][ss][i] << " ";
    //         }
    //         cout << endl;
    //         */
    //     }
    // }


    return 0;
}

