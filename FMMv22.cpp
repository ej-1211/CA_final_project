#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <chrono>
// #include <omp.h>

using namespace std;


double getExecutionTime(const std::chrono::steady_clock::time_point& start, const std::chrono::steady_clock::time_point& end) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000000000.0;
}

int getparentID(int l, int ibox) {
    if (ibox > pow(4, l)) {
        std::cout << "getparentID: boxID and level not consistent!" << std::endl;
        exit(0);
    }
    
    int power = pow(2, l);
    int ibox_x = ((ibox) % power)+1;
    if (ibox_x == 1){
        ibox_x = power+1;
    }
    int ibox_x_p = floor(ibox_x / 2);
    int ibox_y = ceil(ibox / pow(2, l));
    int ibox_y_p = ceil(ibox_y * 0.5);
    int pID = (ibox_y_p - 1) * pow(2, l - 1) + ibox_x_p;
    // if (l==1 && ibox==1){
    //     cout << "power " << power << endl;
    //     cout << "ibox_x " << ibox_x << endl;
    //     cout << "ibox_y " << ibox_y << endl;
    //     cout << "ibox_x_p " << ibox_x_p << endl;
    //     cout << "ibox_y_p " << ibox_y_p << endl;
    //     cout << "pID " << pID << endl;
    //     cout << "ceil(0.5)= " << ceil(0.5) << endl;
    // }
    return pID;
}

std::vector<int> getchildrenID(int l, int ibox) {
    if (ibox > pow(4, l)) {
        std::cout << "getchildrenID: boxID and level not consistent!" << std::endl;
        exit(0);
    }
    int ibox_x = (ibox - 1) % static_cast<int>(pow(2, l)) + 1;
    int ibox_y = ceil(ibox / pow(2, l));
    int ibox_x1 = (ibox_x - 1) * 2 + 1;
    int ibox_x2 = (ibox_x - 1) * 2 + 2;
    int ibox_y1 = (ibox_y - 1) * 2 + 1;
    int ibox_y2 = (ibox_y - 1) * 2 + 2;
    std::vector<int> cID;
    cID.push_back((ibox_y1 - 1) * pow(2, l + 1) + ibox_x1);
    cID.push_back((ibox_y1 - 1) * pow(2, l + 1) + ibox_x2);
    cID.push_back((ibox_y2 - 1) * pow(2, l + 1) + ibox_x1);
    cID.push_back((ibox_y2 - 1) * pow(2, l + 1) + ibox_x2);
    return cID;
}

std::complex<double> getboxcenter(int l, int ibox) {
    if (ibox > pow(4, l)) {
        std::cout << "getboxcenter: boxID and level not consistent!" << std::endl;
        exit(0);
    }
    int ibox_x = (ibox - 1) % static_cast<int>(pow(2, l)) + 1;
    int ibox_y = ceil(ibox / pow(2, l));
    double xC_real = (ibox_x - 0.5) / pow(2, l);
    double xC_imag = (ibox_y - 0.5) / pow(2, l);
    return std::complex<double>(xC_real, xC_imag);
}

std::vector<int> getneighbors(int l, int ibox) {
    std::vector<int> nlist;
    int ibox_x = (ibox - 1) % static_cast<int>(pow(2, l)) + 1;
    int ibox_y = ceil(ibox / pow(2, l));

    if (ibox_y > 1) {
        if (ibox_x > 1) {
            int nID = (ibox_y - 2) * pow(2, l) + ibox_x - 1;    // bottom-left neighbor
            nlist.push_back(nID);
        }

        int nID = (ibox_y - 2) * pow(2, l) + ibox_x;    // bottom neighbor
        nlist.push_back(nID);

        if (ibox_x < pow(2, l)) {
            int nID = (ibox_y - 2) * pow(2, l) + ibox_x + 1;    // bottom-right neighbor
            nlist.push_back(nID);
        }
    }

    if (ibox_x > 1) {
        int nID = (ibox_y - 1) * pow(2, l) + ibox_x - 1;    // left neighbor
        nlist.push_back(nID);
    }

    if (ibox_x < pow(2, l)) {
        int nID = (ibox_y - 1) * pow(2, l) + ibox_x + 1;    // right neighbor
        nlist.push_back(nID);
    }

    if (ibox_y < pow(2, l)) {
        if (ibox_x > 1) {
            int nID = ibox_y * pow(2, l) + ibox_x - 1;    // top-left neighbor
            nlist.push_back(nID);
        }

        int nID = ibox_y * pow(2, l) + ibox_x;    // top neighbor
        nlist.push_back(nID);

        if (ibox_x < pow(2, l)) {
            int nID = ibox_y * pow(2, l) + ibox_x + 1;    // top-right neighbor
            nlist.push_back(nID);
        }
    }

    return nlist;
}

std::vector<int> getinteractionlist(int l, int ibox) {
    std::vector<int> ilist;
    int pID = getparentID(l, ibox);
    std::vector<int> p_nlist = getneighbors(l - 1, pID);
    for (int i = 0; i < p_nlist.size(); ++i) {
        std::vector<int> cID = getchildrenID(l - 1, p_nlist[i]);
        ilist.insert(ilist.end(), cID.begin(), cID.end());
    }
    std::vector<int> nlist = getneighbors(l, ibox);
    for (int i = 0; i < nlist.size(); ++i) {
        ilist.erase(std::remove(ilist.begin(), ilist.end(), nlist[i]), ilist.end());
    }
    // showlist(l, ibox, ilist);
    return ilist;
}


int main() {
    // omp_set_num_threads(1);
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    // Initialize
    int N = 25600;  // number of source/target points
    std::vector<std::complex<double> > z(N);
    // Generating random complex numbers for z
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // # pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        double realPart = dis(gen);
        double imagPart = dis(gen);
        z[i] = std::complex<double>(realPart, imagPart);
    }

    std::vector<double> m(N, 1.0);  // strength of charges
    // for (int i = 0; i < N; ++i) {
    //     cout << m[i] << endl;
    // }
    int p = 3;  // multipole terms
    int n = 4;  // level of refinement

    bool VAL = true;  // validation flag
    struct Box {
                int level;
                int ID;
                std::complex<double> center;
                int pID;
                std::vector<int> cID;
                std::vector<int> nlist;
                std::vector<int> ilist;
                std::vector<std::complex<double> > a;
                std::vector<std::complex<double> > b;
                std::vector<std::complex<double> > c;
                std::vector<std::complex<double> > d;
                std::vector<int> particlelist;
            };
    std::vector<Box> Boxes;
    std::cout << "Initialization" << std::endl;
    std::vector<int> ibox_global;
    
    int base_ID_global = -1;
    for (int l = 1; l <= n; ++l) {
        base_ID_global += std::pow(4, l - 1);
        for (int ibox = 1; ibox <= std::pow(4, l); ++ibox) {
            int ibox_globall = base_ID_global + ibox;
            ibox_global.push_back(ibox_globall-1);
            // Create and initialize Box struct
            
            Box box;
            box.level = l;
            box.ID = ibox;
            box.center = getboxcenter(l, ibox);
            box.pID = getparentID(l, ibox);
            box.cID = getchildrenID(l, ibox);
            box.nlist = getneighbors(l, ibox);
            box.ilist = getinteractionlist(l, ibox);
            box.a.resize(p + 1, 0.0);
            box.b.resize(p + 1, 0.0);
            box.c.resize(p + 1, 0.0);
            box.d.resize(p + 1, 0.0);
            box.particlelist = std::vector<int>();  // this step may need to be performed at every time step

            // Add Box to the Box vector
            
            Boxes.push_back(box);
        }
    }
    // for (int i=0 ; i<Boxes.size() ; i++)
    // {
    //     for (int j=0 ; j<Boxes[i].ilist.size() ; j++)
    //     {
    //         std::cout << Boxes[i].ilist[j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    
    // for (int t=0 ; t<ibox_global.size() ; t++)
    // {
    //     std::cout << "Particl " << t << " Coor= "<<z[ibox_global[t]] << endl;
    // }

    // // Preparation
    // Construct a matrix of useful binomial coefficients
    int numRows = 2 * p + 1;
    int numCols = p;
    std::vector<std::vector<int> > nCk(numRows, std::vector<int>(numCols, 0));
    for (int i_n = 0; i_n < numRows; i_n++) {
        for (int i_k = 0; i_k < std::min(i_n + 1, p + 1); i_k++) {
            nCk[i_n][i_k] = std::tgamma(i_n + 1) / (std::tgamma(i_k + 1) * std::tgamma(i_n - i_k + 1));
        }
    }


    // for (int i=0;i<nCk.size();i++){
    //     for (int j=0;j<nCk[i].size();j++){
    //         std::cout << nCk[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Step 1: S2M
    std::cout << "1.S2M" << std::endl;
    for (int level_1=0;level_1<4;level_1++){
        double boxsize = 0.25;
        std::complex<double> zC = Boxes[level_1].center;
        std::vector<int> idx;
        idx.clear();

        // # pragma omp parallel for
        for (int i = 0; i < N; ++i) {
        // cout << std::abs(std::imag(z[i] - zC)) << endl;
        if (std::abs(std::real(z[i] - zC)) < boxsize && std::abs(std::imag(z[i] - zC)) < boxsize){
            idx.push_back(i);
        }
        }

        // 1~4
        Boxes[level_1].particlelist = idx;
        // for (int i=0;i<idx.size();i++){
        //     std::cout << idx[i] << " ";
        // }
        // cout << endl;


        // if (level_1==0){
        //     for (int i=0 ; i<Boxes[level_1].particlelist.size();i++){
        //     std::cout << Boxes[level_1].particlelist[i] << " ";
        // } 
        // cout << endl;
        
    }

    cout << "first level done" << endl;
    for (int rest_level=2; rest_level<=n;++rest_level){

        for (int ibox = 1; ibox<=pow(4,rest_level);++ibox){
            // calculate global ibox
            int ibox_globalll = 0;
            for (int t=1;t<rest_level;t++){
                ibox_globalll += pow(4,t);
            }
            
            ibox_globalll += ibox;
            ibox_globalll -= 1;
            // cout << "ibox_globalll "<< ibox_globalll << endl; 
            // cout << "ibox_globalll center "<< Boxes[ibox_globalll].center << endl; 
            //calculate global PID
            int pID_global = Boxes[ibox_globalll].pID;
            // cout << "PID "<< pID_global << endl; 
            for (int t=1;t<rest_level-1;t++){
                pID_global += pow(4,t);
            }
            

            // select the particle candidates
            std::vector<int> particlecandidates;
            particlecandidates.clear();
            particlecandidates = Boxes[pID_global-1].particlelist;
            // cout << "parcan didates";
            // for (int i=0;i<particlecandidates.size();i++){
            //     cout << particlecandidates[i] << " ";
            // }
            // cout << endl;
            // choose the box size
            double box_size = pow(0.5,(rest_level+1));

            std::complex<double> zC = Boxes[ibox_globalll].center;
            std::vector<int> idx;

            idx.clear();
            for (int i = 0; i < particlecandidates.size(); ++i) {
            // cout << std::abs(std::imag(z[i] - zC)) << endl;
            if (std::abs(std::real(z[particlecandidates[i]] - zC)) < box_size && std::abs(std::imag(z[particlecandidates[i]] - zC)) < box_size){
             idx.push_back(particlecandidates[i]);
            //  cout << particlecandidates[i] << " ";
            }
            // cout << endl;
            

            Boxes[ibox_globalll].particlelist = idx;

            // if goes to the finest level
            if (rest_level==n){
                std::complex<double> sum_m = (0.0,0.0);
                // cout << "m[idx[i]]= " ; 
                for (int i=0;i<idx.size();++i){
                    // cout  << m[idx[i]] << " ";
                    sum_m += static_cast<double> (m[idx[i]]);
                }
                // cout << endl;
                // cout << "sum= " << sum_m << endl; 
                Boxes[ibox_globalll].a[0] = sum_m;

                for (int k=1;k<=p;k++){
                    std::complex<double> sum_m2;
                    for (int i=0;i<idx.size();i++){
                        std::complex<double> pow_z_minus_zC = std::pow(z[idx[i]] - zC, k);
                        sum_m2 += pow_z_minus_zC * m[idx[i]] /static_cast<double>(k);
                    }
                    Boxes[ibox_globalll].a[k] -= sum_m2;
                }

            }

        }




        }


    }
    cout << "2 level done" << endl;
    // cout << "A" << endl;
    // for (int i=0 ; i<Boxes.size() ; i++)
    // {
    //     for (int j=0 ; j<Boxes[i].ilist.size() ; j++)
    //     {
    //         std::cout << Boxes[i].ilist[j] << " ";
    //     }
    //     std::cout << std::endl;
    //     // for (int j=0 ; j<Boxes[i].particlelist.size() ; j++)
    //     // {
    //     //     std::cout << Boxes[i].particlelist[j] << " ";
    //     // }
    //     // std::cout << std::endl;
    // }



    // Step 2 : M2M
    std::cout << "2.M2M" << std::endl;
    // cout<< "n= "<< n <<endl;
    for (int level=n-1;level>=1;level--){
        
        for (int ibox = 1; ibox<=pow(4,level);++ibox){
            // calculate global ibox
            int ibox_globalll = 0;
            int sum = 0;
            for (int t=1;t<level;t++){
                sum += pow(4,t);
            }
            

            ibox_globalll += sum;
            // cout << "sum " << sum << endl;
            ibox_globalll += ibox;
            // cout << "ibox_globalll "<< ibox_globalll << endl;
            ibox_globalll -= 1;
            // cout << "ibox_globalll "<< ibox_globalll << endl;
            // MOM
            std::complex<double> zM = Boxes[ibox_globalll].center;
            vector<int>  cID = getchildrenID(level,ibox);
            // cout << "cIDb= " << cID[0] << " " << cID[1] << " " << cID[2] << " " << cID[3] << endl;
            int sum2 = 0;
            for (int t=1;t<=level;t++){
                sum2 += pow(4,t);
            }
            //  cout << "sum2 " << sum2 << endl;
            for (int i=0;i<cID.size();i++){
                cID[i] += sum2;
                // cID[i] -= 1;
            }
            
            

            // cout << "cID= " << cID[0] << " " << cID[1] << " " << cID[2] << " " << cID[3] << endl;
            for (int ichild=1;ichild<=4;++ichild){

                // cout << cID[ichild-1]-1 << " a " << Boxes[cID[ichild-1]-1].a[0] << endl;
                Boxes[ibox_globalll].b[0] += Boxes[cID[ichild-1]-1].a[0]; 
                // cout << ibox_globalll << " b " << Boxes[ibox_globalll].b[0] << endl;
                std::complex<double> zC = Boxes[cID[ichild-1]-1].center;
                
                // cout<< "onnnn"<<endl;
                for (int expan=1;expan<=p;++expan){

                    Boxes[ibox_globalll].b[expan] -= Boxes[cID[ichild-1]-1].a[0]*pow((zC-zM),expan)/static_cast<double>(expan);
                    // cout<< "onnnn"<<endl;
                    for (int k=1;k<=expan;++k){
                        Boxes[ibox_globalll].b[expan] += Boxes[cID[ichild-1]-1].a[k]*pow((zC-zM),expan-k)*static_cast<double>(nCk[expan-1][k-1]);
                        // cout << "expan= " << expan << " k= " << k << " nCk= " << nCk[expan-1][k-1] << endl;
                    }
                    
                }



            }
            // cout<< "innnn"<<endl;
            Boxes[ibox_globalll].a = Boxes[ibox_globalll].b;
        }
    }

    // for (int i=0 ; i<Boxes.size() ; i++)
    // {
    //     for (int j=0 ; j<Boxes[i].b.size() ; j++)
    //     {
    //         std::cout << Boxes[i].b[j] << " ";
    //     }
    //     std::cout << std::endl;
    //     // for (int j=0 ; j<Boxes[i].particlelist.size() ; j++)
    //     // {
    //     //     std::cout << Boxes[i].particlelist[j] << " ";
    //     // }
    //     // std::cout << std::endl;
    // }



    // step 3 :Local expansion
    std::cout << "3. Local expansion" << std::endl;
    for (int l=1;l<n;l++){
        
        
        for (int ibox = 1;ibox<=pow(4,l);ibox++){
            // calculate global ibox
            int ibox_globalll = 0;
            
            for (int t=1;t<l;t++){
                ibox_globalll += pow(4,t);
               
            }
            
            ibox_globalll += ibox;
            ibox_globalll -= 1;
        // cout << "ibox_globalll " << ibox_globalll << endl;
        Boxes[ibox_globalll].c = Boxes[ibox_globalll].d;
        
        std::vector<int> iboxlist = Boxes[ibox_globalll].ilist;
        // cout << "ibox_globalll " << ibox_globalll << " iboxlist.size() " << iboxlist.size() << endl;
        if (!iboxlist.empty()){
            
            std::complex<double> zL = Boxes[ibox_globalll].center;
            for (int i=0;i<iboxlist.size();i++){

           
                int ibox_global_m = 0;
                for (int t=1;t<l;t++){
             
                    ibox_global_m += pow(4,t);
                }


                ibox_global_m = ibox_global_m + iboxlist[i]-1;
                // cout << "innn" << endl;
                std::complex<double> zM = Boxes[ibox_global_m].center;
                // cout << "zL " << zL << " zM " << zM << endl;
                // for (int i=0;i<iboxlist.size();i++){
                //     cout <<  iboxlist[i] << " ";
                // }
                // cout << endl;
                // cout << "ibox_globallll " << ibox_globalll << endl;
                // cout << "ibox_global_m " << ibox_global_m << endl;
                Boxes[ibox_globalll].c[0] += Boxes[ibox_global_m].b[0]*log(zL-zM);
                // cout << "onnn" << endl; 
                for (int expan=1;expan<=p;++expan){
                    Boxes[ibox_globalll].c[0] += Boxes[ibox_global_m].b[expan]/pow((zL-zM),expan)*pow(-1,expan);
                }
                // cout << "innn" << endl;
                std::vector<std::complex<double> > pref(p);
                
                for (int kk = 1; kk <= p; kk++) {
                pref[kk-1] = Boxes[ibox_global_m].b[kk] / pow((zM - zL), kk) * pow(-1, kk);
                for (int ll = 1; ll <= p; ll++) {
                    std::complex<double> sum = (0.0,0.0);
                    for (int pref_long=0;pref_long<pref.size();pref_long++){
                        sum = pref[pref_long]*static_cast<double>(nCk[pref_long+ll-1][ll])/pow(zM-zL,ll);
                    }

                Boxes[ibox_globalll].c[ll] +=  sum  - Boxes[ibox_global_m].b[0]/ static_cast<double>(ll) / pow((zM - zL), ll);
                }            
                }
                // cout << "onnn" << endl;
                
            }
        }
        }




        
        
        
        for (int ibox = 1;ibox<=pow(4,l);ibox++){
            // calculate global ibox
            int ibox_globallll = 0;
            for (int t=1;t<l;t++){
                ibox_globallll += pow(4,t);
            }
            
            ibox_globallll += ibox;
            
            ibox_globallll -= 1; 
            // cout << "idx_globallll " << ibox_globallll << endl;
            // cout << "innn" << endl;
            // cout << "idx_globallll " << ibox_globallll << endl;
            std::complex<double> zL = Boxes[ibox_globallll].center;
            vector<int>  cID = getchildrenID(l,ibox);
            // cout << "cIDb= " << cID[0] << " " << cID[1] << " " << cID[2] << " " << cID[3] << endl;
            int sum2 = 0;
            for (int t=1;t<l;t++){
                sum2 += pow(4,t);
            }

            int sum3 = 0;
            for (int t=1;t<l+1;++t){
                sum3 += pow(4,t);
            }
            //  cout << "sum2 " << sum2 << endl;
            for (int i=0;i<cID.size();i++){
                cID[i] += sum3;
                cID[i] -= 1;
            }
            // cout << "cIDa= " << cID[0] << " " << cID[1] << " " << cID[2] << " " << cID[3] << endl;
            for (int i = 0;i<4;i++){
                int ibox_global_ci = cID[i];
                
                std::complex<double> zT = Boxes[ibox_global_ci].center;
                // cout << "ibox_global_ci" << ibox_global_ci << endl;
                for (int ll=0;ll<p+1;ll++){
                    for(int k=ll;k<p+1;k++){
                        Boxes[ibox_global_ci].d[ll] += Boxes[ibox_globallll].c[k]*static_cast<double>(nCk[k][ll])*pow(zT-zL,k-ll);
                    }
                }

            }
        }

            


    }

    
    // step4:compute interactions at finesest level
    cout << "step4" << endl;
    for (int ibox = 1;ibox<=pow(4,n);ibox++){
            // calculate global ibox
            int ibox_globalll = 0;
            
            for (int t=1;t<n;t++){
                ibox_globalll += pow(4,t);
               
            }
            
            ibox_globalll += ibox;
            ibox_globalll -= 1;
        // cout << "ibox_globalll " << ibox_globalll << endl;
        Boxes[ibox_globalll].c = Boxes[ibox_globalll].d;
        
        std::vector<int> iboxlist = Boxes[ibox_globalll].ilist;
        // cout << "ibox_globalll " << ibox_globalll << " iboxlist.size() " << iboxlist.size() << endl;
        if (!iboxlist.empty()){
            
            std::complex<double> zL = Boxes[ibox_globalll].center;
            for (int i=0;i<iboxlist.size();i++){

           
                int ibox_global_m = 0;
                for (int t=1;t<n;t++){
             
                    ibox_global_m += pow(4,t);
                }


                ibox_global_m = ibox_global_m + iboxlist[i]-1;
                // cout << "innn" << endl;
                std::complex<double> zM = Boxes[ibox_global_m].center;
                // cout << "zL " << zL << " zM " << zM << endl;
                // for (int i=0;i<iboxlist.size();i++){
                //     cout <<  iboxlist[i] << " ";
                // }
                // cout << endl;
                // cout << "ibox_globallll " << ibox_globalll << endl;
                // cout << "ibox_global_m " << ibox_global_m << endl;
                Boxes[ibox_globalll].c[0] += Boxes[ibox_global_m].a[0]*log(zL-zM);
                // cout << "onnn" << endl; 
                for (int expan=1;expan<=p;++expan){
                    Boxes[ibox_globalll].c[0] += Boxes[ibox_global_m].a[expan]/pow((zL-zM),expan)*pow(-1,expan);
                }
                // cout << "innn" << endl;
                std::vector<std::complex<double> > pref(p);
                
                for (int kk = 1; kk <= p; kk++) {
                pref[kk-1] = Boxes[ibox_global_m].a[kk] / pow((zM - zL), kk) * pow(-1, kk);
                for (int ll = 1; ll <= p; ll++) {
                    std::complex<double> sum = (0.0,0.0);
                    for (int pref_long=0;pref_long<pref.size();pref_long++){
                        sum = pref[pref_long]*static_cast<double>(nCk[pref_long+ll-1][ll])/pow(zM-zL,ll);
                    }

                    Boxes[ibox_globalll].c[ll] +=  sum  - Boxes[ibox_global_m].a[0]/ static_cast<double>(ll) / pow((zM - zL), ll);
                }            
                }
                // cout << "onnn" << endl;
                
            }
        }
        }



    // Step 5: Evaluate local expansions at particles positions
    cout << "step5:Evaluate local expansions at particles positions" << endl;
    std::vector<std::complex<double> > Phi(N);

    for (int ibox=1;ibox<=pow(4,n);ibox++){
        
        int ibox_globalll = 0;
            
            for (int t=1;t<n;t++){
                ibox_globalll += pow(4,t);
               
            }
            
        ibox_globalll += ibox;
        ibox_globalll -= 1;

        std::complex<double> zC = Boxes[ibox_globalll].center;
        for (int i : Boxes[ibox_globalll].particlelist){
             for (int ll = 0; ll <= p; ll++) {
                Phi[i]+=Boxes[ibox_globalll].c[ll]*pow((z[i-1]-zC),ll);
             }
        }
    
    }

    // for (int i=0;i<Phi.size();i++){
    //     cout << "Phi" << "[" << i << "] " << Phi[i] << endl;
    // }



    // Step 6: Evaluate potential due to nearest neighbors directly
    cout << "step6:Evaluate potential due to nearest neighbors directly" << endl;

    std::vector<std::complex<double> > Phi_direct0(N);
    std::vector<std::complex<double> > Phi_direct1(N);

    for (int ibox=1;ibox<=pow(4,n);ibox++){
        
        int ibox_globalll = 0;
            
            for (int t=1;t<n;t++){
                ibox_globalll += pow(4,t);
               
            }
            
        ibox_globalll += ibox;
        ibox_globalll -= 1;

        vector<int> neilist = Boxes[ibox_globalll].nlist;
        std::vector<int> particles_inbox = Boxes[ibox_globalll].particlelist;
        std::vector<int> particles_neighbors;
        particles_neighbors.clear();
        for (int i:neilist){
            int ibox_global_n = 0;
            
            for (int t=1;t<n;t++){
                ibox_global_n += pow(4,t);
            }
            
            ibox_global_n += i;
            ibox_global_n -= 1;
            particles_neighbors.insert(particles_neighbors.end(), Boxes[ibox_global_n].particlelist.begin(), Boxes[ibox_global_n].particlelist.end());
        }
        
        // cout << "inn" << endl;
        // cout << "particles_inbox.size() " << particles_inbox.size() << endl;
        std::vector<std::vector<std::complex<double> > > matrix(particles_inbox.size(), std::vector<std::complex<double> >(particles_inbox.size()));
        for (int i = 0;i<particles_inbox.size();i++){
            for (int j = 0;j<particles_inbox.size();j++){

                matrix[i][j] = log(z[particles_inbox[j]-1]-z[particles_inbox[i]-1]);
                if (i==j){
                    matrix[i][j] = 0;
                }
            }
        }
        // cout << "nnnn" << endl;

        
        std::vector<std::complex<double> >sum(particles_inbox.size());
        for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[i].size(); j++) {
            sum[j] += m[particles_inbox[j]-1] * matrix[i][j];
            if (j==particles_inbox.size()-1){
                Phi_direct0[particles_inbox[i]-1] = sum[j];
            }
            }
        }
        
        std::vector<std::vector<std::complex<double> > > matrix2(particles_neighbors.size(), std::vector<std::complex<double> >(particles_inbox.size()));
        for (int i = 0;i<particles_neighbors.size();i++){
            for (int j = 0;j<particles_inbox.size();j++){
                matrix2[i][j] = log(z[particles_inbox[j]-1]-z[particles_neighbors[i]-1]);
            }
        }
        
        std::vector<std::complex<double> >sum2(particles_inbox.size());


        // cout << "matrix2.size() " << matrix2.size() << " particles_neghbors.size() " << particles_neighbors.size()<< endl;
        // cout << "matrix2[0].size() " << matrix2[0].size() << " particles_inbox.size() " << particles_inbox.size()<< endl;
        for (size_t i = 0; i < matrix2.size(); i++) {
            
        for (size_t j = 0; j < matrix2[i].size(); j++) {
            
            sum2[j] += m[particles_neighbors[j]-1] * matrix2[i][j];
            

            // if (j==particles_inbox.size()-1){
            //     // cout << "i="<< i << " j=" << j << " particles_inbox.size()-1= " << particles_inbox.size()-1 << endl;
            //     // cout << "sum2 " << sum2[j] << endl;
            //     // cout << "particles_inbox[j]-1 = " << particles_inbox[j]-1 << endl;
                
            //     // for (int i=0;i<Phi_direct1.size();i++){
            //     //     cout << "Phi_direct1" << "[" << i << "] " << Phi_direct1[i] << endl;
            //     // }
            //     // cout << "Phi_direct1[particles_inbox[i]-1] " << Phi_direct1[particles_inbox[j]-1] << endl;
            //     Phi_direct1[particles_inbox[j]-1] += sum2[j];
            //     sum2[j]=0;
            //     }
            // cout << "i " << i << " j " << j << " end-------- " << endl;
            }
            
        }
        
    }
    std::vector<std::complex<double> > Phi_direct(N);
    // # pragma omp parallel for
    for (int i=0;i<N;i++){
        Phi_direct[i] = Phi_direct0[i] + Phi_direct1[i];
    }


    // Step7:Evaluate Combine direct and far-field terms together

    std::vector<double> Potential(N);
    // # pragma omp parallel for
    for (int i=0;i<N;i++){
        Potential[i] = real(Phi[i]) + real(Phi_direct[i]);
    }



    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double elapsedSeconds = getExecutionTime(start, end);
    std::cout << "Simulation time: " << elapsedSeconds << " seconds" << std::endl;


    // Step8 : Calculate directly
    std::chrono::steady_clock::time_point start2 = std::chrono::steady_clock::now();
    std::vector<std::complex<double> > Phi_direct_n_body(N);
    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            if (i!=j){
                Phi_direct_n_body[i] += m[j]*log(z[i]-z[j]);
            }
        }
    }
    std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
    double elapsedSeconds_Nbody = getExecutionTime(start2, end2);
    std::cout << "Simulation time(Direct N-Body): " << elapsedSeconds_Nbody << " seconds" << std::endl;




    // for (int i=0;i<N;i++){
    //     std::cout << "Phi_direct_n_body[" << i << "] = " << real(Phi_direct_n_body[i]) << "\tPotential[" << i << "] = " << Potential[i] << std::endl;
    // }

    // //     for (int i=0 ; i<Boxes.size() ; i++)
    // // {
    // //     for (int j=0 ; j<Boxes[i].d.size() ; j++)
    // //     {
    // //         std::cout << Boxes[i].d[j] << " ";
    // //     }
    // //     std::cout << std::endl;
    // //     // for (int j=0 ; j<Boxes[i].particlelist.size() ; j++)
    // //     // {
    // //     //     std::cout << Boxes[i].particlelist[j] << " ";
    // //     // }
    // //     // std::cout << std::endl;
    // // }

    return 0;
}
