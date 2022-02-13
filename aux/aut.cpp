#include <bits/stdc++.h>


typedef long long ll;

using namespace std;


void generate(vector<int>& rho, vector<int>& tau, vector<int> perm, set<vector<int>>& seen, int depth) {
    if (seen.size() % 1000000 == 0) {
        cout << seen.size() << endl;
    }
    if (seen.find(perm) != seen.end()) {
        return;
    }
    if (depth > 1000) {
        return;
    }
    seen.insert(perm);
    vector<int> perm_tau;
    vector<int> perm_rho;
    for (int i = 0; i < perm.size(); i++) {
        perm_tau.push_back(perm[tau[i]]);
        perm_rho.push_back(perm[rho[i]]);
    }
    generate(rho, tau, perm_tau, seen, depth+1);
    generate(rho, tau, perm_rho, seen, depth+1);
}

int main() {
    vector<int> rho = {2,0,1,5,3,4,8,6,7,11,9,10,14,12,13,17,15,16,20,18,19,23,12,22};
    vector<int> tau = {19, 0, 1, 2, 3, 23, 17, 6, 15, 4, 8, 18, 20, 14, 10, 9, 21, 16, 12, 13, 5, 11, 7, 22};
    vector<int> id  = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

    set<vector<int>> seen;
    
    generate(rho,tau,id,seen,0);

    cout << seen.size() << endl;
}
