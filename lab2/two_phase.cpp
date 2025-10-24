#include <bits/stdc++.h>
using namespace std;

static const double EPS = 1e-9;

struct LP {
    int m, n;
    vector<double> c;
    vector<vector<double>> A;
    vector<char> sign;
    vector<double> b;
};

struct Tableau {
    int m, n;
    vector<vector<double>> t;
    vector<int> basis;
};

int bland_entering(const Tableau& T) {
    int m=T.m, n=T.n;
    int enter=-1;
    for (int j=0;j<n-1;j++) if (T.t[m][j] < -EPS) { enter=j; break; }
    return enter;
}

int bland_leaving(const Tableau& T, int l) {
    int m=T.m, n=T.n;
    int k=-1;
    double best = numeric_limits<double>::infinity();
    for (int i=0;i<m;i++) if (T.t[i][l] > EPS) {
        double r = T.t[i][n-1] / T.t[i][l];
        if (r < best - EPS) { best=r; k=i; }
        else if (fabs(r - best) <= EPS) {
            if (k==-1 || T.basis[i] < T.basis[k]) k=i;
        }
    }
    return k;
}

bool simplex(Tableau& T) {
    int m=T.m, n=T.n;
    while (true) {
        int l = bland_entering(T);
        if (l==-1) break;
        int k = bland_leaving(T, l);
        if (k==-1) return false;
        double piv = T.t[k][l];
        for (int j=0;j<n;j++) T.t[k][j] /= piv;
        for (int i=0;i<=m;i++) if (i!=k) {
            double f = T.t[i][l];
            if (fabs(f) > EPS) for (int j=0;j<n;j++) T.t[i][j] -= f*T.t[k][j];
        }
        T.basis[k] = l;
    }
    return true;
}

Tableau build_phase1(const LP& lp, vector<int>& orig_idx, vector<int>& art_idx) {
    int m=lp.m, n=lp.n;
    int slack=0, surplus=0, art=0;
    for (int i=0;i<m;i++) {
        if (lp.sign[i]=='<') slack++;
        else if (lp.sign[i]=='>') { surplus++; art++; }
        else { art++; }
    }
    int cols = n + slack + surplus + art + 1;
    Tableau T;
    T.m = m; T.n = cols;
    T.t.assign(m+1, vector<double>(cols, 0.0));
    T.basis.assign(m, -1);

    orig_idx.resize(n);
    for (int j=0;j<n;j++) orig_idx[j]=j;

    int s_pos = n;
    int r_pos = s_pos + slack;
    int a_pos = r_pos + surplus;

    int cur_s=0, cur_r=0, cur_a=0;
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) T.t[i][j]=lp.A[i][j];
        if (lp.sign[i]=='<') {
            T.t[i][s_pos+cur_s]=1.0;
            T.basis[i]=s_pos+cur_s;
            cur_s++;
        } else if (lp.sign[i]=='>') {
            T.t[i][r_pos+cur_r]=-1.0;
            T.t[i][a_pos+cur_a]=1.0;
            T.basis[i]=a_pos+cur_a;
            cur_r++; cur_a++;
        } else {
            T.t[i][a_pos+cur_a]=1.0;
            T.basis[i]=a_pos+cur_a;
            cur_a++;
        }
        T.t[i][cols-1]=lp.b[i];
    }

    for (int j=0;j<cols;j++) T.t[m][j]=0.0;
    for (int j=0;j<art;j++) T.t[m][a_pos+j]=1.0;
    for (int i=0;i<m;i++) {
        int bj=T.basis[i];
        if (bj>=a_pos && bj<a_pos+art) {
            for (int j=0;j<cols;j++) T.t[m][j]-=T.t[i][j];
        }
    }

    art_idx.clear();
    for (int j=0;j<art;j++) art_idx.push_back(a_pos+j);
    return T;
}

void drop_artificials(Tableau& T, const vector<int>& art_idx) {
    vector<int> keep;
    for (int j=0;j<T.n-1;j++) {
        bool is_art=false;
        for (int a: art_idx) if (j==a) { is_art=true; break; }
        if (!is_art) keep.push_back(j);
    }
    int new_n = (int)keep.size() + 1;
    vector<vector<double>> nt(T.m+1, vector<double>(new_n,0.0));
    for (int i=0;i<=T.m;i++) {
        for (int c=0;c<(int)keep.size();c++) nt[i][c]=T.t[i][keep[c]];
        nt[i][new_n-1]=T.t[i][T.n-1];
    }
    vector<int> map_idx(T.n, -1);
    for (int c=0;c<(int)keep.size();c++) map_idx[keep[c]]=c;
    for (int i=0;i<T.m;i++) {
        if (T.basis[i]>=0) {
            int nb = map_idx[T.basis[i]];
            if (nb==-1) {
                int enter=-1;
                for (int c=0;c<(int)keep.size();c++) if (fabs(nt[i][c])>EPS) { enter=c; break; }
                if (enter!=-1) {
                    double piv=nt[i][enter];
                    for (int j=0;j<new_n;j++) nt[i][j]/=piv;
                    for (int r=0;r<=T.m;r++) if (r!=i) {
                        double f=nt[r][enter];
                        if (fabs(f)>EPS) for (int j=0;j<new_n;j++) nt[r][j]-=f*nt[i][j];
                    }
                    nb=enter;
                } else nb=-1;
            }
            T.basis[i]=nb;
        }
    }
    T.t.swap(nt);
    T.n=new_n;
}

void set_phase2_objective(Tableau& T, const LP& lp, const vector<int>& orig_idx) {
    int m=T.m, n=T.n;
    vector<double> obj(n,0.0);
    for (int j=0;j<(int)orig_idx.size();j++) obj[j]=lp.c[orig_idx[j]];
    for (int i=0;i<=m;i++) T.t[m][i]=0.0;
    for (int j=0;j<n;j++) T.t[m][j]=0.0;
    for (int j=0;j<(int)orig_idx.size();j++) T.t[m][j]=obj[j];
    T.t[m][n-1]=0.0;
    for (int i=0;i<m;i++) {
        int bj = T.basis[i];
        if (bj>=0 && bj<(int)orig_idx.size()) {
            double coeff = T.t[m][bj];
            if (fabs(coeff)>EPS) {
                for (int j=0;j<n;j++) T.t[m][j]-=coeff*T.t[i][j];
            }
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string which;
    if (!(cin>>which)) return 0;
    int m,n;
    cin>>m>>n;
    LP lp;
    lp.m=m; lp.n=n;
    lp.c.assign(n,0.0);
    for (int j=0;j<n;j++) cin>>lp.c[j];
    lp.A.assign(m, vector<double>(n,0.0));
    lp.sign.assign(m,'=');
    lp.b.assign(m,0.0);
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) cin>>lp.A[i][j];
        string s; cin>>s;
        lp.sign[i] = (s=="<=" ? '<' : (s==">=" ? '>' : '='));
        cin>>lp.b[i];
    }

    vector<int> orig_idx, art_idx;
    Tableau T1 = build_phase1(lp, orig_idx, art_idx);
    bool ok1 = simplex(T1);
    if (!ok1) { cout<<"infeasible\n"; return 0; }
    if (fabs(T1.t[T1.m][T1.n-1]) > 1e-6) { cout<<"infeasible\n"; return 0; }
    drop_artificials(T1, art_idx);
    set_phase2_objective(T1, lp, orig_idx);
    bool ok2 = simplex(T1);
    if (!ok2) { cout<<"unbounded\n"; return 0; }

    vector<double> sol(lp.n,0.0);
    for (int i=0;i<T1.m;i++) {
        int bj=T1.basis[i];
        if (bj>=0 && bj<lp.n) sol[bj]=T1.t[i][T1.n-1];
    }
    double val = T1.t[T1.m][T1.n-1];
    if (val>0) val = -val;

    cout.setf(std::ios::fixed); cout<<setprecision(6);
    cout<<"solution:\n";
    for (int j=0;j<lp.n;j++) cout<<"x"<<(j+1)<<"="<<sol[j]<<"\n";
    cout<<"z_min="<<val<<"\n";
    return 0;
}
