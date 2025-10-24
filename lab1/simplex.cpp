#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>
using namespace std;

void run_simplex(const string &infile, const string &outfile) {
    ifstream fin(infile);
    ofstream fout(outfile);
    int m, n;
    fin >> m >> n;
    vector<vector<double>> t(m + 1, vector<double>(n));
    for (int i = 0; i <= m; i++)
        for (int j = 0; j < n; j++)
            fin >> t[i][j];
    fin.close();

    vector<int> basis(m);
    for (int i = 0; i < m; i++) basis[i] = n - m + i - 1;

    while (true) {
        int l = -1;
        for (int j = 0; j < n - 1; j++)
            if (t[m][j] < 0) { l = j; break; }  // prima coloană negativă — Bland
        if (l == -1) break;

        bool unbounded = true;
        for (int i = 0; i < m; i++)
            if (t[i][l] > 0) { unbounded = false; break; }
        if (unbounded) {
            fout << "Problema este neacotata\n";
            fout.close();
            return;
        }

        int k = -1;
        double minRatio = numeric_limits<double>::infinity();
        for (int i = 0; i < m; i++) {
            if (t[i][l] > 0) {
                double ratio = t[i][n - 1] / t[i][l];
                if (ratio < minRatio || (ratio == minRatio && (k == -1 || basis[i] < basis[k]))) {
                    minRatio = ratio;
                    k = i;
                }
            }
        }

        for (int i = 0; i <= m; i++)
            if (i != k)
                for (int j = 0; j < n; j++)
                    if (j != l)
                        t[i][j] = (t[i][j] * t[k][l] - t[i][l] * t[k][j]) / t[k][l];
        for (int i = 0; i <= m; i++)
            if (i != k)
                t[i][l] = 0;
        for (int j = 0; j < n; j++)
            if (j != l)
                t[k][j] /= t[k][l];
        t[k][l] = 1;
        basis[k] = l;
    }

    fout << fixed << setprecision(3);
    fout << "Coeficientii functiei de minimizat:\n";
    for (int j = 0; j < n - 1; j++) fout << t[m][j] << " ";
    fout << "\n\nSolutie de baza optima:\n";

    vector<double> sol(n - 1, 0.0);
    for (int i = 0; i < m; i++)
        if (basis[i] < n - 1)
            sol[basis[i]] = t[i][n - 1];
    for (int j = 0; j < n - 1; j++)
        fout << "x" << j + 1 << " = " << sol[j] << "\n";

    double val = t[m][n - 1];
    if (val > 0) val = -val;
    fout << "\nValoarea minima: " << val << "\n";
    fout.close();
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Utilizare: ./simplex_exec <fisier_intrare> <fisier_iesire>\n";
        return 1;
    }
    run_simplex(argv[1], argv[2]);
    return 0;
}
