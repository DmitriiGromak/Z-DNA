#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <unordered_map>
#include <string>
#include <fstream>
#include <array>
#include <algorithm>


double dienerg[4][16] = {
/*         AA    AT    AG    AC    TA    TT    TG    TC    GA    GT    GG    GC    CA    CT    CG    CC */
/* AS-AS */
        {4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66, 2.40},
/* SA-AS */
        {6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00},
/* AS-SA */
        {6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00},
/* SA-SA */
        {4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20, 2.40}
};
double RT = 0.59004, a = 0.357, b = 0.4, sig = exp(-10 / RT), K = 1100 * RT / 4363;
int popul = 80000;
std::unordered_map<char, int> nucl;

void calculate_best_conf(std::vector<std::array<double, 2>> &best_conf, std::vector<double> &best_en, std::string &DNA, std::string &conf, int i,
                         int win) {
    best_conf[0] = {0, 0};
    for (int j = 0; j < win; ++j) {
        best_conf[j + 1][0] = std::min(best_conf[j][0] + dienerg[0][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]],
                                  best_conf[j][1] + dienerg[1][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]]);
        best_conf[j + 1][1] = std::min(best_conf[j][0] + dienerg[2][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]],
                                  best_conf[j][1] + dienerg[3][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]]);
    }
    int ind = 0;
    if (best_conf[win][1] < best_conf[win][0]) {
        ind = 1;
    }
    std::array<std::string, 2> con = {"SA", "AS"};
    for (int j = win - 1; j >= 0; --j) {
        conf += con[ind];
        if (best_conf[j][0] + dienerg[ind * 2][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]] ==
            best_conf[j + 1][ind]) {
            best_en[j] = dienerg[ind * 2][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]];
            ind = 0;
        } else {
            best_en[j] = dienerg[ind * 2 + 1][nucl[DNA[i + 2 * j]] * 4 + nucl[DNA[i + 2 * j + 1]]];
            ind = 1;
        }
    }
    reverse(conf.begin(), conf.end());
}

double calculate_tension(std::vector<double> &best_en, int win) {
    double l = -20, r = 40, eps = 1e-4;
    while (r - l > eps) {
        double mid = (l + r) / 2;
        double q = exp(-K * mid * mid / RT), tw = 0;
        for (int j = 1; j <= win; ++j) {
            double sm = 0, prod = 1;
            for (int k = 0; k < win; ++k) {
                prod *= exp(-best_en[k] / RT);
                if (k >= j - 1) {
                    sm += prod;
                    prod /= exp(-best_en[k - j + 1] / RT);
                }
            }
            q += (win - j + 1) * sig * sm * exp(-K * (mid - a * j - 2 * b) * (mid - a * j - 2 * b) / RT);
            tw += (a * j + 2 * b) * (win - j + 1) * sig * sm *
                  exp(-K * (mid - a * j - 2 * b) * (mid - a * j - 2 * b) / RT);
        }
        tw /= q;
        if (tw > 1) r = mid;
        else l = mid;
    }
    return r;
}

int main() {
    std::string DNA;
    int mn_window, mx_window, length, top;
    std::string name;
    std::cin >> mn_window >> mx_window >> top >> name;
    nucl['t'] = 1, nucl['g'] = 2, nucl['c'] = 3;
    std::vector<double> pop_av, dev;
    for (int win = mn_window; win <= mx_window; ++win) {
        pop_av.push_back(0), dev.push_back(0);
        std::vector<double> sc;
        for (int i = 0; i < popul; ++i) {
            std::string seq;
            std::array<char, 4> nu = {'a', 't', 'g', 'c'};
            for (int j = 0; j < 2 * win; ++j) {
                seq += nu[rand() % 4];
            }
            std::vector<std::array<double, 2>> best_conf(win + 1);
            std::vector<double> best_en(win);
            std::string conf;
            calculate_best_conf(best_conf, best_en, seq, conf, 0, win);
            sc.push_back(calculate_tension(best_en, win));
            pop_av.back() += sc.back();
        }
        pop_av.back() /= popul;
        for (int i = 0; i < popul; ++i) {
            dev.back() += (sc[i] - pop_av.back()) * (sc[i] - pop_av.back());
        }
        dev.back() = sqrt(dev.back() / popul);
    }
    std::ifstream fin(name);
    std::string row;
    getline(fin, row);
    while (getline(fin, row)) DNA += row;
    length = DNA.size();
    for (int i = 0; i < length; ++i) {
        DNA[i] = tolower(DNA[i]);
    }
    std::vector<std::pair<double, int>> z_sc;
    int cnt = 0;
    std::vector<std::pair<int, int>> pos;
    std::vector<std::string> all_conf;
    for (int win = mn_window; win <= mx_window; ++win) {
        std::vector<std::array<double, 2>> best_conf(win + 1);
        std::vector<double> best_en(win);
        for (int i = 0; i <= length - 2 * win; ++i) {
            std::string conf;
            calculate_best_conf(best_conf, best_en, DNA, conf, i, win);
            double res = calculate_tension(best_en, win);
            double score = (res - pop_av[win - mn_window]) / dev[win - mn_window];
            double erf_z = erf(score / sqrt(2));
            double probability = 0.5 * (1 + erf_z);
            z_sc.push_back({2 * win / probability, cnt});
            pos.push_back({i, i + 2 * win - 1});
            all_conf.push_back(conf);
            cnt++;
        }
    }
    sort(z_sc.rbegin(), z_sc.rend());
    std::set<std::pair<int, int>> reg;
    for (int i = 0; top && i < z_sc.size(); ++i) {
        int ind = z_sc[i].second;
        auto it = reg.lower_bound(std::make_pair(pos[ind].first, 0));
        if ((it == reg.end() || it->first > pos[ind].second) &&
            (it == reg.begin() || (--it)->second < pos[ind].first)) {
            std::cout << z_sc[i].first << " " << pos[ind].first << " " << pos[ind].second << " " << all_conf[ind] << '\n';
            reg.insert(pos[ind]);
            top--;
        }
    }
}
