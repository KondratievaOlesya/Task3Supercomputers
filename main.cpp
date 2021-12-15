/*
    1 Вариант
    Кондратьева Олеся, 628 группа, МГУ
*/
#include <math.h>
#include <vector>
#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <float.h>
#include <string>
#include <sstream>

#define M_PI 3.14159265358979323846 /* pi */

#define debug false

// Border indexies
#define RIGHT 0
#define LEFT 1
#define UPPER 2
#define BOTTOM 3
#define FRONT 4
#define BACK 5

#define MPI_Isend_TAG 1
using namespace std;

struct block
{
    int i_start;
    int i_end;
    int j_start;
    int j_end;
    int k_start;
    int k_end;
};

struct variables {
    // Boundaries
    double Lx;
    double Ly;
    double Lz;
    // Steps
    double hx, hz, hy;
    // Time step
    double tao;
    // Amount of nodes
    int N;

    vector<int> factor;
};

int find_distance(vector<int> &three){
    return abs(three[0] - three[1]) + abs(three[0] - three[2]) + abs(three[1] - three[2]);
}

// Split num on three factors with minimum distance
vector<int> split_by_three_factors(int num) {
    vector<int> res(3);
    int dist = num * 3 + 1;
    for (int i = 1; i <= num; i++) 
        for (int j = 1; j <= num / i; j++)
        {
            if (num % (i * j) == 0) {
                vector<int> comp;
                comp.push_back(i);
                comp.push_back(j);
                comp.push_back(num / (i * j));
                int d = find_distance(comp);
                if (d < dist) {
                    res = comp;
                    dist = d;
                }
            }

        }
    return res;
}

block make_block(int size, int rank, variables& L) {
    int blocks_length_i = ceil(L.N * 1.0/ L.factor[0]); 
    int blocks_length_j = ceil(L.N * 1.0/ L.factor[1]); 
    int blocks_length_k = ceil(L.N * 1.0/ L.factor[2]); 

    int level = rank / (L.factor[2] * L.factor[1]); 
    int row = (rank - level * (L.factor[2] * L.factor[1])) / L.factor[2]; 
    int col = rank - level * (L.factor[2] * L.factor[1]) - row * L.factor[2];
    block b;

    b.i_start = col * blocks_length_k;
    b.i_end = min(L.N, b.i_start + blocks_length_k);

    b.j_start = row * blocks_length_j;
    b.j_end = min(L.N, b.j_start + blocks_length_j);

    b.k_start = level * blocks_length_i;
    b.k_end = min(L.N, b.k_start + blocks_length_i);
    return b;
}

int volume(const block &b)
{
    return ((b.i_end - b.i_start) * (b.j_end - b.j_start) * (b.k_end - b.k_start));
}

// u (x, y, z, t)
double u_analytical(double x, double y, double z, double t, variables &L)
{
    double at = M_PI * sqrt(1.0 / (L.Lx * L.Lx) + 1.0 / (L.Ly * L.Ly) + 1.0 / (L.Lz * L.Lz));
    return sin(M_PI / L.Lx * x) * sin(M_PI / L.Ly * y) * sin(M_PI / L.Lz * z) * cos(at * t);
}

bool in_same_square(int rank1, int rank2, vector<int> &factor) {
    // L.factor[2] - i; L.factor[1] - j; L.factor[0] - k

    int level1 = rank1 / (factor[2] * factor[1]);
    int level2 = rank2 / (factor[2] * factor[1]);
    return level1 == level2;
}
// Translate i, j, k to x, y, z
double u_ijk(int i, int j, int k, double t, variables & L)
{
    double x = i * L.hx;
    double y = j * L.hy;
    double z = k * L.hz;
    return u_analytical(x, y, z, t, L);
};


// Return index in u_calculated for (i, j, k)
int index(int i, int j, int k, block &b)
{
    int idx = (i - b.i_start) * (b.j_end - b.j_start) * (b.k_end - b.k_start) + (j - b.j_start) * (b.k_end - b.k_start) + (k - b.k_start);
    return idx;
}

// Finds any value including borders
double find_border_value(int i, int j, int k, block &b, vector<double> &u, vector<vector<double> > &borders, variables &L)
{
    if (b.i_start <= i && i < b.i_end && b.j_start <= j && j < b.j_end && b.k_start <= k && k < b.k_end)
    {
        // It is not a border
        return u[index(i, j, k, b) ];
    }

    bool is_edge_of_grid = i == 0 || j == 0 || k == 0;
    is_edge_of_grid = is_edge_of_grid || (i == (L.N - 1)) || (j == (L.N - 1)) || (k == (L.N - 1));    
    if (is_edge_of_grid) {
        return 0.0;
    }


    if (i < b.i_start)
    {
        int idx = (b.k_end - b.k_start) * (j - b.j_start) + (k - b.k_start);
        return borders[LEFT][idx];
    }
    if (i >= b.i_end)
    {
        int idx = (b.k_end - b.k_start) * (j - b.j_start) + (k - b.k_start);
        return borders[RIGHT][idx];
    }

    if (k < b.k_start)
    {
        int idx = (b.j_end - b.j_start) * (i - b.i_start) + (j - b.j_start);
        return borders[BOTTOM][idx];
    }

    if (k >= b.k_end)
    {
        int idx = (b.j_end - b.j_start) * (i - b.i_start) + (j - b.j_start);
        return borders[UPPER][idx];
    }

    if (j < b.j_start)
    {
        int idx = (b.k_end - b.k_start) * (i - b.i_start) + (k - b.k_start);
        return borders[FRONT][idx];
    }

    if (j >= b.j_end)
    {
        int idx = (b.k_end - b.k_start) * (i - b.i_start) + (k - b.k_start);
        return borders[BACK][idx];
    }
    return -9999.0;
}

// Laplas operator on u_calculated
double laplas(int i, int j, int k, block &b, vector<double> &u_calculated, vector<vector<double> > &borders, variables & L)
{
    double sum_i = (find_border_value(i - 1, j, k, b, u_calculated, borders, L) - 2.0 * u_calculated[index(i, j, k, b)] + find_border_value(i + 1, j, k, b, u_calculated, borders, L)) / (L.hx * L.hx);
    double sum_j = (find_border_value(i, j - 1, k, b, u_calculated, borders, L) - 2.0 * u_calculated[index(i, j, k, b)] + find_border_value(i, j + 1, k, b, u_calculated, borders, L)) / (L.hy * L.hy);
    double sum_k = (find_border_value(i, j, k - 1, b, u_calculated, borders, L) - 2.0 * u_calculated[index(i, j, k, b)] + find_border_value(i, j, k + 1, b, u_calculated, borders, L)) / (L.hz * L.hz);
    return sum_i + sum_j + sum_k;
};

void fill_borders(vector<vector<double> > &borders, block &b, vector<double> &layer, int rank, int size, variables & L)
{
    // I need to send 6 blocks
    // Right block wants i_end

    if (rank + 1 < size)
    {
        int i = b.i_end - 1;
        for (int j = b.j_start; j < b.j_end; j++)
            for (int k = b.k_start; k < b.k_end; k++)
                borders[RIGHT][(j - b.j_start) * (b.k_end - b.k_start) + (k - b.k_start)] = layer[index(i, j, k, b)];
    }
    // Left block wants i_start 
    if (rank - 1 >= 0)
    {
        int i = b.i_start;

        for (int j = b.j_start; j < b.j_end; j++)
            for (int k = b.k_start; k < b.k_end; k++)
                borders[LEFT][(j - b.j_start) * (b.k_end - b.k_start) + (k - b.k_start)] = layer[index(i, j, k, b)];
    }
    // Back block wants j_end
    if (rank + L.factor[2] < size && in_same_square(rank, rank + L.factor[2], L.factor))
    {
        int j = b.j_end - 1;
        for (int i = b.i_start; i < b.i_end; i++)
            for (int k = b.k_start; k < b.k_end; k++)
                borders[BACK][(i - b.i_start) * (b.k_end - b.k_start) + (k - b.k_start)] = layer[index(i, j, k, b)];
    }

    // Front block wants j_start
    if (rank - L.factor[2] >= 0 && in_same_square(rank, rank - L.factor[2], L.factor))
    {
        int j = b.j_start;
        for (int i = b.i_start; i < b.i_end; i++)
            for (int k = b.k_start; k < b.k_end; k++)
                borders[FRONT][(i - b.i_start) * (b.k_end - b.k_start) + (k - b.k_start)] = layer[index(i, j, k, b)];
    }
    // Bottom block wants k_start
    if (rank - L.factor[2] * L.factor[1] >= 0)
    {
        int k = b.k_start;
        for (int i = b.i_start; i < b.i_end; i++)
            for (int j = b.j_start; j < b.j_end; j++)
                borders[BOTTOM][(i - b.i_start) * (b.j_end - b.j_start) + (j - b.j_start)] = layer[index(i, j, k, b)];
    }
    // Upper block wants k_end
    if (rank + L.factor[2] * L.factor[1] < size)
    {
        int k = b.k_end - 1;
        for (int i = b.i_start; i < b.i_end; i++)
            for (int j = b.j_start; j < b.j_end; j++)
                borders[UPPER][(i - b.i_start) * (b.j_end - b.j_start) + (j - b.j_start)] = layer[index(i, j, k, b)];
    }
}

void send_recv_borders(vector<vector<double> > &s_borders, vector<vector<double> >& r_borders, int rank, int size, block &b, variables& L, vector<double>& layer)
{
    fill_borders(s_borders, b, layer, rank, size, L);
    // L.factor[2] - i; L.factor[1] - j; L.factor[0] - k

    vector<MPI_Request> request(12);
    vector<MPI_Status> status(12);

    // Right block
    int c = 0;
    if (rank + 1 < size)
    {

        MPI_Isend(s_borders[RIGHT].data(), s_borders[RIGHT].size(), MPI_DOUBLE, rank + 1, MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[RIGHT].data(), r_borders[RIGHT].size(), MPI_DOUBLE, rank + 1, MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
    }
    // Left block
    if (rank - 1 >= 0)
    {

        MPI_Isend(s_borders[LEFT].data(), s_borders[LEFT].size(), MPI_DOUBLE, rank - 1, MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[LEFT].data(), r_borders[LEFT].size(), MPI_DOUBLE, rank - 1, MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;

        
    }
    // Upper block
    if (rank + L.factor[2] * L.factor[1] < size)
    {

        MPI_Isend(s_borders[UPPER].data(), s_borders[UPPER].size(), MPI_DOUBLE, rank + L.factor[2] * L.factor[1], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[UPPER].data(), r_borders[UPPER].size(), MPI_DOUBLE, rank + L.factor[2] * L.factor[1], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        
    }
    // Bottom block
    if (rank - L.factor[2] * L.factor[1] >= 0)
    {

        MPI_Isend(s_borders[BOTTOM].data(), s_borders[BOTTOM].size(), MPI_DOUBLE, rank - L.factor[2] * L.factor[1], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[BOTTOM].data(), r_borders[BOTTOM].size(), MPI_DOUBLE, rank - L.factor[2] * L.factor[1], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        
    }
    // Block in front
    if (rank - L.factor[2] >= 0 && in_same_square(rank, rank - L.factor[2], L.factor))
    {

        MPI_Isend(s_borders[FRONT].data(), s_borders[FRONT].size(), MPI_DOUBLE, rank - L.factor[2], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[FRONT].data(), r_borders[FRONT].size(), MPI_DOUBLE, rank - L.factor[2], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        
    }
    // Block behind
    if (rank + L.factor[2] < size && in_same_square(rank, rank + L.factor[2], L.factor))
    {

        MPI_Isend(s_borders[BACK].data(), s_borders[BACK].size(), MPI_DOUBLE, rank + L.factor[2], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        MPI_Irecv(r_borders[BACK].data(), r_borders[BACK].size(), MPI_DOUBLE, rank + L.factor[2], MPI_Isend_TAG, MPI_COMM_WORLD, &request[c]);
        c += 1;
        
    }
    if (c != 0)
        MPI_Waitall(c, request.data(), status.data());
}


double calculate_max_error(vector<double> &u, block &b, double t, variables &L)
{
    //double res = -DBL_MAX;
    double res = 0;
// #pragma omp parallel for collapse(3) reduction(max: res)
    for (int i = b.i_start; i < b.i_end; i++)
        for (int j = b.j_start; j < b.j_end; j++)
            for (int k = b.k_start; k < b.k_end; k++)
            {
                bool is_edge_of_grid = i == 0 || j == 0 || k == 0;
                is_edge_of_grid = is_edge_of_grid || (i == (L.N - 1)) || (j == (L.N - 1)) || (k == (L.N - 1));
                double error;
                if (is_edge_of_grid) {
                    error = fabs(u[index(i, j, k, b)] - 0.0);
                }
                else 
                    error = fabs(u[index(i, j, k, b)] - u_ijk(i, j, k, t, L));
                res = max(error, res);
            }
    return res;
}

ostream &operator<<(ostream &out, const block &b)
{
    out << "Block:" << endl;
    out << "i: from " << b.i_start << " to " << b.i_end << endl;
    out << "j: from " << b.j_start << " to " << b.j_end << endl;
    out << "k: from " << b.k_start << " to " << b.k_end << endl;
    out << "volume = " << volume(b);
    return out;
}


int main(int argc, char *argv[])
{
    // Check the number of parameters
    if (argc < 7)
    {
        // Tell the user how to run the program
        cerr << "Usage: " << argv[0] << " L.Lx L.Ly L.Lz T K L.N [OUTPUT]" << endl;
        cerr << "    L.Lx, L.Ly, L.Lz     Input parameters of analytical function" << endl;
        cerr << "    T              Time end range" << endl;
        cerr << "    K              Time steps total" << endl;
        cerr << "    L.N              Number of points in grid" << endl;
        return 1;
    }
    variables L;
    // Input region boundaries
    L.Lx = atof(argv[1]);
    L.Ly = atof(argv[2]);
    L.Lz = atof(argv[3]);

    // Time parameters
    int T, K;
    T = atoi(argv[4]);
    K = atoi(argv[5]);
    L.N = atoi(argv[6]);

    // Time step
    L.tao = T * 1.0 / K;

    L.hx = L.Lx * 1.0 / L.N;
    L.hy = L.Ly * 1.0 / L.N;
    L.hz = L.Lz * 1.0 / L.N;

    // File to write result to
    const char *out_file = "result.txt";
    if (argc == 8)
    {
        out_file = argv[7];
    }
    // Final error of each layer
    vector<double> total_error(K);
    double time_res;
    MPI_Init(&argc, &argv);
    // MPI parameters
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // Running time
    double start = MPI_Wtime();
    L.factor = split_by_three_factors(size);
    
    ofstream tlog;
    if (debug) {
        // Because BG/P don't have to_string
        ostringstream ss;
        ss << "log" << rank << ".txt";
        tlog.open((ss.str()).c_str());
    }
    if (debug)
        tlog << L.factor[0] << ", " << L.factor[1] << ", " << L.factor[2] << endl;

    // Block edge start and end
    block b = make_block(size, rank, L);
    if (debug)
        tlog << rank << ": " << b << endl;

    // Errors of each layer
    vector<double> error(min(K, 21));

    vector<vector<double> > u_calculated(3, vector<double>(volume(b))); // u[n][i, j, k] > >:)
    if (debug)
        tlog << rank << ": " << "u_calculated size = " << u_calculated[0].size() << endl;

    // 0 layer
    for (int i = b.i_start; i < b.i_end; i++)
        for (int j = b.j_start; j < b.j_end; j++)
            for (int k = b.k_start; k < b.k_end; k++)
            {
                // It is allowed to use analytical value, right?
                // phi (xi, yj, zk)
                u_calculated[0][index(i, j, k, b)] = u_ijk(i, j, k, 0, L);
            }
    if (debug)
        tlog << rank << ": " << "u[0] layer calculated" << endl;

    // It needs border values of u_calculated[0]
    // There is six of edges and I need to store it somehow and use
    // Laplas operator will try to adress -1 value and +1 from size in each direction i, j, k
    vector<vector<double> > s_borders(6); 
    vector<vector<double> > r_borders(6);

    s_borders[RIGHT].resize((b.j_end - b.j_start) * (b.k_end - b.k_start));
    s_borders[LEFT].resize((b.j_end - b.j_start) * (b.k_end - b.k_start));
    s_borders[BACK].resize((b.i_end - b.i_start) * (b.k_end - b.k_start));
    s_borders[FRONT].resize((b.i_end - b.i_start) * (b.k_end - b.k_start));
    s_borders[BOTTOM].resize((b.j_end - b.j_start) * (b.i_end - b.i_start));
    s_borders[UPPER].resize((b.j_end - b.j_start) * (b.i_end - b.i_start));

    for (int i = 0; i < 6; i++) {
        r_borders[i].resize(s_borders[i].size());
    }
    send_recv_borders(s_borders, r_borders, rank, size, b, L, u_calculated[0]);
    if (debug)
        tlog << rank << ": " << "Exchange finished" << endl;
    // 1 layer
    for (int i = b.i_start; i < b.i_end; i++)
        for (int j = b.j_start; j < b.j_end; j++)
            for (int k = b.k_start; k < b.k_end; k++)
            {
                bool is_edge_of_grid = i == 0 || j == 0 || k == 0;
                is_edge_of_grid = is_edge_of_grid || (i == (L.N - 1)) || (j == (L.N - 1)) || (k == (L.N - 1));
                if (!is_edge_of_grid)
                {
                    u_calculated[1][index(i, j, k, b)] = u_calculated[0][index(i, j, k, b)] + L.tao * L.tao / 2.0 * laplas(i, j, k, b, u_calculated[0], r_borders, L);
                }
                else
                {
                    // First-order conditions (all)
                    u_calculated[1][index(i, j, k, b)] = 0.0;
                }
                // It needs border values of 0 layer for laplas operation
            }

    send_recv_borders(s_borders, r_borders, rank, size, b, L, u_calculated[1]);

    // Caculate max error on layer
    error[1] = calculate_max_error(u_calculated[1], b, L.tao, L);
    if (debug)
        tlog << rank << ": " << "u[1] layer calculated with error = " << error[1] << endl;

    // Rest of layers
    for (int n = 2; n < min(21, K); n++)
    {
        for (int i = b.i_start; i < b.i_end; i++)
            for (int j = b.j_start; j < b.j_end; j++)
                for (int k = b.k_start; k < b.k_end; k++)
                {
                    bool is_edge_of_grid = i == 0 || j == 0 || k == 0;
                    is_edge_of_grid = is_edge_of_grid || (i == (L.N - 1)) || (j == (L.N - 1)) || (k == (L.N - 1));
                    if (!is_edge_of_grid)
                    {
                        u_calculated[n % 3][index(i, j, k, b)] = laplas(i, j, k, b, u_calculated[(n - 1) % 3], r_borders, L) * L.tao * L.tao + 2.0 * u_calculated[(n - 1) % 3][index(i, j, k, b)] - u_calculated[(n - 2) % 3][index(i, j, k, b)];
                    }
                    else
                    {
                        // First-order conditions (all)
                        u_calculated[n % 3][index(i, j, k, b)] = 0.0;
                    }
                }
        // u_calculated[0] = u_calculated[1];
        // u_calculated[1] = u_calculated[2];
        if (debug)
            tlog << rank << ": " << "u[" << n << "] layer calculated "<< endl;
        // Caculate max error on layer
        error[n] = calculate_max_error(u_calculated[n % 3], b, n * L.tao, L);
        if (debug)
            tlog << rank << ": " << "u[" << n << "] layer calculated with error = " << error[n] << endl;
        // Prepare for next layer
        send_recv_borders(s_borders, r_borders, rank, size, b, L, u_calculated[n % 3]);
        if (debug)
            tlog << rank << ": " << "Exchange finished" << endl;
    }
    double end = MPI_Wtime();
    double curr_time_res = end - start;
    MPI_Reduce(&curr_time_res, &time_res, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(error.data(), total_error.data(), K, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0){
        ofstream fout(out_file);
        fout.precision(10);
        for(int i = 1; i < min(21, K); i++){
            fout << "Error of layer " << i << " = " << total_error[i] << endl;
        }
        fout << "Running time = " << time_res << endl;
        fout.close();
    }
    MPI_Finalize();
}