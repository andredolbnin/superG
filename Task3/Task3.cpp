#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <string>

using namespace std;

int M; // по x
int N; // по y

double h1;
double h2;

int x_amount;
int x_offset;
int y_amount;
int y_offset;

int debug = 0;

double u(double x, double y) {
    return sqrt(4.0 + x * y);
}

double q(double x, double y) {
    return x + y;
}

double k(double x, double y) {
    return 4.0 + x + y;
}

double F(double x, double y) {
    return -(1.0 / (4.0 * (4.0 + x * y) * sqrt(4.0 + x * y))) * (y * ((x - 4.0) * y - y * y + 8.0) + x * ((y - 4.0) * x - x * x + 8.0)) + q(x, y) * u(x, y);
}

double psi_L(double y) {
    return -(4.0 + y) * y / 4.0 + 2.0;
}

double psi_R(double y) {
    return (8.0 + y) * y / (4.0 * sqrt(1.0 + y)) + 2.0 * sqrt(1.0 + y);
}

double psi_B(double x) {
    return -(4.0 + x) * x / 4.0;
}

double psi_T(double x) {
    return (7.0 + x) * x / (2.0 * sqrt(4.0 + 3.0 * x));
}

double w_x(double* w, int i, int j) {
    return (w[(i + 1) * (y_amount + 2) + j] - w[i * (y_amount + 2) + j]) / h1;
}

double w_x_back(double* w, int i, int j) {
    return (w[i * (y_amount + 2) + j] - w[(i - 1) * (y_amount + 2) + j]) / h1;
}

double w_y(double* w, int i, int j) {
    return (w[i * (y_amount + 2) + j + 1] - w[i * (y_amount + 2) + j]) / h2;
}

double w_y_back(double* w, int i, int j) {
    return (w[i * (y_amount + 2) + j] - w[i * (y_amount + 2) + j - 1]) / h2;
}

void fill(double* A, double* w, double* q_arr, double* k_arr) {   
    double coef1 = 2.0 / h1;
    double coef2 = 2.0 / h2;
    for (int i = 1; i <= x_amount; ++i) {
        int real_i = x_offset + i - 1;
        int i1 = i * (y_amount + 2);
        for (int j = 1; j <= y_amount; ++j) {

            int real_j = y_offset + j - 1;

            if (real_i != 0 && real_i != M && real_j != 0 && real_j != N) { // внутренние точки
                double a = k_arr[4 * (i1 + j)] * w_x(w, i, j) - k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j) + 2] * w_y(w, i, j) - k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = -a / h1 - b / h2 + q_arr[i1 + j] * w[i1 + j];

            } else if (real_i == 0 && real_j == 0) { // (A1, B1)
                double a = k_arr[4 * ((i + 1) * (y_amount + 2) + j) + 1] * w_x_back(w, i + 1, j);
                double b = k_arr[4 * (i1 + j + 1) + 3] * w_y_back(w, i, j + 1);
                A[i1 + j] = - coef1 * a - coef2 * b + (q_arr[i1 + j] + coef1) * w[i1 + j];

            }
            else if (real_i == M && real_j == 0) { // (A2, B1)
                double a = k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j + 1) + 3] * w_y_back(w, i, j + 1);
                A[i1 + j] = coef1 * a - coef2 * b + (q_arr[i1 + j] + coef1) * w[i1 + j];

            }
            else if (real_i == M && real_j == N) { // (A2, B2)
                double a = k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = coef1 * a + coef2 * b + (q_arr[i1 + j] + coef1) * w[i1 + j];

            }
            else if (real_i == 0 && real_j == N) { // (A1, B2)
                double a = k_arr[4 * ((i + 1) * (y_amount + 2) + j) + 1] * w_x_back(w, i + 1, j);
                double b = k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = - coef1 * a + coef2 * b + (q_arr[i1 + j] + coef1) * w[i1 + j];

            }
            else if (real_i == 0) { // левая граница
                double a = k_arr[4 * ((i + 1) * (y_amount + 2) + j) + 1] * w_x_back(w, i + 1, j);
                double b = k_arr[4 * (i1 + j) + 2] * w_y(w, i, j) - k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = - coef1 * a + (q_arr[i1 + j] + coef1) * w[i1 + j] - b / h2;

            }
            else if (real_i == M) { // правая граница
                double a = k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j) + 2] * w_y(w, i, j) - k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = coef1 * a + (q_arr[i1 + j] + coef1) * w[i1 + j] - b / h2;

            }
            else if (real_j == 0) { // нижняя граница
                double a = k_arr[4 * (i1 + j)]  * w_x(w, i, j) - k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j + 1) + 3] * w_y_back(w, i, j + 1);
                A[i1 + j] = - coef2 * b + q_arr[i1 + j] * w[i1 + j] - a / h1;

            }
            else if (real_j == N) { // верхняя граница
                double a = k_arr[4 * (i1 + j)] * w_x(w, i, j) - k_arr[4 * (i1 + j) + 1] * w_x_back(w, i, j);
                double b = k_arr[4 * (i1 + j) + 3] * w_y_back(w, i, j);
                A[i1 + j] = coef2 * b + q_arr[i1 + j] * w[i1 + j] - a / h1;

            }


            if (debug == 1) {
                if (real_i == N - 1 - real_j) {
                    A[i * (y_amount + 2) + j] += 10.0;
                }
            }
        }
    }
}

void fill_B(double* B) {
    double coef1 = 2.0 / h1;
    double coef2 = 2.0 / h2;
    for (int i = 1; i <= x_amount; ++i) {
        int real_i = x_offset + i - 1;
        for (int j = 1; j <= y_amount; ++j) {

            int real_j = y_offset + j - 1;

            if (real_i != 0 && real_i != M && real_j != 0 && real_j != N) { // внутренние точки
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2);

            }
            else if (real_i == 0 && real_j == 0) { // (A1, B1)
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + (coef1 + coef2) * (psi_L(real_j * h2) + psi_B(real_i * h1)) / 2.0;

            }
            else if (real_i == M && real_j == 0) { // (A2, B1)
                 B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + (coef1 + coef2) * (psi_B(real_i * h1) + psi_R(real_j * h2)) / 2.0;

            }
            else if (real_i == M && real_j == N) { // (A2, B2)
                 B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + (coef1 + coef2) * (psi_R(real_j * h2) + psi_T(real_i * h1)) / 2.0;

            }
            else if (real_i == 0 && real_j == N) { // (A1, B2)
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + (coef1 + coef2) * (psi_T(real_i * h1) + psi_L(real_j * h2)) / 2.0;

            }
            else if (real_i == 0) { // левая граница
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + coef1 * psi_L(real_j * h2);

            }
            else if (real_i == M) { // правая граница
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + coef1 * psi_R(real_j * h2);

            }
            else if (real_j == 0) { // нижняя граница
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + coef2 * psi_B(real_i * h1);

            }
            else if (real_j == N) { // верхняя граница
                B[i * (y_amount + 2) + j] = F(real_i * h1, real_j * h2) + coef2 * psi_T(real_i * h1);

            }
        }
    }
}

void fill_q(double* q_arr) {
    for (int i = 1; i <= x_amount; ++i) {
        int real_i = x_offset + i - 1;
        int i1 = i * (y_amount + 2);
        for (int j = 1; j <= y_amount; ++j) {
            int real_j = y_offset + j - 1;
            q_arr[i1 + j] = q(real_i * h1, real_j * h2);
        }
    }
}
void fill_k(double* k_arr) {
    for (int i = 1; i <= x_amount; ++i) {
        int real_i = x_offset + i - 1;
        int i1 = i * (y_amount + 2);
        for (int j = 1; j <= y_amount; ++j) {
            int real_j = y_offset + j - 1;
            k_arr[4 * (i1 + j)] = k((real_i + 0.5) * h1, real_j * h2);
            k_arr[4 * (i1 + j) + 1] = k((real_i - 0.5) * h1, real_j * h2);
            k_arr[4 * (i1 + j) + 2] = k(real_i * h1, (real_j + 0.5) * h2);
            k_arr[4 * (i1 + j) + 3] = k(real_i * h1, (real_j - 0.5) * h2);
        }
    }
}

void transfer(double* w, int* dims, int* coords, MPI_Comm grid_comm) {
    int tmp;
    int other_rank[8];
    MPI_Request reqs[8];
    MPI_Status stats[8];
    double* buf_x_t_send = new double[x_amount];
    double* buf_x_t_recv = new double[x_amount];
    double* buf_x_b_send = new double[x_amount];
    double* buf_x_b_recv = new double[x_amount];
    double* buf_y_l_send = new double[y_amount];
    double* buf_y_l_recv = new double[y_amount];
    double* buf_y_r_send = new double[y_amount];
    double* buf_y_r_recv = new double[y_amount];

    for (int i = 0; i < x_amount; ++i) {
        buf_x_t_send[i] = 0.0;
        buf_x_t_recv[i] = 0.0;
        buf_x_b_send[i] = 0.0;
        buf_x_b_recv[i] = 0.0;
    }
    for (int i = 0; i < y_amount; ++i) {
        buf_y_l_send[i] = 0.0;
        buf_y_l_recv[i] = 0.0;
        buf_y_r_send[i] = 0.0;
        buf_y_r_recv[i] = 0.0;
    }

    //--------------------отправляем--------------------------------------
    int count = 0;

    // вверх
    if (dims[1] != 1 && coords[1] != dims[1] - 1)
    {
        MPI_Cart_shift(grid_comm, 1, 1, &tmp, &other_rank[count]);
        for (int i = 1; i <= x_amount; ++i)
        {
            buf_x_t_send[i - 1] = w[i * (y_amount + 2) + y_amount];
        }
        MPI_Isend(buf_x_t_send, x_amount, MPI_DOUBLE, other_rank[count], 0, grid_comm, &reqs[count]);
        ++count;
    }

    // вниз
    if (dims[1] != 1 && coords[1] != 0)
    {
        MPI_Cart_shift(grid_comm, 1, -1, &tmp, &other_rank[count]);
        for (int i = 1; i <= x_amount; ++i)
        {
            buf_x_b_send[i - 1] = w[i * (y_amount + 2) + 1];
        }
        MPI_Isend(buf_x_b_send, x_amount, MPI_DOUBLE, other_rank[count], 1, grid_comm, &reqs[count]);
        ++count;
    }

    // влево
    if (dims[0] != 1 && coords[0] != 0)
    {
        MPI_Cart_shift(grid_comm, 0, -1, &tmp, &other_rank[count]);
        for (int j = 1; j <= y_amount; ++j)
        {
            buf_y_l_send[j - 1] = w[1 * (y_amount + 2) + j];
        }
        MPI_Isend(buf_y_l_send, y_amount, MPI_DOUBLE, other_rank[count], 2, grid_comm, &reqs[count]);
        ++count;
    }

    // вправо
    if (dims[0] != 1 && coords[0] != dims[0] - 1)
    {
        MPI_Cart_shift(grid_comm, 0, 1, &tmp, &other_rank[count]);
        for (int j = 1; j <= y_amount; ++j)
        {
            buf_y_r_send[j - 1] = w[x_amount * (y_amount + 2) + j];
        }
        MPI_Isend(buf_y_r_send, y_amount, MPI_DOUBLE, other_rank[count], 3, grid_comm, &reqs[count]);
        ++count;
    }

    //---------------принимаем--------------------------------------------------------------

    // вверх
    if (dims[1] != 1 && coords[1] != 0) {
        MPI_Cart_shift(grid_comm, 1, -1, &tmp, &other_rank[count]);
        MPI_Irecv(buf_x_t_recv, x_amount, MPI_DOUBLE, other_rank[count], 0, grid_comm, &reqs[count]);
        ++count;
    }

    // вниз
    if (dims[1] != 1 && coords[1] != dims[1] - 1) {
        MPI_Cart_shift(grid_comm, 1, 1, &tmp, &other_rank[count]);
        MPI_Irecv(buf_x_b_recv, x_amount, MPI_DOUBLE, other_rank[count], 1, grid_comm, &reqs[count]);
        ++count;
    }

    // влево
    if (dims[0] != 1 && coords[0] != dims[0] - 1) {
        MPI_Cart_shift(grid_comm, 0, 1, &tmp, &other_rank[count]);
        MPI_Irecv(buf_y_l_recv, y_amount, MPI_DOUBLE, other_rank[count], 2, grid_comm, &reqs[count]);
        ++count;
    }

    // вправо
    if (dims[0] != 1 && coords[0] != 0) {
        MPI_Cart_shift(grid_comm, 0, -1, &tmp, &other_rank[count]);
        MPI_Irecv(buf_y_r_recv, y_amount, MPI_DOUBLE, other_rank[count], 3, grid_comm, &reqs[count]);
        ++count;
    }

    MPI_Waitall(count, reqs, stats);

    for (int i = 1; i <= x_amount; ++i)
    {
        w[i * (y_amount + 2) + 0] = buf_x_t_recv[i - 1];
        w[i * (y_amount + 2) + y_amount + 1] = buf_x_b_recv[i - 1];
    }
    for (int j = 1; j <= y_amount; ++j)
    {
        w[(x_amount + 1) * (y_amount + 2) + j] = buf_y_l_recv[j - 1];
        w[0 * (y_amount + 2) + j] = buf_y_r_recv[j - 1];
    }

    delete[] buf_x_t_send;
    delete[] buf_x_t_recv;
    delete[] buf_x_b_send;
    delete[] buf_x_b_recv;
    delete[] buf_y_l_send;
    delete[] buf_y_l_recv;
    delete[] buf_y_r_send;
    delete[] buf_y_r_recv;
}

void fill_rho(double* rho_arr) {
    for (int i = 1; i <= x_amount; ++i) {
        int real_i = x_offset + i - 1;
        double rho1 = (real_i == 0 || real_i == M) ? 0.5 : 1.0;
        for (int j = 1; j <= y_amount; ++j) {
            int real_j = y_offset + j - 1;
            double rho2 = (real_j == 0 || real_j == N) ? 0.5 : 1.0;
            rho_arr[i * (y_amount + 2) + j] = rho1 * rho2;
        }
    }
}

double scal(double* a, double* b, double* rho_arr) {
    double res = 0.0;
    double joint_res = 0.0;
    for (int i = 1; i <= x_amount; ++i) {
        for (int j = 1; j <= y_amount; ++j) {
            res += rho_arr[i * (y_amount + 2) + j] * a[i * (y_amount + 2) + j] * b[i * (y_amount + 2) + j];
        }
    }

    MPI_Allreduce(&res, &joint_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return h1 * h2 * joint_res;
}

double scal_full(double* a, double* b) {
    double res = 0.0;
    for (int i = 0; i <= M; ++i) {
        double rho1 = (i == 0 || i == M) ? 0.5 : 1.0;
        for (int j = 0; j <= N; ++j) {
            double rho2 = (j == 0 || j == N) ? 0.5 : 1.0;
            res += rho1 * rho2 * a[i * (N + 1) + j] * b[i * (N + 1) + j];
        }
    }

    return h1 * h2 * res;
}

double find_tau(double* r, double* Ar, double* rho_arr) {
    double res[] = { 0.0, 0.0 };
    double joint_res[] = { 0.0, 0.0 };
    for (int i = 1; i <= x_amount; ++i) {
        for (int j = 1; j <= y_amount; ++j) {
            double Ar_elem = Ar[i * (y_amount + 2) + j];
            res[0] += rho_arr[i * (y_amount + 2) + j] * r[i * (y_amount + 2) + j] * Ar_elem;
            res[1] += rho_arr[i * (y_amount + 2) + j] * Ar_elem * Ar_elem;
        }
    }

    MPI_Allreduce(res, joint_res, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return joint_res[0] / joint_res[1];
}

void sub(double* res, double* a, double* b) {
    for (int i = 0; i < (x_amount + 2) * (y_amount + 2); ++i) {
        res[i] = a[i] - b[i];
    }
}

void sub_full(double* res, double* a, double* b) {
    for (int i = 0; i < (M + 1) * (N + 1); ++i) {
        res[i] = a[i] - b[i];
    }
}

int main(int argc, char* argv[])
{
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    double needed_eps = atof(argv[3]);
    int is_res_needed = atoi(argv[4]);

    h1 = 4.0 / M;
    h2 = 3.0 / N;

    int rank, size;

    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double max_time = 0.0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[] = {0, 0};
    MPI_Dims_create(size, 2, dims);
    int tmp = dims[0];
    dims[0] = dims[1];
    dims[1] = tmp;

    int periods[] = { 0, 0 };
    MPI_Comm grid_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, true, &grid_comm);

    int coords[2];
    MPI_Cart_coords(grid_comm, rank, 2, coords);

    //cout << coords[0] << coords[1] << endl;

    int x_quotient = (M + 1)/ dims[0];
    int x_remainder = (M + 1) % dims[0];
    x_amount = x_remainder > coords[0] ? x_quotient + 1 : x_quotient;
    x_offset = x_remainder > coords[0]
        ? coords[0] * (x_quotient + 1)
        : x_remainder * (x_quotient + 1) + (coords[0] - x_remainder) * x_quotient;

    int y_quotient = (N + 1) / dims[1];
    int y_remainder = (N + 1) % dims[1];
    y_amount = y_remainder > coords[1] ? y_quotient + 1 : y_quotient;
    y_offset = y_remainder > coords[1]
        ? coords[1] * (y_quotient + 1)
        : y_remainder * (y_quotient + 1) + (coords[1] - y_remainder) * y_quotient;

    cout << rank << " x amount, offset " << x_amount << " " << x_offset << " y amount, offset " << y_amount << " " << y_offset << endl;

    int subarea_amount = (x_amount + 2) * (y_amount + 2);
    double* Aw = new double[subarea_amount];
    double* Ar = new double[subarea_amount];
    double* w = new double[subarea_amount];
    double* next_w = new double[subarea_amount];
    double* r = new double[subarea_amount];
    double* B = new double[subarea_amount];
    double* q_arr = new double[subarea_amount];
    double* k_arr = new double[subarea_amount * 4];
    double* rho_arr = new double[subarea_amount];

    for (int i = 0; i < subarea_amount; ++i) {
        Aw[i] = 0.0;
        Ar[i] = 0.0;
        w[i] = 0.0;
        next_w[i] = 0.0;
        r[i] = 0.0;
        B[i] = 0.0;
        q_arr[i] = 0.0;
        k_arr[4 * i] = 0.0;
        k_arr[4 * i + 1] = 0.0;
        k_arr[4 * i + 2] = 0.0;
        k_arr[4 * i + 3] = 0.0;
        rho_arr[i] = 0.0;
    }

    fill_B(B);
    fill_q(q_arr);
    fill_k(k_arr);
    fill_rho(rho_arr);

    double eps;
    int c = 0;
    do {
        for (int i = 1; i <= x_amount; ++i) {
            int i1 = i * (y_amount + 2);
            for (int j = 1; j <= y_amount; ++j) {
                w[i1 + j] = next_w[i1 + j];
            }
        }

        transfer(w, dims, coords, grid_comm);
        fill(Aw, w, q_arr, k_arr);
        sub(r, Aw, B);
        transfer(r, dims, coords, grid_comm);
        fill(Ar, r, q_arr, k_arr);

        /*double tmp1 = scal(Ar, r, rho_arr);
        double tmp2 = scal(Ar, Ar, rho_arr);
        double tau = tmp1 / tmp2;*/
        
        double tau = find_tau(r, Ar, rho_arr);

        //if (debug == 1) tau *= 0.8;
        for (int i = 1; i <= x_amount; ++i) {
            int i1 = i * (y_amount + 2);
            for (int j = 1; j <= y_amount; ++j) {
                next_w[i1 + j] = w[i1 + j] - tau * r[i1 + j];
            }
        }

        sub(w, next_w, w);

        //double tmp3 = scal(w, w, rho_arr);
        //eps = sqrt(tmp3);

        double max = 0.0;
        for (int i = 1; i <= x_amount; ++i) {
            for (int j = 1; j <= y_amount; ++j) {
                if (abs(w[i * (y_amount + 2) + j]) > max) max = abs(w[i * (y_amount + 2) + j]);
            }
        }
        MPI_Allreduce(&max, &eps, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        ++c;

    } while (eps > needed_eps);

    double time = MPI_Wtime() - start_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        //cout << h1 << " " << h2 << endl;
        cout << "Time = " << max_time << endl;
        cout << "Number of iterations = " << c << endl;
        cout << "Eps = " << eps << endl;
    }

    double* fin = new double[(M + 1) * (N + 1)];
    if (rank != 0) {
        MPI_Send(&x_amount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&x_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&y_amount, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&y_offset, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(next_w, subarea_amount, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // 
    }
    else {
        for (int i = 1; i <= x_amount; ++i) {
            for (int j = 1; j <= y_amount; ++j) {
                int real_i = x_offset + i - 1;
                int real_j = y_offset + j - 1;
                fin[real_i * (N + 1) + real_j] = next_w[i * (y_amount + 2) + j]; //
            }
        }
        for (int i = 1; i < size; ++i) {
            MPI_Recv(&x_amount, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&x_offset, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&y_amount, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&y_offset, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            double* buf = new double[(x_amount + 2) * (y_amount + 2)];
            MPI_Recv(buf, (x_amount + 2) * (y_amount + 2), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 1; i <= x_amount; ++i) {
                for (int j = 1; j <= y_amount; ++j) {
                    int real_i = x_offset + i - 1;
                    int real_j = y_offset + j - 1;
                    fin[real_i * (N + 1) + real_j] = buf[i * (y_amount + 2) + j];
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        double* exact_solution = new double[(M + 1) * (N + 1)];
        for (int i = 0; i <= M; ++i) {
            for (int j = 0; j <= N; ++j) {
                exact_solution[i * (N + 1) + j] = u(i * h1, j * h2);
            }
        }
        sub_full(exact_solution, fin, exact_solution);
        cout << "Norm = " << sqrt(scal_full(exact_solution, exact_solution)) << endl;

        if (is_res_needed) {
            //ofstream out("C:/Users/andre/Desktop/test.csv");
            cout << "res: " << endl;
            for (int i = 0; i <= M; ++i) {
                for (int j = 0; j <= N; ++j) {
                    //out << fin[i * (N + 1) + j];
                    //if (j != N) out << ",";
                    cout << fin[i * (N + 1) + j];
                    if (j != N) cout << ",";
                }
                //out << '\n';
                cout << '\n';
            }
        }

        delete[] exact_solution;
    }

    MPI_Finalize();

    delete[] Aw;
    delete[] Ar;
    delete[] w;
    delete[] next_w;
    delete[] r;
    delete[] B;
    delete[] q_arr;
    delete[] k_arr;
    delete[] rho_arr;
    delete[] fin;

    return 0;
}


