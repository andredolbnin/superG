#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <string>

using namespace std;

double F(long double x, long double y, long double z) {
    return y <= x && z <= x * y ? x * y * y * z * z * z : 0;
}

/*void weighting_func(long double x, long double y, long double z, int len) {
    for (int i = 0; i < len; i++) {
        long double res = sin(x) * cos(y) * z;
    }
    return;
}*/

int main(int argc, char* argv[])
{
    long double eps = atof(argv[1]);

    int rank, size;

    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const long double exact_solution = 1.0L / 364;
    const int root_rank = 0;

    long double error = 1.0L;
    long double joint_sum = 0.0L;
    long double sum = 0.0L;
    double max_time = 0.0;

    int dots_count = 31 * 15 * 250;
    int sub_dots_count = dots_count / (size - 1);
    long double* sub_dots = new long double[3 * sub_dots_count];
    long double* dots = new long double[3 * dots_count];

    long long int k = 0;
    int isWorking = 1;

    if (rank == root_rank) {
        while (error > eps) {
            if (error > 1.0) {
                cerr << "Info: error > 1.0 !" << endl;
                MPI_Abort(MPI_COMM_WORLD, 0);
                return 1;
            }

            for (int i = 0; i < 3 * dots_count; i++) {
                dots[i] = (long double)rand() / RAND_MAX;
            }

            for (int i = 0; i < size - 1; i++) {
                MPI_Send(&dots[i * 3 * sub_dots_count], 3 * sub_dots_count, MPI_LONG_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
            }

            MPI_Reduce(&sum, &joint_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
            k++;
            error = abs(joint_sum / (dots_count * k) - exact_solution);

            isWorking = error > eps ? 1 : 0;
            MPI_Bcast(&isWorking, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
        }
    }
    else {
        do {
            MPI_Recv(sub_dots, 3 * sub_dots_count, MPI_LONG_DOUBLE, root_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (int i = 0; i < 3 * sub_dots_count; i += 3) {
                sum += F(sub_dots[i], sub_dots[i + 1], sub_dots[i + 2]);
                //weighting_func(sub_dots[i], sub_dots[i + 1], sub_dots[i + 2], 1000);
            }

            MPI_Reduce(&sum, &joint_sum, 1, MPI_LONG_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
            MPI_Bcast(&isWorking, 1, MPI_INT, root_rank, MPI_COMM_WORLD);
        
        } while (isWorking == 1);
    }

    double time = MPI_Wtime() - start_time;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, root_rank, MPI_COMM_WORLD);

    if (rank == root_rank) {
        cout << "I = " << joint_sum / (dots_count * k) << endl;
        cout << "Error = " << error << endl;
        cout << "n = " << dots_count * k << endl;
        cout << "Time = " << max_time << endl;
    }

    MPI_Finalize();

    delete[] dots;
    delete[] sub_dots;

    return 0;
}
