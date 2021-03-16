/**
 * Serial implementation of subspace iteration.
 * First implementation December 2020 by Emil Loevbak
 * emil.loevbak@kuleuven.be
 */

#include <iostream>
#include <cassert>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <chrono>

namespace ublas = boost::numeric::ublas;

// Fills the matrix A with random values
void randomMatrix(ublas::matrix<double, ublas::column_major> &A)
{
    for (size_t j = 0; j < A.size2(); ++j)
    {
        for (size_t i = 0; i < A.size1(); ++i)
        {
            A(i, j) = static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
        }
    }
}

// Computes the product of the matrices A and B and stores the result in C.
void matrixProduct(ublas::matrix<double, ublas::row_major> &A, ublas::matrix<double, ublas::column_major> &B, ublas::matrix<double, ublas::column_major> &C)
{
    assert(A.size2() == B.size1());
    assert(A.size1() == C.size1());
    assert(B.size2() == C.size2());

#pragma omp parallel for shared(A,B,C) schedule(static) collapse(2)
    for (size_t j = 0; j < C.size2(); ++j)
    {
        for (size_t i = 0; i < C.size1(); ++i)
        {
            C(i, j) = 0.0;
            for (size_t k = 0; k < A.size2(); ++k)
            {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }
}

// Computes the inner product of v and w.
double innerProduct(ublas::matrix_column<ublas::matrix<double, ublas::column_major>> &v, ublas::matrix_column<ublas::matrix<double, ublas::column_major>> &w)
{
    assert(v.size() == w.size());

    double result = 0.0;
    // #pragma omp for schedule(static) reduction(+: result)
    for (size_t i = 0; i < v.size(); ++i)
    {
        result += v(i) * w(i);
    }
    return result;
}

// Performs a Gram-Schmidt orthogonalization on Q. Only the diagonal elements of R are returned, as these converge to eigenvalues.
void GramSchmidt(ublas::matrix<double, ublas::column_major> &Q, ublas::vector<double> &R_diag)
{
    // We assume a full decomposition of an NxN matrix.
    assert(Q.size2() == R_diag.size());
    size_t subspaceSize = R_diag.size();

// #pragma omp parallel shared(Q, R_diag) //for older omp versions
    for (size_t i = 0; i < subspaceSize; ++i)
    {
        ublas::matrix_column<ublas::matrix<double, ublas::column_major>> Q_i(Q, i);
        R_diag(i) = sqrt(innerProduct(Q_i, Q_i));
        Q_i /= R_diag(i);
        {
#pragma omp parallel for schedule(dynamic) shared(Q_i, Q, i)
            for (size_t j = i + 1; j < subspaceSize; ++j)
            {
                ublas::matrix_column<ublas::matrix<double, ublas::column_major>> Q_j(Q, j);
                // Re-orthogonalize for stabilty
                for (size_t u = 0; u < 2; ++u)
                {
                    double R_ij = innerProduct(Q_i, Q_j);
                    noalias(Q_j) -= R_ij * Q_i;
                }
            }
        }
    }
}

// Subspace iteration method for computing eigenvalues
ublas::vector<double> subspaceIteration(ublas::matrix<double, ublas::row_major> &A, size_t subspaceSize, size_t iterations)
{
    assert(A.size1() == A.size2());

    auto Q = std::make_unique<ublas::matrix<double, ublas::column_major>>(A.size1(), subspaceSize);
    auto Q_tmp = std::make_unique<ublas::matrix<double, ublas::column_major>>(A.size1(), subspaceSize);
    ublas::vector<double> R_diag(subspaceSize);
    randomMatrix(*Q);

    for (size_t k = 0; k < iterations; ++k)
    {
        matrixProduct(A, *Q, *Q_tmp);
        std::swap(Q, Q_tmp);
        GramSchmidt(*Q, R_diag);
    }

    return R_diag;
}

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cout << "Expected exactly three command line arguments, the matrix dimension, the number of eigenvalues and the number of iterations!" << std::endl;
        return EXIT_FAILURE;
    }

    size_t N = strtoul(argv[1], NULL, 0);
    size_t subspaceSize = strtoul(argv[2], NULL, 0);
    size_t iterations = strtoul(argv[3], NULL, 0);
    std::srand(123456789);


    // A is a diagonal matrix, so we know the eigenvalues.
    ublas::matrix<double, ublas::row_major> A(N, N);
    for (size_t i = 0; i < N; ++i)
    {
        A(i, i) = i + 1;
    }

    const int skips=1; 
    int avg_time = 0;
    std::cout << "Performing 10 computations with timings..." << std::endl;
    for (size_t i = 0; i < 10; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();
        ublas::vector<double> eigenvalues = subspaceIteration(A, subspaceSize, iterations);
        auto end = std::chrono::high_resolution_clock::now();
        auto timing = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Computed eigenvalues: " << eigenvalues << std::endl;
        std::cout << "Time needed: " << timing.count() << " milliseconds" << std::endl
                  << std::endl;
        if (i>skips-1) avg_time += timing.count();
    }
    avg_time /= 10-skips;
    std::cout << "Average time: " << avg_time << " milliseconds" << std::endl
              << std::endl;
    return EXIT_SUCCESS;
}
