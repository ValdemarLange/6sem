#include <iostream>
#include "nr3.h"
#include "utilities.h"
#include <vector>
#include <cmath>

// Computes the dot product of two VecDoub vectors
double dotProduct(const VecDoub &a, const VecDoub &b)
{
    if (a.size() != b.size())
    {
        throw std::runtime_error("Vector sizes do not match for dot product.");
    }

    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

// Function to compute the Euclidean norm of a vector
double VectorLength(const VecDoub &in)
{
    double sum = 0;
    for (int i = 0; i < in.size(); i++)
    {
        sum += in[i] * in[i];
    }
    return std::sqrt(sum); // Added return statement
}

// Gram-Schmidt orthonormalization
std::vector<VecDoub> gramSchmidt(const std::vector<VecDoub> &x)
{
    int k = x.size();
    std::vector<VecDoub> e(k, VecDoub(x[0].size(), 0.0)); // Initialize all vectors

    // Normalize the first vector
    double norm = VectorLength(x[0]);
    for (int i = 0; i < x[0].size(); i++)
    {
        e[0][i] = x[0][i] / norm;
    }

    // Process the remaining vectors
    for (int i = 1; i < k; i++)
    {
        VecDoub sum(x[i].size(), 0.0);

        for (int j = 0; j < i; j++)
        { // Orthogonalize against previous vectors
            double dot_product = dotProduct(x[i], e[j]);
            for (int m = 0; m < x[i].size(); m++)
            {
                sum[m] += dot_product * e[j][m];
            }
        }

        // Compute orthogonalized vector
        for (int m = 0; m < x[i].size(); m++)
        {
            e[i][m] = x[i][m] - sum[m];
        }

        // Normalize
        double norm = VectorLength(e[i]);
        if (norm > 1e-10) { // Avoid division by zero
            for (int m = 0; m < x[i].size(); m++) {
                e[i][m] /= norm;
            }
        }
    }

    return e;
}

int main() {
    // Define 3 linearly independent vectors using VecDoub
    VecDoub v1(5);
    v1[0] = 2.0; v1[1] = 8.0; v1[2] = 4.0; v1[3] =  2.0; v1[4] = 1.0;

    VecDoub v2(5);
    v2[0] = 1.0; v2[1] = 1.0; v2[2] = 5.0; v2[3] = 7.0; v2[4] = 8.0;

    VecDoub v3(5);
    v3[0] = 4.0; v3[1] = -5.0; v3[2] = 1.0; v3[3] = -4.0; v3[4] = 3.0;

    // Store them in a vector
    std::vector<VecDoub> x = {v1, v2, v3};

    // Perform Gram-Schmidt orthonormalization
    std::vector<VecDoub> orthonormalBasis = gramSchmidt(x);

    // Print the orthonormal vectors
    std::cout << "Orthonormal basis:\n";
    for (size_t i = 0; i < orthonormalBasis.size(); i++) {
        std::cout << "e" << i + 1 << " = (";
        for (size_t j = 0; j < orthonormalBasis[i].size(); j++) {
            std::cout << orthonormalBasis[i][j];
            if (j < orthonormalBasis[i].size() - 1) std::cout << ", ";
        }
        std::cout << ")\n";
    }

    return 0;
}


