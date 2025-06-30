#include <mat.h>

//Util funcs
static void luDecompose(const mat& m, mat& L, mat& U, std::vector<uint32_t>& pivot) {
    if (m.n != m.m) throw std::invalid_argument("LU decomposition requires a square matrix.");
    uint32_t n = m.n;
    L = mat(n);
    U = mat(n);
    pivot.resize(n);
    for (uint32_t i = 0; i < n; ++i) pivot[i] = i;

    mat A(m);

    for (uint32_t k = 0; k < n; ++k) {
        double maxVal = fabs(A.matrix[k][k]);
        uint32_t maxRow = k;
        for (uint32_t i = k + 1; i < n; ++i) {
            if (fabs(A.matrix[i][k]) > maxVal) {
                maxVal = fabs(A.matrix[i][k]);
                maxRow = i;
            }
        }
        if (maxVal < 1e-12) throw std::runtime_error("Matrix is singular.");

        std::swap(pivot[k], pivot[maxRow]);
        for (uint32_t j = 0; j < n; ++j)
            std::swap(A.matrix[k][j], A.matrix[maxRow][j]);

        for (uint32_t i = k + 1; i < n; ++i) {
            A.matrix[i][k] /= A.matrix[k][k];
            for (uint32_t j = k + 1; j < n; ++j)
                A.matrix[i][j] -= A.matrix[i][k] * A.matrix[k][j];
        }
    }

    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < n; ++j) {
            if (i > j) L.matrix[i][j] = A.matrix[i][j];
            else if (i == j) {
                L.matrix[i][j] = 1.0;
                U.matrix[i][j] = A.matrix[i][j];
            } else {
                U.matrix[i][j] = A.matrix[i][j];
            }
        }
    }
}
//

mat::mat(uint32_t size) : n(size), m(size) {
    matrix = new double*[n];
    for (uint32_t i = 0; i < n; ++i) {
        matrix[i] = new double[m]{};
    }
}

mat::mat(uint32_t n, uint32_t m) : n(n), m(m) {
    matrix = new double*[n];
    for (uint32_t i = 0; i < n; ++i) {
        matrix[i] = new double[m]{};
    }
}

mat::mat(const mat& other) : n(other.n), m(other.m) {
    matrix = new double*[n];
    for (uint32_t i = 0; i < n; ++i) {
        matrix[i] = new double[m];
        for (uint32_t j = 0; j < m; ++j) {
            matrix[i][j] = other.matrix[i][j];
        }
    }
}

mat::~mat() {
    for (uint32_t i = 0; i < n; ++i)
        delete[] matrix[i];
    delete[] matrix;
}

void mat::build(vector<double> values) {
    if (values.size() != n * m)
        throw invalid_argument("Incorrect number of elements to build matrix");

    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < m; ++j)
            matrix[i][j] = values[i * m + j];
}

void mat::print() const {
    for (uint32_t i = 0; i < n; ++i) {
        for (uint32_t j = 0; j < m; ++j) {
            double value = matrix[i][j];
            if (fabs(value) < 0.00000001) value = 0;
            cout << setw(8) << value << " ";
        }
        cout << '\n';
    }
}

double mat::norm() const{
    double sum = 0.0f;
    for (uint32_t i = 0; i < n; ++i)
        for (uint32_t j = 0; j < m; ++j)
            sum += matrix[i][j] * matrix[i][j];
    return std::sqrt(sum);
}

mat mat::add(const mat& m1, const mat& m2) {
    if (m1.n != m2.n || m1.m != m2.m)
        throw invalid_argument("Matrix size mismatch");

    mat result(m1.n, m1.m);
    for (uint32_t i = 0; i < m1.n; ++i)
        for (uint32_t j = 0; j < m1.m; ++j)
            result.matrix[i][j] = m1.matrix[i][j] + m2.matrix[i][j];
    return result;
}

mat mat::sub(const mat& m1, const mat& m2) {
    if (m1.n != m2.n || m1.m != m2.m)
        throw invalid_argument("Matrix size mismatch");

    mat result(m1.n, m1.m);
    for (uint32_t i = 0; i < m1.n; ++i)
        for (uint32_t j = 0; j < m1.m; ++j)
            result.matrix[i][j] = m1.matrix[i][j] - m2.matrix[i][j];
    return result;
}

mat mat::mult(const mat& m1, const mat& m2) {
    if (m1.m != m2.n)
        throw invalid_argument("Matrix dimension mismatch for multiplication");

    mat result(m1.n, m2.m);
    for (uint32_t i = 0; i < m1.n; ++i)
        for (uint32_t j = 0; j < m2.m; ++j)
            for (uint32_t k = 0; k < m1.m; ++k)
                result.matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
    return result;
}

mat mat::mult(const mat& m, double scalar) {
    mat result(m.n, m.m);
    for (uint32_t i = 0; i < m.n; ++i)
        for (uint32_t j = 0; j < m.m; ++j)
            result.matrix[i][j] = m.matrix[i][j] * scalar;
    return result;
}

mat mat::div(const mat& m, double scalar) {
    if (scalar == 0.0f)
        throw invalid_argument("Division by zero");

    mat result(m.n, m.m);
    for (uint32_t i = 0; i < m.n; ++i)
        for (uint32_t j = 0; j < m.m; ++j)
            result.matrix[i][j] = m.matrix[i][j] / scalar;
    return result;
}

mat mat::normalize(const mat& m) {
    double norm_val = m.norm();
  
    if (norm_val == 0.0f) {
        throw std::runtime_error("Cannot normalize a zero-norm matrix.");
    }

    mat result(m.n, m.m);
    for (uint32_t i = 0; i < m.n; ++i) {
        for (uint32_t j = 0; j < m.m; ++j) {
            result.matrix[i][j] = m.matrix[i][j] / norm_val;
        }
    }
    return result;
}

mat mat::transpose(const mat& m) {
    mat result(m.m, m.n);
    for (uint32_t i = 0; i < m.n; ++i)
        for (uint32_t j = 0; j < m.m; ++j)
            result.matrix[j][i] = m.matrix[i][j];
    return result;
}

double mat::det(const mat& m) {
    if (m.n != m.m) throw std::invalid_argument("Matrix must be square.");
    mat L(m.n), U(m.n);
    std::vector<uint32_t> pivot;
    luDecompose(m, L, U, pivot);

    double det = 1.0;
    for (uint32_t i = 0; i < m.n; ++i)
        det *= U.matrix[i][i];

    uint32_t swaps = 0;
    for (uint32_t i = 0; i < pivot.size(); ++i)
        if (pivot[i] != i) ++swaps;

    if (swaps % 2 != 0)
        det = -det;

    return det;
}

mat mat::inverse(const mat& m) {
    if (m.n != m.m) throw std::invalid_argument("Matrix must be square.");
    uint32_t n = m.n;
    mat L(n), U(n);
    std::vector<uint32_t> pivot;
    luDecompose(m, L, U, pivot);

    mat inv(n);

    for (uint32_t col = 0; col < n; ++col) {
        std::vector<double> y(n, 0.0);
        for (uint32_t i = 0; i < n; ++i) {
            y[i] = (pivot[i] == col) ? 1.0 : 0.0;
            for (uint32_t j = 0; j < i; ++j)
                y[i] -= L.matrix[i][j] * y[j];
        }

        std::vector<double> x(n, 0.0);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (uint32_t j = i + 1; j < n; ++j)
                x[i] -= U.matrix[i][j] * x[j];
            x[i] /= U.matrix[i][i];
        }

        for (uint32_t i = 0; i < n; ++i)
            inv.matrix[i][col] = x[i];
    }

    return inv;
}

mat mat::identity(uint32_t size){
    mat identity(size, size);

    for (int i = 0; i < size; ++i){
        identity.matrix[i][i] = 1.0;
    }

    return identity;
}

bool mat::isSymmetric(const mat& m) {
    if (m.n != m.m) return false;
    for (uint32_t i = 0; i < m.n; ++i) {
        for (uint32_t j = i + 1; j < m.m; ++j) {
            if (std::abs(m.matrix[i][j] - m.matrix[j][i]) > 0.000001f)
                return false;
        }
    }
    return true;
}

mat mat::operator+(const mat& other) const {
    return add(*this, other);
}

mat mat::operator-(const mat& other) const {
    return sub(*this, other);
}

mat mat::operator*(const mat& other) const {
    return mult(*this, other);
}

mat mat::operator*(double scalar) const {
    return mult(*this, scalar);
}

mat mat::operator/(double scalar) const {
    return div(*this, scalar);
}

mat& mat::operator=(const mat& other) {
    if (this == &other)
        return *this;

    for (uint32_t i = 0; i < n; ++i)
        delete[] matrix[i];
    delete[] matrix;

    n = other.n;
    m = other.m;

    matrix = new double*[n];
    for (uint32_t i = 0; i < n; ++i) {
        matrix[i] = new double[m];
        for (uint32_t j = 0; j < m; ++j) {
            matrix[i][j] = other.matrix[i][j];
        }
    }

    return *this;
}