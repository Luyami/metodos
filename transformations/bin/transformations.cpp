#include <transformations.h>

mat Transformations::getHouseholderForSymmetricMat(const mat& A, int i) {
    if (A.n < 3 || A.n != A.m) throw std::invalid_argument("Must be square matrix and size >= 3!");
    if (!mat::isSymmetric(A)) throw std::invalid_argument("Must be symmetric!");

    int n = A.n;
    int k = i + 1;

    int len = n - k;
    mat x(len, 1);
    for (int r = 0; r < len; ++r) x.matrix[r][0] = A.matrix[k + r][i];

    mat e(len, 1);
    for (int r = 0; r < len; ++r) e.matrix[r][0] = 0.0;
    e.matrix[0][0] = 1.0;

    double norm_x = x.norm();
    double sign = (x.matrix[0][0] >= 0) ? 1.0 : -1.0;
    mat u = mat::add(x, mat::mult(e, sign * norm_x));
    double norm_u = u.norm();
    if (std::fabs(norm_u) < 0.00000001) throw std::runtime_error("Cannot have zero vector!");
    u = mat::div(u, norm_u);

    mat I_hat(len);
    for (int r = 0; r < len; ++r) I_hat.matrix[r][r] = 1.0;
    mat uuT = mat::mult(u, mat::transpose(u));
    mat H_hat = mat::sub(I_hat, mat::mult(uuT, 2.0));

    mat H(n);
    for (int r = 0; r < n; ++r) H.matrix[r][r] = 1.0;
    for (int r = 0; r < len; ++r) {
        for (int c = 0; c < len; ++c) {
            double val = H_hat.matrix[r][c];
            H.matrix[k + r][k + c] = val;
        }
    }

    return H;
}

mat Transformations::tridiagonalize(const mat& A) {
    if (A.n < 3 || A.n != A.m) throw std::invalid_argument("Must be square matrix and size >= 3!");
    if (!mat::isSymmetric(A)) throw std::invalid_argument("Must be symmetric!");    

    int n = A.n;
    mat B = A;

    for (int i = 0; i < n - 2; ++i) {
        mat H = getHouseholderForSymmetricMat(B, i);
        B = H * B * mat::transpose(H);
    }

    return B;
}

std::pair<mat, mat> Transformations::qrDecompose(const mat& A) {
    if (A.n != A.m) throw std::invalid_argument("Must be square matrix!");

    int n = A.n;
    mat R = A;
    mat Q = mat::identity(n);

    for (int i = 0; i < n - 1; ++i) {
        int len = n - i;
        mat x(len, 1);
        for (int r = 0; r < len; ++r) x.matrix[r][0] = R.matrix[i + r][i];

        mat e(len, 1);
        e.matrix[0][0] = 1.0;

        double norm_x = x.norm();
        double sign = (x.matrix[0][0] >= 0) ? 1.0 : -1.0;
        mat u = mat::add(x, mat::mult(e, sign * norm_x));
        double norm_u = u.norm();
        if (std::abs(norm_u) < 1e-12) continue;
        u = mat::div(u, norm_u);

        mat I_hat(len);
        for (int r = 0; r < len; ++r) I_hat.matrix[r][r] = 1.0;
        mat H_hat = mat::sub(I_hat, mat::mult(mat::mult(u, mat::transpose(u)), 2.0));

        mat H = mat::identity(n);
        for (int r = 0; r < len; ++r)
            for (int c = 0; c < len; ++c)
                H.matrix[i + r][i + c] = H_hat.matrix[r][c];

        R = H * R;
        Q = Q * mat::transpose(H);
    }

    return std::make_pair(Q, R);
}

std::pair<mat, mat> Transformations::qrEigen(const mat& A){
    if (A.n != A.m) throw std::invalid_argument("Must be square matrix!");
    if (!mat::isSymmetric(A)) throw std::invalid_argument("Matrix must be symmetric!");

    int n = A.n;
    mat Ak = A;
    mat V = mat::identity(n);

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        auto [Q, R] = qrDecompose(Ak);
        Ak = R * Q;
        V = V * Q;

        double off = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j) off += std::abs(Ak.matrix[i][j]);
        if (off < 0.00000001) break;
    }

    return std::make_pair(Ak, V);
}