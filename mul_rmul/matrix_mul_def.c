#include <stdio.h>

// Iterative multipy C = A * B.
void mul(const int n, int A[n][n], int B[n][n], int C[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][0] * B[0][j];
            for (int k = 1; k < n; k++) {
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
            }
        }
    }
}

// Iterative sum C = A + B.
void sum(const int n, int A[n][n], int B[n][n], int C[n][n]) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

// Recursive multipy matrices C = A * B.
void rmul(const int n, int A[n][n], int B[n][n], int C[n][n]) {
    if (n == 1) {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }

    int A11[n/2][n/2] = {};
    int A12[n/2][n/2] = {};
    int A21[n/2][n/2] = {};
    int A22[n/2][n/2] = {};

    int B11[n/2][n/2] = {};
    int B12[n/2][n/2] = {};
    int B21[n/2][n/2] = {};
    int B22[n/2][n/2] = {};

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + n/2];
            A21[i][j] = A[i + n/2][j];
            A22[i][j] = A[i+n/2][j+n/2];

            B11[i][j] = B[i][j];
            B12[i][j] = B[i][j + n/2];
            B21[i][j] = B[i + n/2][j];
            B22[i][j] = B[i+n/2][j+n/2];
        }
    }

    int C11[n/2][n/2] = {};
    int C12[n/2][n/2] = {};
    int C21[n/2][n/2] = {};
    int C22[n/2][n/2] = {};

    int tmp[n/2][n/2] = {};

    rmul(n/2, A11, B11, C11);
    rmul(n/2, A12, B21, tmp);
    sum(n/2, C11, tmp, C11);

    rmul(n/2, A11, B12, C12);
    rmul(n/2, A12, B22, tmp);
    sum(n/2, C12, tmp, C12);

    rmul(n/2, A21, B11, C21);
    rmul(n/2, A22, B21, tmp);
    sum(n/2, C21, tmp, C21);

    rmul(n/2, A21, B12, C22);
    rmul(n/2, A22, B22, tmp);
    sum(n/2, C22, tmp, C22);

    for (int i = 0; i < n/2; i++) {
        for (int j = 0; j < n/2; j++) {
            C[i][j] = C11[i][j];
            C[i][j + n/2] = C12[i][j];
            C[i + n/2][j] = C21[i][j];
            C[i + n/2][j + n/2] = C22[i][j];
        }
    }
}

int main() {
    int A[8][8] = {
        {1,1,1,1,1,1,1,1},
        {2,2,2,2,2,2,2,2},
        {3,3,3,3,3,3,3,3},
        {4,4,4,4,4,4,4,4},
        {5,5,5,5,5,5,5,5},
        {6,6,6,6,6,6,6,6},
        {7,7,7,7,7,7,7,7},
        {8,8,8,8,8,8,8,8},
    };

    int B[8][8] = {
        {1,1,1,1,1,1,1,1},
        {2,2,2,2,2,2,2,2},
        {3,3,3,3,3,3,3,3},
        {4,4,4,4,4,4,4,4},
        {5,5,5,5,5,5,5,5},
        {6,6,6,6,6,6,6,6},
        {7,7,7,7,7,7,7,7},
        {8,8,8,8,8,8,8,8},
    };

    int C[8][8] = {
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0}
    };

    rmul(8, A, B, C);

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%d ", C[i][j]);
        }
        printf("\n");
    }
}
