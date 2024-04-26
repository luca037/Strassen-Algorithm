#include <stdio.h>
#include <time.h>

#define N 512

struct m_block {
    int col;
    int row;
    int dim_b;
    int dim_m;
    int* matrix;
};

// Iterative sum.
void sum(struct m_block* a, struct m_block* b, struct m_block* c) {
    for (int ia = a->row, ib = b->row, ic = c->row; (ia < a->row + a->dim_b) && (ib < b->row + b->dim_b) && (ic < c->row + c->dim_b); ia++, ib++, ic++) {
        for (int ja = a->col, jb = b->col, jc = c->col; (ja < a->col + a->dim_b) && (jb < b->col + b->dim_b) && (jc < c->col + c->dim_b); ja++, jb++, jc++) {
            c->matrix[ic * c->dim_m + jc] = a->matrix[ia * a->dim_m + ja] + b->matrix[ib * b->dim_m + jb];
        }
    }
}

// Recursive multipy matrices C = A * B.
void rmul(struct m_block* A, struct m_block* B, struct m_block* C) {
    if (A->dim_b == 1 && B->dim_b == 1 && C->dim_b == 1) {
        C->matrix[C->row * C->dim_m + C->col] = A->matrix[A->row * A->dim_m + A->col] * B->matrix[B->row * B->dim_m + B->col];
        return;
    }

    struct m_block a11, a12, a21, a22;

    a11.matrix = A->matrix; a11.dim_m = A->dim_m;
    a11.dim_b = A->dim_b / 2;
    a11.row = 0 + A->row; a11.col = 0 + A->col;

    a12.matrix = A->matrix; a12.dim_m = A->dim_m;
    a12.dim_b = A->dim_b / 2;
    a12.row = 0 + A->row; a12.col = A->dim_b / 2 + A->col;

    a21.matrix = A->matrix; a21.dim_m = A->dim_m;
    a21.dim_b = A->dim_b / 2;
    a21.row = A->dim_b / 2 + A->row; a21.col = 0 + A->col;

    a22.matrix = A->matrix; a22.dim_m = A->dim_m;
    a22.dim_b = A->dim_b / 2;
    a22.row = A->dim_b / 2 + A->row; a22.col = A->dim_b / 2 + A->col;

    struct m_block b11, b12, b21, b22;

    b11.matrix = B->matrix; b11.dim_m = B->dim_m;
    b11.dim_b = B->dim_b / 2;
    b11.row = 0 + B->row; b11.col = 0 + B->col;

    b12.matrix = B->matrix; b12.dim_m = B->dim_m;
    b12.dim_b = B->dim_b / 2;
    b12.row = 0 + B->row; b12.col = B->dim_b / 2 + B->col;

    b21.matrix = B->matrix; b21.dim_m = B->dim_m;
    b21.dim_b = B->dim_b / 2;
    b21.row = B->dim_b / 2 + B->row; b21.col = 0 + B->col;

    b22.matrix = B->matrix; b22.dim_m = B->dim_m;
    b22.dim_b = B->dim_b / 2;
    b22.row = B->dim_b / 2 + B->row; b22.col = B->dim_b / 2 + B->col;

    struct m_block c11, c12, c21, c22;

    c11.matrix = C->matrix; c11.dim_m = C->dim_m;
    c11.dim_b = C->dim_b / 2;
    c11.row = 0 + C->row; c11.col = 0 + C->col;

    c12.matrix = C->matrix; c12.dim_m = C->dim_m;
    c12.dim_b = C->dim_b / 2;
    c12.row = 0 + C->row; c12.col = C->dim_b / 2 + C->col;

    c21.matrix = C->matrix; c21.dim_m = C->dim_m;
    c21.dim_b = C->dim_b / 2;
    c21.row = C->dim_b / 2 + C->row; c21.col = 0 + C->col;

    c22.matrix = C->matrix; c22.dim_m = C->dim_m;
    c22.dim_b = C->dim_b / 2;
    c22.row = C->dim_b / 2 + C->row; c22.col = C->dim_b / 2 + C->col;

    struct m_block tmp1;
    tmp1.dim_b = C->dim_b / 2;
    tmp1.row = 0;
    tmp1.col = 0;
    tmp1.dim_m = C->dim_b / 2;

    const int n = tmp1.dim_b;
    int tmp_m[n][n] = {};
    tmp1.matrix = (int*) tmp_m;

    rmul(&a11, &b11, &c11);
    rmul(&a12, &b21, &tmp1);
    sum(&c11, &tmp1, &c11);

    rmul(&a11, &b12, &c12);
    rmul(&a12, &b22, &tmp1);
    sum(&c12, &tmp1, &c12);

    rmul(&a21, &b11, &c21);
    rmul(&a22, &b21, &tmp1);
    sum(&c21, &tmp1, &c21);

    rmul(&a21, &b12, &c22);
    rmul(&a22, &b22, &tmp1);
    sum(&c22, &tmp1, &c22);
}

int main() {
    int A[N][N] = {};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = i + 1;
        }
    }

    int B[N][N] = {};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i][j] = i + 1;
        }
    }

    int C[N][N] = {};

    struct m_block a;
    a.dim_m = N; a.dim_b = N;
    a.col = 0; a.row = 0;
    a.matrix = (int*) A;

    //printf("Matrice A\n");
    //for (int i = 0; i < N; i++) {
    //    for (int j = 0; j < N; j++) {
    //        *(a.matrix + i * N + j) = i + 1;
    //        printf("%d ", *(a.matrix + i * N + j));
    //    }
    //    printf("\n");
    //}

    struct m_block b;
    b.dim_m = N; b.dim_b = N;
    b.col = 0; b.row = 0;
    b.matrix = (int*) B;

    //printf("Matrice B\n");
    //for (int i = 0; i < N; i++) {
    //    for (int j = 0; j < N; j++) {
    //        *(b.matrix + i * N + j) = i + 1;
    //        printf("%d ", *(b.matrix + i * N + j));
    //    }
    //    printf("\n");
    //}

    struct m_block c;
    c.dim_m = N; c.dim_b = N;
    c.col = 0; c.row = 0;
    c.matrix = (int*) C;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            *(c.matrix + i * N + j) = 0;
        }
    }
    float start_time = (float)clock() / CLOCKS_PER_SEC;
    rmul(&a, &b, &c);
    float end_time = (float)clock()/CLOCKS_PER_SEC;

    //printf("\nMatrice C\n");
    //for (int i = 0; i < N; i++) {
    //    for (int j = 0; j < N; j++) {
    //        printf("%d ", *(c.matrix + i * N + j));
    //    }
    //    printf("\n");
    //}

    printf("\n time = %f\n", end_time - start_time);

    return 0;
}
// -O1
// -O2
// -O3
// -Ofast
// -Os
