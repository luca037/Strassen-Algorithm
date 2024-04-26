#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 4096

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

// Iterative sub.
void sub(struct m_block* a, struct m_block* b, struct m_block* c) {
    for (int ia = a->row, ib = b->row, ic = c->row; (ia < a->row + a->dim_b) && (ib < b->row + b->dim_b) && (ic < c->row + c->dim_b); ia++, ib++, ic++) {
        for (int ja = a->col, jb = b->col, jc = c->col; (ja < a->col + a->dim_b) && (jb < b->col + b->dim_b) && (jc < c->col + c->dim_b); ja++, jb++, jc++) {
            c->matrix[ic * c->dim_m + jc] = a->matrix[ia * a->dim_m + ja] - b->matrix[ib * b->dim_m + jb];
        }
    }
}

// Recursive multipy matrices C = A * B.
void smul(struct m_block* A, struct m_block* B, struct m_block* C) {
    if (A->dim_b == 1 && B->dim_b == 1 && C->dim_b == 1) {
        C->matrix[C->row * C->dim_m + C->col] = A->matrix[A->row * A->dim_m + A->col] * B->matrix[B->row * B->dim_m + B->col];
        return;
    }

    // Blocks.
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

    // Strassen's blocks.
    struct m_block a1, a2, a3, a4, a7;

    a1.dim_b = a1.dim_m = a12.dim_b;
    a1.col = a1.row = 0;
    a1.matrix = (int*) malloc(a12.dim_b * a12.dim_b * sizeof(int));

    a2.dim_b = a2.dim_m = a12.dim_b;
    a2.col = a2.row = 0;
    a2.matrix = (int*) malloc(a12.dim_b * a12.dim_b * sizeof(int));

    a3.dim_b = a3.dim_m = a12.dim_b;
    a3.col = a3.row = 0;
    a3.matrix = (int*) malloc(a12.dim_b * a12.dim_b * sizeof(int));

    a4.dim_b = a4.dim_m = a12.dim_b;
    a4.col = a4.row = 0;
    a4.matrix = (int*) malloc(a12.dim_b * a12.dim_b * sizeof(int));

    a7.dim_b = a7.dim_m = a12.dim_b;
    a7.col = a7.row = 0;
    a7.matrix = (int*) malloc(a12.dim_b * a12.dim_b * sizeof(int));

    sub(&a12, &a22, &a1);
    sum(&a11, &a22, &a2);
    sub(&a11, &a21, &a3);
    sum(&a11, &a12, &a4);
    sum(&a21, &a22, &a7);

    struct m_block b1, b2, b3, b5, b6;

    b1.dim_b = b1.dim_m = b12.dim_b;
    b1.col = b1.row = 0;
    b1.matrix = (int*) malloc(b12.dim_b * b12.dim_b * sizeof(int));

    b2.dim_b = b2.dim_m = b12.dim_b;
    b2.col = b2.row = 0;
    b2.matrix = (int*) malloc(b12.dim_b * b12.dim_b * sizeof(int));

    b3.dim_b = b3.dim_m = b12.dim_b;
    b3.col = b3.row = 0;
    b3.matrix = (int*) malloc(b12.dim_b * b12.dim_b * sizeof(int));

    b5.dim_b = b5.dim_m = b12.dim_b;
    b5.col = b5.row = 0;
    b5.matrix = (int*) malloc(b12.dim_b * b12.dim_b * sizeof(int));

    b6.dim_b = b6.dim_m = b12.dim_b;
    b6.col = b6.row = 0;
    b6.matrix = (int*) malloc(b12.dim_b * b12.dim_b * sizeof(int));

    sum(&b21, &b22, &b1);
    sum(&b11, &b22, &b2);
    sub(&b11, &b12, &b3);
    sub(&b12, &b22, &b5);
    sub(&b21, &b11, &b6);

    struct m_block p1, p2, p3, p4, p5, p6, p7;
    p1.dim_b = p1.dim_m = a1.dim_b; p1.col = p1.row = 0;
    p2.dim_b = p2.dim_m = a1.dim_b; p2.col = p2.row = 0;
    p3.dim_b = p3.dim_m = a1.dim_b; p3.col = p3.row = 0;
    p4.dim_b = p4.dim_m = a1.dim_b; p4.col = p4.row = 0;
    p5.dim_b = p5.dim_m = a1.dim_b; p5.col = p5.row = 0;
    p6.dim_b = p6.dim_m = a1.dim_b; p6.col = p6.row = 0;
    p7.dim_b = p7.dim_m = a1.dim_b; p7.col = p7.row = 0;

    p1.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p2.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p3.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p4.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p5.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p6.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));
    p7.matrix = (int*) malloc((a1.dim_b * a1.dim_b) * sizeof(int));

    smul(&a1,  &b1,  &p1);
    smul(&a2,  &b2,  &p2);
    smul(&a3,  &b3,  &p3);
    smul(&a4,  &b22, &p4);
    smul(&a11, &b5,  &p5);
    smul(&a22, &b6,  &p6);
    smul(&a7,  &b11, &p7);

    // C11
    sum(&p1, &p2, &c11);
    sub(&c11, &p4, &c11);
    sum(&c11, &p6, &c11);

    // C12
    sum(&p4, &p5, &c12);

    // C21
    sum(&p6, &p7, &c21);

    // C22
    sub(&p2, &p3, &c22);
    sum(&c22, &p5, &c22);
    sub(&c22, &p7, &c22);

    free(a1.matrix);
    free(a2.matrix);
    free(a3.matrix);
    free(a4.matrix);
    free(a7.matrix);

    free(b1.matrix);
    free(b2.matrix);
    free(b3.matrix);
    free(b5.matrix);
    free(b6.matrix);

    free(p1.matrix);
    free(p2.matrix);
    free(p3.matrix);
    free(p4.matrix);
    free(p5.matrix);
    free(p6.matrix);
    free(p7.matrix);
}

int main() {
    struct m_block a;
    a.dim_m = N; a.dim_b = N;
    a.col = 0; a.row = 0;
    a.matrix = (int*) malloc((N * N) * sizeof(int));

    //printf("Matrice A\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            a.matrix[i * N + j] = 1;
            //printf("%d ", a.matrix[i * N + j]);
        }
        //printf("\n");
    }

    struct m_block b;
    b.dim_m = N; b.dim_b = N;
    b.col = 0; b.row = 0;
    b.matrix = (int*) malloc((N * N) *  sizeof(int));

    //printf("Matrice B\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            b.matrix[i * N + j] = 1;
            //printf("%d ", b.matrix[i][j]);
        }
        //printf("\n");
    }

    struct m_block c;
    c.dim_m = N; c.dim_b = N;
    c.col = 0; c.row = 0;
    c.matrix = (int*) malloc((N * N) * sizeof(int));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            c.matrix[i * N + j] = 0;
        }
    }

    float start_time = (float)clock() / CLOCKS_PER_SEC;
    smul(&a, &b, &c);
    float end_time = (float)clock() / CLOCKS_PER_SEC;

    smul(&a, &b, &c);

    printf("\ntime = %f\n", end_time - start_time);

    //printf("\nMatrice C\n");
    //for (int i = 0; i < N; i++) {
    //    for (int j = 0; j < N; j++) {
    //        printf("%d ", c.matrix[i * N + j]);
    //    }
    //    printf("\n");
    //}

    // Dealloc
    free(a.matrix);

    free(b.matrix);

    free(c.matrix);
}
