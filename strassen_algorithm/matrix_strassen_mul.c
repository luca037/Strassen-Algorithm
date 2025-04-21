#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// Sub-matrix block.
struct m_block {
    int col;
    int row;
    int dim_b;
    int dim_m;
    int* matrix;
};

// Allocate a matrix.
struct m_block make_alloc_block(int dim_b) {
    struct m_block block;
    block.row = 0;
    block.col = 0;
    block.dim_b = block.dim_m = dim_b;
    block.matrix = (int*) malloc(dim_b * dim_b * sizeof(int));
    return block;
}

// Crate a sub-matrix block.
struct m_block make_block(struct m_block* parent, int row_offset, int col_offset, int dim_b) {
    struct m_block block;
    block.matrix = parent->matrix;
    block.dim_m = parent->dim_m;
    block.dim_b = dim_b;
    block.row = parent->row + row_offset;
    block.col = parent->col + col_offset;
    return block;
}

// Iterative sum.
void sum(struct m_block* restrict a, struct m_block* restrict b, struct m_block* restrict c) {
    int* pa = &a->matrix[a->row * a->dim_m + a->col];
    int* pb = &b->matrix[b->row * b->dim_m + b->col];
    int* pc = &c->matrix[c->row * c->dim_m + c->col];
    for (int i = 0; i < a->dim_b; ++i) {
        for (int j = 0; j < a->dim_b; ++j) {
            *pc = *pa + *pb;
            ++pa; ++pb; ++pc;
        }
        pa += a->dim_m - a->dim_b;
        pb += b->dim_m - b->dim_b;
        pc += c->dim_m - c->dim_b;
    }
}

// Iterative sub.
void sub(struct m_block* restrict a, struct m_block* restrict b, struct m_block* restrict c) {
    int* pa = &a->matrix[a->row * a->dim_m + a->col];
    int* pb = &b->matrix[b->row * b->dim_m + b->col];
    int* pc = &c->matrix[c->row * c->dim_m + c->col];
    for (int i = 0; i < a->dim_b; ++i) {
        for (int j = 0; j < a->dim_b; ++j) {
            *pc = *pa - *pb;
            ++pa; ++pb; ++pc;
        }
        pa += a->dim_m - a->dim_b;
        pb += b->dim_m - b->dim_b;
        pc += c->dim_m - c->dim_b;
    }
}

// Recursive multipy matrices C = A * B.
void smul(struct m_block* restrict A, struct m_block* restrict B, struct m_block* restrict C) {
    if (A->dim_b == 1) {
        int a_val = A->matrix[A->row * A->dim_m + A->col];
        int b_val = B->matrix[B->row * B->dim_m + B->col];
        C->matrix[C->row * C->dim_m + C->col] = a_val * b_val;
        return;
    }

    // Blocks size.
    int half = A->dim_b / 2;

    // Blocks of matrix A.
    struct m_block a11 = make_block(A, 0,      0,     half);
    struct m_block a12 = make_block(A, 0,      half,  half);
    struct m_block a21 = make_block(A, half,   0,     half);
    struct m_block a22 = make_block(A, half,   half,  half);

    // Blocks of matrix B.
    struct m_block b11 = make_block(B, 0,      0,    half);
    struct m_block b12 = make_block(B, 0,      half, half);
    struct m_block b21 = make_block(B, half,   0,    half);
    struct m_block b22 = make_block(B, half,   half, half);

    // Blocks of matrix C.
    struct m_block c11 = make_block(C, 0,      0,    half);
    struct m_block c12 = make_block(C, 0,      half, half);
    struct m_block c21 = make_block(C, half,   0,    half);
    struct m_block c22 = make_block(C, half,   half, half);

    // Tmp blocks to compute the Strassen's operations.
    struct m_block a = make_alloc_block(half);
    struct m_block b = make_alloc_block(half);
    struct m_block p = make_alloc_block(half);

    // Clean up C.
    memset(C->matrix, 0, C->dim_m * C->dim_m * sizeof(int));

    // Compute p1 ad sum it to c11.
    sub(&a12, &a22, &a); // a = a12 - a22
    sum(&b21, &b22, &b); // b = b21 + b22
    smul(&a,  &b,  &p);  // p = a * b
    sum(&c11, &p, &c11); // c11 += p

    // Compute p2 and sum it to c11 and c22.
    sum(&a11, &a22, &a); // a = a11 + a22
    sum(&b11, &b22, &b); // b = b11 + b22
    smul(&a, &b, &p);    // p = a * b
    sum(&c11, &p, &c11); // c11 += p
    sum(&c22, &p, &c22); // c22 += p

    // Compute p3 and sum it to c22.
    sub(&a11, &a21, &a); // a = a11 - a21
    sum(&b11, &b12, &b); // b = b11 + b12
    smul(&a, &b, &p);    // p = a * b
    sub(&c22, &p, &c22); // c22 -= p

    // Compute p4 and sum it to c11 and c12.
    sum(&a11, &a12, &a); // a = a11 + a12
    smul(&a, &b22, &p);  // p = a * b22
    sub(&c11, &p, &c11); // c11 -= p
    sum(&c12, &p, &c12); // c12 += p

    // Compute p5 and sum it to c12 and c22.
    sub(&b12, &b22, &b); // b = b12 - b22
    smul(&a11, &b, &p);  // p = a11 * b
    sum(&c12, &p, &c12); // c12 += p
    sum(&c22, &p, &c22); // c22 += p

    // Compute p6 and sum it to c11 and c21.
    sub(&b21, &b11, &b); // b = b21 - b11
    smul(&a22, &b, &p);  // p = a22 * b
    sum(&c11, &p, &c11); // c11 += p
    sum(&c21, &p, &c21); // c21 += p

    // Compute p7 and add it to c21 and c22.
    sum(&a21, &a22, &a); // a = a21 + a22
    smul(&a, &b11, &p);  // p = a * b11
    sum(&c21, &p, &c21); // c21 += p
    sub(&c22, &p, &c22); // c22 -= p

    free(a.matrix);
    free(b.matrix);
    free(p.matrix);
}

// Execute Strassen algorithm, print CPU time.
void exec_strassen(const int N) {
    struct m_block a = make_alloc_block(N);
    struct m_block b = make_alloc_block(N);

    // Init matrices A, B with ones.
    for (int i = 0; i < N * N; i++) {
        a.matrix[i] = 1;
        b.matrix[i] = 1;
    }

    struct m_block c = make_alloc_block(N);

    // Matrix C is cleaned up in smul.

    float start_time = (float)clock() / CLOCKS_PER_SEC;
    smul(&a, &b, &c);
    float end_time = (float)clock() / CLOCKS_PER_SEC;

    printf("%d, %f\n", N, end_time - start_time);

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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("You need to specify the size (N).\n");
        return 1;
    }

    exec_strassen(atoi(argv[1]));

    return 0;
}
