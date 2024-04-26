#include <stdio.h>
#include <stdlib.h>

struct m_block {
    int col;
    int row;
    int dim_b;
    int dim_m;
    int** matrix;
};

struct m_block2 {
    int col;
    int row;
    int dim_b;
    int dim_m;
    int* matrix;
};
// Iterative sum.
void sum2(struct m_block2* a, struct m_block2* b, struct m_block2* c) {
    for (int ia = a->row, ib = b->row, ic = c->row; (ia < a->row + a->dim_b) && (ib < b->row + b->dim_b) && (ic < c->row + c->dim_b); ia++, ib++, ic++) {
        for (int ja = a->col, jb = b->col, jc = c->col; (ja < a->col + a->dim_b) && (jb < b->col + b->dim_b) && (jc < c->col + c->dim_b); ja++, jb++, jc++) {
            c->matrix[ic * c->dim_m + jc] = a->matrix[ia * a->dim_m + ja] + b->matrix[ib * b->dim_m + jb];
        }
    }
}

// Iterative sum.
void sum(struct m_block* a, struct m_block* b, struct m_block* c) {
    for (int ia = a->row, ib = b->row, ic = c->row; (ia < a->row + a->dim_b) && (ib < b->row + b->dim_b) && (ic < c->row + c->dim_b); ia++, ib++, ic++) {
        for (int ja = a->col, jb = b->col, jc = c->col; (ja < a->col + a->dim_b) && (jb < b->col + b->dim_b) && (jc < c->col + c->dim_b); ja++, jb++, jc++) {
            c->matrix[ic][jc] = a->matrix[ia][ja] + b->matrix[ib][jb];
        }
    }
}

void test_sum() {
    int** a;
    a = (int**) malloc(8 * sizeof(int*));
    for (int i = 0; i < 8; i++){
        a[i] = (int*) malloc(8 * sizeof(int));
    }

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            a[i][j] = i + 1;
            printf("%d ", a[i][j]);
        }
        printf("\n");
    }

    int** c;
    c = (int**) malloc(8 * sizeof(int*));
    for (int i = 0; i < 8; i++){
        c[i] = (int*) malloc(8 * sizeof(int));
    }

    printf("\n");
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            c[i][j] = 0;
            printf("%d ", c[i][j]);
        }
        printf("\n");
    }

    struct m_block a11, a12, a21, a22;

    a11.matrix = a12.matrix = a21.matrix = a22.matrix = a;

    a11.dim_m = a12.dim_m = a21.dim_m = a22.dim_m = 8;
    a11.dim_b = a12.dim_b = a21.dim_b = a22.dim_b = 4;

    a11.col = a11.row = 0;
    a12.col = 4; a12.row = 0;
    a21.col = 0; a21.row = 4;
    a22.col = 4; a22.row = 4;

    struct m_block c11;
    c11.matrix = c;
    c11.dim_m = 8;
    c11.dim_b = 4;
    c11.row = 0; c11.row = 0;

    printf("\nPrimo blocco:\n");
    for (int i = a11.row; i < a11.row + a11.dim_b; i++) {
        for (int j = a11.col; j < a11.col + a11.dim_b; j++) {
            printf("%d ", a11.matrix[i][j]);
        }
        printf("\n");
    }

    printf("\nSecondo blocco:\n");
    for (int i = a12.row; i < a12.row + a12.dim_b; i++) {
        for (int j = a12.col; j < a12.col + a12.dim_b; j++) {
            printf("%d ", a12.matrix[i][j]);
        }
        printf("\n");
    }

    printf("\nSomma:\n");
    sum(&a11, &a12, &c11);
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%d ", c[i][j]);
        }
        printf("\n");
    }

    for (int i = 0; i < 8; i++){
        free(a[i]);
    }
    free(a);

    for (int i = 0; i < 8; i++){
        free(c[i]);
    }
    free(c);
}

void test_sum_2() {

}

int main(void) {
    int A[8][8] = {
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
    };

    int B[8][8] = {
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1},
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

    struct m_block2 a11, a12, a21, a22;

    a11.matrix = a12.matrix = a21.matrix = a22.matrix = (int*) A;

    a11.dim_m = a12.dim_m = a21.dim_m = a22.dim_m = 8;
    a11.dim_b = a12.dim_b = a21.dim_b = a22.dim_b = 4;

    a11.col = a11.row = 0;
    a12.col = 4; a12.row = 0;
    a21.col = 0; a21.row = 4;
    a22.col = 4; a22.row = 4;

    struct m_block2 c11;
    c11.matrix = (int*) C;
    c11.dim_m = 8;
    c11.dim_b = 4;
    c11.row = 0; c11.col = 0;

    printf("\nPrimo blocco:\n");
    for (int i = a11.row; i < a11.row + a11.dim_b; i++) {
        for (int j = a11.col; j < a11.col + a11.dim_b; j++) {
            printf("%d ", a11.matrix[i * a11.dim_b + j]);
        }
        printf("\n");
    }

    printf("\nSecondo blocco:\n");
    for (int i = a12.row; i < a12.row + a12.dim_b; i++) {
        for (int j = a12.col; j < a12.col + a12.dim_b; j++) {
            printf("%d ", a12.matrix[i * a12.dim_b + j]);
        }
        printf("\n");
    }

    printf("\nSomma:\n");
    sum2(&a11, &a12, &c11);
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%d ", C[i][j]);
        }
        printf("\n");
    }

}
