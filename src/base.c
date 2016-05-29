#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>


void rowmatch_double(double* x, int* _xrn, int* _xcn,
                     double* tab, int* _tabrn,
                     int* r) {
    /*
     * in R:
     * - check shapes
     * - rowsort x, tab
     */
    int xrn = *_xrn;
    int xcn = *_xcn;
    int tabrn = *_tabrn;
    int i1, i2, j, j1, j2;
    double a, b;
    bool nextX = FALSE;
    bool nextTab = FALSE;

    i2 = 0;
    for (i1 = 0; i1 < xrn; i1++) { // for each row in x
        r[i1] = R_NaReal;
        for (; i2 < tabrn; i2++) { // find row in tab
            j1 = i1;
            j2 = i2;
            for (j = 0; j < xcn; j++) { // for each column
                a = x[j1];
                b = tab[j2];

                if (a < b) {
                    nextX = TRUE;
                    break;
                } else if (b < a) {
                    nextTab = TRUE;
                    break;
                }

                j1 += xrn;
                j2 += tabrn;
            }

            if (nextX) {
                nextX = FALSE;
                break;
            }

            if (nextTab) {
                nextTab = FALSE;
                continue;
            }

            r[i1] = i2;
            break;
        }
    }
}
