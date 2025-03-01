#include <stdbool.h>

/**
 * @brief Calculates the midpoint between two integers.
 *
 * This inline function computes the midpoint between two integer values \c a and \c b.
 * It uses bitwise operations to avoid potential overflow that might occur with large integers.
 *
 * @param a The first integer.
 * @param b The second integer.
 * @return The midpoint between \c a and \c b.
 *
 * @note This function assumes that the midpoint calculation does not overflow.
 *       For very large integers, consider using a different method to prevent overflow.
 */
static inline int midpoint(int a, int b) {
    return a + ((b - a) >> 1);
}

/**
 * @brief Searches for the interval containing a specific value within a sorted array.
 *
 * The \c hunt function performs an efficient search to find the index \c jlo such that:
 * - If the array \c X is sorted in ascending order:
 *     \f$ X[jlo] \leq XP < X[jlo + 1] \f$
 * - If the array \c X is sorted in descending order:
 *     \f$ X[jlo + 1] \leq XP < X[jlo] \f$
 *
 * The search starts from an initial guess provided through \c jlo and dynamically adjusts
 * the search bounds using a combination of exponential search (hunting phase) and
 * binary search (bisection phase) to locate the correct interval efficiently.
 *
 * @param X Pointer to the first element of a sorted array of double values.
 * @param N The number of elements in the array \c X.
 * @param XP The value to search for within the array \c X.
 * @param jlo Pointer to an integer representing the initial guess index.
 *            On return, it holds the lower index of the interval containing \c XP.
 *
 * @return void
 *
 * @note 
 * - If \c XP is outside the bounds of \c X, \c jlo is set to \c -1 or \c N accordingly.
 * - The function handles both ascending and descending sorted arrays.
 * - The initial value of \c jlo should be a reasonable guess to optimize search performance.
 */
void hunt(const double *X, int N, double XP, int *jlo) {
    bool ascend = (X[N - 1] > X[0]);
    int jhi;
    int inc = 1;

    // Validate initial guess
    if (*jlo < 0 || *jlo >= N) {
        *jlo = -1;
        jhi = N;
    } else {
        // Determine direction to hunt
        if ((XP >= X[*jlo]) == ascend) { // Hunt up
            jhi = *jlo + inc;
            while (jhi < N && ((XP >= X[jhi]) == ascend)) {
                *jlo = jhi;
                inc <<= 1; // Equivalent to inc *= 2
                jhi = *jlo + inc;
            }
            if (jhi >= N) {
                jhi = N;
            }
        } else { // Hunt down
            jhi = *jlo;
            while (true) {
                *jlo = jhi - inc;
                if (*jlo < 0) {
                    *jlo = -1;
                    break;
                }
                if (((XP < X[*jlo]) == ascend)) {
                    jhi = *jlo;
                    inc <<= 1; // Equivalent to inc *= 2
                } else {
                    break;
                }
            }
        }
    }

    // Bisection phase
    while (jhi - *jlo > 1) {
        int jm = midpoint(*jlo, jhi);
        if ((XP > X[jm]) == ascend) {
            *jlo = jm;
        } else {
            jhi = jm;
        }
    }
}
