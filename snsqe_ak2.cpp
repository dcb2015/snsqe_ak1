// snsqe_ak2 - A program for finding a zero of a system of
// N non-linear functions in N variables by a modification of
// the Powell hybrid method.
//
// This program is a translation of the FORTRAN routine SNSQE, by K. L. Hiebert(SNLA), which is posted on the SLATEC site:
// https://www.netlib.org/slatec/src/snsqe.f
//
// To distinguish this program from others, an '_ak1' suffix has been appended to its name.
// 
// ***            IMPORTANT!         ***
// A major change to this version of the program, compared
// to the original FORTRAN version, is that the Jacobian matrix, fjac, is
// the transpose of the corresponding FORTRAN matrix FJAC.
// 
// In the original FORTRAN version, FJAC is [J(i,j)] = d(fi)/dxj.
// Given a small change in the j-th variable, 
// how is the i-th function affected?
// FJAC is oriented column-wise.
// In the present version fjac is [J(i,j)]^T
// This format makes it more convenient in C++ and some of the routines
// can be more concise.
// fjac is oriented row-wise and all the sub-routines in the program have
// been edited to take into account this change.
//
// 27 January 2019
//
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
//

#define _SCL_SECURE_NO_WARNINGS

#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

#define N 2	// Number of dimensions and functions in this program
#define LR 3                                /* Equals (N*(N + 1))/2 */

double Euclid2Norm_ak1(double v[], int lim);
void F(double in[N], double av[N][8], double out[N]);
void J_ak1(double in[N], double coeffVec[N][8], double out[N][N]);
void fdjac1_ak1(double X[N], double coeffVec[N][8], double FVEC[N], double a[N][N]);
void qrfac_ak1(double A[N][N], double sigma[N]);
void qform_ak1(double q[N][N], double wa[N]);
void dogleg(double r[LR], double qtb[N], double del, double x[N]);
void r1updt(double s[LR], double u[N], double v[N], double w[N]);
void r1mpyq(int m, double a[], double v[N], double w[N]);
void snsqe_ak1(double X[N], double A[N][8], double FVEC[N], int iopt, double tol, int *info);

double Euclid2Norm_ak1(double v[], int lim){

	// This subroutine is based on the LINPACK routine SNRM2, written 25 October 1982, modified
	// on 14 October 1993 by Sven Hammarling of NAG Ltd.
	// I have modified it for use with a row of an array rather than a column vector.
	//
	// v - Matrix
	// lim - col index, col to right of which (inclusive) norm of row is to be calculated

	double absxi, dummy, scale = 0.0, ssq = 1.0;

	if (v[0] != 0) scale = fabs(v[0]);

	for (int i = 1; i < lim; ++i){
		if (v[i] != 0){
			absxi = fabs(v[i]);
			dummy = scale / absxi;
			ssq = 1.0 + ssq*dummy*dummy;
			scale = absxi;
		}//End if (v[i] != 0)
	}//End for i

	return scale*sqrt(ssq);

} //End Euclid2Norm_ak1

void F(double in[N], double av[N][8], double out[N]){
	// Evaluates the N functions of the N variables, outputs results in out.

	double f1xsum, f1ysum;

	f1xsum = ((av[0][2] + av[0][3] * in[1])*in[0] + av[0][1] * in[1] + av[0][0])*in[0];
	f1ysum = ((av[0][5] + av[0][4] * in[0])*in[1] + av[0][6])*in[1];
	out[0] = f1xsum + f1ysum + av[0][7];

	f1xsum = ((av[1][2] + av[1][3] * in[1])*in[0] + av[1][1] * in[1] + av[1][0])*in[0];
	f1ysum = ((av[1][5] + av[1][4] * in[0])*in[1] + av[1][6])*in[1];
	out[1] = f1xsum + f1ysum + av[1][7];

	return;   // End F
}

void J_ak1(double in[N], double coeffVec[N][8], double out[N][N]){

	// User-supplied routine to create the Jacobian matrix.
	// Note that it creates [J(i,j)]^T, the transpose of [J(i,j)] = d(fi)/dxj
	// (Given a small change in the j-th variable, how is the i-th function affected?)
	// This format makes it more convenient in C++ for some of the other sub-
	// routines in this program.

	out[0][0] = coeffVec[0][0] + 2.0 * in[0] * (coeffVec[0][2] + coeffVec[0][3] * in[1]) + in[1] * (coeffVec[0][1] + coeffVec[0][4] * in[1]);
	out[1][0] = in[0] * (coeffVec[0][1] + coeffVec[0][3] * in[0]) + 2.0 * in[1] * (coeffVec[0][4] * in[0] + coeffVec[0][5]) + coeffVec[0][6];
	out[0][1] = coeffVec[1][0] + 2.0 * in[0] * (coeffVec[1][2] + coeffVec[1][3] * in[1]) + in[1] * (coeffVec[1][1] + coeffVec[1][4] * in[1]);
	out[1][1] = in[0] * (coeffVec[1][1] + coeffVec[1][3] * in[0]) + 2.0 * in[1] * (coeffVec[1][4] * in[0] + coeffVec[1][5]) + coeffVec[1][6];

	return;  // End J_ak1
}

void fdjac1_ak1(double X[N], double coeffVec[N][8], double FVEC[N], double a[N][N]){

	// Computes the Jacobian matrix by Forward Difference approximation.
	// Note that it creates [J(i,j)]^T, the transpose of [J(i,j)] = d(fi)/dxj
	// (Given a small change in the j-th variable, how is the i-th function affected?)
	// This format makes it more convenient in C++ for some of the other sub-
	// routines in this program.

	double h, temp, w1[N];
	int i, j;
	static const double EPS = sqrt(DBL_EPSILON);

	for (j = 0; j < N; ++j)	{
		temp = X[j];
		h = EPS*fabs(temp);
		if (h == 0.0) h = EPS;
		X[j] = temp + h;
		F(X, coeffVec, w1);
		X[j] = temp;
		for (i = 0; i < N; ++i)	a[j][i] = (w1[i] - FVEC[i]) / h;
	}  // End for j
	return;
} // End fdjac1_ak1

void qrfac_ak1(double A[N][N], double sigma[N]){
	double ajnorm, sum;
	int i, j, k;

	for (j = 0; j < N; ++j)	{

		ajnorm = Euclid2Norm_ak1(&A[j][j], N - j);
		if (ajnorm != 0.0) {
			if (A[j][j] < 0.0)	ajnorm = -ajnorm;
			for (i = j; i < N; ++i) A[j][i] /= ajnorm;
			A[j][j] += 1.0;
			for (k = j + 1; k < N; ++k)	{
				sum = 0.0;
				for (i = j; i < N; ++i) sum += A[j][i] * A[k][i];
				sum /= A[j][j];
				for (i = j; i < N; ++i) A[k][i] -= sum*A[j][i];
			}  // End for k
		}  // End if ajnorm
		sigma[j] = -ajnorm;

	}  // End for j
	return;
}  // End qrfac_ak1

void qform_ak1(double q[N][N], double wa[N]) {
	double sum;
	int i, j, k, kp1 = N + 1;

	//Elements below diagonal become zero
	// This loop has been moved out of this routine and into snsqe, in the assignment loop for r.
	//for (i = 1; i < N; ++i){
	//	for (j = 0; j < i; ++j) q[i][j] = 0.0;
	//}

	// Accumulate q from its factored form

	for (k = (N - 1); k >= 0; --k)	{

		--kp1;  // kp1 = k + 1

		for (i = k; i < N; ++i)	{
			wa[i] = q[k][i];
			q[k][i] = 0.0;
		} //End for i
		q[k][k] = 1.0;

		if (wa[k] != 0.0){

			// Take j = k iteration out of for j loop
			q[k][k] = -wa[k] + 1.0;
			for (i = kp1; i < N; ++i) q[k][i] = -wa[i];

			for (j = kp1; j < N; ++j)	{
				sum = 0.0;
				for (i = kp1; i < N; ++i) sum += q[j][i] * wa[i];
				sum /= wa[k];
				q[j][k] = -sum*wa[k];
				for (i = kp1; i < N; ++i) q[j][i] -= sum*wa[i];
			} // End for j

		} // End if (wa[k] != 0)

	} // End for k

	return;
} // End qform_ak1

void dogleg(double r[LR], double qtb[N], double del, double x[N]){
	double alpha, bn, deldivqn, deldqn2, gn, qn, sg, sgdivdel, sgddel2, sum, temp, w1[N], w2[N];
	int i, j, jj = LR, l, NM1 = N - 1;

	for (j = NM1; j >= 0; --j) {
		jj -= N - j;
		l = jj + 1;
		w1[j] = sum = 0.0;
		for (i = (j + 1); i < N; ++i) sum += r[l++] * x[i];
		temp = r[jj];
		if (temp == 0.0) {
			l = j;
			for (i = 0; i <= j; ++i) {
				bn = fabs(r[l]);  // Use bn as a temp variable
				if (bn > temp) temp = bn;
				l += NM1 - i;
			} // End for i
			temp *= DBL_EPSILON;
			if (temp == 0.0) temp = DBL_EPSILON;
		} // End if (temp == 0.0)
		x[j] = (qtb[j] - sum) / temp;
	} // End for j

	qn = Euclid2Norm_ak1(x, N);
	if (qn <= del) return;

	l = 0;
	for (j = 0; j < N; ++j)	{
		temp = qtb[j];
		for (i = j; i < N; ++i)	w1[i] += r[l++] * temp;
	} // End for j

	gn = Euclid2Norm_ak1(w1, N);
	sg = 0.0;
	alpha = del / qn;

	if (gn != 0.0){

		for (j = 0; j < N; ++j)	w1[j] /= gn;

		l = 0;
		for (j = 0; j < N; ++j) {
			sum = 0.0;
			for (i = j; i < N; ++i)	sum += r[l++] * w1[i];
			w2[j] = sum;
		} // End for j

		temp = Euclid2Norm_ak1(w2, N);
		sg = gn / temp / temp;
		alpha = 0.0;

		if (sg < del){

			sgdivdel = sg / del;
			sgddel2 = sgdivdel * sgdivdel;
			deldivqn = del / qn;
			deldqn2 = deldivqn * deldivqn;
			bn = Euclid2Norm_ak1(qtb, N);
			temp = (bn / gn) * (bn / qn) * sgdivdel;
			sum = (temp - deldivqn) * (temp - deldivqn);
			sum = sqrt(sum + (1.0 - deldqn2) * (1.0 - sgddel2));
			temp += sum - deldivqn * sgddel2;
			alpha = (deldivqn * (1.0 - sgddel2)) / temp;

		} // End if (sg < del)

	} // End if (gn != 0.0)

	//twenty (or 120 in FORTRAN)

	// Form appropriate convex combination of the Gauss-Newton direction
	// and the scaled gradient direction

	temp = sg;
	if (del < temp) temp = del;
	temp *= 1.0 - alpha;
	for (j = 0; j < N; ++j)	x[j] = temp * w1[j] + alpha * x[j];

	return;
} // End dogleg

void r1updt(double s[LR], double u[N], double v[N], double w[N]){
	double co, cot, si, ta, tau, temp;
	int i, j, jj, l, NM1 = N - 1;

	jj = LR - 1;
	w[NM1] = s[jj];

	for (j = NM1 - 1; j >= 0; --j) {

		jj -= N - j;
		w[j] = 0.0;

		if (v[j] == 0.0) continue;

		if (fabs(v[NM1]) < fabs(v[j])) {

			cot = v[NM1] / v[j];
			si = 0.5 / sqrt(0.25 + 0.25 * cot * cot);
			co = si * cot;
			tau = 1.0;

			if (fabs(co)*DBL_MAX > 1.0) tau = 1.0 / co;

		} // End if fabs(v[NM1]) < fabs(v[j])

		else {

			ta = v[j] / v[NM1];
			co = 0.5 / sqrt(0.25 + 0.25 * ta * ta);
			tau = si = co * ta;

		} // End else

		v[NM1] = si * v[j] + co * v[NM1];
		v[j] = tau;

		l = jj;

		for (i = j; i < N; ++i) {
			temp = co * s[l] - si * w[i];
			w[i] = si * s[l] + co * w[i];
			s[l++] = temp;
		} // End for i

	} // End for j

	for (i = 0; i < N; ++i) w[i] += v[NM1] * u[i];

	for (j = 0; j < NM1; ++j) {

		if (w[j] != 0.0){

			if (fabs(s[jj]) < fabs(w[j])) {
				cot = s[jj] / w[j];
				si = 0.5 / sqrt(0.25 + 0.25 * cot * cot);
				co = si * cot;
				tau = 1.0;
				if (fabs(co)*DBL_MAX > 1.0)	tau = 1.0 / co;
			} // End if fabs
			else {
				ta = w[j] / s[jj];
				co = 0.5 / sqrt(0.25 + 0.25 * ta * ta);
				tau = si = co * ta;
			} // End else

			l = jj;

			for (i = j; i < N; ++i)	{
				temp = co * s[l] + si * w[i];
				w[i] = co * w[i] - si * s[l];
				s[l++] = temp;
			} // End for i

			w[j] = tau;

		} // End if (w[j] != 0.0)

		jj += N - j;

	} // End for j

	s[jj] = w[NM1];

	return;
} // End r1updt

void r1mpyq(int m, double a[], double v[N], double w[N]) {
	double co, si, temp;
	int NM1 = N - 1, i, j, k1, k2 = NM1 * m, k4;

	for (j = NM1 - 1; j >= 0; --j) {

		k1 = m * j;
		if (fabs(v[j]) > 1.0) {
			co = 1.0 / v[j];
			si = sqrt(1.0 - co * co);
		} // End if
		else {
			si = v[j];
			co = sqrt(1.0 - si * si);
		} // End else

		k4 = k2;
		for (i = 0; i < m; ++i) {
			temp = co*a[k1] - si*a[k4];
			a[k4] = si*a[k1] + co*a[k4++];
			a[k1++] = temp;
		} // End for i

	} // End for j

	k1 = 0;
	for (j = 0; j < NM1; ++j) {

		if (fabs(w[j]) > 1.0) {
			co = 1.0 / w[j];
			si = sqrt(1.0 - co*co);
		} // End if
		else {
			si = w[j];
			co = sqrt(1.0 - si*si);
		} // End else

		k4 = k2;
		for (i = 0; i < m; ++i)	{
			temp = co*a[k1] + si*a[k4];
			a[k4] = co*a[k4++] - si*a[k1];
			a[k1++] = temp;
		} // End for i

	} // End for j

	return;
} // End r1mpyq

void snsqe_ak1(double X[N], double A[N][8], double FVEC[N], int iopt, double tol, int *info) {
	/* This function finds the zero of a system of N non-linear
	functions in N variables by a modification of the Powell hybrid method.

	Parameters
	iopt = 1  Program employs user-supplied subroutine J to determine Jacobian matrix.
	iopt = 2  Program calculates Jacobian with a forward-difference approximation.
	tol       Specifies desired error between solution, X and actual zero.
	Unless high precision is required, the recommended value for
	tol is the square root of the machine precision.
	info < 0  User terminated execution.
	info = 0  Input parameters or LR in err.
	info = 1  Success. Relative error between X and zero is at most	tol.
	info = 2  Number of function calls has exceeded 100*(N + 1) for iopt = 1
	or 200*(N + 1) for iopt = 2. If iopt = 2 may solve the problem
	by supplying J, which may also allow tolerance to be tightened
	for more accuracy.
	info = 3  tol is too small. No further improvement in the approximate solution X is possible.
	info = 4  iteration is not making good progress, as measured by the improvement
	from the last five Jacobian evaluations. Maybe try another starting point.
	info = 5  iteration is not making good progress, as measured by the improvement
	from the last ten iterations. Maybe try another starting point.

	Author:    Hiebert, K. L. (SNLA)
	*/

	double actred, delta, fjac[N][N], fnorm, fnorm1, pnorm, prered, qtf[N],
		r[LR], ratio, sum, wa1[N], wa2[N], wa3[N], wa4[N], xnorm;
	int i, iter, j, l, maxfev, ncsuc, ncfail, nfev, NM1 = N - 1, nslow1, nslow2;
	bool jeval, x_val_Changed;

	*info = nslow2 = nslow1 = ncfail = ncsuc = 0;

	if ((iopt < 1) || (iopt > 2) || (N <= 0) || (tol <= 0.0) || (LR < (N*(N + 1)) / 2)){
		cout << "Invalid input parameter." << endl;
		return;
	}

	maxfev = 100 * (N + 1);
	if (iopt == 2)  maxfev *= 2;

	F(X, A, FVEC);
	iter = nfev = 1;

//	cout << "FVEC[0] = " << FVEC[0] << endl;
//	cout << "FVEC[1] = " << FVEC[1] << endl << endl;

	fnorm = Euclid2Norm_ak1(FVEC, N);

	// Calculate the norm of X and initialize the step bound delta

	xnorm = Euclid2Norm_ak1(X, N);
	delta = 100.0 * xnorm;
	if (delta == 0.0)  delta = 100.0;

	// 30 ***  Beginning of outer loop  ***
	for (;;){
		
		// Calculate the Jacobian matrix

		jeval = 1;

		if (iopt == 1) { // User-supplied Jacobian
			J_ak1(X, A, fjac);
		}
		else {  // Compute Jacobian by Forward Difference approximation
			fdjac1_ak1(X, A, FVEC, fjac);
			nfev += N;
		}

		// Compute the QR Factorization of the Jacobian

		qrfac_ak1(fjac, wa1);

		// Form Q^T * FVEC and store in qtf

		//for (i = 0; i < N; ++i) qtf[i] = FVEC[i];
		copy(FVEC, FVEC + N, qtf);		// Use the built-in copy function

		for (j = 0; j < N; ++j) {
			if (fjac[j][j] != 0.0) {
				sum = 0.0;
				for (i = j; i < N; ++i) sum += fjac[j][i] * qtf[i];
				sum /= -fjac[j][j];
				for (i = j; i < N; ++i) qtf[i] += fjac[j][i] * sum;
			}  // End if fjac
		}  // End for j

		// Copy the triangular factor of the QR factorization into r

		// Take first iteration out of loop
		r[0] = wa1[0];

		for (j = 1; j < N; ++j)	{
			l = j;
			for (i = 0; i < j; ++i)	{
				r[l] = fjac[j][i];
				fjac[j][i] = 0.0;  // This line is from the first loop of qform.
				l += NM1 - i;
			} // End for i
			r[l] = wa1[j];
		} // End for j

		// Accumulate the orthogonal factor in fjac
		qform_ak1(fjac, wa1);

		// ***  180 Beginning of inner loop  ***

		x_val_Changed = 0;

		for (;;){
			
			// Determine the direction P

			dogleg(r, qtf, delta, wa1);

			// Store the direction p and x + p. Calculate the norm of p

			pnorm = Euclid2Norm_ak1(wa1, N);

			for (j = 0; j < N; ++j) {
				wa1[j] = -wa1[j];
				wa2[j] = X[j] + wa1[j];
			} // End for j

			if ((iter == 1) && (pnorm < delta)) delta = pnorm;
			
			// Evaluate the function at x + p and calculate its norm

			F(wa2, A, wa4);
			++nfev;
			fnorm1 = Euclid2Norm_ak1(wa4, N);

			// Compute the actual reduction

			actred = -1.0;
			if (fnorm1 < fnorm)	actred = 1.0 - (fnorm1 / fnorm)*(fnorm1 / fnorm);

			// Compute the predicted reduction

			l = 0;
			for (i = 0; i < N; ++i) {
				sum = 0.0;
				for (j = i; j < N; ++j)	sum += r[l++] * wa1[j];
				wa3[i] = qtf[i] + sum;
			} // End for i

			sum = Euclid2Norm_ak1(wa3, N);
			ratio = prered = 0.0;
			if (sum < fnorm) prered = 1.0 - (sum / fnorm)*(sum / fnorm);

			// Compute the ratio of the actual to the predicted reduction

			if (prered > 0.0) ratio = actred / prered;

			// Update the step bound

			if (ratio < 0.1) {
				ncsuc = 0;
				++ncfail;
				delta *= 0.5;
			} // End if (ratio < 0.1)
			else {
				ncfail = 0;
				++ncsuc;
				if ((ratio >= 0.5) || (ncsuc > 1)) delta = (double)max(delta, 2 * pnorm);
				if (fabs(ratio - 1.0) <= 0.1) delta = 2.0 * pnorm;
			} // End else

			// Test for successful iteration

			if (ratio >= 0.0001) {

				// Successful iteration. Update X, FVEC, and their norms

				//for (j = 0; j < N; ++j) {
					//X[j] = wa2[j];
					copy(wa2, wa2 + N, X);		// Use the built-in copy function
					//FVEC[j] = wa4[j];
					copy(wa4, wa4 + N, FVEC);		// Use the built-in copy function
				//} // End for j

				xnorm = Euclid2Norm_ak1(wa2, N);
				fnorm = fnorm1;
				++iter;
				x_val_Changed = 1;

			} // End if ratio

			// Determine the progress of the iteration

			++nslow1;
			if (actred >= 0.001) nslow1 = 0;
			if (jeval) ++nslow2;
			if (actred >= 0.1) nslow2 = 0;

			// Test for convergence
			// ***  The following statements test exit conditions for the outer loop  ***

			if ((delta <= tol*xnorm) || (fnorm == 0.0))	{
				*info = 1;
				return;  // Break out of outer loop, leave sub-routine
			}

			// Tests for termination and stringent tolerances

			if (nfev >= maxfev)	*info = 2;

			if (0.1*(double)max(0.1*delta, pnorm) <= DBL_EPSILON*xnorm) *info = 3;

			if (nslow2 == 5) *info = 4;

			if (nslow1 == 10) *info = 5;

			if (*info != 0) return;  // Break out of outer loop, leave sub-routine

			// Criterion for recalculating Jacobian

			if ((x_val_Changed) && (ncfail == 2)) break; // Break out of inner loop

			// Calculate the rank-one modification to the Jacobian and
			// update QTF if necessary

			for (j = 0; j < N; ++j)	{
				sum = 0.0;
				for (i = 0; i < N; ++i)	sum += fjac[j][i] * wa4[i];
				wa2[j] = (sum - wa3[j]) / pnorm;
				wa1[j] /= pnorm;
				if (ratio >= 0.0001) qtf[j] = sum;
			} // End for j

			// Compute the QR factorization of the updated Jacobian

			r1updt(r, wa1, wa2, wa3);
			r1mpyq(N, &fjac[0][0], wa2, wa3);
			r1mpyq(1, qtf, wa2, wa3);

			jeval = 0;

		} // End for (;;) ***  End of inner loop  ***

	} // End for (;;) ***  End of outer loop  ***

	return;
} // End snsqe_ak1

int main() {
	char rflag;			//Readiness flag

	cout << "                       snsqe_ak2 (27 January 2019)" << endl;
	cout << "======================================================================" << endl;
	cout << "\nThis program finds a zero of a system of N non-linear functions in N variables:" << endl << endl;
	cout << "f1 = a00*x + a01*x*y + a02*x^2 + a03*x^2*y + a04*x*y^2 + a05*y^2 + a06*y + c1" << endl << endl;
	cout << "f2 = a10*x + a11*x*y + a12*x^2 + a13*x^2*y + a14*x*y^2 + a15*y^2 + a16*y + c2" << endl << endl;
	cout << "The results are calculated to double precision--" << DBL_DIG << " decimal places." << endl << endl;
	cout << "\nEverything ready? If yes, press y." << endl;
	cout << "Otherwise, press any other key." << endl;
	cin >> rflag;
	if (toupper(rflag) == 'Y') {

		double X[N], FVEC[N], A[N][8], tol;
		int info, iopt;

		A[0][0] = 0.1;
		A[0][1] = -0.2;
		A[0][2] = 1.4;
		A[0][3] = -2.2;
		A[0][4] = -1.4;
		A[0][5] = 0.7;
		A[0][6] = -0.5;
		A[0][7] = 0.633;

		A[1][0] = -0.333;
		A[1][1] = 1.2;
		A[1][2] = 0.637;
		A[1][3] = 0.618;
		A[1][4] = -0.736;
		A[1][5] = -.824;
		A[1][6] = 2.1;
		A[1][7] = 0.75;

		// Echo the array to the console
		for (unsigned int i = 0; i < N; ++i) {
			for (unsigned int j = 0; j < 8; ++j) {
				cout << A[i][j] << " ";
			} // End for j
			cout << "\n";
		} // End for i
		cout << "\n";

		X[0] = 0.5;
		X[1] = 0.5;                           /* Initial guess for zero. */
		tol = sqrt(DBL_EPSILON);

		// iopt = 1: User-supplied function to compute Jacobian
		// iopt = 2: Jacobian computed by Forward Difference approximation
		iopt = 1;

		snsqe_ak1(X, A, FVEC, iopt, tol, &info);

		cout << endl;
		cout << "Just returned from snsqe." << endl;
		cout << "The error code, info, is " << info << endl;
		cout << endl;
		cout << "The roots are as follows:" << endl;
		cout << endl;
		cout << "X[0] is " << X[0] << endl;
		cout << "X[1] is " << X[1] << endl;
		cout << endl;
		cout << "The function values are:" << endl;
		cout << endl;
		cout << "FVEC[0] is " << FVEC[0] << endl;
		cout << "FVEC[1] is " << FVEC[1] << endl;

	}	//End if ready
	else cout << "Not ready. Try again when ready with information.\n";

	cout << "\nEnter any key to continue." << endl;
	cin >> rflag;

	return 0;

} // End main