#include "matrix.h"

matrix22::matrix22() : a(1.0f), b(0.0f), c(0.0f), d(1.0f) { }
matrix22::matrix22(double a, double b, double c, double d) : a(a), b(b), c(c), d(d) { }

matrix22 matrix22::operator*(const matrix22& m) {
	return matrix22(this->a * m.a + this->b * m.c,
			this->a * m.b + this->b * m.d,
			this->c * m.a + this->d * m.c,
			this->c * m.b + this->d * m.d);
}

vector2 matrix22::operator*(const vector2& v) {
	return vector2(this->a * v.x + this->b * v.y, this->c * v.x + this->d * v.y);
}

matrix22 matrix22::operator+(const matrix22& m) {
	return matrix22(this->a + m.a, this->b + m.b, this->c + m.c, this->d + m.d);
}

matrix22 matrix22::operator-(const matrix22& m) {
	return matrix22(this->a - m.a, this->b - m.b, this->c - m.c, this->d - m.d);
}

matrix22 matrix22::operator*(const double s) {
	return matrix22(this->a*s, this->b*s, this->c*s, this->d*s);
}

matrix22& matrix22::operator=(const matrix22& m) {
	this->a = m.a; this->b = m.b; this->c = m.c; this->d = m.d;
	return *this;
}

double matrix22::det() {
	return this->a * this->d - this->b * this->c;
}

matrix22 matrix22::inv() {
	return matrix22(this->d, -this->b, -this->c, this->a) * (1.0 / this->det());
}

double matrix22::trace() {
	return this->a + this->d;
}

void matrix22::eigenvalues(double& lambda0, double& lambda1) {
	double a0 = (this->a + this->d) / 2.0;
	double a1 =  this->a - this->d;
	double a2 = sqrt(4*this->b*this->c + a1*a1)/2.0;
	lambda0 = a0 + a2;
	lambda1 = a0 - a2;
}

double* gaussJordanElimination(const double m_[], int n, bool& singular) {
	double *m = new double[n * n];
	memcpy(m, m_, n * n * sizeof(double));

	double *mInv = new double[n * n];
	double epsilon = 0.000001;
	
	// Identity;
	for (int i = 0; i < n; i++) { // column
		for (int j = 0; j < n; j++) { // row
			mInv[i * n + j] = i == j ? 1.0 : 0.0;
		}
	}

	// Row additions
	for (int i = 0; i < n - 1; i++) { // column
		int pivot = -1;
		for (int j = 0; j < n; j++) { // row
			if (fabs(m[i * n + j]) > epsilon) {
				pivot = j;
				break;
			}
			// singular
		}
		if (pivot < 0) {
			singular = true;
			break;
		}
		for (int j = i + 1; j < n; j++) { // column
			if (fabs(m[j * n + pivot]) > fabs(m[i * n + pivot])) {
				double t;
				for (int l = 0; l < n; l++) {
					t = m[i * n + l];
					m[i * n + l] = m[j * n + l];
					m[j * n + l] = t;
					t = mInv[i * n + l];
					mInv[i * n + l] = mInv[j * n + l];
					mInv[j * n + l] = t;
				}
			}
		}
		for (int j = i + 1; j < n; j++) { // column
			double s = m[j * n + pivot]/m[i * n + pivot];
			for (int k = 0; k < n; k++) { // row
				m[j * n + k] -= s * m[i * n + k];
				mInv[j * n + k] -= s * mInv[i * n + k];
			}
		}
	}

	// Row swaps
	for (int i = 0; i < n - 1 && !singular; i++) { // column
		int pivot = -1;
		for (int j = 0; j < n; j++) { // row
			if (fabs(m[i * n + j]) > epsilon) {
				pivot = j;
				break;
			}
			// singular
		}
		if (pivot < 0) {
			singular = true;
			break;
		}
		for (int j = i + 1; j < n; j++) { // column
			int pivot1 = -1;
			for (int k = 0; k < n; k++) { // row
				if (fabs(m[j * n + k]) > epsilon) {
					pivot1 = k;
					break;
				}
				// singular
			}
			if (pivot1 < 0) {
				singular = true;
				break;
			}
			if (pivot > pivot1) { // swap rows
				double t;
				for (int l = 0; l < n; l++) {
					t = m[i * n + l];
					m[i * n + l] = m[j * n + l];
					m[j * n + l] = t;
					t = mInv[i * n + l];
					mInv[i * n + l] = mInv[j * n + l];
					mInv[j * n + l] = t;
				}
			}
		}
	}

	// Row multiplications
	for (int i = 0; i < n && !singular; i++) { // column
		double s = m[i * n + i];
		for (int j = 0; j < n; j++) { // row
			m[i * n + j] /= s;
			mInv[i * n + j] /= s;
		}
	}

	// M is in reduced row echelon form

	// Final row additions to find inverse
	for (int i = n - 1; i > 0 && !singular; i--) { // column and row
		for (int j = i - 1; j >= 0; j--) { // column
			double s = m[j * n + i];
			for (int k = 0; k < n; k++) { // row
				m[j * n + k] -= s * m[i * n + k];
				mInv[j * n + k] -= s * mInv[i * n + k];
			}
		}
	}

	delete [] m;
	return mInv;
}



void cMatrix::output() {
	for (int j = 0; j < m; j++) {
		for (int i = 0; i < n; i++) {
			std::cout << std::setw(16) << A[j * n + i] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

cMatrix::cMatrix(const int m, const int n) : m(m), n(n) {
	A = new double[m * n];
	this->eye();
}

cMatrix::cMatrix(const int m, const int n, const double A[]) : m(m), n(n) {
	this->A = new double[m * n];
	for (int i = 0; i < m * n; i++) this->A[i] = A[i];
}

cMatrix::~cMatrix() {
	delete [] A;
}

cMatrix& cMatrix::set(const double A[]) {
	for (int i = 0; i < m * n; i++) this->A[i] = A[i];
	return *this;
}

cMatrix& cMatrix::eye() {
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			A[j * n + i] = i == j ? 1.0 : 0.0;
	return *this;
}

cMatrix& cMatrix::zero() {
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			A[j * n + i] = 0.0;
	return *this;
}

cMatrix cMatrix::transpose() {
	cMatrix R(n, m);
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			R.A[i * m + j] = A[j * n + i];
	return R;
}

cMatrix cMatrix::operator*(const cMatrix& B) {
	cMatrix R(m, B.n);
	for (int j = 0; j < m; j++)
		for (int i = 0; i < B.n; i++) {
			R.A[j * B.n + i] = 0.0;
			for (int k = 0; k < n; k++)
				R.A[j * B.n + i] += A[j * n + k] * B.A[k * B.n + i];
		}
	return R;
}

cMatrix cMatrix::operator*(const double s) {
	cMatrix R(m, n);
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			R.A[j * n + i] = A[j * n + i] * s;
	return R;
}

cMatrix cMatrix::operator+(const cMatrix& B) {
	cMatrix R(m, n);
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			R.A[j * n + i] = A[j * n + i] + B.A[j * n + i];
	return R;
}

cMatrix cMatrix::operator-(const cMatrix& B) {
	cMatrix R(m, n);
	for (int j = 0; j < m; j++)
		for (int i = 0; i < n; i++)
			R.A[j * n + i] = A[j * n + i] - B.A[j * n + i];
	return R;
}

cMatrix& cMatrix::operator=(const cMatrix& B) {
	if (m * n != B.m * B.n) {
		delete [] A;
		A = new double[B.m * B.n];
	}
	m = B.m;
	n = B.n;

	for (int i = 0; i < m * n; i++) A[i] = B.A[i];
	return *this;
}

void cMatrix::householderDecomposition(cMatrix& Q, cMatrix& R) {
	double mag, alpha;
	cMatrix u(m, 1), v(m, 1);
	cMatrix P(m, m), I(m, m);

	Q = cMatrix(m, m);
	R = *this;

	for (int i = 0; i < n; i++) {
		u.zero(); v.zero();
		
		mag = 0.0;
		for (int j = i; j < m; j++) {
			u.A[j] = R.A[j * n + i];
			mag += u.A[j] * u.A[j];
		}
		mag = sqrt(mag);
		
		alpha = u.A[i] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i; j < m; j++) {
			v.A[j] = j == i ? u.A[j] + alpha : u.A[j];
			mag += v.A[j] * v.A[j];
		}
		mag = sqrt(mag);

		if (mag < 0.0000000001) continue;

		for (int j = i; j < m; j++) v.A[j] /= mag;

		P = I - (v * v.transpose()) * 2.0;

		R = P * R;
		Q = Q * P;
	}
}

void cMatrix::householderBidiagonalization(cMatrix& Q, cMatrix& R, cMatrix& S) {
	double mag, alpha;

	cMatrix u(m, 1),  v(m, 1),
		u_(n, 1), v_(n, 1);

	cMatrix P(m, m),  I(m, m),
		P_(n, n), I_(n, n);

	Q = cMatrix(m, m);
	R = *this;
	S = cMatrix(n, n);

	for (int i = 0; i < n; i++) {
		u.zero(); v.zero();
		
		mag = 0.0;
		for (int j = i; j < m; j++) {
			u.A[j] = R.A[j * n + i];
			mag += u.A[j] * u.A[j];
		}
		mag = sqrt(mag);
		
		alpha = u.A[i] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i; j < m; j++) {
			v.A[j] = j == i ? u.A[j] + alpha : u.A[j];
			mag += v.A[j] * v.A[j];
		}
		mag = sqrt(mag);

		if (mag > 0.0000000001) {
			for (int j = i; j < m; j++) v.A[j] /= mag;

			P = I - (v * v.transpose()) * 2.0;

			R = P * R;
			Q = Q * P;
		}
/////////////////////////
		u_.zero(); v_.zero();

		mag = 0.0;
		for (int j = i + 1; j < n; j++) {
			u_.A[j] = R.A[i * n + j];
			mag += u_.A[j] * u_.A[j];
		}

		mag = sqrt(mag);

		alpha = u_.A[i + 1] < 0 ? mag : -mag;

		mag = 0.0;
		for (int j = i + 1; j < n; j++) {
			v_.A[j] = j == i + 1 ? u_.A[j] + alpha : u_.A[j];
			mag += v_.A[j] * v_.A[j];
		}
		mag = sqrt(mag);

		if (mag > 0.0000000001) {
			for (int j = i + 1; j < n; j++) v_.A[j] /= mag;

			P_ = I_ - (v_ * v_.transpose()) * 2.0;

			R = R * P_;
			S = P_ * S;
		}
	}
}

vector3 cMatrix::mult(const vector3& v) {
	vector3 r;
	if (m == 4 && n == 4) {
		r.x = this->A[0] * v.x + this->A[1] * v.y + this->A[2]  * v.z + this->A[3]  * 1;
		r.y = this->A[4] * v.x + this->A[5] * v.y + this->A[6]  * v.z + this->A[7]  * 1;
		r.z = this->A[8] * v.x + this->A[9] * v.y + this->A[10] * v.z + this->A[11] * 1;
	}
	return r;
}

vector3 cMatrix::mult_(const vector3& v) {
	vector3 r;
	if (m == 4 && n == 4) {
		r.x = this->A[0] * v.x + this->A[1] * v.y + this->A[2];
		r.y = this->A[4] * v.x + this->A[5] * v.y + this->A[6];
		r.z = this->A[8] * v.x + this->A[9] * v.y + this->A[10];
	}
	return r;
}
