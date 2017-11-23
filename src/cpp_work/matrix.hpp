#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <random>

using Vec = std::vector<std::vector<double>>;

void round_number(double input, int digit, double* output)
{
	if(input != 0.0) {
		double d = pow(10.0, (double)digit);
		input *= d;
		*output = (int)(input + 0.5);
		*output /= d;
	} else {
		*output = 0.0;
	}
}

class Matrix
{
	private:
		Vec matrix;
	public:
		Matrix(int row, int column);
		void zeroMatrix(int row, int column);
		void zeroMatrix(int row, int column, Vec *v);
		void unitMatrix(int row, int column);
		void setMatrix(Vec v);
		Vec getMatrix();
		double getElement(int row, int column);
		void emptyMatrix();
		void randomCreateMatrix();
		
		void addMatrix(Vec v);
		void subtractMatrix(Vec v);
		void productMatrix(Vec v);
		void productMatrix(double num);
		
		void rotateMatrix(double theta);
		
		void disp();
		void disp(Vec v);
		
		void LU_decomposition(Vec *l, Vec *u);
		void determinantTri(Vec v, double *d);
		
		void invertMatrix();
		void invertMatrix(Vec l, Vec u, Vec *v);
};

Matrix::Matrix(int row, int column)
{
	Matrix::zeroMatrix(row, column, &matrix);
}

void Matrix::zeroMatrix(int row, int column)
{
	std::vector<double> _v;
	Vec().swap(matrix);
	
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < column; j++) {
			_v.push_back(0.0);
		}
		matrix.push_back(_v);
		_v.clear();
	}
}

void Matrix::zeroMatrix(int row, int column, Vec *v)
{
	std::vector<double> _v;
	Vec().swap(*v);
	
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < column; j++) {
			_v.push_back(0.0);
		}
		(*v).push_back(_v);
		_v.clear();
	}
}

void Matrix::unitMatrix(int row, int column)
{
	if(row != column) {
		std::cout << "(row != column) can't make unit matrix..." << std::endl;
		exit(1);
	}
	
	Matrix::zeroMatrix(row, column, &matrix);
	
	for(int i = 0; i < row; i++) {
		matrix[i][i] = 1.0;
	}
}

void Matrix::setMatrix(Vec v)
{
	// substitute v for matrix
	if((v.size() == matrix.size()) && (v[0].size() == matrix[0].size())) {
		for(int i = 0; i < v.size(); i++) {
			for(int j = 0; j < v[i].size(); j++) {
				matrix[i][j] = v[i][j];
			}
		}
	} else {
		std::cout << "matrix size error" << std::endl;
	}
}

Vec Matrix::getMatrix()
{
	return matrix;
}

double Matrix::getElement(int row, int column)
{
	return matrix[row][column];
}

void Matrix::productMatrix(Vec v)
{
	// prepare answer matrix
	Vec answer;
	int v_column = v[0].size();
	int v_row = matrix.size();

	// check row and column
	if(matrix[0].size() != v.size()) {
		std::cout << "can't product matrix..." << std::endl;
		exit(1);
	}

	Matrix::zeroMatrix(v_row, v_column, &answer);

	// product (matrix * v)
	for(int i = 0; i < answer.size(); i++) {
		for(int j = 0; j < answer[i].size(); j++) {
			for(int k = 0; k < v_row; k++) {
				answer[i][j] += matrix[i][k]*v[k][j];
			}
		}
	}

	// substitute answer for matrix
	Matrix::emptyMatrix();
	Matrix::zeroMatrix(v_row, v_column, &matrix);

	matrix = answer;
}

void Matrix::productMatrix(double num)
{
	for(int i = 0; i < matrix.size(); i++) {
		for(int j = 0; j < matrix[i].size(); j++) {
			matrix[i][j] *= num;
		}
	}
}

void Matrix::addMatrix(Vec v)
{
	for(int i = 0; i < v.size(); i++) {
		for(int j = 0; j < v[i].size(); j++) {
			matrix[i][j] += v[i][j];
		}
	}
}

void Matrix::subtractMatrix(Vec v)
{
	for(int i = 0; i < v.size(); i++) {
		for(int j = 0; j < v[i].size(); j++) {
			matrix[i][j] -= v[i][j];
		}
	}
}

void Matrix::rotateMatrix(double theta)
{
	Matrix m(2, 2);
	m.setMatrix(matrix);
	m.productMatrix({{cos(theta), -1*sin(theta)}, {sin(theta), cos(theta)}});
	
	matrix = m.getMatrix();
}

void Matrix::emptyMatrix()
{
	Vec().swap(matrix);
}

void Matrix::disp()
{
	double d;
	std::string str_matrix = "";
	
	for(int i = 0; i < matrix.size(); i++) {
		str_matrix += "[ ";
		for(int j = 0; j < matrix[i].size(); j++) {
			round_number(matrix[i][j], 6, &d);
			if(d == 0.0) str_matrix += "0        ";
			else {
				str_matrix += std::to_string(d) + "  ";
			}
		}
		str_matrix += "]\n\n";
	}
	
	std::cout << str_matrix << std::endl;
}

void Matrix::disp(Vec v)
{
	double d;
	std::string str_matrix = "";
	
	for(int i = 0; i < v.size(); i++) {
		str_matrix += "[ ";
		for(int j = 0; j < v[i].size(); j++) {
			round_number(v[i][j], 6, &d);
			if(d == 0.0) str_matrix += "0        ";
			else {
				str_matrix += std::to_string(d) + "  ";
			}
		}
		str_matrix += "]\n\n";
	}
	
	std::cout << str_matrix << std::endl;
}

void Matrix::LU_decomposition(Vec *l, Vec *u)
{
	int row = matrix.size();
	int column = matrix[0].size();
	
	if(row != column) {
		std::cout << "LU decomposition failed..." << std::endl;
		exit(1);
	}
	
	// set L and U matrix
	Vec().swap(*l);
	Vec().swap(*u);
	
	Matrix::zeroMatrix(row, column, l);
	Matrix::zeroMatrix(row, column, u);

	// set U_00 ~ U_0n and L_00, L_11, ...
	for(int i = 0; i < column; i++) {
		(*l)[i][i] = 1;
		(*u)[0][i] = matrix[0][i];
	}
	
	// calc L and U matrix
	for(int i = 1; i < row; i++) {
		for(int j = i; j < column; j++) {
			// L
			(*l)[j][i - 1] = matrix[j][i - 1];
			for(int k = 0; k < i - 1; k++) {
				(*l)[j][i - 1] -= (*l)[j][k]*(*u)[k][i - 1];
			}
			(*l)[j][i - 1] /= (*u)[i - 1][i - 1];
			// U
			(*u)[i][j] = matrix[i][j];
			for(int k = 0; k < i; k++) {
				(*u)[i][j] -= (*l)[i][k]*(*u)[k][j];
			}
		}
	}
	/*
	std::cout << "A=" << std::endl;
	Matrix::disp(matrix);
	std::cout << "L=" << std::endl;
	Matrix::disp(*l);
	std::cout << "U=" << std::endl;
	Matrix::disp(*u);
	*/
}

void Matrix::randomCreateMatrix()
{
	int row = matrix.size();
	int column = matrix[0].size();

	std::random_device rnd;
	
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < column; j++) {
			matrix[i][j] = rnd()%10 + 1.0;
		}
	}
}

void Matrix::determinantTri(Vec v, double *d)
{
	*d = 1.0;
	
	for(int i = 0; i < v.size(); i++) {
		(*d) *= v[i][i];
	}
}

void Matrix::invertMatrix()
{
	Vec l_mat, u_mat;
	Matrix m(matrix.size(), matrix[0].size());
	m.setMatrix(matrix);
	m.LU_decomposition(&l_mat, &u_mat);
	
	Matrix::invertMatrix(l_mat, u_mat, &matrix);
}

void Matrix::invertMatrix(Vec l, Vec u, Vec *x)
{
	Vec e, y;
	int row = l.size();
	int column = l[0].size();

	// check regular matrix
	double d;
	Matrix::determinantTri(u, &d);
	if(d == 0) {
		std::cout << "Matrix U is not regular." << std::endl;
		exit(1);
	}
	
	// set matrix X, Y and E
	Matrix::zeroMatrix(row, column, x);
	Matrix::zeroMatrix(row, column, &e);
	Matrix::zeroMatrix(row, column, &y);
	
	for(int i = 0; i < row; i++) {
		e[i][i] = 1.0;
	}
	
	// calculate Y (L*Y = E)
	for(int c = 0; c < column; c++) {
		for(int r = 0; r < row; r++) {
			y[r][c] = e[r][c];
			for(int i = 0; i < r; i++) {
				y[r][c] -= y[i][c]*l[r][i];
			}
		}
	}
	
	// calculate X (U*X = Y)
	for(int c = 0; c < column; c++) {
		for(int r = row - 1; r >= 0; r--) {
			(*x)[r][c] = y[r][c];
			for(int i = row - r - 1; i > 0; i--) {
				(*x)[r][c] -= (*x)[row - i][c]*u[r][row - i];
			}
			(*x)[r][c] /= u[r][r];
		}
	}
}

/*
int main(void)
{
	Vec l, u, m;
	double l_det, u_det;
	
	Matrix matrix(2, 2);
	matrix.randomCreateMatrix();
	matrix.LU_decomposition(&l, &u);
	
	matrix.invertMatrix(l, u, &m);
	matrix.productMatrix(m);
	matrix.disp();

	return 0;
}
*/


