#include <iostream>
#include <vector>
#include <math.h>

using Vec = std::vector<std::vector<double>>;

class Matrix
{
	private:
		Vec matrix;
	public:
		Matrix(int row, int column);
		void setMatrix(Vec v);
		void getMatrix(Vec *v);
		void emptyMatrix();
		void productMatrix(Vec v);
		void rotateMatrix(double theta);
		void disp();
};

Matrix::Matrix(int row, int column)
{
	std::vector<double> v;
	Matrix::emptyMatrix();
	
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < column; j++) {
			v.push_back(0.0);
		}
		matrix.push_back(v);
		std::vector<double>().swap(v);
	}
	std::cout << "initMatrix success!" << std::endl;
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
		std::cout << "setMatrix success!" << std::endl;
	} else {
		std::cout << "matrix size error" << std::endl;
	}
}

void Matrix::getMatrix(Vec *v)
{
	*v = matrix;
	std::cout << "getMatrix success!" << std::endl;
}

void Matrix::productMatrix(Vec v)
{
	// prepare answer matrix
	Vec answer;
	std::vector<double> tmp;
	int v_column = v[0].size();
	int v_row = matrix.size();

	// set answer
	for(int i = 0; i < v_row; i++) {
		for(int j = 0; j < v_column; j++) {
			tmp.push_back(0.0);
		}
		answer.push_back(tmp);
		std::vector<double>().swap(tmp);
	}

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
	for(int i = 0; i < answer.size(); i++) {
		for(int j = 0; j < answer[i].size(); j++) {
			tmp.push_back(answer[i][j]);
		}
		matrix.push_back(tmp);
		std::vector<double>().swap(tmp);
	}
	
	std::cout << "productMatrix success!" << std::endl;
}

void Matrix::rotateMatrix(double theta)
{
	Matrix::productMatrix({{sin(theta), -1*cos(theta)}, {cos(theta), sin(theta)}});
}

void Matrix::emptyMatrix()
{
	Vec().swap(matrix);
}

void Matrix::disp()
{
	for(int i = 0; i < matrix.size(); i++) {
		for(int j = 0; j < matrix[i].size(); j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

int main(void)
{
	Vec m;
	Matrix mat(2, 2);
	
	mat.setMatrix({{1, 2}, {3, 4}});
	mat.disp();
	mat.productMatrix({{1, 2}, {3, 4}});
	mat.disp();
	mat.setMatrix({{1, 2}, {3, 4}});
	mat.productMatrix({{1},{2}});
	mat.disp();

	return 0;
}


