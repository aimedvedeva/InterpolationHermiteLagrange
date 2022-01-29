#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
#define pi 3.1415926536

double basicPolynomial(double value, int i, int n, double table[]);
double auxilaryHermiteCoefficient(int j, int n, double table[]);

double ContinuousFun(double value) {
	return sinl(value * pi);
}

double PiecewiseFun(double value) {
	if (value < -1.25) {
		return 3 * sqrt(abs(value - 1)) - 3.25;
	}
	else if (-1.25 <= value && value < 1.25) {
		return -value;
	}
	else if (value >= -1.25) {
		return 3.35 - 3 * sqrt(value + 1);
	}
}

double PiecewiseFunDerivative(double value) {
	if (value < -1.25) {
		return ( 3 * (value - 1) ) / ( 2 * sqrt( powl(abs(value - 1), 3) ) );
	}
	else if (-1.25 <= value && value > 1.25) {
		return -1;
	}
	else if (value >= -1.25) {
		return -3 / (2 * sqrt(value + 1));
	}
}

double HermitePilynomial(double value, int n, double table[], double (*DerFun) (double), double(*Fun) (double)) {
	double result = 0;
	for (int j = 0; j < n; j++) {
		double tmpValue = basicPolynomial(value, j, n, table);
		tmpValue *= tmpValue;

		double auxilaryVal = auxilaryHermiteCoefficient(j, n, table);
		double phy = (1 + auxilaryVal * (value - table[j]));
		double psy = (value - table[j]);
		double derivative = DerFun(table[j]);
		double funValue = Fun(table[j]);
		result += tmpValue * ( derivative * psy + funValue * phy );
	}
	return result;
}

double auxilaryHermiteCoefficient(int j, int n, double table[]) {
	double result = 0;
	for (int i = 0; i < n; i++) {
		if (i != j) {
			result += 1 / (table[j] - table[i]);
		}
	}
	return -2 * result;
}

double LagrangePolynomial(double value, int n, double table[], double(*Fun) (double)) {
	double result = 0;
	for (int i = 0; i < n; i++) {
		double tmpValue = basicPolynomial(value, i, n, table);
		double funValue = Fun(table[i]);
		result += funValue * tmpValue;
	}
	return result;
}

double basicPolynomial(double value, int i, int n, double table[]) {
	//i - num of value in table
	double result = 1;
	for (int j = 0; j < n; j++) {
		if (j != i) {
			result *= value - table[j];
			result /= table[i] - table[j];
		}
	}
	return result;
}

double *GetChebyshevGreed(double a, double b, int n) {

	double *table = new double[n];
	for (int i = 0; i < n; i++) {
		double sum = b + a;
		double subtract = b - a;
		double val = (2 * (double)i + 1) / (2 * ((double)n + 1)) * pi;
		double cosVal = cos(val);
		table[i] = sum / 2 + (subtract / 2) * cosVal;
	}
	return table;
}

double ErrDataLagrange(int numOfDots, double table[], double (*Fun) (double)) {
	double *midX = new double[numOfDots - 1];

	for (int i = 0; i < numOfDots - 1; i++) {
		midX[i] = (table[i] + table[i + 1]) / 2;
	}

	double maxVal = 0;
	for (int i = 0; i < numOfDots - 1; i++) {
		double curVal = abs(LagrangePolynomial(midX[i], numOfDots, table, ContinuousFun) - ContinuousFun(midX[i]));
		if (curVal > maxVal) {
			maxVal = curVal;
		}
	}

	delete[] midX;
	return maxVal;
}

double ErrDataHermite(int numOfDots, double table[], double(*Fun) (double)) {
	double *midX = new double[numOfDots - 1];

	for (int i = 0; i < numOfDots - 1; i++) {
		midX[i] = (table[i] + table[i + 1]) / 2;
	}

	double maxVal = 0;
	for (int i = 0; i < numOfDots - 1; i++) {
		double curVal = abs(HermitePilynomial(midX[i], numOfDots, table, PiecewiseFunDerivative, PiecewiseFun) - PiecewiseFun(midX[i]));
		if (curVal > maxVal) {
			maxVal = curVal;
		}
	}

	delete[] midX;
	return maxVal;
}
double der(double x){
	return cos(x);
}

double *GetValuesForGraph(double a, double b, int numOfValues) {
	double *values = new double[numOfValues];
	double delta = (b - a) / (double)numOfValues;
	values[0] = a + delta;
	for (int i = 1; i < numOfValues - 1; i++) {
		values[i] = values[i - 1] + delta;
	}
	values[numOfValues - 1] = b;
	return values;
}

int main(void) {
	ofstream fileX, fileY;
	fileY.open("Y.txt");
	fileX.open("X.txt");
    
	
	
	double a = -5;
	double b = 2;
	int numOfValues = 100;
	double *Values = GetValuesForGraph(a, b, numOfValues);
    int numOfDots = 3;
	double *ChebyshovTable = GetChebyshevGreed(a, b, numOfDots);

	double result = HermitePilynomial(-3, numOfDots, ChebyshovTable, PiecewiseFunDerivative, PiecewiseFun);
	for (int i = 0; i < numOfValues; i++) {
		fileX << Values[i] << endl;
		//fileY << LagrangePolynomial(Values[i], numOfDots, ChebyshovTable, ContinuousFun) << endl;
		//fileY << HermitePilynomial(Values[i], numOfDots, ChebyshovTable, der, ContinuousFun) << endl;
		//fileY << LagrangePolynomial(Values[i], numOfDots, ChebyshovTable, PiecewiseFun) << "   ";
		
		fileY << HermitePilynomial(Values[i], numOfDots, ChebyshovTable, PiecewiseFunDerivative, PiecewiseFun) << endl;
	}


	/*for (int i = 0; i < numOfDots; i++) {
		cout << ChebyshovTable[i] << endl;
		cout << LagrangePolynomial(ChebyshovTable[i], numOfDots, ChebyshovTable, PiecewiseFun) << endl;
		cout << PiecewiseFun(ChebyshovTable[i]) << endl;
		cout << endl;
	}*/

	fileX.close();
	fileY.close();

	//cout << ErrDataHermite(numOfDots, unevenTable, PiecewiseFun) << endl;*/
}