//-------------------------------------------------------------------------------------------
// Simplex.h
// Description: Класс, реализующий алгоритм двойственного симплекс-метода
// Date: 28.02.2023
// Authors: Хьюго М.А.
// Ⓒ Sibsutis university
//-------------------------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <iomanip>
#include <map>

#include "../SimpleFraction/SimpleFraction.h"
#include "../JordanGauss/JordanGauss.h"

typedef std::map<std::string, std::vector<Fraction::SimpleFraction>> TSimplexTable;

namespace LinearEquations {
class Simplex {
public:
	Simplex();

	void Calculate();

	void PrintMatrix(FractionMatrix& matrix);
	void PrintSimplexTable();

	FractionMatrix ReadMatrixFromFile(std::string fileName = "Simplex\\Matrix.txt");

protected:
	int m_variablesCount = 0;
	int m_solutionsCount = 0;
	int m_nextSolutionZInd = -1;

	bool m_tableIsValid = false;
	bool m_solutionIsOptimal= false;
	bool m_existNextSolution = false;
	bool m_haveNoValidStr = false;

	Fraction::SimpleFraction m_solutionValue;

	JordanGauss m_jordanGauss;

	FractionMatrix m_matrix;
	TBasicVariables m_basicVariables;

	TSimplexTable m_simplexTable;
	TSimplexTable m_prevSimplexTable;
	std::vector<Fraction::SimpleFraction> m_simplexDivisions;

	std::vector<Fraction::SimpleFraction> m_function;
	std::map<int, Fraction::SimpleFraction> m_expressedFunction;
	std::vector<std::vector<Fraction::SimpleFraction>> m_solutions;

	void FunctionExpressing();
	void CreateSourceTable();
	void FindResolvingElement(std::string& row, int& column);
	void FindNextSolutionResElem(std::string& row, int& column);
	void RecalculateTable();
	void FillNullCells(std::string& mainRow, int& mainColumn);
	void Finishing();

	bool TableIsValid();
	bool SolutionIsOptimal();

};
}