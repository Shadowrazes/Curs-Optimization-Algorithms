#include "JordanGauss.h"

namespace LinearEquations {
//-----------------------------------------------------------------------------------
JordanGauss::JordanGauss() {

}
//-----------------------------------------------------------------------------------
void JordanGauss::Process(FractionMatrix& matrix) {
	m_result.clear();
	PrintMatrix(matrix);

	std::vector<int> readyRows;
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[i].size(); j++) {
			if (matrix[i][j].IsOne()) {
				int zeroCount = 0;
				for (int k = 0; k < matrix.size(); k++) {
					if (k != i && matrix[k][j].IsZero()) {
						zeroCount++;
					}
				}

				if (zeroCount == matrix.size() - 1) {
					readyRows.push_back(i);
					break;
				}
			}
		}
	}

	for (int i = 0, j = 0; i < matrix.size() && j < m_variablesCount; i++, j++) {
		if (CheckFalseEqualities(matrix)) {
			std::cout << "The system have no solution\n\n";
			m_result.clear();
			return;
		}

		if (std::find(readyRows.begin(), readyRows.end(), i) != readyRows.end()) continue;

		if (FindAndSwapMaxColumnElem(matrix, i, j)) {
			PrintMatrix(matrix);
		}

		if (matrix[i][j].IsZero()) {
			i--;
			continue;
		}

		if (!matrix[i][j].IsOne()) {
			LineDivision(matrix, i, j);
			PrintMatrix(matrix);
		}

		bool addicted = false;
		for (int k = 0; k < matrix.size(); k++) {
			if (k != i && !matrix[k][j].IsZero()) {
				LineAdditional(matrix, i, k, j);
				addicted = true;
			}
		}
		
		if(addicted) PrintMatrix(matrix);
	}

	if (CheckFalseEqualities(matrix)) {
		std::cout << "The system have no solution\n\n";
		m_result.clear();
		return;
	}

	if (DeleteZeroLines(matrix)) PrintMatrix(matrix);
	Finishing(matrix);
}
//-----------------------------------------------------------------------------------
void JordanGauss::LineAdditional(FractionMatrix& matrix, int firstLine, int secondLine, int column) {
	Fraction::SimpleFraction coef = matrix[secondLine][column] * Fraction::SimpleFraction(-1, 1);

	for (int i = 0; i < m_variablesCount + 1; i++) {
		matrix[secondLine][i] += matrix[firstLine][i] * coef;
	}

}
//-----------------------------------------------------------------------------------
void JordanGauss::LineDivision(FractionMatrix& matrix, int lineNum, int dividerInd) {
	Fraction::SimpleFraction divider = matrix[lineNum][dividerInd];
	for (auto& elem : matrix[lineNum]) {
		elem /= divider;
	}
}
//-----------------------------------------------------------------------------------
bool JordanGauss::DeleteZeroLines(FractionMatrix& matrix) {
	std::vector<int> deleteInds;
	for (int i = 0; i < matrix.size(); i++) {
		int zeroesCount = std::count_if(matrix[i].begin(), matrix[i].end(),
			[](Fraction::SimpleFraction elem) {return elem.IsZero(); });

		if (zeroesCount == m_variablesCount + 1) {
			deleteInds.push_back(i);
		}
	}

	if (deleteInds.empty())
		return false;

	std::reverse(deleteInds.begin(), deleteInds.end());

	for (const auto& ind : deleteInds) {
		matrix.erase(matrix.begin() + ind);
	}

	return true;
}
//-----------------------------------------------------------------------------------
bool JordanGauss::CheckFalseEqualities(FractionMatrix& matrix) {
	int zeroesCount = 0;
	for (const auto& line : matrix) {
		for (int i = 0; i < line.size(); i++) {
			line[i].IsZero() ? zeroesCount++ : zeroesCount;
		}
		
		if (zeroesCount == m_variablesCount && !line.back().IsZero()) {
			return true;
		}
		zeroesCount = 0;
	}
	return false;
}
//-----------------------------------------------------------------------------------
void JordanGauss::Finishing(FractionMatrix& matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		int zeroesCount = std::count_if(matrix[i].begin(), matrix[i].end(),
			[](Fraction::SimpleFraction elem) {return elem.IsZero(); });

		if (zeroesCount == m_variablesCount - 1 && !matrix[i].back().IsZero() ||
			zeroesCount == m_variablesCount && matrix[i].back().IsZero()) {
			int variableNum = std::distance(matrix[i].begin(), std::find_if(matrix[i].begin(), matrix[i].end(),
				[](Fraction::SimpleFraction elem) {return elem.IsOne(); }));
			std::cout << "x" << variableNum << " = " << matrix[i].back();
			m_result.push_back(matrix[i].back());
		}
		else {
			int oneInd = -1;
			std::vector<std::pair<int, Fraction::SimpleFraction>> elems;

			for (int j = 0; j < m_variablesCount + 1; j++) {
				if (!matrix[i][j].IsZero()) {
					elems.push_back(std::make_pair(j + 1, matrix[i][j]));

					if (matrix[i][j].IsOne() && j != m_variablesCount) {
						bool isBasic = true;

						for (int k = 0; k < matrix.size(); k++) {
							if (!matrix[k][j].IsZero() && k != i) isBasic = false;
						}

						if (isBasic) {
							oneInd = j + 1;
						}
					}
				}
			}

			std::map<int, Fraction::SimpleFraction> basicVariable;
			basicVariable[0] = Fraction::SimpleFraction(0, 1);

			if (!matrix[i].back().IsZero()) {
				std::cout << "x" << oneInd << " = " << elems.back().second;
				basicVariable[0] = elems.back().second;
			}
			else {
				std::cout << "x" << oneInd << " = " << elems.back().second * Fraction::SimpleFraction(-1, 1) << "x" << elems.back().first;
				basicVariable[elems.back().first] = elems.back().second * Fraction::SimpleFraction(-1, 1);
			}

			for (int j = 0; j < elems.size() - 1; j++) {
				if (elems[j].first == oneInd) continue;
				std::cout << " + " << elems[j].second * Fraction::SimpleFraction(-1, 1) << "x" << elems[j].first;
				basicVariable[elems[j].first] = elems[j].second * Fraction::SimpleFraction(-1, 1);
			}

			m_basicVariables[oneInd] = basicVariable;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
//-----------------------------------------------------------------------------------
bool JordanGauss::FindAndSwapMaxColumnElem(FractionMatrix& matrix, int lineNum, int columnNum) {
	int swapLineInd = -1;
	Fraction::SimpleFraction max = matrix[lineNum][columnNum].ToAbs();

	for (int i = lineNum + 1; i < matrix.size(); i++) {
		if (max < matrix[i][columnNum].ToAbs()) {
			max = matrix[i][columnNum].ToAbs();
			swapLineInd = i;
		}
	}

	if (swapLineInd != -1) {
		std::swap(matrix[lineNum], matrix[swapLineInd]);
		return true;
	}

	return false;
}
//-----------------------------------------------------------------------------------
FractionMatrix JordanGauss::ReadMatrixFromFile(std::string fileName) {
	std::fstream file(fileName);

	if (!file.is_open()) {
		std::cout << "Error opening the file" << std::endl;
	}

	std::string line = "";

	while (getline(file, line)) {
		m_matrix.push_back(std::vector<Fraction::SimpleFraction>());
		while (line.size() != 0) {
			int spaceInd = line.find(" ");

			m_matrix.back().push_back(Fraction::SimpleFraction(std::stoi(line.substr(0, spaceInd)), 1));

			line.erase(0, spaceInd + 1);
			if (spaceInd == -1)
				break;
		}
	}

	m_variablesCount = m_matrix[0].size() - 1;
	return m_matrix;
}
//-----------------------------------------------------------------------------------
void JordanGauss::PrintMatrix(FractionMatrix& matrix) {
	for (const auto& line : matrix) {
		for (const auto& column : line) {
			std::cout << column << "\t";
		}
		std::cout << "\n";
	}
	std::cout << std::endl;
}
//-----------------------------------------------------------------------------------
void JordanGauss::FindBasicSolutions() {
	FractionMatrix originalMatrix = m_matrix;
	Process(m_matrix);

	std::vector<int> basicVariablesInds;
	std::vector<std::vector<int>> variableInclusions;
	FractionMatrix result;

	for (int i = 0, j = 0; i < m_matrix.size() && j < m_variablesCount; i++, j++) {
		if (m_matrix[i][j].IsOne()) {
			basicVariablesInds.push_back(j);
			continue;
		}

		i--;
	}

	int variablesCount = m_variablesCount;
	variableInclusions = CombinationsWithoutRepeats(m_variablesCount, m_matrix.size());
	m_variablesCount = m_matrix.size();

	std::cout << "\nCombinations:\n";
	for (auto& inclusion : variableInclusions) {
		std::cout << '(';
		for (int i = 0; i < inclusion.size(); i++) {
			std::cout << inclusion[i];
			if (i < inclusion.size() - 1) std::cout << ", ";
		}
		std::cout << ')' << std::endl;
	}
	std::cout << std::endl;

	for (auto& variables : variableInclusions) {
		FractionMatrix matrix = ExtractMatrixWithIncludedVariables(variables);
		Process(matrix);

		auto a = m_result;
		if (!m_result.empty()) {
			std::vector<Fraction::SimpleFraction> basicSolution;
			for (int i = 0; i < variablesCount; i++) {
				if (std::find(variables.begin(), variables.end(), i) == variables.end()) {
					basicSolution.push_back(Fraction::SimpleFraction(0, 1));
				}
				else {
					basicSolution.push_back(m_result.front());
					m_result.erase(m_result.begin());
				}
			}
			result.push_back(basicSolution);
		}
	}

	for (auto& basicSolution : result) {
		std::cout << '(';
		for (int i = 0; i < basicSolution.size(); i++) {
			std::cout << basicSolution[i];
			if (i < basicSolution.size() - 1) std::cout << ", ";
		}
		std::cout << ')' << std::endl;
	}

	m_variablesCount = variablesCount;
	//SolutionVerification(originalMatrix, result);
}
//-----------------------------------------------------------------------------------
FractionMatrix JordanGauss::ExtractMatrixWithIncludedVariables(std::vector<int>& columnInds) {
	FractionMatrix result;

	for (int i = 0; i < m_matrix.size(); i++) {
		std::vector<Fraction::SimpleFraction> tempLine;
		for (const auto& ind : columnInds) {
			tempLine.push_back(m_matrix[i][ind]);
		}
		tempLine.push_back(m_matrix[i].back());
		result.push_back(tempLine);
	}

	return result;
}
//-----------------------------------------------------------------------------------
std::vector<std::vector<int>> JordanGauss::CombinationsWithoutRepeats(int n, int k) {
	std::string bitmask(k, 1);
	bitmask.resize(n, 0);

	std::vector<int> combination;
	std::vector<std::vector<int>> result;

	do {
		for (int i = 0; i < n; i++) {
			if (bitmask[i]) combination.push_back(i);
		}
		result.push_back(combination);
		combination.clear();
	} while (prev_permutation(bitmask.begin(), bitmask.end()));

	return result;
}
//-----------------------------------------------------------------------------------
void JordanGauss::SolutionVerification(FractionMatrix& matrix, FractionMatrix& solutions) {
	std::cout << "\nVerification\n\n";
	for (auto& solution : solutions) {
		std::cout << '(';
		for (int i = 0; i < solution.size(); i++) {
			std::cout << solution[i];
			if (i < solution.size() - 1) std::cout << ", ";
		}
		std::cout << ')' << std::endl;

		for (int i = 0; i < matrix.size(); i++) {
			Fraction::SimpleFraction res = Fraction::SimpleFraction(0, 1);
			for (int j = 0; j < m_variablesCount; j++) {
				res += matrix[i][j] * solution[j];
			}
			std::cout << res << " = " << matrix[i].back() << (res == matrix[i].back() ? " Passed" : " Failed") << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
//-----------------------------------------------------------------------------------
}