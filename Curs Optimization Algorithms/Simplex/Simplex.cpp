#include "Simplex.h"

namespace LinearEquations {
//-----------------------------------------------------------------------------------
Simplex::Simplex() {

}
//-----------------------------------------------------------------------------------
void Simplex::Calculate() {
	ReadMatrixFromFile();
	m_jordanGauss.SetVariablesCount(m_variablesCount);
	m_jordanGauss.Process(m_matrix);

	TBasicVariables basic = m_jordanGauss.GetBasicVariables();
	std::map<int, Fraction::SimpleFraction> emptyVariable;

	for (int i = 0; i < m_variablesCount + 1; i++) {
		emptyVariable[i] = Fraction::SimpleFraction(0, 1);
	}

	for (const auto& basicVariable : basic) {
		m_basicVariables[basicVariable.first] = emptyVariable;
		m_basicVariables[basicVariable.first][basicVariable.first] = Fraction::SimpleFraction(-1, 1);

		for (const auto& freeVariable : basicVariable.second) {
			m_basicVariables[basicVariable.first][freeVariable.first] = freeVariable.second;
		}
	}

	FunctionExpressing();
	CreateSourceTable();
	PrintSimplexTable();

	while (1) {
		SolutionIsOptimal();
		TableIsValid();

		if (!m_solutionIsOptimal) {
			std::cout << "\nSolution is not optimal.";
		}
		else if(m_existNextSolution && m_solutionsCount < 2) {
			std::cout << "\nSolution is optimal, but not one.";
			Finishing();
		}
		else if (m_solutionsCount == 2) {
			std::vector<std::pair<Fraction::SimpleFraction, Fraction::SimpleFraction>> lambdaResult;

			for (int i = 0; i < m_solutions[0].size(); i++) {
				lambdaResult.push_back(std::make_pair(
					m_solutions[0][i], m_solutions[0][i] * Fraction::SimpleFraction(-1, 1) + m_solutions[1][i]
				)
				);
			}

			std::cout << "\nZmax = Z(";
			for (int i = 0; i < lambdaResult.size(); i++) {
				std::cout << lambdaResult[i].first << " + " << lambdaResult[i].second << "L";
				if (i != lambdaResult.size() - 1) std::cout << ", ";
			}
			std::cout << ") = " << m_solutionValue << std::endl;
			std::cout << "0 <= L <= 1" << std::endl;

			break;
		}
		else {
			std::cout << "\nSolution is optimal, and it is only one.";
			Finishing();
			break;
		}

		if ((!m_tableIsValid || m_haveNoValidStr) && !m_existNextSolution) {
			std::cout << "\nTable is invalid, but solution not found, exit.\n";
			break;
		}
		else {
			std::cout << "\nTable is valid.\n\n";
		}

		RecalculateTable();
		PrintSimplexTable();
	}
}
//-----------------------------------------------------------------------------------
FractionMatrix Simplex::ReadMatrixFromFile(std::string fileName) {
	std::fstream file(fileName);

	if (!file.is_open()) {
		std::cout << "Error opening the file" << std::endl;
	}

	std::string line = "";

	m_function.push_back(Fraction::SimpleFraction(0, 1));
	getline(file, line);
	while (line.size() != 0) {
		int spaceInd = line.find(" ");

		m_function.push_back(Fraction::SimpleFraction(std::stoi(line.substr(0, spaceInd)), 1));

		line.erase(0, spaceInd + 1);
		if (spaceInd == -1)
			break;
	}

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
	m_simplexDivisions = std::vector<Fraction::SimpleFraction>(m_variablesCount + 1, Fraction::SimpleFraction(LLONG_MIN, 1));

	return m_matrix;
}
//-----------------------------------------------------------------------------------
void Simplex::PrintMatrix(FractionMatrix& matrix) {
	for (const auto& line : matrix) {
		for (const auto& column : line) {
			std::cout << column << "\t";
		}
		std::cout << "\n";
	}
	std::cout << std::endl;
}
//-----------------------------------------------------------------------------------
void Simplex::PrintSimplexTable() {
	std::string tableHeader("B.V.\t1\t");
	for (int i = 0; i < m_variablesCount; i++) {
		tableHeader += "x" + std::to_string(i + 1) + "\t";
	}

	std::cout << tableHeader << std::endl;
	for (const auto& stroke : m_simplexTable) {
		std::cout << stroke.first << "\t";
		for (const auto& variable : stroke.second) {
			std::cout << variable << "\t";
		}
		std::cout << std::endl;
	}

	/*std::cout << "S.D.\t";
	for (auto& simplexDivision : m_simplexDivisions) {
		std::cout << (simplexDivision.IsMin() ? "-" : simplexDivision.ToString()) << "\t";
	}
	std::cout << std::endl;*/
}
//-----------------------------------------------------------------------------------
void Simplex::FunctionExpressing() {
	std::cout << "Z = ";
	int notZeroCount = std::count_if(m_function.begin(), m_function.end(), 
		[&](Fraction::SimpleFraction& elem) {return elem != Fraction::SimpleFraction(0, 1); });

	for (int i = 1; i < m_function.size(); i++) {
		if (!m_function[i].IsZero()) {
			notZeroCount--;
			std::cout << m_function[i] << "x" << i << (notZeroCount != 0 ? " + " : "");
		}
	}
	std::cout << std::endl;

	for (int i = 0; i < m_variablesCount + 1; i++) {
		m_expressedFunction[i] = m_basicVariables.find(i) != m_basicVariables.end() ? Fraction::SimpleFraction(0, 1) : m_function[i];
	}

	for (auto& it : m_jordanGauss.GetBasicVariables()) {
		for (auto& freeVariable : it.second) {
			//m_expressedFunction[freeVariable.first] += freeVariable.second * (freeVariable.first == 0 ? Fraction::SimpleFraction() : m_function[it.first]);
			m_expressedFunction[freeVariable.first] += freeVariable.second * m_function[it.first];
		}
	}

	std::cout << "Exp. Z = ";
	if (!m_expressedFunction[0].IsZero()) {
		std::cout << m_expressedFunction[0] << " + ";
	}
	
	using pairType = decltype(m_expressedFunction)::value_type;
	notZeroCount = std::count_if(m_expressedFunction.begin(), m_expressedFunction.end(),
		[&](pairType& elem) {return elem.second != Fraction::SimpleFraction(0, 1); });
	for (auto& variable : m_expressedFunction) {
		if (variable.first == 0) continue;

		if (!variable.second.IsZero()) {
			notZeroCount--;
			std::cout << variable.second << "x" << variable.first << (notZeroCount != 0 ? " + " : "");
		}
	}
	std::cout << std::endl << std::endl;
}
//-----------------------------------------------------------------------------------
void Simplex::CreateSourceTable() {
	for (auto& basicVariable : m_basicVariables) {
		std::string xName = "x" + std::to_string(basicVariable.first);
		m_simplexTable[xName] = std::vector<Fraction::SimpleFraction>();

		for (auto& freeVariable : basicVariable.second) {
			auto coef = freeVariable.first != 0 ? freeVariable.second * Fraction::SimpleFraction(-1, 1)  : freeVariable.second;
			m_simplexTable[xName].push_back(coef);
		}
	}

	m_simplexTable["z"] = std::vector<Fraction::SimpleFraction>();
	for (auto& variable : m_expressedFunction) {
		m_simplexTable["z"].push_back(variable.first == 0 ? variable.second : variable.second * Fraction::SimpleFraction(-1, 1));
	}

	m_prevSimplexTable = m_simplexTable;
}
//-----------------------------------------------------------------------------------
void Simplex::FindResolvingElement(std::string& row, int& column) {
	using pairType = decltype(m_simplexTable)::value_type;

	auto lastStroke = --m_simplexTable.end();
	std::vector<Fraction::SimpleFraction> oneColumn;

	for (auto& it : m_simplexTable) {
		if (it == *lastStroke) break;

		oneColumn.push_back(it.second[0]);
	}

	std::sort(oneColumn.begin(), oneColumn.end());
	Fraction::SimpleFraction strElem;

	for (auto& it : oneColumn) {
		if ((it < Fraction::SimpleFraction(0, 1) || it.IsZero())) {
			auto stroke =
				*std::find_if(m_simplexTable.begin(), --m_simplexTable.end(), [&](pairType p) { return p.second[0] == it; });
			if (1 < std::count_if(stroke.second.begin(), stroke.second.end(), [&](Fraction::SimpleFraction fr) {return fr < Fraction::SimpleFraction(0, 1); })) {
				row = stroke.first;
				break;
			}
		}
	}

	if (row == "") return;

	/*row = std::min_element
	(
		m_simplexTable.begin(), 
		--m_simplexTable.end(),
		[&](pairType& p1, pairType& p2) { return p1.second[0] < p2.second[0]; }
	)->first;*/
	//&& 1 < std::count_if(p1.second.begin(), p1.second.end(), [&](Fraction::SimpleFraction fr) {return fr < Fraction::SimpleFraction(0, 1); });

	auto maxElem = Fraction::SimpleFraction(LLONG_MIN, 1);
	for (int i = 1; i < m_simplexTable["z"].size(); i++) {
		if (m_simplexTable[row][i] > Fraction::SimpleFraction(0, 1) || m_simplexTable[row][i].IsZero() || m_simplexTable["z"][i].IsZero()) {
			m_simplexDivisions[i] = Fraction::SimpleFraction(LLONG_MIN, 1);
			continue;
		}

		m_simplexDivisions[i] = m_simplexTable["z"][i] / m_simplexTable[row][i];
		maxElem = m_simplexDivisions[i];
	}

	for (auto& elem : m_simplexDivisions) {
		if (elem.IsMin()) continue;
		if (maxElem < elem) maxElem = elem;
	}

	for (column = 1; column < m_simplexTable[row].size(); column++) {
		if (m_simplexDivisions[column] == maxElem) {
			break;
		}
	}

	std::cout << "\nS.D.\t";
	for (auto& simplexDivision : m_simplexDivisions) {
		std::cout << (simplexDivision.IsMin() ? "-" : simplexDivision.ToString()) << "\t";
	}
	std::cout << std::endl;

	std::cout << "Resolving element: row - " << row << " element - " << m_simplexTable[row][column] << std::endl << std::endl;
}
//-----------------------------------------------------------------------------------
void Simplex::RecalculateTable() {
	std::string row = "";
	int column = -1;
	m_prevSimplexTable = m_simplexTable;

	if (!m_existNextSolution)
		FindResolvingElement(row, column);
	else
		FindNextSolutionResElem(row, column);

	if (row == "" || column == -1) { 
		m_haveNoValidStr = true;
		return; 
	}

	m_simplexTable.erase(row);
	std::vector<Fraction::SimpleFraction> emptyVector(m_prevSimplexTable["z"].size(), Fraction::SimpleFraction(LLONG_MAX, 1));

	for (auto& stroke : m_simplexTable) {
		stroke.second = emptyVector;
	}

	std::string newBasicVarName = "x" + std::to_string(column);
	m_simplexTable[newBasicVarName] = m_prevSimplexTable[row];

	for (auto& var : m_simplexTable[newBasicVarName]) {
		var /= m_prevSimplexTable[row][column];
	}

	std::vector<int> basicInds;
	for (auto& stroke : m_simplexTable) {
		if (stroke.first != "z") {
			basicInds.push_back(std::stoi(stroke.first.substr(1)));
		}
	}

	for (const auto& basicInd : basicInds) {
		m_simplexTable["z"][basicInd] = Fraction::SimpleFraction(0, 1);
		for (auto& stroke : m_simplexTable) {
			if (stroke.first == "z") continue;
			stroke.second[basicInd] = std::stoi(stroke.first.substr(1)) == basicInd ? Fraction::SimpleFraction() : Fraction::SimpleFraction(0, 1);
		}
	}

	FillNullCells(row, column);
}
//-----------------------------------------------------------------------------------
void Simplex::FillNullCells(std::string& mainRow, int& mainColumn) {
	for (auto& stroke : m_simplexTable) {
		for (int i = 0; i < stroke.second.size(); i++) {
			if (stroke.second[i].IsMax()) {
				stroke.second[i] = m_prevSimplexTable[stroke.first][i] - m_prevSimplexTable[stroke.first][mainColumn] *
					m_prevSimplexTable[mainRow][i] / m_prevSimplexTable[mainRow][mainColumn];
			}
		}
	}
}
//-----------------------------------------------------------------------------------
bool Simplex::TableIsValid() {
	auto lastStroke = --m_simplexTable.end();
	auto zFreeCoef = (*lastStroke).second.begin();
	bool freeVarsIsValid = false;

	for (auto& stroke : m_simplexTable) {
		if (stroke == *lastStroke) break;

		if (stroke.second[0] < Fraction::SimpleFraction(0, 1)) {
			freeVarsIsValid = true;
			break;
		}
	}

	if (m_existNextSolution) freeVarsIsValid = true;
	m_tableIsValid = freeVarsIsValid;
	if (!freeVarsIsValid) return false;

	bool zCoefsIsValid = true;
	for (auto& zCoef : m_simplexTable["z"]) {
		if (zCoef == *zFreeCoef) continue;

		if (zCoef < Fraction::SimpleFraction(0, 1)) {
			zCoefsIsValid = false;
			break;
		}
	}

	m_tableIsValid &= zCoefsIsValid;
	return zCoefsIsValid;
}
//-----------------------------------------------------------------------------------
bool Simplex::SolutionIsOptimal() {
	auto lastStroke = --m_simplexTable.end();
	auto zFreeCoef = (*lastStroke).second.begin();
	bool freeVarsIsOptimal = true;
	bool zCoefsIsOptimal = true;

	for (auto& stroke : m_simplexTable) {
		if (stroke == *lastStroke) break;

		if (stroke.second[0] < Fraction::SimpleFraction(0, 1)) {
			freeVarsIsOptimal = false;
			break;
		}
	}

	m_solutionIsOptimal = freeVarsIsOptimal;
	if (!freeVarsIsOptimal) return false;

	std::vector<int> basicInds;
	for (auto& stroke : m_simplexTable) {
		if (stroke.first != "z") {
			basicInds.push_back(std::stoi(stroke.first.substr(1)));
		}
	}

	for (int i = 1; i < m_simplexTable["z"].size(); i++) {
		if (m_simplexTable["z"][i].IsZero() && std::find(basicInds.begin(), basicInds.end(), i) == basicInds.end()) {
			m_existNextSolution = true;
			m_nextSolutionZInd = i;
		}

		if (m_simplexTable["z"][i] < Fraction::SimpleFraction(0, 1)) {
			zCoefsIsOptimal = false;
			break;
		}
	}

	m_solutionIsOptimal &= zCoefsIsOptimal;
	m_existNextSolution &= zCoefsIsOptimal;

	return zCoefsIsOptimal;
}
//-----------------------------------------------------------------------------------
void Simplex::Finishing() {
	std::vector<Fraction::SimpleFraction> basicVars(m_variablesCount, Fraction::SimpleFraction(0, 1));
	for (auto& stroke : m_simplexTable) {
		if (stroke.first != "z") {
			basicVars[std::stoi(stroke.first.substr(1)) - 1] = stroke.second[0];
		}
	}

	std::cout << "\nZmax = Z(";
	for (int i = 0; i < basicVars.size(); i++) {
		std::cout << basicVars[i];
		if (i != basicVars.size() - 1) std::cout << ", ";
	}
	std::cout << ") = ";

	Fraction::SimpleFraction sum = m_expressedFunction[0];
	for (int i = 1; i < m_expressedFunction.size(); i++) {
		sum += m_expressedFunction[i] * basicVars[i - 1];
	}
	std::cout << sum << std::endl;

	m_solutionValue = sum;
	m_solutions.push_back(basicVars);
	m_solutionsCount++;
}
//-----------------------------------------------------------------------------------
void Simplex::FindNextSolutionResElem(std::string& row, int& column) {
	column = m_nextSolutionZInd;
	Fraction::SimpleFraction minDiv(LLONG_MAX, 1);
	m_simplexDivisions.clear();
	auto lastStroke = --m_simplexTable.end();

	for (auto& stroke : m_simplexTable) {
		if (stroke == *lastStroke) break;

		m_simplexDivisions.push_back(
			stroke.second[m_nextSolutionZInd] > Fraction::SimpleFraction(0, 1) ? 
			stroke.second[0] / stroke.second[m_nextSolutionZInd] : Fraction::SimpleFraction(LLONG_MAX, 1)
		);
	}

	auto minElem = std::min_element(m_simplexDivisions.begin(), m_simplexDivisions.end());
	int minDivInd = minElem - m_simplexDivisions.begin();

	auto it = m_simplexTable.begin();

	int i = 0;
	for (auto& stroke : m_simplexTable) {
		if (i == minDivInd) {
			row = stroke.first;
			break;
		}
		i++;
	}
}
//-----------------------------------------------------------------------------------
}