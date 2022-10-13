//MathFunc.cpp: Calculator类的数学函数
#include "Calculator.h"
#include "BigDecimal.h"

bool Calculator::containsMathFunc(string expression, string funcName) {
  return expression.find(funcName) != -1;
}

string Calculator::replaceMathFunc(string expression) {
  expression = getSqrt(expression);
  expression = getMax(expression);
  expression = getMin(expression);
  expression = getRandom(expression);
  expression = getAbs(expression);
  expression = getSin(expression);
  expression = getCos(expression);
  expression = getTan(expression);
  expression = getCeil(expression);
  expression = getFloor(expression);
  expression = getLog(expression);
  return expression;
}

// 处理表达式中的sqrt数学函数
string Calculator::getSqrt(string expression) {
  while (containsMathFunc(expression, "sqrt")) {
    int indexOfSqrt = expression.find("sqrt");
    int indexOfLeftBracket = expression.find('(', indexOfSqrt);
    int indexOfRightBracket = expression.find(')', indexOfSqrt);
    string sqrtStr = expression.substr(indexOfSqrt, indexOfRightBracket + 1);
    string sqrtValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double sqrtResult = sqrt(stod(sqrtValueStr));
    expression = replaceAll(expression, sqrtStr, to_string(sqrtResult));
  }
  return expression;
}

// 处理表达式中的max函数
string Calculator::getMax(string expression) {
  while (containsMathFunc(expression, "max")) {
    int indexOfMax = expression.find("max");
    int indexOfLeftBracket = expression.find('(', indexOfMax);
    int indexOfRightBracket = expression.find(')', indexOfMax);
    string maxStr = expression.substr(indexOfMax, indexOfRightBracket + 1);
    string maxValuesStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    vector<string> maxValues = split(maxValuesStr, ',');
    double maxResult = stod(maxValues[0]);
    for (int i = 1; i < maxValues.size(); ++i) {
      if (stod(maxValues[i]) > maxResult) {
        maxResult = stod(maxValues[i]);
      }
    }
    expression = replaceAll(expression, maxStr, to_string(maxResult));
  }
  return expression;
}
// 处理表达式中的min函数
string Calculator::getMin(string expression) {
  while (containsMathFunc(expression, "min")) {
    int indexOfMin = expression.find("min");
    int indexOfLeftBracket = expression.find('(', indexOfMin);
    int indexOfRightBracket = expression.find(')', indexOfMin);
    string minStr = expression.substr(indexOfMin, indexOfRightBracket + 1);
    string minValuesStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    vector<string> minValues = split(minValuesStr, ',');
    double minResult = stod(minValues[0]);
    for (int i = 1; i < minValues.size(); ++i) {
      if (stod(minValues[i]) < minResult) {
        minResult = stod(minValues[i]);
      }
    }
    expression = replaceAll(expression, minStr, to_string(minResult));
  }
  return expression;
}
// 处理表达式中的random函数,使用rand()函数
string Calculator::getRandom(string expression) {
  while (containsMathFunc(expression, "random")) {
    int indexOfRandom = expression.find("random");
    int indexOfLeftBracket = expression.find('(', indexOfRandom);
    int indexOfRightBracket = expression.find(')', indexOfRandom);
    string randomStr =
        expression.substr(indexOfRandom, indexOfRightBracket + 1);
    string randomValuesStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    vector<string> randomValues = split(randomValuesStr, ',');
    int randomResult =
        rand() % (stoi(randomValues[1]) - stoi(randomValues[0]) + 1) +
        stoi(randomValues[0]);
    expression = replaceAll(expression, randomStr, to_string(randomResult));
  }
  return expression;
}
// 处理表达式中的abs函数
string Calculator::getAbs(string expression){
  while (containsMathFunc(expression, "abs")) {
    int indexOfAbs = expression.find("abs");
    int indexOfLeftBracket = expression.find('(', indexOfAbs);
    int indexOfRightBracket = expression.find(')', indexOfAbs);
    string absStr = expression.substr(indexOfAbs, indexOfRightBracket + 1);
    string absValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double absResult = abs(stod(absValueStr));
    expression = replaceAll(expression, absStr, to_string(absResult));
  }
  return expression;
}
// 处理表达式中的exp函数
string Calculator::getExp(string expression) {
  while (containsMathFunc(expression, "exp")) {
    int indexOfExp = expression.find("exp");
    int indexOfLeftBracket = expression.find('(', indexOfExp);
    int indexOfRightBracket = expression.find(')', indexOfExp);
    string expStr = expression.substr(indexOfExp, indexOfRightBracket + 1);
    string expValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double expResult = exp(stod(expValueStr));
    expression = replaceAll(expression, expStr, to_string(expResult));
  }
  return expression;
}
// 处理表达式中的sin函数
string Calculator::getSin(string expression) {
  while (containsMathFunc(expression, "sin")) {
    int indexOfSin = expression.find("sin");
    int indexOfLeftBracket = expression.find('(', indexOfSin);
    int indexOfRightBracket = expression.find(')', indexOfSin);
    string sinStr = expression.substr(indexOfSin, indexOfRightBracket + 1);
    string sinValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double sinResult = sin(stod(sinValueStr));
    expression = replaceAll(expression, sinStr, to_string(sinResult));
  }
  return expression;
} 
// 处理表达式中的cos函数
string Calculator::getCos(string expression) {
  while (containsMathFunc(expression, "cos")) {
    int indexOfCos = expression.find("cos");
    int indexOfLeftBracket = expression.find('(', indexOfCos);
    int indexOfRightBracket = expression.find(')', indexOfCos);
    string cosStr = expression.substr(indexOfCos, indexOfRightBracket + 1);
    string cosValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double cosResult = cos(stod(cosValueStr));
    expression = replaceAll(expression, cosStr, to_string(cosResult));
  }
  return expression;
} 
// 处理表达式中的tan函数
string Calculator::getTan(string expression) {
  while (containsMathFunc(expression, "tan")) {
    int indexOfTan = expression.find("tan");
    int indexOfLeftBracket = expression.find('(', indexOfTan);
    int indexOfRightBracket = expression.find(')', indexOfTan);
    string tanStr = expression.substr(indexOfTan, indexOfRightBracket + 1);
    string tanValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double tanResult = tan(stod(tanValueStr));
    expression = replaceAll(expression, tanStr, to_string(tanResult));
  }
  return expression;
} 
// 处理表达式中的ceil函数
string Calculator::getCeil(string expression) {
  while (containsMathFunc(expression, "ceil")) {
    int indexOfCeil = expression.find("ceil");
    int indexOfLeftBracket = expression.find('(', indexOfCeil);
    int indexOfRightBracket = expression.find(')', indexOfCeil);
    string ceilStr = expression.substr(indexOfCeil, indexOfRightBracket + 1);
    string ceilValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double ceilResult = ceil(stod(ceilValueStr));
    expression = replaceAll(expression, ceilStr, to_string(ceilResult));
  }
  return expression;
} 
// 处理表达式中的floor函数
string Calculator::getFloor(string expression) {
  while (containsMathFunc(expression, "floor")) {
    int indexOfFloor = expression.find("floor");
    int indexOfLeftBracket = expression.find('(', indexOfFloor);
    int indexOfRightBracket = expression.find(')', indexOfFloor);
    string floorStr = expression.substr(indexOfFloor, indexOfRightBracket + 1);
    string floorValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double floorResult = floor(stod(floorValueStr));
    expression = replaceAll(expression, floorStr, to_string(floorResult));
  }
  return expression;
} 
// 处理表达式中的log函数
string Calculator::getLog(string expression) {
  while (containsMathFunc(expression, "log")) {
    int indexOfLog = expression.find("log");
    int indexOfLeftBracket = expression.find('(', indexOfLog);
    int indexOfRightBracket = expression.find(')', indexOfLog);
    string logStr = expression.substr(indexOfLog, indexOfRightBracket + 1);
    string logValueStr = expression.substr(
        indexOfLeftBracket + 1, indexOfRightBracket - indexOfLeftBracket - 1);
    double logResult = log(stod(logValueStr));
    expression = replaceAll(expression, logStr, to_string(logResult));
  }
  return expression;
} 