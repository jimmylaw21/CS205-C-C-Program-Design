//MathFunc.cpp: Calculator类的变量表达式实现文件

#include "Calculator.h"
#include "BigDecimal.h"
#include <cstring>

// 变量表达式计算
// is there any unknown variable in the expression. If , then cerr << a message.
void Calculator::checkUnknownVar(string expression) {
  int indexOfAssign = expression.find('=');
  string varValue;
  if (indexOfAssign != -1) {
    varValue = expression.substr(indexOfAssign + 1);
  } else {
    varValue = expression;
  }
  for (int i = 0; i < varValue.length(); ++i) {
    if (!(isdigit(varValue[i]) || varValue[i] == '.' || varValue[i] == '+' ||
          varValue[i] == '-' || varValue[i] == '*' || varValue[i] == '/' ||
          varValue[i] == '^' || varValue[i] == '!' || varValue[i] == '(' ||
          varValue[i] == ')')) {
      cerr << "The variable does not exist!" << endl;
    }
  }
}

// to check if the expression is for assigning
bool Calculator::isAssigning(string expression) {
  return expression.find('=') != string::npos;
}

bool Calculator::isVariableNameValid(string expression) {
  int indexOfAssign = expression.find('=');
  string varName = expression.substr(0, indexOfAssign);
  if (!(isalpha(varName[0]) || varName[0] == '_')) {
    return false;
  }
  for (int i = 1; i < varName.length(); ++i) {
    if (!(isdigit(varName[i]) || isalpha(varName[i]) || varName[i] == '_')) {
      return false;
      break;
    }
  }
  if(varName == "PI" || varName == "E" || varName == "ans"){
    return false;
  }
  return true;
}

// to add some variables Map.
void Calculator::addVariable(string expression) {
  int indexOfAssign = expression.find('=');
  string varName = expression.substr(0, indexOfAssign);
  string varValue = expression.substr(indexOfAssign + 1);
  // valStr=replaceAllMathFunctions(valStr);
  if (isVariableNameValid(expression)) {
    variables[varName] = varValue;
  } else {
    cerr << "Invalid variable name." << endl;
  }
}

// to replace the variables in the expression with their values.
string Calculator::replaceVar(string expression) {
  int indexOfAssign = expression.find('=');
  if (indexOfAssign != -1) {
    string varName = expression.substr(0, indexOfAssign);
    string varValue = expression.substr(indexOfAssign + 1);
    for (auto it = variables.begin(); it != variables.end(); ++it) {
      varValue = replaceAll(varValue, it->first, it->second);
    }
    expression = varName + "=" + varValue;
  } else {
    for (auto it = variables.begin(); it != variables.end(); ++it) {
      expression = replaceAll(expression, it->first, it->second);
    }
  }
  return expression;
}

string Calculator::replaceAll(string str, string oldStr, string newStr) {
  int pos = 0;
  while ((pos = str.find(oldStr, pos)) != string::npos) {
    str.replace(pos, oldStr.length(), newStr);
    pos += newStr.length();
  }
  return str;
}
void Calculator::pretreat(string expression) {
  // 实现正负数
  // for (int i = 0; i < expression.length(); i++) {
  // //string下标调用运算符时可能会导致类型溢出
  for (size_t i = 0; i < expression.size();
       i++) { // string.size()返回size_type类型，避免下标运算时的类型溢出
    if (expression[i] == '-' ||
        expression[i] == '+') { //-x转换为0-x，+x转化为0+x
      if (i == 0) {
        expression.insert(0, 1, '0');
      } else if (stdInfix[i - 1] == '(') {
        expression.insert(i, 1, '0');
      }
    }
  }
}
void Calculator::Assign() {
  if (isAssigning(infix)) {
    addVariable(infix);
  }
  pretreat(infix);
  // infix = replaceMathFunctions(infix);
  infix = replaceVar(infix);
  pretreat(infix);
  // checkUnknownVar(infix);
}



vector<string> Calculator::split(const std::string &strIn, char delim) {
  char *str = const_cast<char *>(strIn.c_str());
  std::string s;
  s.append(1, delim);
  std::vector<std::string> elems;
  char *splitted = strtok(str, s.c_str());
  while (splitted != NULL) {
    elems.push_back(std::string(splitted));
    splitted = strtok(NULL, s.c_str());
  }
  return elems;
}