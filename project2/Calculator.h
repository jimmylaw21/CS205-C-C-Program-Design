// Calculator.h: 头文件
#include "BigDecimal.h"
#include <cmath>
#include <map>
#include <stack>
#include <string>
#include <vector>
using namespace std;

// 计算器类
class Calculator {
public:
  Calculator(); // 构造函数

  // 基本表达式计算
  void getFormat();       // 表达式自定义标准格式化
  int getPrior(char c);   // 获取算术符号优先级
  void getPostfix();      // 后缀表达式转换
  void calResult();       // 计算后缀表达式
  void calculate();       // 计算方法
  BigDecimal getResult(); // 获取结果

  // 变量表达式计算
  bool isAssigning(string expression);         // 检查是否为赋值表达式
  bool isVariableNameValid(string expression); // 检查变量名是否合法
  void addVariable(string expression);         // 添加变量
  string replaceVar(string expression);        // 替换变量
  void checkUnknownVar(string expression); // 检查表达式中是否有未知变量
  void pretreat(string expression);        // 预处理表达式
  void Assign();                           // 赋值方法

  // 数学函数计算
  bool containsMathFunc(string expression,
                        string funcName); // 检查表达式中是否包含某函数
  string replaceMathFunc(string expression); // 替换表达式中所有数学函数

  string getSqrt(string expression);         // 处理sqrt函数
  string getMax(string expression);          // 处理max函数
  string getMin(string expression);          // 处理min函数
  string getRandom(string expression);       // 处理random函数
  string getExp(string expression);          // 处理exp函数
  string getAbs(string expression);          // 处理abs函数
  string getSin(string expression);          // 处理sin函数
  string getCos(string expression);          // 处理cos函数
  string getTan(string expression);          // 处理tan函数
  string getCeil(string expression);         // 处理ceil函数
  string getFloor(string expression);        // 处理floor函数
  string getLog(string expression);          // 处理log函数

  // 工具函数
  vector<string> split(const std::string &strIn, char delim); // 字符串分割
  string replaceAll(string str, string oldStr,
                    string newStr); // 替换str中所有相同的oldStr为newStr

  // 小助手
  void helper();  //帮助系统
  void welcome(); // 欢迎界面
  void help();     // 帮助信息
  void about();    // 关于信息
  void exit();     // 退出程序
  void setLanguage(); // 设置语言
  void setscale(); // 设置精度
  void clear();    // 清空屏幕

  void printWelcomeEnglish(); // 打印欢迎界面英文版
  void printWelcomeChinese(); // 打印欢迎界面中文版
  void printHelpEnglish();    // 打印帮助信息英文版
  void printHelpChinese();    // 打印帮助信息中文版
  void printAboutEnglish();   // 打印关于信息英文版
  void printAboutChinese();   // 打印关于信息中文版


  // 数据成员
  string operatorSym; // 运算符号
  string infix;       // 表达式缓存
  bool isExit;  // 是否退出
  bool isHelp;  // 是否帮助

private:
  int language; // 语言选择
  vector<string> postfix;        // 后缀表达式向量
  stack<char> symStack;          // 符号栈
  stack<BigDecimal> figStack;    // 数字栈
  map<string, string> variables; // 变量表
  string stdInfix;               // 自定义标准格式化表达式
  BigDecimal result;             // 最终计算结果
};
