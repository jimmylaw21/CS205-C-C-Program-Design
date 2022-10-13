#include "BigDecimal.h"
#include "Calculator.h"
#include <iostream>

// 帮助系统
void Calculator::helper() {
  clearEmptySpace();
  checkInfix();
  exit();
  setLanguage();
  welcome();
  help();
  about();
  clear();
  setscale();
}
// 设置精度
void Calculator::setscale() {
  if (infix.find("scale") != -1) {
    std::string scaleStr = infix.substr(infix.find("scale") + 6);
    BigDecimal::setscale(std::stoi(scaleStr));
    isHelp = true;
  }
}
// 设置语言
void Calculator::setLanguage() {
  if (infix == "Chinese") {
    language = 0;
    cout << "中文设置成功！" << endl;
    isHelp = true;
  } else if (infix == "English") {
    cout << "English is set successfully!" << endl;
    language = 1;
    isHelp = true;
  }
}
// 欢迎界面
void Calculator::welcome() {
  if (infix == "welcome") {
    isHelp = true;
    if (language == 0) {
      printWelcomeChinese();
    } else if (language == 1) {
      printWelcomeEnglish();
    }
  }
}
// 帮助信息
void Calculator::help() {
  if (infix == "help") {
    isHelp = true;
    if (language == 0) {
      printHelpChinese();
    } else if (language == 1) {
      printHelpEnglish();
    }
  }
}
// 关于信息
void Calculator::about() {
  if (infix == "about") {
    isHelp = true;
    if (language == 0) {
      printAboutChinese();
    } else if (language == 1) {
      printAboutEnglish();
    }
  }
}
// 退出程序
void Calculator::exit() {
  if (infix == "exit") {
    isHelp = true;
    isExit = true;
  }
}
// 清空屏幕
void Calculator::clear() {
  if (infix == "clear") {
    isHelp = true;
    system("clear");
  }
}

void Calculator::checkInfix() {
  if (infix.length() <= 0) {
    if (language == 0) {
      cout << "表达式为空！" << endl;
    } else if (language == 1) {
      cout << "The expression is empty!" << endl;
    }
    isHelp = true;
  }
}

void Calculator::clearEmptySpace(){
    trim(infix);
}

void Calculator::trim(string &s)
{
	int index = 0;
	if(!s.empty())
	{
		while( (index = s.find(' ',index)) != -1)
		{
			s.erase(index,1);
		}
	}
}

void Calculator::printWelcomeEnglish() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout << "Please enjoy your calculation!" << endl;
  cout << "If you need help, please print help to get some information."
       << endl;
  cout << "If you need welcome, please print welcome to get this information "
          "again."
       << endl;
  cout << "If you need exit, please print exit to get exit information."
       << endl;
  cout << "如果需要显示中文，请输入Chinese" << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}

void Calculator::printWelcomeChinese() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout << "欢迎使用计算器！" << endl;
  cout << "如果需要帮助，请输入help" << endl;
  cout << "如果需要显示英文，请输入English" << endl;
  cout << "如果需要退出，请输入exit" << endl;
  cout << "If you need to display English,please print English" << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}

void Calculator::printAboutEnglish() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout
      << "This calculator is mainly for the Assignment2 in CS205(C/C++) course."
      << endl;
  cout << "This calculator supports basic operations including exponents and "
          "factorials, "
       << endl;
  cout << "as well as variable assignment calculations and mathematical "
          "functions"
       << endl;
  cout << "Course Instructor: Shiqi Yu" << endl;
  cout << "Author: Jimmy Luo" << endl;
  cout << "Version: 1.00" << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}
void Calculator::printAboutChinese() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout << "本计算器主要是用于南方科技大学CS205(C/C++程序设计)课程的第二次项目。"
       << endl;
  cout << "本计算器支持包括指数和阶乘的基本运算，以及变量赋值计算和数学函数"
       << endl;
  cout << "课程教师：于仕琪" << endl;
  cout << "作者：罗皓予" << endl;
  cout << "版本号：1.00" << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}
void Calculator::printHelpEnglish() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout << "0. The constant value in this calculator:" << endl;
  cout << "E=2.718281828459" << endl;
  cout << "PI=3.1415926535898" << endl;
  cout << endl;

  cout << "Supported functions:" << endl;
  cout << "1. Basic mathematical calculation: + - * / % ^ ! ||" << endl;
  cout << "example-1.1:" << endl;
  cout << "3*4" << endl;
  cout << "12" << endl;
  cout << "example-1.2:" << endl;
  cout << "3.45*4.32" << endl;
  cout << "14.904" << endl;
  cout << endl;

  cout << "2. calculation with brackets: ( ) [ ] { }" << endl;
  cout << "example-2.1:" << endl;
  cout << "(2+5)*9.9" << endl;
  cout << "69.3" << endl;
  cout << "example-2.2:" << endl;
  cout << "(4.5+6.8)*1.2" << endl;
  cout << "13.56" << endl;
  cout << endl;

  cout << "3. Definition with single char variable: x, y, z" << endl;
  cout << "example-3.1:" << endl;
  cout << "x=3" << endl;
  cout << "y=6" << endl;
  cout << "x+2*y" << endl;
  cout << "21" << endl;
  cout << "example-3.2:" << endl;
  cout << "x=1" << endl;
  cout << "y=x*2" << endl;
  cout << "x=3" << endl;
  cout << "x+x*y" << endl;
  cout << "9" << endl;
  cout << "example-3.3:" << endl;
  cout << "x=3" << endl;
  cout << "y=4" << endl;
  cout << "x=x+y" << endl;
  cout << "x" << endl;
  cout << "7" << endl;
  cout << endl;

  cout << "4. Some single math functions:" << endl;
  cout << "example-4.1:" << endl;
  cout << "sqrt(3)" << endl;
  cout << "1.732051" << endl;
  cout << "example-4.2:" << endl;
  cout << "sin(2)" << endl;
  cout << "0.909297" << endl;
  cout << "sqrt(x), pow(a,b), max(a,b), abs(x), exp(x),\n"
          "min(a,b), random(), sin(x), cos(x), tan(x), ceil(x), floor(x), "
          "log(x),\n"
          "is supported"
       << endl;
  cout << "However, the compound calculation is not supported." << endl;
  cout << "If input a compound calculation, the calculator will exist." << endl;
  cout << "such as: pow(pow(2,3),4)" << endl;
  cout << endl;

  cout << "5. Arbitrary Precision:" << endl;
  cout << "Any length of integer or decimal value are supported in this "
          "calculator."
       << endl;
  cout << "example-5.1:" << endl;
  cout << "5.44444444444444444444444444444444444444444*2+1."
          "1111111111111111111111111111111111111"
       << endl;
  cout << "11.99999999999999999999999999999999999998888" << endl;
  cout << "example-5.2:" << endl;
  cout << "10.8888888888888888888888888888888888888888888888888888888/2"
       << endl;
  cout << "5.4444444444444444444444444444444444444444444444444444444" << endl;
  cout << endl;

  //   cout << "6. Add Single line annotation or multiple line annotations:" <<
  //   endl; cout << "single line: start with // " << endl; cout << "multiple
  //   lines: start with /*, end with */" << endl; cout << "The input in the
  //   annotation will not be detected by the calculator"
  //        << endl;
  //   cout << endl;

  //   cout << "This calculator is still being improved. Due to time limit, the
  //   "
  //           "calculator supports these functions."
  //        << endl;
  //   cout <<
  //   "-----------------------------------------------------------------"
  //        << endl;

  cout << "The calculator is still being improved. \n"
          "Due to time constraints, "
          "only so many functions are currently"
          "supported."
       << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}

void Calculator::printHelpChinese() {
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
  cout << "0.计算器中的一些常数:" << endl;
  cout << "E=2.718281828459" << endl;
  cout << "PI=3.1415926535898" << endl;
  cout << endl;
  cout << "计算器中支持的功能" << endl;
  cout << "1.基本数学运算：+ - * / ^ ! ||" << endl;
  cout << "示例-1.1:" << endl;
  cout << "3*4" << endl;
  cout << "12" << endl;
  cout << "示例-1.2:" << endl;
  cout << "3.45*4.32" << endl;
  cout << "14.904" << endl;
  cout << endl;

  cout << "2.带括号的运算：( ) [ ] { }" << endl;
  cout << "示例-2.1:" << endl;
  cout << "(2+5)*9.9" << endl;
  cout << "69.3" << endl;
  cout << "示例-2.2:" << endl;
  cout << "(4.5+6.8)*1.2" << endl;
  cout << "13.56" << endl;
  cout << endl;

  cout << "3.单个字符的自变量: x, y, z" << endl;
  cout << "示例-3.1:" << endl;
  cout << "x=3" << endl;
  cout << "y=6" << endl;
  cout << "x+2*y" << endl;
  cout << "15" << endl;
  cout << "示例-3.2:" << endl;
  cout << "x=1" << endl;
  cout << "y=x*2" << endl;
  cout << "x=3" << endl;
  cout << "x+x*y" << endl;
  cout << "9" << endl;
  cout << "示例-3.3:" << endl;
  cout << "x=3" << endl;
  cout << "y=4" << endl;
  cout << "x=x+y" << endl;
  cout << "x" << endl;
  cout << "7" << endl;
  cout << endl;

  cout << "4.一些简单的数学函数:" << endl;
  cout << "示例-4.1:" << endl;
  cout << "sqrt(3)" << endl;
  cout << "1.732051" << endl;
  cout << "示例-4.2:" << endl;
  cout << "sin(2)" << endl;
  cout << "0.909297" << endl;
  cout << "sqrt(a), pow(a,b), max(a,b), abs(a), exp(a),\n"
          "min(a,b), random(), sin(x), cos(x), tan(x), ceil(x), floor(x), "
          "log(x)\n"
          "已支持"
       << endl;
  cout << "该计算器尚未完全支持复合数学函数的运算" << endl;
  cout << "如果输入复合数学函数的运算，则该计算器有可能异常退出" << endl;
  cout << "例如：pow(pow(2,3),4)" << endl;
  cout << endl;

  cout << "5.高精度运算：" << endl;
  cout << "该计算器支持任意长度的浮点数运算" << endl;
  cout << "示例-5.1:" << endl;
  cout << "5.44444444444444444444444444444444444444444*2+1."
          "1111111111111111111111111111111111111"
       << endl;
  cout << "11.99999999999999999999999999999999999998888" << endl;
  cout << "示例-5.2:" << endl;
  cout << "10.8888888888888888888888888888888888888888888888888888888/2"
       << endl;
  cout << "5.4444444444444444444444444444444444444444444444444444444" << endl;
  cout << endl;

  //   cout << "6.单行注释、多行注释实现" << endl;
  //   cout << "单行注释：某一行以//开始 " << endl;
  //   cout << "多行注释：以/*开始，以*/结束" << endl;
  //   cout << "注释内的内容不会被计算器所识别" << endl;
  //   cout << endl;

  cout << "该计算器仍然处于改进当中。\n由于时间的限制，现暂时只支持这么多功能。"
       << endl;
  cout << "--------------------------------------------------------------------"
          "-----------------"
       << endl;
}
