#include <iostream>
#include "Calculator.h"
#include "BigDecimal.h"
using namespace std;

//const int MAX_EXP_LEN = 1000;			//最大表达式长度	防止内存溢出

//main函数
int main()
{
	Calculator cal;
	cal.printWelcomeEnglish();
	// system("chcp 65001");

	while (true) {
		if (cal.isExit) {
			break;
		}
		getline(cin, cal.infix);
		/*
		if (cal.infix.length() > MAX_EXP_LEN) {
			cout << "超出最大长度！" << endl;
			system("pause");
		}
		else {
			cal.calculate();
		}
		*/
		//cal.Assign();
		cal.helper();
		if (cal.isHelp) {
			cal.isHelp = false;
			continue;
		}
		cal.calculate();
	}
	return 0;
}
