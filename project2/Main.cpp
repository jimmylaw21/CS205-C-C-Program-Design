#include "BigDecimal.h"
#include "Calculator.h"
#include <iostream>
using namespace std;

// main函数
int main() {
  Calculator cal;
  cal.printWelcomeEnglish();
  // system("chcp 65001");

  while (true) {
    if (cal.isExit) {
      break;
    }
    getline(cin, cal.infix);
    cal.helper();
    if (cal.isHelp) {
      cal.isHelp = false;
      continue;
    }
    cal.calculate();
  }
  return 0;
}
