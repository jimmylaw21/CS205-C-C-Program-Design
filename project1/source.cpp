#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

using namespace std;

string processString(char a[]);
void mul(string a, string b);

int main(int argc, char *argv[]) {
  string a = processString(argv[1]);
  string b = processString(argv[2]);

  if (a == "error") {
    cout << "The input cannot be interpret as numbers!" << endl;
    return 0;
  }

  if (b == "error") {
    cout << "The input cannot be interpret as numbers!" << endl;
    return 0;
  }
  mul(a, b);
  return 0;
}

// 炒饭大厨，同时能处理科学计数法
string processString(char str[]) {
  // 检查炒饭，拆分底数和指数
  char *nums = (char *)calloc(100000, sizeof(char *));
  char *ans = (char *)calloc(100000, sizeof(char *));
  char *exps = (char *)calloc(100, sizeof(char *));

  int len = strlen(str);
  int start = 0;
  bool isNeg = false;

  if (str[0] == '-') {
    start++;
    isNeg = true;
  }
  // 拆分
  bool flag = false; // false表示底数
  bool hasPoint = false;
  bool hasSign = false;
  int lenNums = 0, lenExps = 0;
  for (int i = start; i < len; i++) {

    if (str[i] == 'e' && flag == false) {
      flag = true;
    } else if (flag == false) // 底数阶段
    {
      if (!((str[i] == '.' && !hasPoint) || (str[i] >= '0' && str[i] <= '9'))) {
        cout << "底数部分输入有误" << endl;
        return "error";
      }
      if (str[i] == '.') {
        hasPoint = true;
      }
      nums[lenNums++] = str[i]; // 含小数点
    } else if (flag == true)    // 指数阶段，带正负号
    {
      if (!((str[i] == '+' || str[i] == '-') && !hasSign) &&
          !(str[i] >= '0' && str[i] <= '9')) {
        cout << "指数部分输入有误" << endl;
        return "error";
      }

      if (str[i] == '+' || str[i] == '-') {
        hasSign = true;
      }

      exps[lenExps++] = str[i];
    }
  }

  int index = 0;            // 记录ans数组
  int ex = 0;               // 指数的值
  int actCnt = lenNums - 2; // 实际小数点后的位数
  sscanf(exps, "%d", &ex);

  // 用（负号），底数和指数组装一个string返回
  if (isNeg) {
    ans[0] = '-';
    index++;
  }

  if (ex == 0) // 指数是0次幂
  {

    for (int i = 0; i < lenNums; i++) {
      ans[index++] = nums[i]; // 若为底数为正数，nums[0]是'+'
    }
  } else if (exps[0] == '+' ||
             exps[0] <= '9' &&
                 exps[0] >= '1') // 指数为正 || exps[0] <= '9' && exps[0] >= '1'
  {

    // 对底数操作
    ans[index++] = nums[0]; // 处理小数点前的一位
    if (ex >= actCnt)       // 指数大于等于小数点后的位数
    {

      for (int i = 2; i < actCnt + 2; i++) // i从小数点后的一位开始
      {
        ans[index++] = nums[i];
      }
      for (int i = 0; i < ex - actCnt; i++) // 多出来的补零
      {

        ans[index++] = '0';
      }
    } else // 指数小于后面的位数，要点小数点
    {

      for (int i = 2; i < ex + 2; i++) // 先写数
      {
        ans[index++] = nums[i];
      }
      ans[index++] = '.';
      for (int i = ex + 2; i < lenNums; i++) {
        ans[index++] = nums[i];
      }
    }
  } else if (exps[0] == '-') // 指数为负
  {

    ex = -ex; // 取相反数
    ans[index++] = '0';
    ans[index++] = '.';
    if (ex > 1) {
      for (int j = 0; j < ex - 1; j++) {
        ans[index++] = '0';
      }
    }
    for (int i = 0; i < lenNums; i++) {
      if (nums[i] != '+' && nums[i] != '.') {
        ans[index++] = nums[i];
      }
    }
  }

  // 输出
  string ansStr = ans;

  free(nums);
  free(ans);
  free(exps);

  return ansStr;
}

// 大数乘法方法组，使用string处理自然进位，可以处理正负的整数和小数

// 读取数字
int toint(char n) {
  if (n - '0' >= 0 && n - '0' <= 9) {
    return n - '0';
  } else {
    return 999999;
  }
}

// 运算时消除小数点，返回小数点的位置
int deletep(string &a) {
  string::iterator it;
  int i = 0;
  for (it = a.end() - 1; it != a.begin(); it--, i++) {
    if (*it == '.') {
      a.erase(it);
      return i;
    }
  }
  return 0;
}

// 去掉小数末尾的零
void deletebackzero(string &a) {
  string::iterator it;
  for (it = a.end() - 1; it != a.begin(); it--) {
    if (*it != '0') {
      if (*it == '.') {
        a.erase(it);
        return;
      }
      return;
    }
    a.erase(it);
  }
}

// 去掉小数前面的零,当数字小于1大于等于0时，保留一个零
void deletefowardzero(string &a) {
  for (int i = 0;;) {
    if (a[0] == '0') {
      a.erase(0, 1);
    } else {
      if (a[0] == '.') {
        a.insert(0, "0");
      }
      break;
    }
  }
}

// 核心乘法算法
void mul(string a, string b) {

  // 测试计算时长
  clock_t start = clock();
  bool isNeg = false;

  // 处理负数
  if (a[0] == '-') {
    a.erase(a.begin());
    isNeg = !isNeg;
  }
  if (b[0] == '-') {
    b.erase(b.begin());
    isNeg = !isNeg;
  }

  // 处理小数点
  int mia = deletep(a);
  int mib = deletep(b);
  int mi = mia + mib;

  // 初始化string长度和进位值
  int s = 0;
  int x = a.size() - 1;
  int y = b.size() - 1;

  // 设置内存占用量，溢出时增大
  string result[100000];

  // 乘法运算
  for (int i = y, mul = 0; i >= 0; i--, mul++) { // 从低位开始乘
    int num1 = toint(b[i]);
    string str = "";
    for (int j = x; j >= 0; j--) {
      int num2 = toint(a[j]);
      int temp = num1 * num2 + s;
      int g = temp % 10;
      s = temp / 10;
      str = to_string(g) + str;
      if (j == 0 && s != 0) {
        str = to_string(s) + str;
      }
    }
    result[mul] = str; // 保存第二个数的一位与第一个数每一位的相乘结果
    s = 0;             // 进位值归0
  }

  for (int i = y, num = 0; i >= 0; i--, num++) // 末尾加零以对位
  {
    for (int j = num; j > 0; j--) {
      result[y - i] = result[y - i] + '0';
    }
  }

  // 加法运算
  string anslen = result[y];
  x = anslen.size() - 1; // 确定结果长度
  string fresult;        // 初始化结果
  s = 0;                 // 进位值归零

  for (int i = x, count = 0; i >= 0; i--, count++) { // 从低位开始加
    int temp = 0;
    for (int j = 0; j <= y; j++) {
      string lim = result[j];     // 从保存的结果中取出一位
      int limit = lim.size() - 1; // 确定该位的长度
      if (i < x - limit) {        // 该位长度不够时，补零
        continue;
      }
      temp += toint(result[j][i - x + limit]); // 将该位的每一位相加
    }
    temp += s;
    int g = temp % 10;
    s = temp / 10;
    fresult = to_string(g) + fresult;
  }

  // 处理小数点
  if (mi != 0) {
    fresult.insert(fresult.size() - mi, ".");

    deletebackzero(fresult);
  }

  deletefowardzero(fresult);

  // 处理负数
  if (isNeg) {
    fresult = "-" + fresult;
  }

  // 输出结果
  cout << "result: " << fresult << endl << endl;

  // 测试计算时长
  clock_t end = clock();
  double endtime = (double)(end - start) / CLOCKS_PER_SEC;
  cout << "Total time:" << endtime * 1000 << "ms" << endl; // ms为单位
}