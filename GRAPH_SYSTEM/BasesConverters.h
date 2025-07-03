#pragma once
#include <iostream>
#include <sstream>
#include <stack>
#include <string.h>

static int FromBaseToBaseTen(const std::string& number , int CurrentBase)
{
	int sum = 0;
	char ch;
	int value;
	for (int i = 0; i < number.length(); i++)
	{
		ch = number[i];
		if (ch >= '0' && ch <= '9')value = ch - '0';
		else if (ch >= 'A' && ch <= 'Z')value = 10 + ch - 'A';
		else return -1;//Error
		if (value >= CurrentBase)return -1;
		sum = value + (sum * CurrentBase);
		
	}
	return sum;
}



static std::string BaseTenToBase(int number , int BASE , int NUMBERS_OF_DIGITS = 0)
{
	if (number == 0)return std::string(std::max(NUMBERS_OF_DIGITS, 1), '0');
	const char digits[] = { '0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F' };
	if (BASE > 16 || BASE <= 1)return "[ERROR] Invalid base. Supported range: 2 to 16";
	std::stack<char> stack;
	std::stringstream result;
	while (number > 0)
	{
		stack.push(digits[number % BASE]);
		number /= BASE;
	}
	int Padding = NUMBERS_OF_DIGITS - stack.size();
	while (Padding > 0)
	{
		result << '0';
		Padding--;
	}

	while (!stack.empty())
	{
		result << stack.top();
		stack.pop();
	}
	return result.str();
}

