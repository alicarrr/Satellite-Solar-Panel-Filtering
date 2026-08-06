#include "stuff.hpp"
#include <iostream>

void doStuff(int numberOfStuffs)
{
    for (int i = 0; i < numberOfStuffs; i++)
    {
        std::cout << "Hi " << i + 1 << std::endl;
    }
}