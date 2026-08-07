// Your main file doesn't have to be called main, but it must contain the main function.
#include "stuff.hpp"
#include <lcm/lcm-cpp.hpp>
#include <vector>
#include "simMessages/sim_data_t.hpp"
#include <iostream>
#include <chrono>
#include <thread>
#include "spectra.hpp"


/*
This function is your executable function, it is the function that gets called when you run your code.
*/
int main()
{
    std::this_thread::sleep_for(std::chrono::milliseconds(3000));

    lcm::LCM lcm;

    if(!lcm.good())
    {
        return 1;
    }
    
    // std::cout<<getSunSpectrum("hi")<<std::endl;

    for (int i=0; i<3; i++)
    {
        simMessages::sim_data_t test;
        test.t = 0;
        test.pos[0] = 0;
        test.pos[1] = 1;
        test.pos[2] = -2;

        lcm.publish("TEST", &test);
        std::cout<<"Published"<<std::endl;
    }
    


    

    // leave at bottom of main - shows no errors
    return 0;
    
}