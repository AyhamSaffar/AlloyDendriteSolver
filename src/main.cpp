#include <iostream>
#include "alloy.h"
#include "approximators.h"
#include "solver.h"
#include "logger.h"


int main()
{
    alloy::Alloy A{alloy::SnAg};
    const double maxError{1e-5};
    const int maxIterations{1000};
    logger::FileLogger logger{"result_file.csv"};
    solver::LGKSolver solver{const A&, approx::getTipRadius(), approx::get};
    
    while (int iteration{0}; iteration<maxIterations; ++iteration)
    {
        logger << solver.step();
        std::cout << logger.params << '\n';

        if (logger.error <= maxError)
            break;
    }
}
