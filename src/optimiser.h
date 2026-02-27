#ifndef OPTIMISER_H
#define OPTIMISER_H

#include <vector>
#include "alloy.h"

namespace optimiser
{
    class BaseOptimiser
    {
        struct Jacobian{double df1dV{}; double df1dR{}; double df2dV{}; double df2dR{};};
        
        private:
            Jacobian getDerivatives();
        
        protected:
            alloy::Alloy m_A{};
            double m_V{};
            double m_R{};
            double m_f1Error{};
            double m_f2Error{};
            virtual double calculatef1() = 0;
            virtual double calculatef2() = 0;

        public:
            BaseOptimiser(const alloy::Alloy& A, double V0, double R0)
                :m_A{A}, m_V{V0}, m_R{R0} {}
            std::vector<int> step();
            friend Logger& operator<<(Logger& log, const BaseOptimiser& optimiser);
    
    };  
}

#endif
