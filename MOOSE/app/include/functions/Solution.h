#pragma once
#include "Function.h"

class SolutionFunc : public Function {
    public:
    static InputParameters validParams();
    SolutionFunc(const InputParameters &params);
    virtual Real value(Real t, const Point &p) const;
    private:
    Real _alpha;
    Real _m;
    Real _epsilon;
};

