#pragma once
#include "Function.h"

class RHSFunc : public Function {
    public:
    static InputParameters validParams();
    RHSFunc(const InputParameters &params);
    virtual Real value(Real t, const Point &p) const;
    private:
    Real _alpha;
    Real _m;
    Real _epsilon;
};

