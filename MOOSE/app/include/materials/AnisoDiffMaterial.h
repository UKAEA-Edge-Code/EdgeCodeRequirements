#pragma once

#include "Material.h"
//#include "RankTwoTensor.h"

class AnisoDiffMaterial : public Material {
public:
    AnisoDiffMaterial(InputParameters const & parameters);

    static InputParameters validParams();

    protected:
    virtual void computeQpProperties() override;
    private:
    MaterialProperty<RealTensorValue> &_A;
    Real _alpha;
    Real _m;
    Real _epsilon;
};

class SimpleDiffMaterial : public Material {
public:
    SimpleDiffMaterial(InputParameters const & parameters);

    static InputParameters validParams();

    protected:
    virtual void computeQpProperties() override;
    private:
    MaterialProperty<RealTensorValue> &_A;
    Real _alpha;
    Real _m;
    Real _epsilon;
};
