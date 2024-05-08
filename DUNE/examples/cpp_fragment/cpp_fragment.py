from dune.fem.function import gridFunction
from dune.fem.plotting import plotComponents
from dune.grid import structuredGrid as leafGridView
from matplotlib import pyplot as plt

# Set pyplot color map
plt.set_cmap("jet")

mins = [0, 0]
maxes = [1, 1]
ranges = [4, 4]
gridView = leafGridView(mins, maxes, ranges)

code = """
#include <cmath>
#include <dune/common/fvector.hh>
template <class GridView>
auto aTimesExact(double a) {
  return [a](const auto& en,const auto& xLocal) -> auto {
    auto x = en.geometry().global(xLocal);
    return a*(1./2.*(std::pow(x[0],2)+std::pow(x[1],2)) - 1./3.*(std::pow(x[0],3) - std::pow(x[1],3)) + 1.);
  };
}
"""
exactCpp = gridFunction(code, gridView, name="aTimesExact", order=2, args=[2.0])
exactCpp.plot()
