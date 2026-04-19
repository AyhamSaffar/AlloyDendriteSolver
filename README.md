# Alloy Dendrite Solver

A fast, modular C++ tool for numerically calculating alloy solidification parameters from thermodynamic constants.

The solidification parameters of interest are:

    dT	Growth undercooling
    V	Dendrite tip velocity
    R	Dendrite tip radius
    C0	Bulk alloy solute concentration.

Where typically dT & C0 are fixed so that V and R can solved iteratively.

## Method

This library is fully modular, offering many different options for each key step. See [*docs*](docs) for more detailed
explanations of supported techniques for each module.

An *alloy* is a container for a binary alloy's thermodnamic constants. A **model** is a pair of coupled equations that
take V, R, and an alloy and return f1 and f2 respectively. Correct values for V and R would give zero f1 and f2, so the
goal is to iteratively update V and R to reduce the model's outputs.

Every model must be differentiable so that the **differentials** module can automatically calculate their
derivatives at compile time. An **optimiser** than takes these and the model outputs to calculate an update to V and R
that reduce model outputs.

The model, differentiation, and optimisation steps are then repeated until the error in V and R are acceptable.

The **approximators** module uses strong assumptions to give a good first guess for V and R. This increases the chance
that the opimiser converges. The **updators** module (work in progress) modifies the optimiser outputs in order to
further reduce the chance of divergence and increase convergence speed.

## Installation

Currently the library only supports building from source. This requires correctly configured basic C++ tools (CMake and
Make) on the command line. It is therefore recommended that Windows users use [WSL](http://ubuntu.com/desktop/wsl) (EG
through [VSCode](https://code.visualstudio.com/docs/remote/wsl)) which installs these tools automatically. 

Install the automatic differentiation compiler (Enzyme) using [Brew](https://brew.sh/)
```
brew install ezyme lld
```

Clone repository with supporting libraries (via git submodules)
```
git clone --recurse-submodules https://github.com/AyhamSaffar/AlloyDendriteSolver.git
```

Navigate into the repository
```
cd AlloyDendriteSolver
```

Now the library can be used following the standard CMake process. Code editors often contain functionality that automate
this (EG [VSCode CMake Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cmake-tools)) but the below
details the manual process.

First the build system (series of scripts specifying how to compile the code) must be generated
```
cmake -B build -DCMAKE_BUILD_TYPE=Debug
```
Setting build type to *Release* significantly speeds up the code but disables debugging at runtime.

Navigate into the build system directory
```
cd build
```

Compile (convert source code to executables) the available files
```
make
```
This can take a couple of minutes the first time as libraries must be compiled, but will be much faster on subsequent
calls.

All available executables will then be present in the directory. On Mac and Linux for example, unit tests can be run
with the following
```
./tests
```

## Usage

All key experiment types exist in [*scripts*](scripts) and are automatically found by CMake. On compilation, the
corresponding executables are created with the same, minus the suffix.

For example, take [*scripts/Minimal_Example.cpp*](scripts/Minimal_Example.cpp):

```C++
// Script used to demonstrate simplest possible workflow

#include <iostream>
#include <tuple>
#include "alloy.h"
#include "approximators.h"
#include "differentials.h"
#include "models.h"
#include "optimiser.h"


int main()
{
    // initialise variables
    const alloy::Alloy A{alloy::SnAg}; // common solder material
    double f1{}, f2{}, dV{}, dR{}, dT{10.0}, C0{5.0}; // only non SI unit is concentration (wt.%)
    double V{approx::getTipVelocity(dT, C0, A)}, R{approx::getTipRadius(dT, C0, A)};
    diff::Jacobian J{};

    // iteratively solve for V and R
    for (int step{0}; step<100; ++step)
    {
        std::tie(f1, f2) = models::LGK(V, R, dT, C0, A);
        J = diff::calculateGrads<models::LGK>(V, R, dT, C0, A);
        std::tie(dV, dR) = optimisers::newtonRaphson(f1, f2, J);
        V += 0.1 * dV; // smaller steps improve convergence
        R += 0.1 * dR;
    }

    // print result
    std::cout << "R = " << R << " m, V = " << V << " m/s\n";
    return 0;
}
```

This experiment can be run on Mac & Linux with the following

```
❯ ./Minimal_Example 
R = 3.25994e-07 m, V = 0.00526955 m/s
```

Any experiments that require data logging can dump to the [*data*](data) directory using the following macro
```C++
#include <fstream>
#include <string>

std::string dataPath{DATA_PATH};
std::ofstream outf{dataPath + "/my_results.csv"};

outf << "V,R,dT,C0\n" // logging column headers
```

[*Data*](data) also contains a Python [uv](https://docs.astral.sh/uv/) environment for data analysis. Each script also
has its own directory there with an interactive Python notebook. This plots the key outputs of any number of experiments
run using the corresponding script.

Custom scripts should be saved to [*scripts*](scripts) so that they can be automatically detected and compiled using the
CMake configuration.

## Support

The best place to report bugs or request features would be
[the project's issues page](https://github.com/AyhamSaffar/AlloyDendriteSolver/issues). This ensures that any current
maintainers are notified by email.

The maintainers will endeavour to respond to new issues as soon as possible.

## Roadmap

This library is still in early and active development. Planned additions include:

- Temperature dependant alloys to better accomodate for high undercoolings. Solute diffusivity could follow an an
Arrhenius model while solidus and liquidus slopes could be modelled by polynomials. Note changing some alloy
thermodynamic constants with temperature may violate a given model's assumptions.

- An [LKT-BCT](https://doi.org/10.1016/0001-6160(87)90174-X) model which maintains accuracy at higher undercoolings.

- Support for higher order gradients and their accompanying optimisers. One example is
[Halley's Method](https://en.wikipedia.org/wiki/Halley%27s_method), which takes advantage of second order tangents for
faster and more reliable convergence.

- An *updators* module that improves convergence through clamping / scaling / adapting optimiser outputs. This could
look like the following:
    ```C++
    #include <tuple>
    #include "updators.h"

    const double momentum{0.6};
    updators::AdaptiveUpdator updator{momentum};
    
    ...
        std::tie(V, R) = updator.apply(V, R, dV, dR); // in minimisation loop

        if (updator.hasConverged() || updator.hasDiverged() || updator.hasStalled())
    ...
    ```

- A Python / Matlab embedding to make common experiment types more accessible.

## Authors and Acknowledgment

This library was written by me, Ayham Saffar, under the supervision of Dr Chris Gourlay. It was inspired by previous
work in his research group.

No AI was used whatsoever in the development of this project.
