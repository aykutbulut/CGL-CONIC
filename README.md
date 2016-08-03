# CGL-Conic [![Build Status](https://travis-ci.org/aykutbulut/CGL-CONIC.svg?branch=master)](https://travis-ci.org/aykutbulut/CGL-CONIC)

CGL-Conic is library for generation of conic cuts for Mixed Integer Second
Order Conic Optimization (MISOCO) problems. CGL-Conic can be considered as a
generalization of COIN-OR's Cut Generation Library (CGL).

CGL-Conic depends on CGL and other COIN-OR projects. CGL-Conic implements the
following cut procedures,

* Disjunctive Cuts given by Belotti et. al.,
see [related paper][Belotti et. al. 2015]
* Conic mixed-integer rounding cuts by Atamturk and Narayanan,
see [related paper][Atamturk and Narayanan 2008]
* Outer approximation cuts, related paper to be published soon.
* Interior Point Method (IPM) approximation cuts, related paper to be published soon.

Cut procedures that are considered for implementation.

* Two-term disjunctions on the second-order cone by Karzan and Yildiz,
see [related paper][Karzan Yildiz 2015]
* Cuts for 0-1 mixed convex programming by Stubbs and Mehrotra,
see [related paper][Stubbs and Mehrotra 1999].
* Cuts for mixed 0-1 conic programming by Cezik and Iyengar,
see [related paper][Cezik and Iyengar 2005]

[Belotti et. al. 2015]: http://link.springer.com/chapter/10.1007/978-3-319-17689-5_1
[Atamturk and Narayanan 2008]: http://link.springer.com/article/10.1007/s10107-008-0239-4
[Karzan and Yildiz 2015]: http://link.springer.com/article/10.1007/s10107-015-0903-4
[Stubbs and Mehrotra 1999]: http://link.springer.com/article/10.1007/s10107-015-0903-4
[Cezik and Iyengar 2005]: http://link.springer.com/article/10.1007/s10107-005-0578-3
[6]: https://github.com/aykutbulut/DisCO

CGL-Conic is used by [DisCO][6] to solve MISOCO problems.

# Installing CGL-Conic

## Basic installation

CGL-Conic is tested/works in Linux environment only for now. You can use
COIN-OR BuildTools fetch and build utility for building and installing
CGL-Conic. After cloning CGL-Conic, typing following commands should work.

```shell
git clone --branch=stable/0.8 https://github.com/coin-or-tools/BuildTools
bash BuildTools/get.dependencies.sh fetch
bash BuildTools/get.dependencies.sh build
```

First command clones BuildTools, second fetches CGL-Conic dependencies and third builds/installs GGL-Conic.

## Advanced Users

CGL-Conic can find dependencies using ```pkg-config```. If your ```.pc``` files
are at your ```PKG_CONFIG_PATH``` CGL-Conic configure script will find
them. Running ```make install``` should work once the configure script is
successful.

## Specifying an IPM Solver

Some cut implementations in CGL-Conic depends on solving SOCO problems using
IPM method. You can specify the solver to use for this purpose. CGL-Conic uses
Ipopt by default. You can use Cplex[1] or Mosek[2] through their OsiConic
interface.

[1]: https://github.com/aykutbulut/OsiCplex
[2]: https://github.com/aykutbulut/OSI-MOSEK

To specify the solver, you need to give ```--with-ipm-solver``` flag to
configure script. For example, following command configures OSI-Conic with
Mosek.

```shell
./configure --with-ipm-solver=mosek
```

Similarly you can use ```cplex``` or ```ipopt``` instead of ```mosek```. If no
ipm solver is specified CGL-Conic will use Ipopt.
