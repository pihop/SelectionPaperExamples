# "Feedback between stochastic gene networks and population dynamics enables cellular decision making"

This repository provides the implementation of the models and parameter fitting.
The stochastic simulation and the finite state projection algorithm for numerical results is maintained in a separate
repository [AgentBasedFSP](https://github.com/pihop/AgentBasedCells.jl) and added as a dependency here.

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/).
It is authored by Paul Piho.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data is not included in the git-history and needs to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

 
