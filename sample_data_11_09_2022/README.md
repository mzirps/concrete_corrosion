Sample corrosion data from one COMSOL similation run, over 9 timesteps, represented in files Corrosion_simulation_1_timeStep_k.txt. Each file contains a header with metadata, followed by two columns: x and Height. x represents the point along the rebar lengthwise, and height represents the amount of corrosion.

The full dataset contains ~10k simulation runs, each with 9 timesteps.

For each full COMSOL simulation dataset (such as the one example here), there is a single corresponding output file (in this sample, named output_1_1.mat) containing properties of the concrete, as well as a binary label indicating whether there was a surface crack in the FEM model.
