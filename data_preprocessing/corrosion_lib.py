import argparse
import glob
import re
import zipfile


NUM_TIMESTEPS_PER_SIMULATION = 9

# Unzip specified number of corrosion files into the same directory,
# ignoring any edge files.
def extract_corrosion_output(path, num_simulations = 1):
  corrosion_filenames = get_corrosion_filenames(num_simulations)

  corrosion_path_lst = path.split('/')
  corrosion_dir_base = '/'.join(corrosion_path_lst[:-1])
  corrosion_dir = corrosion_dir_base

  with zipfile.ZipFile(path, 'r') as zip_obj:
    file_names = zip_obj.namelist()
    for file_name in file_names:
      file_name = file_name.split("/")[-1]
      if file_name not in corrosion_filenames:
        continue
      full_file_name = path.split("/")[-1].split(".")[0] + "/" + file_name
      print("extracting " + full_file_name + " to " + corrosion_dir)
      zip_obj.extract(full_file_name, corrosion_dir)

# Returns a list of corrosion simulation filenames, for the first num_simulation
# datapoints. Note that the files are 1-indexed, so this returns filenames from
# simulation 1 ... (num_simulations+1). Each simulation contains 
# NUM_TIMESTEPS_PER_SIMULATION files, one per timestep.
def get_corrosion_filenames(num_simulations = 1):
  filenames = []
  # file names from COMSOL are 1-indexed
  for simulation_idx in range(1, num_simulations + 1):
    for timestep in range(1, NUM_TIMESTEPS_PER_SIMULATION + 1):
      filename = "Corrosion_simulation_%d_timeStep_%d.txt" % (simulation_idx, timestep)
      filenames.append(filename)
  return filenames

# Returns a 1d corrosion map given a single filepath. A 1d corrosion map is
# represented as a python dictionary, with keys representing the location on the
# x-axis along a horizontal rebar, and the values representing corrosion depth
# at that point.
def extract_1d_corrosion_map_from_filepath(filepath):
  with open(filepath, 'r') as f:
    lines = f.readlines()
  corrosion = {}
  for line in lines:
    if line.startswith("%"):
      continue
    spl = re.split(r'\s+', line.strip())
    assert len(spl) == 2, spl
    rebar_location = float(spl[0])
    corrosion_depth = float(spl[1])
    corrosion[rebar_location] = corrosion_depth
  return corrosion

# Returns a list of pairs, each containing the filename of the simulation, and
# a 1-d corrosion map.
def extract_1d_corrosion_maps(output_dir, num_simulations = 1):
  corrosion_dir = output_dir + "/corrosion"

  file_and_corrosion_map = []
  for filename in get_corrosion_filenames(num_simulations):
    filepath = corrosion_dir + '/' + filename
    corrosion_map = extract_1d_corrosion_map_from_filepath(filepath)
    file_and_corrosion_map.append((filepath, corrosion_map))
  return file_and_corrosion_map

# Asserts that all corrosion simulations are sampled from the same rebar
# locations- that is, the points along the x-axis of the rebar are the same for
# all samples. If they are not, then we will need to rescale the inputs before
# training.
def verify_rebar_locations(file_and_corrosion_map):
  rebar_locations = [tuple(x[1].keys()) for x in file_and_corrosion_map]
  all(x == rebar_locations[0] for x in rebar_locations)
