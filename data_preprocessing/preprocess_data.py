import argparse
import corrosion_lib
import output_lib

parser = argparse.ArgumentParser(description="Extract and process corrosion and output files")
parser.add_argument('--output_path', help='Path to location of corrosion data. If --extract=True, the extracted files will also be extracted to here.')
parser.add_argument('--num_simulations', type=int, help='Number of simulations to process')
parser.add_argument('--extract', action='store_true', help='If true, unzip and extract files from zipped filepath.')
parser.add_argument('--corrosion_zipped_filename', default='corrosion.zip')
parser.add_argument('--output_zipped_filename', default='Data_outputs.zip')

args = parser.parse_args()

print(args)

def preprocess():
  # Extract the first num_simulations experiments to output_path.
  if (args.extract):
    corrosion_lib.extract_corrosion_output(args.output_path + '/' + args.corrosion_zipped_filename, args.num_simulations)
    output_lib.extract_FEM_output(args.output_path + '/' + args.output_zipped_filename, args.num_simulations)

  corrosion_maps = corrosion_lib.extract_1d_corrosion_maps(args.output_path)
  output_maps = output_lib.extract_concrete_outputs(args.output_path)
  
  print(len(corrosion_maps), len(output_maps))
  import pdb; pdb.set_trace()



if __name__ == "__main__":
   preprocess()
