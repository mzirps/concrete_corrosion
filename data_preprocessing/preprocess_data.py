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
  if (args.extract):
    corrosion_lib.extract_corrosion_output(args.output_path + '/' + args.corrosion_zipped_filename)
    output_lib.extract_FEM_output(args.output_path + '/' + args.output_zipped_filename)

if __name__ == "__main__":
   preprocess()
