import os
import argparse
import subprocess

root = os.path.dirname(os.path.abspath(__file__))


def main() -> None:
    args = parseArgs()

    if args.num_confs == 1:
        cmd = "python {0}/tools/gen_3d_structs.py -i {1} -o {2} -t {3}".format(root, args.in_file, args.out_file, args.max_time)
    else:
        cmd = "python {0}/tools/gen_confs.py -i {1} -o {2}".format(root, args.in_file, args.out_file)

    subprocess.Popen(cmd, shell=True).wait()


def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Generates conformers for a given list of SMILES strings.')

    parser.add_argument('-i',
                        dest='in_file',
                        required=True,
                        metavar='<file>',
                        help='Molecule input file in CSV format. A column named "smiles" is required.')
    parser.add_argument('-o',
                        dest='out_file',
                        required=True,
                        metavar='<file>',
                        help='Output file in SDF SDF format.')
    parser.add_argument('-t',
                        dest='max_time',
                        required=False,
                        metavar='<int>',
                        type=int,
                        default=60,
                        help='Max. allowed molecule processing time (default: 60 sec)')
    parser.add_argument('-q',
                        dest='quiet',
                        required=False,
                        action='store_true',
                        default=False,
                        help='Disable progress output (default: false)')
    parser.add_argument('-n',
                        dest='num_confs',
                        required=False,
                        metavar='<int>',
                        type=int,
                        default=10,
                        help='Number of conformers to generate (default: 10)')
    
    return parser.parse_args()


if __name__ == '__main__':
    main()