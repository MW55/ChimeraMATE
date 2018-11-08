import chimeramate_main
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("otus", type=str,
                    help="input file in FASTA format, with headers formatted in UCHIME style")
parser.add_argument("-soft_k", type=int, default=24,dest='soft_k')
parser.add_argument("-cutoff", type=int, default=10,dest='cutoff')
parser.add_argument("-softmask_file", type=str, default='softmasked.fasta',dest='softmask_file')
parser.add_argument("-chim_k", type=int, default=29,dest='chim_k')
parser.add_argument("-abskew", type=float, default=0.03,dest='abskew')

args = parser.parse_args()
#otus = sys.argv[1]
#soft_k = sys.argv[2]
#cutoff = sys.argv[3]
#softmask_file = sys.argv[4]
#chim_k = sys.argv[5]
#abskew = sys.argv[6]

chimeramate_main.run(args.otus, args.soft_k, args.cutoff, args.softmask_file, args.chim_k, args.abskew)
