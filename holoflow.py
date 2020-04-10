# subprocess
import argparse

#Argument parsing
parser = argparse.ArgumentParser(description='Runs holoflow pipeline.')
parser.add_argument('-f', help="megahit input", dest="input", required=True)
parser.add_argument('-d', help="spades input", dest="path", required=True)
parser.add_argument('-w', help="output directory", dest="workflow", required=True)
parser.add_argument('-c', help="assembler", dest="config", required=True)
args = parser.parse_args()

input=args.input
path=args.path
workflow=args.workflow
config=args.config
