#!/usr/bin/env python3

# consolidate multiple conda yamls into one

# usage python3 consolidate_envs.py <env_1> <env_2> ... <env_n> <output>

import yaml
import sys

def main(args):
	
	yamls = args[1:-1]
	outfile = args[-1]
	
	out = {'name': 'isling', 'channels': ['conda-forge', 'bioconda'], 'dependencies': []}

	# check each yaml for new info
	for y in yamls:
		env = import_yaml(y)
		if env is None:
			continue
		for c in env['channels']:
			if c not in out['channels']:
				out['channels'].append(c)
		for d in env['dependencies']:
			if d not in out['dependencies']:
				out['dependencies'].append(d)
		
	# write output		
	with open(outfile, 'w') as outhandle:
		yaml.dump(out, outhandle)	
		
	print(f"wrote output to {outfile}")


def import_yaml(config):
	# import simulation config yaml
	with open(config, 'r') as stream:
		try:
			in_config = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)	
			
	return in_config

	
if __name__ == "__main__":
	main(sys.argv)
