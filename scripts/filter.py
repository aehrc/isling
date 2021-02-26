#!/usr/bin/env python3

# usage python3 filter.py -i <input_file> -k <keep_file> -e <exclude_file> -c <criteria>

## filter integrations output by perl scripts by user-defined criteria
## user specifies criteria which define integrations to KEEP

## three terms define each criterion: 1) a column to apply to, 2) a comparator, and 3) a value

## allowable compararators are '>', '<', '<=', '>=', '==', '!='
## allowable values depend on the column - for integer columns, the value must be an integer, 
##		and for text columns, the text must be one of the allowed options

### NoAmbiguousBases (integer)
### OverlapType ('none', 'gap', 'overlap', 'discordant')
### Orientation ('hv', 'vh')
### ViralOrientation ('+', '-')
### HostEditDist (integer)
### ViralEditDist (integer)
### TotalEditDist (integer)
### PossibleHostTranslocation ('yes', 'no')
### PossibleVectorRearrangement ('yes', 'no')
### HostPossibleAmbiguous ('yes', 'no')
### ViralPossibleAmbiguous ('yes', 'no')
### Type ('chimeric', 'discordant')
### HostMapQ (integer)
### ViralMapQ (integer)

## critera may be separaterated by AND and OR
## brackets may be used to group criteria
 
from sys import argv
import numpy as np
import argparse
import csv
import pdb

columns = {
	'NoAmbiguousBases': 'integer',
	'OverlapType' : {'none', 'gap', 'overlap', 'discordant'},
	'Orientation': {'hv', 'vh'},
	'ViralOrientation': {'+', '-'},
	'HostEditDist': 'integer',
	'ViralEditDist': 'integer',
	'TotalEditDist': 'integer',
	'PossibleHostTranslocation': {'yes', 'no'},
	'PossibleVectorRearrangement': {'yes', 'no'},
	'HostPossibleAmbiguous': {'yes', 'no'},
	'ViralPossibleAmbiguous': {'yes', 'no'},
	'Type': {'chimeric', 'discordant'},
	'HostMapQ': 'integer',
	'ViralMapQ': 'integer'
}

def main(args):
	parser = argparse.ArgumentParser()
	parser.add_argument('--input', '-i', help='Input integration file to filter', required=True)
	parser.add_argument('--keep', '-k', help='Output file with integrations to keep', default='integrations.keep.tsv')
	parser.add_argument('--exclude', '-e', help='Output file with integrations to exclude', default='integrations.exclude.tsv')	
	parser.add_argument('--criteria', '-c', help='Criteria to use for filering', nargs='+', required=True)
	args = parser.parse_args()

	crit = Criteria(args.criteria, columns)
	print(f"filtering {args.input} using criteria {crit.criteria}")
	
	with open(args.input, 'r') as in_handle, open(args.keep, 'w') as keep_handle, open(args.exclude, 'w') as exclude_handle:
		in_csv = csv.DictReader(in_handle, delimiter = '\t')
		keep = csv.DictWriter(keep_handle, delimiter = '\t', fieldnames = in_csv.fieldnames)
		keep.writeheader()
		exclude = csv.DictWriter(exclude_handle, delimiter = '\t', fieldnames = in_csv.fieldnames)
		exclude.writeheader()
		
		kept, excluded = 0, 0
		
		for row in in_csv:
			if crit.filter_row(row):
				keep.writerow(row)
				kept += 1
			else:
				exclude.writerow(row)	
				excluded += 1

	print(f"found {kept+excluded} candidate integrations: kept {kept} and excuded {excluded}")
	print(f"Saved kept output to {args.keep} and excluded output to {args.exclude}")



class Criteria():
	"""
	Class to store and enforce criteria
	"""
	
	def __init__(self, criteria_list, column_spec=columns):
	
		assert isinstance(criteria_list, list)
		assert isinstance(column_spec, dict)
		for value in column_spec.values():
			assert isinstance(value, set) or value == 'integer'
	
		# assign inputs to self
		self.row_var_name = 'row'
		self.criteria_list = []
		for string in criteria_list:
			self.criteria_list += [i for i in string.split() if i != '']
		self.column_spec = column_spec
		
		# check that criteria are valid
		self.__check_criteria_valid()
		
		
	def filter_row(self, row_original):

		assert isinstance(row_original, dict)
		
		row = dict(row_original)
		
		# when comparing integers, we need to convert strings imported by DictReader into ints
		for int_col in [i for i in self.column_spec.keys() if self.column_spec[i] == 'integer']:
			try:
				row[int_col] = int(row[int_col] )
			except ValueError:
				row[int_col] = np.nan
		
		return eval(self.criteria)		
	
	def __check_criteria_valid(self):
		"""
		Check if the criteria list supplied is valid.  Criteria list is valid if if has at least three
		elements, all elements are strings, the number of open brackets is the same as the number of closed brackets,
		and each element in it is a valid element (according to the definitions for each element type)
		"""
		
		# criteria list can be just 'True' to include all integrations
		if len(self.criteria_list) == 1:
			assert self.criteria_list[0] == "True"
			self.criteria = "True"
			return
		
		self.comparators = {'>', '<', '<=', '>=', '==', '!='}
		self.logical = {'and', 'or', 'not'}
		self.paren = {'(', ')'}
		self.categorical_values = set()
		for values in self.column_spec.values():
			if isinstance(values, set):
				self.categorical_values = self.categorical_values | values
 
		
		# must be at least three elments in criteria list
		if len(self.criteria_list) < 3:
			raise ValueError("Invalid criteria: {" ".join(self.criteria_list)}")
		
		# every element in list must be a string
		if not all([isinstance(i, str) for i in self.criteria_list]):
			raise ValueError("Invalid criteria: {" ".join(self.criteria_list)}")
			
		# split any brackets from the end of terms
		self.__split_brackets()
				
		# check that the number of open brackets equals the number of closed brackets
		open_brackets = sum([1 for i in self.criteria_list if i == '('])
		closed_brackets = sum([1 for i in self.criteria_list if i == ')'])
		if open_brackets != closed_brackets:
			raise ValueError("Invalid criteria: "
					"Number of open brackets must equal number of closed brackets" 
					f"(found {open_brackets} open brackets and {closed_brackets} closed brackets)")		
			
		# remove any redundant brackets
		self.__remove_redundant_brackets()
					
		# each element is valid
		for i in range(len(self.criteria_list)):
			if self.criteria_list[i] == '(':
				self.__check_open_paren(i)

			elif self.criteria_list[i] == ')':
				self.__check_closed_paren(i)

			elif self.criteria_list[i] in self.comparators:
				self.__check_comparator(i)

			elif self.criteria_list[i] == 'not':
				self.__check_not(i)
			
			elif self.criteria_list[i] in {'and', 'or'}:
				self.__check_and_or(i)
			
			elif self.criteria_list[i] in self.column_spec.keys():
				self.__check_column_name(i)
				
			elif self.criteria_list[i] in self.categorical_values or self.__valid_int(self.criteria_list[i]):
				self.__check_value(i)
			
			else:
				raise ValueError("Invalid criteria: "
					f'Disallowed term {self.criteria_list[i]} was found in criteria')
					
		# create string used for filtering
		self.criteria = list(self.criteria_list)
		for i in range(len(self.criteria_list)):
			
			# need to replace row names with row['row_name']
			if self.criteria[i] in self.column_spec.keys():
				self.criteria[i] = f"{self.row_var_name}['{self.criteria_list[i]}']"
			
			# need to quote values that are strings
			if self.criteria_list[i] in self.comparators:
				col = self.criteria_list[i-1]
				val = self.criteria_list[i+1]
				if self.column_spec[col] != 'integer':
					self.criteria[i+1] = f"'{val}'"
					
		self.criteria = " ".join(self.criteria)
		

	def __split_brackets(self):
		"""
		If there are brackets at the start or end of any terms, split these into a separate term
		Only split for open brackets at the start and closed brackets at the end
		"""
		brackets = self.__has_unsplit_brackets()
		while len(brackets) != 0:
			ins = 0
			for i in sorted(brackets.keys()):
				if brackets[i] == 'start':
					self.criteria_list.insert(i+ins, "(")
					self.criteria_list[i+ins+1] = self.criteria_list[i+ins+1][1:]
				else:
					self.criteria_list.insert(i+ins+1, ")")
					self.criteria_list[i+ins] = self.criteria_list[i+ins][:-1]
				ins += 1
			brackets = self.__has_unsplit_brackets()		
		
	def __has_unsplit_brackets(self):
		"""
		Check if there are brackets at the start or end of any terms
		Only check for open brackets at the start and closed brackets at the end
		"""
		brackets = {}
		for i in range(len(self.criteria_list)):
			if len(self.criteria_list[i]) == 1:
				continue
			if i in brackets.keys():
				continue
			if self.criteria_list[i][0] == "(":
				brackets[i] = "start"
			if i in brackets.keys():
				continue
			if self.criteria_list[i][-1] == ")" :
				brackets[i] = "end"
				
		return brackets
		

	def __check_value(self, i):
		"""
		Check if element at index i is a valid value
		A valid value is followed by nothing, a closed parenthesis, 'and', or 'or'
		A valid value is part of a valid criterion
		"""
		assert (self.criteria_list[i] in self.categorical_values) or self.__valid_int(self.criteria_list[i])
		
		# followed by nothing, a closed parenthesis, 'and', or 'or'
		if i + 1 <= (len(self.criteria_list) - 1):
			if self.criteria_list[i+1] not in {')', 'and', 'or'}:
				raise ValueError(f"Invalid criteria: "
					"values must be followed by ')', 'and' or 'or'")
					
		# part of a valid criterion
		self.__check_valid_criterion(i-2)
		
	
	def __check_column_name(self, i):
		"""
		Check if element at index i is a valid column name
		A valid column name is preceded by nothing, and open parenthesis, 'and', 'or' or 'not'
		A valid column name is part of a valid criterion
		"""		
		assert self.criteria_list[i] in self.column_spec.keys()
		
		# if previous term, '(' 'and', 'or' or 'not'
		if (i - 1) > 0:
			if self.criteria_list[i-1] not in self.logical | { "(" }:
				raise ValueError(f"Invalid criteria: "
					"column names must be preceded by '(', 'not', 'and' or 'or'")
		
		# part of valid criterion
		self.__check_valid_criterion(i)

	def __check_not(self, i):
		"""
		Check if element at index i is a valid 'not'
		A valid 'not' is preceded by nothing, and open paranthesis, 'and', or 'or'
		A valid 'not' is followed by an open parenthesis or a valid criterion
		"""	
		assert self.criteria_list[i] == 'not'	
		
		# precceeding element, if it exists, must be 'and' or 'or'
		if (i - 1) > 0:
			if self.criteria_list[i-1] not in {'(' 'and', 'or'}:
				raise ValueError("Invalid criteria: "
					"'not' must be preceded by '(', 'and' or 'or'"
				)
		
		# following element must be open bracket or valid criterion
		if i >= len(self.criteria_list) - 1:
			raise ValueError(f"Invalid criteria: "
				"'not' cannot be last element in criteria")
				
		if self.criteria_list[i+1] != "(":
			try:
				self.__check_valid_criterion(i+1)
			except ValueError:
				raise ValueError("Invalid criteria: "
					"'not' must be followed by '(' or a valid criterion"
				)
				
		
		
	def __check_and_or(self, i):	
		"""
		Check if element at index i is a valid 'and' or 'or'
		A valid 'and' or 'or' is preceeded by a valid criterion or a closed bracket
		A valid 'and' or 'or' is followed by a a valid criterion, 'not' or an open bracket
		"""	
		assert self.criteria_list[i] in {'or', 'and'}
		
		# check not first or last element
		if (i - 1) < 0:
			raise ValueError("Invalid criteria: "
					f"'or' or 'and' is first element")
		if i  >= len(self.criteria_list) - 1:
			raise ValueError("Invalid criteria: "
					f"'or' or 'and' is last element")
					
		# check preceeded by closed bracket or valid criterion
		if self.criteria_list[i-1] != ")":
			try:
				self.__check_valid_criterion(i-3)
			except ValueError:
				raise ValueError(f"Invalid criteria"
					"'and' and 'or' must be preceeded by a closed bracket or a valid criterion")
			
		# check followed by open bracket, 'not' or valid criterion
		if not (self.criteria_list[i+1] == "(" or self.criteria_list[i+1] == 'not'):
			try:
				self.__check_valid_criterion(i+1)
			except ValueError:
				raise ValueError(f"Invalid criteria"
					"'and' and 'or' must be followed by by an open bracket, 'not', or a valid criterion")					
					
					

	def __check_comparator(self, i):
		"""
		Check if element at index i is a valid comparator ('==', '!=', '<', '>', '>=', '<=')
		A valid compartor is preceeded by a valid column name
		A valid comprator is followed by a value valid for the column name
		"""
		assert self.criteria_list[i] in self.comparators
		
		if (i - 1) < 0:
			raise ValueError("Invalid criteria: "
					f"Comparator {self.criteria_list[i]} has nothing preceeding it, " 
					"but it should be preceeded by a column name")
		if i  >= (len(self.criteria_list) - 1):
			raise ValueError("Invalid criteria: "
						f"Comparator {self.criteria_list[i]} has nothing following it, " 
						"but it should be followed by a value")	
		self.__check_valid_criterion(i-1)
				
	
	def __check_closed_paren(self, i):
		"""
		Check if element at index i is a valid closed parenthesis
		A valid open parenthesis must be preceeded by another element
		A valid open paranthesis must be preceeded by a valid criterion or another closed parenthesis
		A valid open paranthesis must be followed by nothing, a logical operator or another closed parenthesis
		"""
		assert self.criteria_list[i] == ')'
		
		if i == 0:
			raise ValueError("Invalid criteria: "
				"Found a closed bracket ')' that has nothing preceeding it")
		
		# closed paren not preceeded by a valid criterion or another closed paren is invalid
		prev_term = self.criteria_list[i - 1]
		
		# if previous terms is not ")", then must be previous criterion
		if prev_term != ")":
			
			if (i - 3) < 0:
				raise ValueError("Invalid criteria: "
				"Closed parentheses must be preceded by another closed parenthesis or a valid criterion")
			
			try:
				self.__check_valid_criterion(i-3)
			except ValueError:
				raise ValueError("Invalid criteria: "
				"Closed parentheses must be preceded by another closed parenthesis or a valid criterion")
		
		# closed paren not followed by 'or', 'and', or another closed paren is invalid
		if i + 1 <= (len(self.criteria_list) - 1):
			next_term = self.criteria_list[i + 1]
			if not next_term in (self .logical | { ')'}):
				raise ValueError("Invalid criteria: "
				"Closed parentheses must be followed by 'and', 'or', or another closed parenthesis")
		
		
	def __check_open_paren(self, i):
		"""
		Check if element at index i is a valid open parenthesis
		A valid open parenthesis must be followed by another element
		A valid open paranthesis must be followed by 'not', a column name, or another open parenthesis
		A valid open paranthesis must be preceeded by nothing or a logical operator
		"""
		assert self.criteria_list[i] == '('
		# can't be last term in list
		if i == len(self.criteria_list):
			raise ValueError("Invalid criteria: "
			"Found an open bracket that has nothing following it")	
		
		next_term = self.criteria_list[i + 1]
		
		# open paren not followed by column name or another open paren is invalid
		if not (next_term in self.column_spec.keys() or next_term == '(' or next_term == 'not'):
			raise ValueError("Invalid criteria: "
			"Open parentheses must be followed by 'not', a column name or another open parenthesis")
		
		# open paren not preceeded by 'and', 'or', or 'not' or another open paren is invalid
		if i > 0:
			prev_term = self.criteria_list[i - 1]
			if prev_term not in {'and', 'or', 'not', '('}:
				raise ValueError("Invalid criteria: "
				"Open parentheses must be preceeded by 'and', 'or', 'not', or another open parenthesis")
		
			
	def __contains_redundant_brackets(self):
		"""
		Check if a criteria list contains redundant brackets.
		Redundant brackets are '()'
		"""
		for i in range(len(self.criteria_list)):
			if self.criteria_list[i] != '(':
				continue
			if i + 1 not in range(len(self.criteria_list)):
				continue
			if self.criteria_list[i+1] == ')':
				return True
				
		return False
		
	def __remove_redundant_brackets(self):
		while self.__contains_redundant_brackets():
			remove = []
			for i in range(len(self.criteria_list)):
				if self.criteria_list[i] != '(':
					continue
				if i+1 not in range(len(self.criteria_list)):
					continue
				if self.criteria_list[i+1] != ')':
					continue
				remove += [i, i+1]
			for i in reversed(sorted(remove)):
				self.criteria_list.pop(i)
	
	def __check_valid_criterion(self, i):
		"""
		check if three terms in a list, starting at term i, comprise a valid criterion.
		A valid criterion is a list with 3 , where all elements are strings
		
		A valid three-term criterion contatins a column name, a comparator valid for that column
		name, and a value valid for that column name.  
		
		Valid comparators for categorical columns are '!=' and '=='
		Valid comparators for integer columns are  '!=' and '==', '>', '<', '<=', '>='
		
		Valid values for integer columns are strings that can be converted into integers
		Valid values for categorical columns are the set of values defined in self.column_spec for that column name
		"""
		
		try:
			criteria_sublist = self.criteria_list[i:i+3]
		except IndexError:
			raise ValueError(f'Criteria list "{" ".join(self.criteria_list)}" contains an invalid criterion')
		
		assert all([isinstance(i, str) for i in criteria_sublist])
		
		error_str = f"Criterion '{' '.join(criteria_sublist)}' is not valid: "
		
		# check first element is column name
		if criteria_sublist[0] not in self.column_spec.keys():
			raise ValueError(error_str + f"column '{criteria_sublist[0]}' is invalid")
			
		# check second element is a comparator
		if criteria_sublist[1] not in self.comparators:
			raise ValueError(error_str + f"comparator '{criteria_sublist[1]}' is invalid")
			
		# if value is from a set (ie non-integer), the only valid comprators are '==' and '!='
		if self.column_spec[criteria_sublist[0]] != 'integer':
			if criteria_sublist[1] not in {'==', '!='}:
				raise ValueError(error_str + f"only comparators '==' and '!=' are valid for categorical columns, "
													f"but found comparator {criteria_sublist[1]} for categorical column {criteria_sublist[0]}")
		# check third element is a value valid for the column name
		if not self.__valid_value(criteria_sublist[0], criteria_sublist[2]):
			raise ValueError(error_str + f"value '{criteria_sublist[2]}' is "
						f"invalid for column '{criteria_sublist[0]}'")

		
	def __valid_value(self, column_name, value):
		"""
		check if a value is valid for a particular column
		"""	
		assert column_name in self.column_spec.keys()
		# if this is an integer column
		if self.column_spec[column_name] == "integer":
			if not self.__valid_int(value):
				return False
		# if this is a column with a string from a set of allowed values
		else:
			if value not in self.column_spec[column_name]:
				return False
				
		return True
		
	def __valid_int(self, value):
		"""
		Check if a string can be converted to a valid integer
		"""
		assert isinstance(value, str)
		try:
			int(value)
		except ValueError:
			return False
			
		return True
		
	def __repr__(self):
		return "Criteria: " + " ".join(self.criteria_list)

if __name__ == "__main__":
	main(argv)



