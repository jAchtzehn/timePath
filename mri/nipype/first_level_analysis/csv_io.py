# _*_ coding: utf-8 _*_
'''
Define two simple functions to write and read from csv/text files.
'''
import csv

# define a function to write to a file
def write(data, csvfile, delimiter=','):
	'''Simply write a row of data to csvfile appending
	contents. There is another method to write several
	rows (csv.writerows).
	The delimiter can be changed from , to any other passing
	the delimiter argument. By default is ','.
	'''
	# windows is a nightmare, normally in serious OS
	# you just write 'ab' in the write options
	with open(csvfile, 'a+b') as f:
		writer = csv.writer(f, delimiter=delimiter, quoting=csv.QUOTE_MINIMAL)
		writer.writerow(data)

# define a function to read the data
def read(csvfile):
	# type: (object) -> object
	'''Read contents from the CSV file csvfile.
	returns : a list with data
	'''
	data = []
	add_to_data = data.append
	with open(csvfile, 'r+b') as f:
		reader = csv.reader(f)
		for row in reader:
			add_to_data(row)
	return data