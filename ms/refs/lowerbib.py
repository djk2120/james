#!/usr/bin/python
# Copyright 2010 Google Inc.
# Licensed under the Apache License, Version 2.0
# http://www.apache.org/licenses/LICENSE-2.0

# Google's Python Class
# http://code.google.com/edu/languages/google-python-class/

import sys
import re

"""Baby Names exercise

Define the extract_names() function below and change main()
to call it.

For writing regex, it's nice to include a copy of the target
text for inspiration.

Here's what the html looks like in the baby.html files:
...
<h3 align="center">Popularity in 1990</h3>
....
<tr align="right"><td>1</td><td>Michael</td><td>Jessica</td>
<tr align="right"><td>2</td><td>Christopher</td><td>Ashley</td>
<tr align="right"><td>3</td><td>Matthew</td><td>Brittany</td>
...

Suggested milestones for incremental development:
 -Extract the year and print it
 -Extract the names and rank numbers and just print them
 -Get the names data into a dict and print it
 -Build the [year, 'name rank', ... ] list and print it
 -Fix main() to use the extract_names list
"""

def lower_case(filename):
  f= open(filename,'rU')
  file_str = f.read()
  new_str = file_str
  matches = re.findall(r'author = {(.+)}',file_str)
  for match in matches:
    names = re.findall(r'(\w+)',match)
  for name in names:
    uppers = re.findall(r'([A-Z][A-Z]+)',name)
    for upper in uppers:
      lwr = upper[0]+upper[1:].lower()
      new_str=re.sub(upper,lwr,new_str)
  return new_str

def extract_names(filename,dict):
  """
  Given a file name for baby.html, returns a list starting with the year string
  followed by the name-rank strings in alphabetical order.
  ['2006', 'Aaliyah 91', Aaron 57', 'Abagail 895', ' ...]
  """
  # +++your code here+++
  # open file and read to string
  f= open(filename,'rU')
  file_str= f.read()

  # extract the year
  # looking for: <h3 align="center">Popularity in 1990</h3>

  # <h3 align="center">Popularity in 1990</h3>
  matchyear= re.search(r'>Popularity in (\d\d\d\d)<',file_str)
  namelist= [matchyear.group(1)]
  yr = matchyear.group(1)

  # extract the names and ranks
  # example: <tr align="right"><td>1</td><td>Michael</td><td>Jessica</td>
  matches= re.findall(r'<tr align="right"><td>(\d+)</td><td>(\w+)</td><td>(\w+)</td>',file_str)

  # build a dict from the matches
  for rank in matches:
    if dict.get(rank[1]):
      if dict.get(rank[1])>rank[0]:
        dict[rank[1]]=(rank[0],yr)
    else:
      dict[rank[1]]=(rank[0],yr)
    if dict.get(rank[2]):
      if dict.get(rank[2])>rank[0]:
        dict[rank[2]]=(rank[0],yr)
    else:
      dict[rank[2]]=(rank[0],yr)

  # sort the matches and output as a list, with the year as header
  for item in sorted(dict.items()):
    namelist.append(item[0]+' '+item[1][0])

  f.close()
  return (namelist,dict)


def main():
  # This command-line parsing code is provided.
  # Make a list of command line arguments, omitting the [0] element
  # which is the script itself.
  args = sys.argv[1:]

  if not args:
    print 'usage: file [file ...]'
    sys.exit(1)

  # Notice the summary flag and remove it from args if it is present.
  #summary = False
  #if args[0] == '--summaryfile':
  #  summary = True
  #  del args[0]

  #trackname = False
  #if args[0] == '--trackname':
  #  trackname = args[1]
  #  del args[0]
  #  del args[0]
  #  if trackname:
  #    print trackname

  # +++your code here+++
  # For each filename, get the names, then either print the text output
  # or write it to a summary file
  dict = {}
  

  #if trackname:
  #  for filename in args:
  #    (names,dict) = extract_names(filename,dict)

  for filename in args: 
    print filename
    new_str = lower_case(filename)
    f = open(filename+'.lower','w')
    f.write(new_str)
    f.close()
    print new_str
#    for name in names:
#      print name
  
if __name__ == '__main__':
  main()
