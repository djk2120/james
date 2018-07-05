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
  ct = 1
  for name in names:
    if len(name)==3:
      if 'and' in name:
        ct=0
    if ct==1:
        new_str=re.sub(name,name+',',new_str)
    ct += 1
        
  return new_str


def main():

  args = sys.argv[1:]

  if not args:
    print 'usage: file [file ...]'
    sys.exit(1)


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
