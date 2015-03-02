#!/usr/bin/python
import re
import os
import os.path
import sys
import re
import gzip
import argparse

def read_coords( f ):
    # Skip two lines of header
    f.readline()
    f.readline()
    # Read the epoch
    l = f.readline()
    m = re.search(r'\sEPOCH\:\s(\d\d\d\d\-\d\d\-\d\d\s\d\d\:\d\d\:\d\d)\s*$',l)
    if not m:
        raise ValueError('Invalid epoch')
    epoch=m.group(1).replace(' ','T')
    # Skip to first station
    f.readline()
    f.readline()
    f.readline()
    # Read the station data
    for l in f:
        name=l[5:9].strip()
        if not name:
            continue
        flag=l[66:].strip()
        xyz=l[20:66].split()
        if len(xyz) != 3:
            raise ValueError('Invalid coordinate')
        xyz = '\t'.join(xyz)
        yield name, ('\t'.join((name,epoch,xyz,flag)))+'\n'

def dirkey(x):
    return '0'*(10-len(re.match(r'\d*',x).group(0)))+x

ap = argparse.ArgumentParser(description='Extract coordinate time series from directory .../prefix###.CRD.Z')
ap.add_argument('output_file',type=str,help='Filename can include {code}, replaced with station id')
ap.add_argument('root_dir',type=str,default='.')
ap.add_argument('-p','--prefix',help='File name prefix',type=str,default='h1')

args=ap.parse_args()
prefix = args.prefix.lower()
outputfile=args.output_file

usecode= '{code}' in outputfile
data={}

for parent,dirs,files in os.walk(args.root_dir):
    dirs.sort(key=dirkey)
    files.sort(key=dirkey)
    for f in files:
        fref = f.lower()
        if fref.startswith(prefix) and fref.endswith('crd.gz'):
            filename=os.path.join(parent,f)
            print "Processing",filename
            try:
                gz=gzip.open(filename)
                for name,c in read_coords(gz):
                    name = name if usecode else ''
                    if name not in data:
                        data[name]=[]
                    data[name].append(c)
            except:
                msg = str(sys.exc_info()[1])
                print "Error:",filename,": ",msg

for code in data.keys():
    of = outputfile.replace('{code}',code)
    with open(of,'w') as out:
        out.write("name\tepoch\tx\ty\tz\tflag\n")
        out.writelines(data[code])
                    


