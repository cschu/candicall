#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile

import httplib
import csv
import json

#os.environ['http_proxy'] = ''

import urllib


TESTDATA = {"ref": "aThalianaTAIR10",
 "data": [
    {
        "is_synonymous": "FALSE",
        "change": "-",
        "effect": "INTERGENIC",
        "gene": "-",
        "allele_freq": 0.411764706,
        "reference_base": "C",
        "position": 59214,
        "chromosome": "1",
        "is_ctga": "FALSE",
        "in_cds": "FALSE",
        "alternate_base": "A"
    },
    {
        "is_synonymous": "FALSE",
        "change": "-",
        "effect": "INTERGENIC",
        "gene": "-",
        "allele_freq": 0.272727273,
        "reference_base": "C",
        "position": 78140,
        "chromosome": "1",
        "is_ctga": "TRUE",
        "in_cds": "FALSE",
        "alternate_base": "T"
    },
    {
        "is_synonymous": "FALSE",
        "change": "E/K",
        "effect": "NON_SYNONYMOUS_CODING",
        "gene": "AT1G01320",
        "allele_freq": 0.291666667,
        "reference_base": "C",
        "position": 126483,
        "chromosome": "1",
        "is_ctga": "TRUE",
        "in_cds": "TRUE",
        "alternate_base": "T"
    }]}



def main(argv):

    """<command interpreter="python">candicall.py
                --message="${candiMessage.value}"
                $candicall_tsv
        </command>
    """
    
    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--message', default='Hi Candi!', help='A message to Candi.')
    parser.add_argument('candicall_tsv', type=str, help='The output file.')
    args = parser.parse_args()

    
    fo = open(args.candicall_tsv, 'wb')
    fo.write('%s\n' % args.message)

    
    headers = {
               'Content-Type': 'application/json',
               #'Content-type': 'application/x-www-form-urlencoded',
               'Accept': 'text/plain'
              }
    """params = urllib.urlencode({'@number': 12524, '@type': 'issue', '@action': 'show'})
    """

    inputfile = '/tsl/services/galaxy/dist/galaxy-dist/database/files/024/dataset_24818.dat'
    data = []
    for row in csv.DictReader(open(inputfile), delimiter='\t', quotechar='"'):
        data.append(json.dumps(row))
    # fo.write('\n'.join([string for string in data]) + '\n')
    candicall =  ','.join(data)
    #params = urllib.urlencode({'candicall': candicall})
    # params = urllib.urlencode(json.dumps(TESTDATA))
    params = json.dumps(TESTDATA)
    fo.write(params + '\n')
    
    # url = 'v0518.nbi.ac.uk'
    
    # url = 'darkjade.net' #'galaxy.tsl.ac.uk'
    url = 'candisnp.tsl.ac.uk:8080'

    import urllib2

    req = urllib2.Request("http://"+url, headers = headers,data = params)
    f = urllib2.urlopen(req)
    fo.write('%s\n' % f.read())

    """
    conn = httplib.HTTPConnection(url)
    conn.request('POST', '', params, headers)
    # conn.request('GET', '')
    r1 = conn.getresponse()
    fo.write('%s\n%s\n' % (r1.status, r1.reason))
    data1 = r1.read()
    fo.write('%s\n' % data1)
    fo.write('LEngth = %i\n' % len(candicall))
    conn.close()
    """
    fo.close()
    pass


 
    
# main(sys.argv[1:])

if __name__ == '__main__': main(sys.argv[1:])
