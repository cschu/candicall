#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import struct
import shutil
import tempfile

import urllib2
import csv
import json

    
CANDISNP_SERVER = 'http://candisnp.tsl.ac.uk' #:8080' 

CANDI_TESTDATA = {
    "ref": "athalianaTair10",
    "data": [{
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
               "alternate_base": "T",
              },
            ],
}
"""
'is_ctga', 
'effect', 
'change', 
'gene', 

"""


SNP_DATA_HEADERS = {                    
                    0: 'chromosome', 
                    1: 'position', 
                    3: 'reference_base',  
                    4: 'alternate_base',
                    7: 'info'               
}

contentHeaders = {
    'Content-Type': 'application/json',
    'Accept': 'text/plain'
}

class SNPEffect(object):
    def __init__(self, string):
        effectAndData = string.strip(')').split('(') 
        # sys.stderr.write(str(effectAndData)+'\n')
        self.effect, effectData = effectAndData[0], effectAndData[1]
        effectData = effectData.split('|')
        assert len(effectData) >= 11, 'Invalid effect field %s' % string
        self.impact, self.fClass, changeDist, self.aaChange, self.aaLength, self.gene, self.transcriptBioType, isCoding, self.transcriptID, exonIntronRank, genotypeNumber, warning = (effectData + [''])[:12]
        if self.effect in ('UPSTREAM', 'DOWNSTREAM'):
            self.transcriptDist = changeDist
            self.codonChange = None
        else:
            self.transcriptDist = None
            self.codonChange = changeDist
        self.isCoding = (isCoding == 'CODING')
        self.exonIntronRank = int(exonIntronRank) if exonIntronRank else 'NA'
        self.genotypeNumber = int(genotypeNumber)
        self.warning = warning if warning else None

        self.is_synonymous = self.effect != 'NON_SYNONYMOUS_CODING'
        if self.effect in ('INTRON', 'INTERGENIC'):
            self.in_cds, self.is_synonymous = False, False
        else:
            self.in_cds = True
        try:
            self.change = '%s/%s' % (self.aaChange[0], self.aaChange[-1])
        except:
            self.change = ''


        pass
    def toDict(self):
        def is_ctga(ref, alt):
            return (ref == 'C' and alt == 'T') or (ref == 'G' and alt == 'A')
        assert hasattr(self, 'is_synonymous'), 'Missing attribute: is_synonymous'
        assert hasattr(self, 'allele_frequency'), 'Missing attribute: allele_frequency'
        assert hasattr(self, 'reference_base'), 'Missing attribute: reference_base' 
        assert hasattr(self, 'alternate_base'), 'Missing attribute: alternate_base' 
        return { "is_synonymous": str(self.is_synonymous).upper(),
                 "change": self.change if self.change else 'NA',
                 "effect": self.effect if self.effect else 'NA',
                 "gene": self.gene if self.gene else 'NA',
                 "allele_freq": self.allele_frequency,
                 "reference_base": self.reference_base,
                 "position": self.position,
                 "chromosome": self.chromosome,
                 "is_ctga": str(is_ctga(self.reference_base, self.alternate_base)).upper(),
                 "in_cds": str(self.in_cds).upper(),
                 "alternate_base": self.alternate_base, } 

        
    pass
        
        

class AnnotatedSNPEffectFactory(object):
    def __init__(self, *args, **kwargs):
        assert len(args) >= 7, 'Not enough values in args (%s)' % args
        self.log = kwargs['log']
        self.chromosome = args[0]
        self.position = int(args[1])
        self.reference_base = args[3]
        self.alternate_base = args[4]
     
        # info = dict([field.split('=') for field in args[7].strip().split(';')])
        # info = dict(field.split('=') for field in args[7].strip().split(';'))
        try:
            info = dict(field.split('=') for field in args[7].strip().split(';') if len(field.split('=')) == 2)
        except:
            sys.stderr.write('Found a problem in your input data. Please make sure it is a VCF produced by SnpEff.\n')
            sys.exit(1)
        # self.log.write('***'+str(info.get('EFF', '@@@')+'***\n'))
        
        self.allele_frequency = float(info.get('AF', '0.0'))
        self.effects = [effect 
                        for effect in info.get('EFF', '').split(',')
                        if effect]
        # assert self.effects, 'No effects in info-string %s' % info.get('EF', '')
        pass

    def getEffects(self, notWanted=('UPSTREAM', 'DOWNSTREAM', 'FRAME_SHIFT')):
        for effectString in self.effects:            
            effect = SNPEffect(effectString)
            if effect.effect not in notWanted:
                for attr in ('chromosome', 'position', 'reference_base', 'alternate_base', 'allele_frequency'):
                    setattr(effect, attr, getattr(self, attr))
                yield effect.toDict()
                

        pass
    pass

def extractCandiDataFromSnpEffVCF(snpEffVCF, fo):
    snps = []
    for line in snpEffVCF:
        if not line.startswith('#'):
            #fo.write(line)
            try:
                asf = AnnotatedSNPEffectFactory(*line.strip().split('\t'), log=fo)
                snpEffects = list(asf.getEffects())                        
            except:
                continue
            #fo.write(str(snpEffects)+'\n')
            snps.extend(snpEffects)
    return snps



def main(argv):

    descr = ''
    parser = argparse.ArgumentParser(description=descr)
    parser.add_argument('--ref', help='The snpEff genome reference.')
    parser.add_argument('snpEff_output', type=str, help='The input file.')
    parser.add_argument('candisnp_html', type=str, help='The output file.')
    args = parser.parse_args()

    fo = open(args.candisnp_html, 'wb')
    """
    import subprocess
    import urllib2
    response = urllib2.urlopen('http://ruup.xyz:8080/monitors')
    html = response.read()


    # p = subprocess.Popen(['wget http://ruup.xyz:8080/monitors'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    # fo.write(p.communicate()[0])
    fo.write(html)
    fo.close()

    return None
    """
    """ 
    # This works.
    candiData = CANDI_TESTDATA 
    candiMessage = json.dumps(candiData)
    """
    candiData = extractCandiDataFromSnpEffVCF(open(args.snpEff_output), fo)
    candiMessage = json.dumps({'ref': args.ref, 'data': candiData})
    
    try:
        request = urllib2.Request(CANDISNP_SERVER + ':8080', 
                                  headers=contentHeaders,
                                  data=candiMessage)
    except:
        sys.stderr.write('There was a problem with your request. Sorry.\n')  


    try:
        response = urllib2.urlopen(request)
        candiURL = response.read()
    except urllib2.URLError, e:
        candiURL = ''
        sys.stderr.write(str(e.reason) + '\n')
   
    if candiURL:
        # body = urllib2.urlopen(candiURL)
        # fo.write(body.read())
        #fo.write('<iframe name="galaxy_main" id="galaxy_main" frameborder="0" style="position: absolute; width: 100%; height: 100%;" src="%s"></iframe>' % candiURL)
        #fo.write('<iframe name="galaxy_main" id="galaxy_main" frameborder="0" style="position: absolute; width: 100%; height: 100%;" src="%s"></iframe>' % candiURL)     
        fo.write('<iframe src="%s" frameborder="0" style="position: absolute; width: 100%%; height: 100%%;"></iframe>\n' % candiURL)
        # fo.write('<iframe src="%s"></iframe>\n' % candiURL)
    else:
        fo.write('I am sorry. CandiSNP does not pick up. Maybe (<a href="%s" target="_blank">try it manually?</a>)\n' % CANDISNP_SERVER)
     
    fo.close()
    pass

# main(sys.argv[1:])

if __name__ == '__main__': main(sys.argv[1:])
