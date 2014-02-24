#import re
#import sys
#import os
import numpy as np

def readcol(file,skipline=0l,numline=-1l,format=[1,2],verbose=False, twod=False):
    """Routine to emulate IDL readcol.pro"""
    if verbose:
        print ' > Readcol: being talkative ', verbose
    if verbose:
        print ' > Readcol routine. D.Pietrobon Oct 2013'
    if verbose:
        print ' > Readcol: reading file "%s"' %file
    nlines   = int( get_lineNum(file,verbose=verbose) )
    skipline = int( skipline )
    numline  = int( numline )
    #print nlines, skipline, numline
    if verbose:
        print ' > Readcol: total number of lines: %s' %nlines
    try:
        file_content = open(file, 'r')
        
        if skipline > 0:
            if verbose:
                print ' > Readcol: skipping %s lines...' %skipline
        if numline < 0:
            numline = nlines-skipline
        else:
            if numline > (nlines-skipline):
                print ' > Readcol: WARNING, not enough lines. Reading %s lines:' %nlines-skipline
                numline = nlines-skipline
        # Getting format
        ifields = np.asarray( format, np.int16 )
        nfields = len( format )
        if verbose:
            print ' > readcol: number of fields to be read = %s ' %nfields
        if verbose:
            print ' > readcol: fields read = %s' %ifields
        ifields[:] -= 1
        values = np.zeros( (numline, nfields) )
#        for iline in np.arange(np.max( [nlines,numline] ))+skipline:
        cnt = 0
        for iline,line in enumerate( file_content ):
            if (iline >= skipline) and (iline < skipline+numline):
                if verbose:
                    print ' > line %s' %iline
                    print ' >      %s' %line
                line = line.strip('\n')
                entries = np.asarray( line.split() )
                #print ifields
                #print entries[ ifields[0] ]
                #print entries[ ifields ]
                if verbose:
                    print entries[ [ifields] ]
                #print entries[ifields]
                entries = np.asarray( entries[ifields] )
                if verbose:
                    print entries
                #values.append(entries)
                values[cnt,:] = entries
                if verbose:
                    print values[cnt,:]
                cnt += 1


        file_content.close()
        #print [ values[:,icol] for icol in np.arange(nfields) ]
        if twod:
            return values
        else:
            return [ values[:,icol] for icol in np.arange(nfields) ]

    except IOError:
        print "File not found: %s" %file
        
# ----------------------------------------------------------------------
def get_lineNum(file,verbose=False):
    if verbose:
        print ' > get_lineNum routine. D.Pietrobon Oct 2013'
        nlines = -1
    try:
        f = open(file, 'r')
        nlines=0
        for iline in f:
            nlines+=1
        if verbose:
            print " > get_lineNum: line = %s" %nlines
        return nlines
    except IOError:
        print "File not found: %s" %file


