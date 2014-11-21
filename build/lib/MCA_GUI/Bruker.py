#! /usr/bin/env python
# -*- encoding:utf-8 -*-
#
# Copyright (C) Yu 2011-2013, Tra 2014
#
# convert bruker raw v3 format to a two column ascii format
# get Burker data from uxd file
# usage : as a script
# python bruker.py infile.raw output.uxd
# 
# usage as a python module
# from bruker import convert_raw_to_uxd
# convert_raw_to_uxd('infile.raw', 'output.uxd')

import os
import struct
import copy 
import types
import collections
from cStringIO import StringIO
import numpy as np

#
# this implements a ordered user dictionary , 
# serves as an interface for (read-)accessing _metadata using dictionary syntax
#
class Metadata(object):
	
	def __init__(self):
		# use OrderedDictionary if possible ( python > 2.7 )
		if hasattr( collections, 'OrderedDict' ):
			self._metadata = collections.OrderedDict()
		else:
			self._metadata = { } 

	def __getitem__( self, key ):
		return self._metadata[ key ]

	def __setitem__( self, key, val ):
		self._metadata[ key ] = val

	def __delitem__( self, key ):
		del self._metadata[ key ]
	
	def __getattr__(self, key):
		if self._metadata.has_key(key):
			return self._metadata[ key]
		else:
			raise AttributeError('Key\'%s\' does not exists' % key)
	
	
	#def __setattr__(self, key, val):
		#if self._metadata.has_key(key):
			#self._metadata[ key ] = val
		#else:
			#raise AttributeError("Please don't try to create metadata like this." % name)
	
	

class Scan(Metadata):
	def __init__(self):
		Metadata.__init__(self)
		
		self.data = []# array.array('f')          # holding scan data

	def __len__(self):
		return len(self.data)

	def __repr__(self):
		return "<Scan object containing {0:d} data point(s) at {1:d}".format(len(self), id(self))

	def __str__(self):
		return self.pretty_format()

	def pretty_format(self, print_header = True):
		import StringIO
		out = StringIO.StringIO()
		# meta data
		if print_header :
			for k,v in self._metadata.items() :
				if types.FloatType == type(v) :
					print >>out, "_{0:<15} = {1:.6f}".format(k,v)
				else:
					print >>out, "_{0:<15} = {1}".format(k,v)

		# print scan data
		if self.data:
			for (x,y) in self.data:
				print >> out ,  "{0:.6f}\t{1:.6f}".format(x,y)
		return out.getvalue()
	

		
class Dataset(Metadata):
	''' base dataset class defines contents and access specifications. parsing functions should be in derived class'''

	def __init__(self):
		Metadata.__init__(self)
		self.scans =  []
	
	# len(ds) : return number of scans 
	def __len__(self):
		return len( self.scans )
	
	def __str__(self):
		return self.pretty_format()
	
	def __repr__(self):
		if 0 == len(self) :
			return '<Empty dataset at 0x{1:x}>'
		else:
			return '<Dataset containing {0:d} scan(s) at 0x{1:x}>'.format(len(self), id(self))
			
	def pretty_format(self, print_header = True ):
		'''print dataset's header(optional) and each scan'''
		import StringIO
		out = StringIO.StringIO()

		# print dataset header if told so
		if print_header:
			for k,v in self._metadata.items():
				if types.FloatType == type(v) :
					print >> out,  "_{0} = {1:.6f}".format(k,v) # 6 decimal for float
				else:
					print >> out , "_{0} = {1}".format(k,v)
		
		for i in range( len(self.scans) ) :
			if print_header:
				out.write( "\n; ( Data for Range number {0:d} )\n".format(i) )
			out.write( self.scans[i].pretty_format(print_header) + "\n" )
			pass
		
		
		return out.getvalue()

						

	# check whether selected dataset handler can parse given
						
	def validate(self,f=None):
		return False
	
	def parse(self,fh):
		pass

	def merge( self, foreign ):
		c = copy.deepcopy(self)
		c.scans.extend( foreign.scans )
		return c
	
	
class DatasetDiffractPlusV3( Dataset ):
	'''
	DataSet - container for a series of scans
	DataSetBruker
	'''
	
	FILE_HEADER_LENGTH = 712

	#
	# file_header_desc_tbl file header is the first 712 bytes of the
	# file it contains information common to all the following scan
	# 'ranges' it is followed by the header of the first scan range
	#
	# name - name used in meta data dictionary
	# type - code for understood by 'unpack'
	# start -
	# len  - should correspond to the type declared
	#
	file_header_desc_tbl = [
		# ( name                  , type,   start, len )
		( '_FILE_STATUS_CODE'      , 'I'  ,    8,     4 ), # le tiret signifie qu'il va etre supprime/remplace
		( 'RANGE_CNT'             , 'I'  ,   12,     4 ),
		( 'DATE'                  , 'str',   16,  10 ), 
		( 'TIME'                  , 'str',   26,  10 ), 
		( 'USER'                  , 'str',   36,  72 ), 
		( 'SAMPLE'             , 'str',   326, 60 ), 
		( '+SAMPLE'             , 'str',   386, 160 ),  
		( 'GNONIOMETER_RADIUS'    , 'f',     564,  4 ), 
		( 'ANODE_MATERIAL'        , 'str',   608,  4 ), 
		( 'WL1'                   , 'd',     624,  8 ), 
		( 'WL2'                   , 'd',     632,  8 ), 
		( 'WL_UNIT'               , 'str',     656,  4 ), 
		( 'MEASUREMENT TIME'      , 'f',     664,  4 )
	] # end of file_header_desc_tbl
	
	#
	# _header_desc_tbl
	#
	# range header is the first few bytes of each range block ; this
	# length is variable and its length is given in the first byte
	#
	# todo : how to tell what the scan type is ?, e.g. omg-2th, omg,
	# etc
	#
	range_header_desc_tbl = [
			#(name , type , length )
			( 'HEADER_LENGTH',      'I', 0,     4 ), # 0 , 304
			( 'STEPS',              'I',  4,    4 ), # 4
			( 'OMEGA',        'd',     8, 8 ), # 8
			( 'TWOTHETA',       'd',   16,   8 ), # 16
			( 'KHI',          'd',     24, 8 ), # 24
			( 'PHI',          'd',     32, 8 ), # 32
			( 'X',            'd',      40 ,8 ), # 40
			( 'Y',            'd',      48, 8 ), # 48
			( 'Z',            'd',      56,8 ), # 56
			( 'DETECTOR',      'I',      96, 4 ), # 96 0x60
			( 'HIGH_VOLTAGE',       'f',  100,    4 ), # 100
			( 'AMPLIFIER_GAIN',     'f',   104,   4 ), # 10
			( 'AUX1',     'd',   144,   8 ), # 144
			( 'AUX2',     'd',   152,   8 ), # 144
			( 'AUX3',     'd',   160,   8 ), # 144
			( 'SCAN_MODE',          'I',   168,   4 ), # 176
			( 'STEP_SIZE',          'd',   176,   8 ), # 176
			( 'STEP_SIZE_B',          'd',   184,   8 ), # 176
			( 'STEPTIME',      'f',     192, 4 ), # 192
			( '_STEPPING_DRIVE_CODE',      'I',     196, 4 ), 
			( 'TIMESTARTED',                'f',  204,     4 ), # 204
			( 'TEMP_RATE',                'f',     212, 4 ), # 212
			( 'TEMP_DELAY',                'f',    216,  4 ), # 216
			( 'KV',  'I',   224,   4 ), # 224
			( 'MA',  'I',    228,  4 ), # 228
			( 'RANGE_WL',       'd',      240,  8 ), # 240
			( '_VARYINGPARAMS',       'I',      248,  4 ),
			( '_DATUM_LENGTH',       'I',      252,  4 ),
			( 'SUPPLEMENT_HEADER_SIZE','I',  256, 4 ) # 256
			]   # end of range_header_desc_tbl
			
	tbl_stepping_drives = {
		0 : ( "locked coupled" , "TWOTHETA" ),
		1 : ( "unlocked coupled" , "TWOTHETA" ),
		2 : ( "detector scan" , "TWOTHETA" ),
		3 : ( "rocking curve" , "OMEGA" ),
		4 : ( "khi scan" , "KHI" ),
		5 : ( "phi scan" , "PHI" ),
		6 : ( "x-scan" , "X" ),
		7 : ( "y-scan" , "Y" ),
		8 : ( "z-scan" , "Z" ),
		9 : ( "aux1 scan" , "AUX1" ),
		10 : ( "aux2 scan" , "AUX2" ),
		11 : ( "aux3 scan" , "AUX3" ),
		12 : ( "psi scan" , "TWOTHETA" ),
		13 : ( "hkl scan" , "TWOTHETA" ),
		129 : ( "psd fixed scan" , "TWOTHETA" ),
		130 : ( "psd fast scan" "TWOTHETA" )
		}

	def __init__(self, ifh = None ):
		Dataset.__init__( self )
		self.scans =  [] # a list of scans
		if ifh:
			ifh.seek(0,0)
			self.parse( ifh )
			self.determin_dataset_type()

	#
	# only bruker raw file V3 is served
	#
	def validate(self, fh ):
		pos = fh.tell() # remember where it was
		fh.seek( 0, os.SEEK_SET)
		is_rawfile = ( "RAW1." == fh.read( 5 ) )
		is_v3 = ( "01" == fh.read( 2 ) )
		fh.seek( pos, os.SEEK_SET) # just being nice
		return is_rawfile and is_v3
	

	# determin_dataset_type and set stepping_drive1,2 
	def determin_dataset_type(self):
		
		if len( self.scans ) < 1 :
			raise Exception( "Empty Dataset" )
			
		if len( self.scans ) == 1 : 
			self['TYPE'] = 'SingleScanPlot'
			t = self.scans[0]['_STEPPING_DRIVE_CODE']
			self['SCAN_TYPE'], self['STEPPING_DRIVE1'] = self.tbl_stepping_drives[t]
			return
		
		a = self.scans[0]
		b = self.scans[1]
		
		if not a['TYPE'] == b['TYPE'] :
			raise Exception( " Donno how to deal with this kinds of scan " )
		
		if a['_STEPPING_DRIVE_CODE'] in [13] :
			raise Exception( " Donno how to deal with this kinds of scan " )
			
		# PSD scans
		if a['_STEPPING_DRIVE_CODE'] in [129, 130] :
			self['TYPE'] = 'RSMPlot'
			self['STEPPING_DRIVE1'] = 'OMEGA'
			self['STEPPING_DRIVE2'] = 'TWOTHETA'
			return
			
		# NOTE :
		# it is assumed that only one axis move during each range
		# so that it is safe to say for 2d scan
		self['TYPE'] = 'TwoDPlot'
		for drv in [ 'KHI', 'PHI', 'X', 'Y', 'Z','AUX1','AUX2','AUX3' ]:
			if a[drv] != b[drv]:
				self['STEPPING_DRIVE1'] = drv
				self['STEPPING_DRIVE2'] = a['STEPPING_DRIVE']
		
		
		
		pass
	def parse(self, ifh ):
		# read into seekable buffer
		f = StringIO( ifh.read() )
		f.seek( 0 , os.SEEK_SET )
		
		# valid file type signature
		if not self.validate(f):
			raise Exception("The file format is not of 'diffract plus raw file version 3'.")
		
		
		#(key,type,start,len) in list_table
		for (k,t,s,l) in self.file_header_desc_tbl :
			f.seek(s, os.SEEK_SET)
			buf = f.read(l)
			if 'str' == t:
				self[k] = buf.rstrip('\0')
			elif 'c' == t :
				# TODO rather ackward 
				self[k] = ord( struct.unpack(t,buf)[0] )
			else:
				self[k] = struct.unpack(t,buf)[0]

		if   1 == self['_FILE_STATUS_CODE'] :
			self['FILE_STATUS'] = "done"
		elif 2 == self['_FILE_STATUS_CODE'] :
			self['FILE_STATUS'] = "active"
		elif 3 == self['_FILE_STATUS_CODE'] :
			self['FILE_STATUS'] = "aborted"
		elif 4 == self['_FILE_STATUS_CODE'] :
			self['FILE_STATUS'] = "interrupted"

		
		# beginning of first range
		f.seek(self.FILE_HEADER_LENGTH, os.SEEK_SET)
		
		range_start = self.FILE_HEADER_LENGTH # start of first range
		for i in range( self['RANGE_CNT'] ):
			scn = Scan()
			scn['SEQ'] = i
			
			# read range headers
			for (k,t,s,l) in  self.range_header_desc_tbl :
				f.seek(s + range_start, os.SEEK_SET)
				buf = f.read(l)
				if 'str' == t :
					scn[k] = buf.rstrip('\0')
				elif 'c' == t :
					scn[k] = ord( struct.unpack( t, buf )[0] )
				else:
					scn[k] = struct.unpack(t,buf)[0]
			
			# some constraint described in file-exchange help file
			assert scn['_VARYINGPARAMS'] == 0,"non-conforming file format: more than 1 varying parameters in one range " 
			assert scn['_DATUM_LENGTH'] == 4,"non-conforming file format : datum length more than 4 byte " 
			assert scn['_STEPPING_DRIVE_CODE'] not in [9,10,11] , "non-conforming file format, using AUX* drives" 
			
			# seek to the start of header
			f.seek( range_start + scn['HEADER_LENGTH'], os.SEEK_SET )       
			# skip supplement header if exists
			if scn['SUPPLEMENT_HEADER_SIZE'] > 0 :
				f.read( scn['SUPPLEMENT_HEADER_SIZE'] ) 
			
			scn['TYPE'] = self.tbl_stepping_drives[ scn['_STEPPING_DRIVE_CODE']  ][0]
			scn['STEPPING_DRIVE'] = self.tbl_stepping_drives[ scn['_STEPPING_DRIVE_CODE']  ][1]
			x = scn[scn['STEPPING_DRIVE']] 
			xstep  = scn['STEP_SIZE']
			
			for i in range( scn['STEPS'] ):
				y = struct.unpack('f',f.read(4))[0]
				scn.data.append((x,y))
				#scn.xdata.append( x )
				#scn.xdata.append( y )
				x = x + xstep
			
			self.scans.append(scn)
			range_start = f.tell()
		#end of range
	# end of def parse(self, ifstream )

def convert_raw_to_uxd( ifn, ofn ):
	print ifn, ofn
	ds = DatasetDiffractPlusV3( open( ifn , 'rb') )
	with open(ofn, 'wb') as ofh:
		ofh.write(';(content of file %s)\n' % ifn)  
		ofh.write(ds.pretty_format(print_header = True ))

def get_Bruker(uxd_file):
	#-- Read UXD file, output: theta file, 2Theta file and intensity file
	try:
		data = open(uxd_file,'r')
	except:
		print "No such a Data file! Please check if you entered the correct path in ini file."
	omega = []
	tth_2d   = []
	intensity_2d = []
	tth = []            
	intensity=[]
	for lines in data.readlines():
		if not lines.strip():
			continue
		elif lines.startswith("_"):
			if lines.startswith("_OMEGA"):
				line=lines.split("=")
				theta=line[1]
			else:
				continue
		elif lines.startswith(";"):
			if lines.startswith("; ( Data for Range number") or lines.startswith("; Data for range"):
				tth_2d.append(tth)
				intensity_2d.append(intensity)
				tth = []            
				intensity=[]
		else:
			line=lines.split()
			dTheta=line[0]
			i=line[1]
			omega.append(float(theta))
			tth.append(float(dTheta))
			intensity.append(float(i))
	tth_2d.append(tth)
	intensity_2d.append(intensity)
	data.close()
	
	tth_2d = tth_2d[1:]
	tth_2d = np.asarray(tth_2d)
	intensity_2d = intensity_2d[1:]
	intensity_2d = np.asarray(intensity_2d)
	omega = np.asarray(omega)
	omega = omega.reshape(intensity_2d.shape)
	return {"omega":omega, "tth":tth_2d, "data":intensity_2d}

if '__main__' == __name__:
    import sys
    if not 3 == len( sys.argv ):
        print >> sys.stderr, "this script accept two and only two arguments"
        sys.exit('')
        
    _, ifn, ofn = sys.argv
    
    convert_raw_to_uxd( ifn, ofn ) 