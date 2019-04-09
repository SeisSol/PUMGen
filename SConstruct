#! /usr/bin/python

# @file
#  This file is part of PUMGen
#
#  For conditions of distribution and use, please see the copyright
#  notice in the file 'COPYING' at the root directory of this package
#  and the copyright notice at https://github.com/SeisSol/PUMGen
#
# @copyright 2017 Technical University of Munich
# @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
#

import os
import sys

import libs
import utils.compiler
import utils.variables

# set possible variables
vars = utils.variables.Variables()

# Add build type
vars.AddBuildType()

# Add prefix path
vars.AddPrefixPathVariable()

# Add compiler variables
vars.AddCompilerVariable()

# PUMGen specific variables
vars.AddVariables(
	PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),

	EnumVariable( 'logLevel',
		'logging level. \'debug\' prints all information available, \'info\' prints information at runtime (time step, plot number), \'warning\' prints warnings during runtime, \'error\' is most basic and prints errors only',
		'info',
		allowed_values=('debug', 'info', 'warning', 'error')
	),

	BoolVariable( 'netcdf', 'compile netCDF support (required to read old SeisSol meshes)',
		False,
	),

	BoolVariable( 'BeforeSim11', 'compile using simModeler up to 11 (enables  Mesh_Export and change call in setSurfaceShapeMetric)',
		False,
	),

	BoolVariable( 'simModSuite', 'compile with support for simModSuite from Simmetrix',
		False
	),

	( 'mpiLib', 'MPI library against this is linked (only required for simModSuite)',
		'mpich2',
		None, None
	),

	PathVariable( 'cc',
		'C compiler (default: mpicc)',
		None,
		PathVariable.PathAccept
	),

	PathVariable( 'cxx',
		'C++ compiler (default: mpicxx)',
		None,
		PathVariable.PathAccept
	),

	BoolVariable( 'useEnv',
		'set variables set in the execution environment',
		True
	)
)

vars.ParseVariableFile()

# Create the environment
env = Environment(variables=vars)

# generate help text
vars.SetHelpText(env)

# Check for any unknown (maybe misspelled) variables
vars.CheckUnknownVariables(env)

# Set environment
if env['useEnv']:
	env['ENV'] = os.environ

#
# precompiler, compiler and linker flags
#

# Set compiler
env['CC'] = 'mpicc'
env['CXX'] = 'mpicxx'
vars.SetCompiler(env)

# set level of logger
if env['logLevel'] == 'debug':
	env.Append(CPPDEFINES=['LOG_LEVEL=3'])
elif env['logLevel'] == 'info':
	env.Append(CPPDEFINES=['LOG_LEVEL=2'])
elif env['logLevel'] == 'warning':
	env.Append(CPPDEFINES=['LOG_LEVEL=1'])
elif env['logLevel'] == 'error':
	env.Append(CPPDEFINES=['LOG_LEVEL=0'])
else:
	assert(false)

# compiler flags for generated kernels
env.Append(CXXFLAGS = ['-Wall', '-ansi', '-std=c++0x'])
if utils.compiler.optimizationEnabled(env):
	env.Append(CPPDEFINES=['NDEBUG'])
	env.Append(CXXFLAGS=['-O2'])
else:
	env.Append(CXXFLAGS=['-O0'])
if utils.compiler.debugEnabled(env):
	env.Append(CXXFLAGS=['-g'])

# add pathname to the list of directories which are search for include
env.Append(CPPPATH=['#/src'])

# Enable openmp
env.Append(CXXFLAGS = ['-fopenmp'])
env.Append(LINKFLAGS= ['-fopenmp'])

# Set prefix pathes for libraries
vars.SetPrefixPathes(env)

# utils library
env.Append(CPPPATH=['#/submodules'])

# APF (one of the Zoltan functions is required)
libs.find(env, 'apf', simmetrix=env['simModSuite'], zoltan=True)
env.Tool('Hdf5Tool', required=True, parallel=True)

# netCDF
env['use_netcdf'] = libs.find(env, 'netcdf', required=env['netcdf'], parallel=True)
if env['use_netcdf']:
	env.Append(CPPDEFINES=['USE_NETCDF'])

# SimModSuite
env['use_simmodsuite'] = libs.find(env, 'simmodsuite', required=env['simModSuite'], mpiLib=env['mpiLib'])
if env['use_simmodsuite']:
	env.Append(CPPDEFINES=['USE_SIMMOD'])

if env['BeforeSim11']:
   env.Append(CPPDEFINES=['BeforeSim11'])

# get the source files
env.sourceFiles = []

Export('env')
SConscript('src/SConscript', variant_dir='#/'+env['buildDir']+'/src', src_dir='#/src', duplicate=0)
Import('env')

# build tools
env.Program('#/'+env['buildDir']+'/pumgen', env.sourceFiles)
