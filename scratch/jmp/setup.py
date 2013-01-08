#!/usr/bin/env python

import subprocess

from distutils.errors import CompileError

from distutils.command.build_ext import build_ext as BuildExtCommand
from distutils.command.clean import clean as CleanCommand
from distutils.cmd import Command
from distutils.util import get_platform

class SconsBuildExt(Command):
    description = BuildExtCommand.description

    user_options = [
        ('build-lib=', 'b',
         "directory for compiled extension modules"),
        ('build-temp=', 't',
         "directory for temporary files (build by-products)"),
        ('plat-name=', 'p',
         "platform name to cross-compile for, if supported "
         "(default: %s)" % get_platform())]

    def initialize_options(self):
        self.build_lib = None
        self.plat_name = None
        self.build_temp = None        
    
    def finalize_options(self):
        from distutils import sysconfig
        self.set_undefined_options('build',
                                  ('build_lib', 'build_lib'),
                                  ('build_temp', 'build_temp'),
                                  ('plat_name', 'plat_name'))
    
    def get_source_files(self): 
        return []

    def run(self):
        try:
            subprocess.check_call(['libtbx.scons', 
                                   '--build-base={0}'.format(self.build_lib)])
            #subprocess.check_call(['libtbx.scons',
            #                       '--install',
            #                       '--prefix={0}'.format(self.build_lib)])
        except subprocess.CalledProcessError:
            raise CompileError("Error while building Python Extensions")
        self.extensions=[]

class SconsClean(CleanCommand):
    def run(self):
        CleanCommand.run(self)
        try:
             subprocess.check_call(['libtbx.scons', 
                                    '--build-base={0}'.format(self.build_temp), 
                                    '--remove'])
        except subprocess.CalledProcessError:
            raise CompileError("Error while cleaning Python Extensions")


def setup_package():
    from distutils.core import setup, Extension

    setup(name='dials',
          version='0.1.0',
          description='Dials Integration',
          packages=['dials', 
                    'dials.old',
                    'dials.old.reflection',
                    'dials.util',
                    'dials.io',
                    'dials.equipment',
                    'dials.geometry',
                    'dials.geometry.transform'],
          
          scripts=['dials/bin/generate_spot_positions.py'],
          
          ext_modules=['dials'],
          cmdclass = { 'build_ext' : SconsBuildExt, 'clean': SconsClean }
     )

if __name__ == '__main__':
    setup_package()
