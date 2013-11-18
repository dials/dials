from __future__ import division

def run():
  from dials.util.config import CompletionGenerator
  gen = CompletionGenerator()
  gen.generate()

try:
  run()
except Exception, e:
  pass
