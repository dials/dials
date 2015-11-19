import sys
import cmd

class clui(cmd.Cmd):
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = 'DIALS> '
        self._datablock = None
        self._strong = None
        self._experiments = None
        self._indexed = None
        return

    def do_import(self, args):
        '''Do some importing'''
        print 'Hello from importing!'
        print 'Args were:', args
        self._datablock = 'datablock.json'
        return

    def do_findspots(self, args):
        '''Do some spot finding'''
        print 'Hello from spot finding!'
        print 'Args were:', args
        self._strong = 'strong.pickle'
        return

    def do_index(self, args):
        '''Do some indexing'''
        print 'Hello from indexing!'
        print 'Args were:', args
        self._indexed = 'indexed.pickle'
        self._experiments = 'experiments.json'
        return

    def do_status(self, args):
        '''Show the status'''
        print 'Datablock:     %s' % self._datablock
        print 'Strong spots:  %s' % self._strong
        print 'Experiments:   %s' % self._experiments
        print 'Indexed spots: %s' % self._indexed

    def do_quit(self, args):
        import sys
        sys.exit(0)

if __name__ == '__main__':
    interp = clui()
    interp.cmdloop('DIALS command interpreter')
