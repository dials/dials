
class _TR:
    '''Singleton class to handle template regular expressions.'''

    def __init__(self):
        self._setup = False
        self._patterns = None

    def __setup(self):

        import re
        patterns = [r'([0-9]{2,12})\.(.*)',
                    r'(.*)\.([0-9]{2,12})_(.*)',
                    r'(.*)\.([0-9]{2,12})(.*)']

        self._joiners = ['.', '_', '']
        self._compiled_patterns = [re.compile(pattern) for pattern in patterns]

        self._setup = True

    def __call__(self, filename):
        '''Get the template, # digits.'''

        if not self._setup:
            self.__setup()

        rfilename = filename[::-1]
        template = None

        for j, cp in enumerate(self._compiled_patterns):
            match = cp.match(rfilename)
            if not match:
                continue
            groups = match.groups()

            if len(groups) == 3:
                exten = '.' + groups[0][::-1]
                digits = groups[1][::-1]
                prefix = groups[2][::-1] + self._joiners[j]
            else:
                exten = ''
                digits = groups[0][::-1]
                prefix = groups[1][::-1] + self._joiners[j]

            template = prefix + ''.join(['#' for d in digits]) + exten
            break

        if not template:
            raise RuntimeError, 'cannpt parse filename %s' % filename

        return template, int(digits)

TR = _TR()

def work_TR():
    questions_answers = {
        'foo_bar_001.img':'foo_bar_###.img',
        'foo_bar001.img':'foo_bar###.img',
        'foo_bar_1.8A_001.img':'foo_bar_1.8A_###.img',
        'foo_bar.001':'foo_bar.###',
        'foo_bar_001.img1000':'foo_bar_###.img1000',
        'foo_bar_00001.img':'foo_bar_#####.img'
        }

    for filename in questions_answers:
        answer = TR(filename)
        assert answer[0] == questions_answers[filename]

def image2template(filename):
    return TR(filename)[0]

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        work_TR()
        print 'OK'
    else:
        for arg in sys.argv[1:]:
            print '%s -> %s' % (arg, TR(arg)[0])
