import sys
sys.stdout.write('Hello\n')
sys.stdout.write('World\n')
sys.stdout.flush()
sys.stdout.write('\033[2A\rHey Hey')
sys.stdout.flush()
