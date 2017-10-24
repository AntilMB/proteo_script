import urllib.request
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-user', default='')
parser.add_argument('-db', default='')
parser.add_argument('-task', default='', nargs='*')
args = parser.parse_args()

user = args.user
db = args.db
task = '+'.join(args.task)

url = 'http://mascot.ripcm.com/mascot/x-cgi/ms-review.exe?CalledFromForm=1&lo' + \
          'gfile=..%2Flogs%2Fsearches.log&start=-1&howMany=20000&pathToData=&column=0' + \
          '&s0=1&s1=1&s2=1&s3=1&s4=1&s5=1&s6=1&s7=1&s8=1&s9=1&s10=1&s11=1&s12=1&s' + \
          '13=1&s14=1&f0=&f1=&f2=' + db + '&f3=' + user + '&f4=&f5=' + task + '&f' + \
          '6=&f7=&f8=&f9=&f10=&f11=&f12=&f13=&f14='

wp = urllib.request.urlopen(url)

pw = wp.read()

tt = re.findall('(\.\./data/[0-9]{8}/F[0-9]{6}\.dat)', str(pw))

f = open('export.list', 'w')
tt_pared = set([re.search('\.\./data/(.*)', strs).group(1) for strs in tt])
for i in sorted(tt_pared):
        print(i, file=f)

print(len(tt_pared), 'entries were found')