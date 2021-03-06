#!/usr/bin/env python3
import os
import datetime
import argparse

parser = argparse.ArgumentParser(
    description='Options for generating condor scripts')
parser.add_argument('--beamx',
                    '-x',
                    type=float,
                    nargs='+',
                    required=True,
                    help='List of x values of beam center')
parser.add_argument('--beamz',
                    '-z',
                    type=float,
                    nargs='+',
                    default=[0],
                    help='List of z values of beam center')
parser.add_argument('--TileX',
                    '-X',
                    type=float,
                    nargs='+',
                    default=[50],
                    help='List of x values of Tile X')
parser.add_argument('--TileY',
                    '-Y',
                    type=float,
                    nargs='+',
                    default=[10],
                    help='List of Y values of Tile X')
parser.add_argument('--tilewidth',
                    '-l',
                    type=float,
                    nargs='+',
                    default=[30],
                    help='Tile width')
parser.add_argument('--beamidth',
                    '-w',
                    type=float,
                    nargs='+',
                    default=[1.5],
                    help='beam width parameter')
parser.add_argument('--fiberZ',
                    '-f',
                    type=float,
                    nargs='+',
                    default=[5.2],
                    help='fiber length [m]')
parser.add_argument('--fiberZshift',
                    '-s',
                    type=float,
                    nargs='+',
                    default=[1.7],
                    help='fiber z shift [m]')
parser.add_argument('--absmult',
                    '-a',
                    type=float,
                    nargs='+',
                    default=[1000.],
                    help='List of tile absorption length')
parser.add_argument('--tiledecay',
                    '-td',
                    type=float,
                    nargs='+',
                    default=[2.5],
                    help='tiledecay')
parser.add_argument('--tilerise',
                    '-tr',
                    type=float,
                    nargs='+',
                    default=[0.5],
                    help='tilerise')
parser.add_argument('--Y11abs',
                    '-b',
                    type=float,
                    nargs='+',
                    default=[100.],
                    help='List of tile Y11 absorption length')
parser.add_argument('--LY',
                    '-y',
                    type=float,
                    nargs='+',
                    default=[10.],
                    help='light yield /MeV')
parser.add_argument('--wrapreflect',
                    '-m',
                    type=float,
                    nargs='+',
                    default=[0.985],
                    help='List of wrap reflectivity')
parser.add_argument('--Y11decayTime',
                    '-d',
                    type=float,
                    nargs='+',
                    default=[11.5],
                    help='Y11 time constant')
parser.add_argument('--sipmeff',
                    '-e',
                    type=float,
                    nargs='+',
                    default=[1],
                    help='SiPM eff')
parser.add_argument('--NEvents',
                    '-N',
                    type=int,
                    default=100,
                    help='Number of events to run')
parser.add_argument('--useProton',
                    '-P',
                    type=int,
                    default=1,
                    help='useProton')
parser.add_argument('--handwrap',
                    '-H',
                    type=int,
                    default=1,
                    help='handwrap')
parser.add_argument('--cladlayer',
                    '-c',
                    type=int,
                    default=1,
                    help='cladlayer')
parser.add_argument('--fiberR',
                    '-r',
                    type=float,
                    nargs='+',
                    default=[0.7],
                    help='fiber radius')
parser.add_argument('--holeR',
                    '-R',
                    type=float,
                    nargs='+',
                    default=[1.0],
                    help='hole radius')
parser.add_argument('--spare',
                    '-S',
                    type=float,
                    nargs='+',
                    default=[1.0],
                    help='spare parameter')
parser.add_argument('--prefix',
                    type=str,
                    default='',
                    help='Additional string to place in filename')

args = parser.parse_args()

#BASE_DIR = os.path.abspath(os.environ['CMSSW_BASE'] + '/src/' +
#                           '/HGCalTileSim/condor/')
BASE_DIR = os.path.abspath(os.environ['CMSSW_BASE'] + '/src/' +
                           '/MATHUSLA/condor/')
DATA_DIR = os.path.abspath(BASE_DIR + '/JDL_root/')
Store_DIR = DATA_DIR
#Store_DIR = os.path.abspath('/store/user/yilai/test_T3/FNALNew_trial0329_Holeradius/')

CONDOR_JDL_TEMPLATE = """
universe              = vanilla
Executable            = {0}/condor-LYSquareTrigger_CMSSW.sh
should_transfer_files = NO
Requirements          = TARGET.FileSystemDomain == "privnet" 
#Requirements          = TARGET.Machine =?= "hepcms-hn2.umd.edu"
#Requirements = TARGET.Machine =!= "r510-0-1.privnet" && TARGET.Machine =!= "r510-0-10.privnet" && TARGET.Machine =!= "r510-0-11.privnet" && TARGET.Machine =!= "r510-0-4.privnet" && TARGET.Machine =!= "r510-0-5.privnet" && TARGET.Machine =!= "r510-0-6.privnet" && TARGET.Machine =!= "r510-0-9.privnet" && TARGET.Machine =!= "r540-0-20.privnet" && TARGET.Machine =!= "r540-0-21.privnet" && TARGET.Machine =!= "r720-0-1.privnet" && TARGET.Machine =!= "r720-0-2.privnet" && TARGET.Machine =!= "compute-0-11.privnet" && TARGET.Machine =!= "compute-0-5.privnet" &&  TARGET.Machine =!= "compute-0-7.privnet" &&  TARGET.Machine =!= "hepcms-namenode2.umd.edu"
request_memory        = 1 GB
#Output                = {1}.$(ClusterId).stdout
#Error                 = {1}.$(ClusterId).$(ProcId).stderr
#Log                   = {1}.$(ClusterId).condor
Output                = ../log/out.$(ClusterId).stdout
Error                 = ../log/err.$(ClusterId).$(ProcId).stderr
Log                   = ../log/log.$(ClusterId).condor
Arguments             = {2}
Queue 1
"""
#Requirements = TARGET.Machine =!= "r510-0-1.privnet" && TARGET.Machine =!= "r510-0-10.privnet" && TARGET.Machine =!= "r510-0-11.privnet" && TARGET.Machine =!= "r510-0-4.privnet" && TARGET.Machine =!= "r510-0-5.privnet" && TARGET.Machine =!= "r510-0-6.privnet" && TARGET.Machine =!= "r510-0-9.privnet" && TARGET.Machine =!= "r540-0-20.privnet" && TARGET.Machine =!= "r540-0-21.privnet" && TARGET.Machine =!= "r720-0-1.privnet" && TARGET.Machine =!= "r720-0-2.privnet" && TARGET.Machine =!= "compute-0-11.privnet" && TARGET.Machine =!= "compute-0-5.privnet" &&  TARGET.Machine =!= "compute-0-7.privnet" &&  TARGET.Machine =!= "hepcms-namenode2.umd.edu"

for x, z,X,Y, l, w, f, s, a,td, tr, b, y, m, d,e,r,R,S in [
    (x, z,X,Y, l, w, f, s, a,td,tr, b, y, m, d,e,r,R, S) for x in args.beamx
    for z in args.beamz 
    for X in args.TileX
    for Y in args.TileY
    for l in args.tilewidth 
    for w in args.beamidth 
    for f in args.fiberZ
    for s in args.fiberZshift
    for a in args.absmult
    for td in args.tiledecay
    for tr in args.tilerise
    for b in args.Y11abs
    for y in args.LY
    for m in args.wrapreflect 
    for d in args.Y11decayTime
    for r in args.fiberR
    for R in args.holeR
    for S in args.spare
    for e in args.sipmeff
]:

  def make_str(prefix):
    args_string = '_'.join([
        'bmx{0:.1f}'.format(x), 'bmz{0:.1f}'.format(z), 'bmw{0:.1f}'.format(w), 'bmTp{0:.0f}'.format(args.useProton),
        'TlX{0:.1f}'.format(X),'TlY{0:.1f}'.format(Y),'TlL{0:.1f}'.format(l),'WrapTp{0:.0f}'.format(args.handwrap),
        'FbS{0:.1f}'.format(s),'FbL{0:.1f}'.format(f),'FbR{0:.2f}'.format(r),'HolR{0:.2f}'.format(R),'FbTp{0:.0f}'.format(args.cladlayer),
        'TlAbs{0:.1f}'.format(a), 'FbAbs{0:.1f}'.format(b),
        'TlTR{0:.0f}'.format(tr*100),'TlTD{0:.0f}'.format(td*100),'FbTD{0:.0f}'.format(d*100),
        'WrapRef{0:.0f}'.format(m*100), 'PDE{0:.0f}'.format(e*100),
        'TlLY{0:.0f}'.format(y),'Spare{0:.3f}'.format(S),
    ])
    return prefix + args.prefix + '_' + args_string.replace('.', 'p')

  usehn2=False
  if usehn2:
      save_filename = os.path.abspath('/data3/yihuilai/' + '/' + 
                                      make_str('extruded_') + '.root')
  else:
      #save_filename = os.path.abspath(DATA_DIR + '/root/' + '/' +
      save_filename = os.path.abspath(Store_DIR + '/root/' + '/' +
                                      make_str('extruded_') + '.root')

  condor_args = ' '.join([
      '-x {}'.format(x), '-z {}'.format(z), '-X {}'.format(X),'-Y {}'.format(Y),
      '-l {}'.format(l), '-w {}'.format(w), '-f {}'.format(f), '-s {}'.format(s),
      '-a {}'.format(a), '--tiledecay {}'.format(td),'--tilerise {}'.format(tr),
      '-b {}'.format(b), '-y {}'.format(y), '-m {}'.format(m),  '-d {}'.format(d), '-e {}'.format(e),
      '-P {}'.format(args.useProton), '-H {}'.format(args.handwrap), '-c {}'.format(args.cladlayer),
      '-r {}'.format(r),'-R {}'.format(R),'-S {}'.format(S),
      '-N {}'.format(args.NEvents), '-o {}'.format(os.path.abspath(save_filename)),
  ])

  log_filename = os.path.abspath(Store_DIR + '/log/' + '/' +
                                 make_str('log_tilesim'))
  jdl_filename = os.path.abspath(DATA_DIR + '/condor/' + '/' +
                                 make_str('extruded_') + '.jdl')
  jdl_content = CONDOR_JDL_TEMPLATE.format(BASE_DIR, log_filename, condor_args)

  ## Writing jdl files
  os.makedirs(os.path.dirname(jdl_filename), exist_ok=True)
  with open(jdl_filename, 'w') as file:
    file.write(jdl_content)

  os.makedirs(DATA_DIR + '/log/', exist_ok=True)
  os.makedirs(Store_DIR + '/root/', exist_ok=True)
  os.makedirs(Store_DIR + '/log/', exist_ok=True)

  ## Making directory for output files
  if not usehn2:
      os.makedirs(os.path.dirname(save_filename), exist_ok=True)

  print(jdl_filename)
