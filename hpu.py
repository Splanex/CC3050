import warnings
warnings.filterwarnings("ignore")

import hpf
import sys

class C:
     def avg(args):
          if args == 'help':
               return "avg out_file file1 ... filen\n\t\tcombine files by averaging"
          if len(args) < 2:
               print('No files given!')
               return 
          hpf.save(args[0], hpf.average(args[1:]))
     def cdf(args):
          if args == 'help':
               return "cdf out_file df_file ff_file file1 ... filen\n\t\tperform darkframe correction"
          if len(args) < 4:
               print('No files given!')
               return 
          img = hpf.darkframe_correction_masterlight(  darkframe=args[1], 
                                                       flatfield=args[2],
                                                       image_list=args[3:])
          hpf.save(args[0], img)
     def chp(args):
          if args == 'help':
               return "chp out_file df_file ff_file k file1 ... filen\n\t\tperform hot pixel correction"
          if len(args) < 5:
               print('No files given!')
               return 
          img = hpf.hotpixel_correction_masterlight(   darkframe=args[1], 
                                                       flatfield=args[2],
                                                       image_list=args[4:],
                                                       k=float(args[3]))
          hpf.save(args[0], img)
     def diff(args):
          if args == 'help':
               return "diff diff_file file1 file2\n\t\tcompare two files"
          if len(args) < 3:
               print("Invalid use!", C.diff('help'))
               return 
          hpf.save(args[0], hpf.diff(args[1], args[2]))
     def help(args):
          if args == 'help':
               return 'help\n\t\tdisplay help'
          print(
'''Usage:
\thpu.py command arg1 ... argn\n
Available commands:''')
          for (c,e) in C.__dict__.items():
               if '_' not in c:
                    print('\t', e('help'))
     def hph(args):
          if args == 'help':
               return 'hph out_file src_file k\n\t\tshow hot pixels detected in file'
          if len(args) != 3:
               print('Invalid arguments')
               return 
          hp = hpf.highlight_hotpixels(args[1],k=float(args[2]))
          hpf.save(args[0], hp)
     def info(args):
          if args == 'help':
               return "info file1 ... filen\n\t\tshow FITS file info"
          if len(args) == 0:
               print('No files given!')
               return 
          for file in args:
               hpf.fits_info(file)
     def psnr(args):
          if args == 'help':
               return "psnr base_file file1 ... filen\n\t\tshow PSNR values vs base file"
          if len(args) < 2:
               print('No files given!')
               return  
          base = args[0] 
          base_img = hpf.load(base)
          for file in args[1:]:
               print("PSNR %s,%s = %.3f " % (base, file, hpf.psnr(base,file)))
     def register(args):
          if args == 'help':
               return "register file1 ... filen\n\t\timage alignment"
          if len(args) == 0:
               print('No files given!')
               return  
          hpf.register(args)
     def snr(args):
          if args == 'help':
               return "snr file1 ... filen\n\t\tshow SNR values"
          if len(args) == 0:
               print('No files given!')
               return 
          for file in args:
               print("SNR %s = %.6f " % (file, hpf.snr(file)))


if len(sys.argv) == 1:
     C.help(None)
     sys.exit(0)

try:
     command = C.__dict__[sys.argv[1]]
     command(sys.argv[2:])
except KeyError:
     print('Invalid command:', sys.argv[1], '\nUse "help" for available commands')


