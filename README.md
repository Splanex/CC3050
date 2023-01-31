Usage:
	python3 hpu.py command arg1 ... argn

Available commands:
  1. Combine files by averaging:
     avg out_file file1 ... filen	
	
  2. Perform darkframe correction:
     cdf out_file df_file ff_file file1 ... filen
     
  3. Perform hot pixel correction:
     chp out_file df_file ff_file k file1 ... filen
     
  4. Compare two files:
     diff diff_file file1 file2
     
  5. Display help:
     help
     
  6. Show hot pixels detected in file:
     hph out_file src_file k
     
  7. Show FITS file info:
     info file1 ... filen
     
  8. Show PSNR values vs base file:
     psnr base_file file1 ... filen

  9. Image alignment:
     register file1 ... filen
     
  10. Show SNR values:
      snr file1 ... filen

