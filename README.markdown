<code>
Convert directory of Spec Files to CSV

optional arguments:
  -h, --help            show this help message and exit
  -i IN_DIR, --in-dir IN_DIR
                        the input directory containing the spec files.
  -o OUT_FILE, --out-file OUT_FILE
                        a csv file to contain the merged specs.
  --header              setting this flag will skip headers.
  --min-nm MIN_NM       lowest nm to include. Default is 400 nm.
  --max-nm MAX_NM       highest nm to include. Default is 700 nm.
  --intrp INTRP         interpolate nm increments. Default is 1 nm
  -s, --smooth          add smoothing function. Default is a 100 nm hanning
                        window
  --window-type {flat,hanning,hamming,bartlett,blackman}
  --window-length WINDOW_LENGTH
                        window size for smoothing. Longer is more aggressive.
                        Default is 100.
  -p, --plot            produce interactive plot with matplotlib.
  -v, -verbose          write verbose output (non functional)
  --version             prints version.
</code>


![screenshot](https://github.com/ngcrawford/Coloration/raw/master/web/screenshot_01.png)
