Basic Instructions:
-----------------

        usage: Spec.py [-h] [-i IN_DIR] [-o OUT_FILE] [--header] [--min-nm MIN_NM]
                       [--max-nm MAX_NM] [--intrp INTRP] [-s]
                       [--window-type {flat,hanning,hamming,bartlett,blackman}]
                       [--window-length WINDOW_LENGTH] [-p] [-v] [--version]

        Convert directory of Spec Files to CSV, interpolate nanometers, and smooth and
        plot data.

        optional arguments:
          -h, --help            show this help message and exit
          -i IN_DIR, --input-dir IN_DIR
                                The input directory containing the spec files.
          -o OUT_FILE, --output-file OUT_FILE
                                A csv file to contain the merged specs suitable for
                                opening in excel.
          --header              Setting this flag will skip headers.
          --min-nm MIN_NM       Lowest nm to include. Default is 400 nm.
          --max-nm MAX_NM       Highest nm to include. Default is 700 nm.
          --intrp INTRP         Interpolate nm increments. Default is 1 nm.
          -s, --smooth          Add smoothing function. Default is a 100 nm hanning
                                window.
          --window-type {flat,hanning,hamming,bartlett,blackman}
          --window-length WINDOW_LENGTH
                                Window size for smoothing. Longer is more aggressive.
                                Default is 100.
          -p, --plot            Produce interactive plots with matplotlib.
          -v, -verbose          Write verbose output (non functional).
          --version             Print version.

Screen Shot:
-----------

![screenshot](https://github.com/ngcrawford/Coloration/raw/master/web/screenshot_01.jpg)
