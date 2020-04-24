# bmrsApiUtils
M. Deakin, April 2020

bmrs_run.py is a script to pull data.
bmrs_utils.py is a small set of utility functions.
bmrs_data_.py is a (template) file for specifying api key & save directory

The main output of bmrs_run.py is an object bmd with the following properties:
- tDict: has start, stop, and timestep parameters
- dataOut: is an Nt x Nh numpy array of floats, with nans in missing values
- headings: a list of the names of headings of each column of dataOut.

It also contains:
- dtmsFull: is the full time series, for plotting (suggested not to save this)
- dtms: not used (the datetimes from the downloaded data)

It is intended that bmrs_data_.py is copied and saved as bmrs_data.py (this 
is in gitignore), and then your api key and the output directory for saving
data specified in here.

To run a query not in the list of queries already implemented, then refer
to the API guide. The XML file may also likely need to be inspected too to
populate all of the dicts of api_d, for calling and post processing the data
into a usable format.