[CRAB]

jobtype = cmssw
#scheduler = glite
scheduler = caf
### NOTE: just setting the name of the server (pi, lnl etc etc )
###       crab will submit the jobs to the server...
#server_name = bari

[CMSSW]

### The data you want to access (to be found on DBS)
## runs with 0T
#runselection=69557,69559,69564,69572,69573,69587,69594,69728,69743,70147,70170,70195,70410,70411,70412,70413,70414,70415,70416,70417,70421,70454,70664,70674,70675
##first part of 3.8T
runselection=66612,66615,66637,66644,66657,66662,66668,66676,66692,66703,66711,66714,66716,66720,66722,66739,66740,66746,66748,66752,66757,66783,66887,66893,66904,67126,67128,67139,67147,67173,67225,67539,67541,67544,67548,67557,67573,67645,67647,67810,67818,67838,67977,68000
##second part of 3.8T
#runselection=68021,68087,68094,68100,68129,68141,68264,68273,68276,68279,68286,68288,68500,68665,68926,68949,68958,69044,69046,69330,69332,69333,69335,69343,69351,69357,69365,69382,69396,69491,69750,69788,69797,69800,69850,69892,69902,69912,69997,70036,70088
#runselection=66644,66711,66733,66748,66904,67147,67489,67531,67573,66709,66720,66732,66739,67541,67544,67548,66615,66676,66668,66722,66752,66756,66783,66887,67124,67217,67529,66612,66627,67141,67145,67234,66637,66657,66692,66703,66716,66893,67126,67122,67139,67173,67527,67538,67644,66656,66662,66740,66878,67534,66621,67128,67219,67225,67547,66613,66706,66714,66746,66757,67114,67214,67539

datasetpath=/Cosmics/Commissioning08_CRAFT_ALL_V9_225_ReReco_FromTrackerPointing_v1/RECO
### The ParameterSet you want to use
pset=pixelNtuplizer_RealData_cfg38T.py

### Splitting parameters
total_number_of_events=-1
#total_number_of_events=100
events_per_job = 5000
#number_of_jobs = 20

### The output files (comma separated list)
output_file = Commissioning08_CRAFT_ALL_V9_225_ReReco_FromTrackerPointing_v1_ntupl.root

dls_phedex_url=http://cmsweb.cern.ch/phedex/datasvc/xml/prod/
#dbs_url=cms_dbs_caf_dbs_analysis_01
[USER]

### OUTPUT files Management
##  output back into UI
return_data = 0

### To use a specific name of UI directory where CRAB will create job to submit (with full path).
### the default directory will be "crab_0_data_time"
#ui_working_dir = /full/path/Name_of_Directory

### To specify the UI directory where to store the CMS executable output
### FULL path is mandatory. Default is  <ui_working_dir>/res will be used.
#outputdir= /full/path/yourOutDir

### To specify the UI directory where to store the stderr, stdout and .BrokerInfo of submitted jobs
### FULL path is mandatory. Default is <ui_working_dir>/res will be used.
#logdir= /full/path/yourLogDir

### OUTPUT files INTO A SE
copy_data = 1

### if you want to copy data in a "official CMS site"
### you have to specify the name as written in 
#storage_element = T2_IT_Bari
### the user_remote_dir will be created under the SE mountpoint
### in the case of publication this directory is not considered
#user_remote_dir = name_directory_you_want

### if you want to copy your data at CAF
#storage_element = T2_CH_CAF
### the user_remote_dir will be created under the SE mountpoint
### in the case of publication this directory is not considered
#user_remote_dir = name_directory_you_want

### if you want to copy your data to your area in castor at cern
### or in a "not official CMS site" you have to specify the complete name of SE
storage_element=srm-cms.cern.ch
### this directory is the mountpoin of SE 
storage_pool=None 
storage_path=/castor/cern.ch
### directory or tree of directory under the mounpoint 
user_remote_dir=/user/t/trommers/test


### To publish produced output in a local istance of DBS set publish_data = 1
publish_data=0
### Specify the dataset name. The full path will be <primarydataset>/<publish_data_name>/USER
publish_data_name = name_you_prefer
### Specify the URL of DBS istance where CRAB has to publish the output files
#dbs_url_for_publication = https://cmsdbsprod.cern.ch:8443/cms_dbs_caf_analysis_01_writer/servlet/DBSServlet 

#if server
#thresholdLevel = 100
#eMail = your@Email.address

[EDG]

## RB/WMS management:
rb = CERN

##  Black and White Lists management:
## By Storage
#se_black_list = T0,T1
#se_white_list =

## By ComputingElement
#ce_black_list = ucsd
#ce_white_list = cern

[CONDORG]

# Set this to condor to override the batchsystem defined in gridcat.
#batchsystem = condor

# Specify addition condor_g requirments
# use this requirment to run on a cms dedicated hardare
# globus_rsl = (condor_submit=(requirements 'ClusterName == \"CMS\" && (Arch == \"INTEL\" || Arch == \"X86_64\")'))
# use this requirement to run on the new hardware
#globus_rsl = (condor_submit=(requirements 'regexp(\"cms-*\",Machine)'))

