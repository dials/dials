distl/dials apache server
==============================

Introduction
------------

Follow these instructions to install an apache webserver configured to run the distl and dials spotfinders. For more information on the distl spotfinder server, see http://cci.lbl.gov/labelit/html/client_server.html.  The version described here allows the beamline operator to use either distl or the dials spotfinder to assess the quality of single images, typically during a raster scan of a crystal cluster. Use of the apache server allows automatic load balancing and scaling of many parallel requests, typical for a high throughput environment.

These instructions were valid on 06/07/16 on CentOS 7 using apache version 2.4.20, mod_python version 3.5.0 and the apache dependencies 'apr', version 1.5.2 and 'apr-util', version 1.5.4.

Build instructions
------------------

These instructions include downloading apache directly from their webpage.  The apache maintainers ask that users use a mirror.  Your closest mirror can be found by visiting this web page: http://www.apache.org/dyn/closer.cgi/httpd

Detailed instructions:

* Create a new empty folder, call it $project_root
* cd $project_root
* mkdir -p apache/project_src
* cd $project_root/apache
* wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/spotfinder/servers/apache_install_dials.csh
* cd $project_root/apache/project_src
* wget https://www.apache.org/dist/httpd/httpd-2.4.20.tar.gz
* wget http://dist.modpython.org/dist/mod_python-3.5.0.tgz
* wget http://www.motorlogy.com/apache//apr/apr-1.5.2.tar.gz
* wget http://www.motorlogy.com/apache//apr/apr-util-1.5.4.tar.gz
* cd $project_root

The following command will download, configure and compile cctbx, distl and dials, then configure, build and install apache, its dependencies, and mod_python.  This takes about an hour.

* ./apache/apache_install_dials.csh

After the build is complete, instructions for starting the server and testing it will be displayed.  You may also want to modify the processor defaults stored in the mpm_prefork_module key in apache/httpd/conf/extra/httpd-mpm.conf to better suit your server capabilities.

Testing the server
------------------

These tests assume the files are accessible by the machine hosting the apache server processes.

Test the distl server with a cbf image:
curl "http://localhost:8125/spotfinder/distl.signal_strength?distl.image=/net/cci/dials/from_sunbird/sauter/rawdata/pilatus/ssrl_P6/all/I3_1_0001.cbf&distl.res.outer=3.3&distl.bins.verbose=True"

Test the dials server with a cbf image:
curl "http://localhost:8125/spotfinder/dials.find_spots?file_name=/net/cci/dials/from_sunbird/sauter/rawdata/pilatus/ssrl_P6/all/I3_1_0001.cbf&stats=True&spotfinder.filter.d_min=3.3"

Test the dials server with an Eiger image:
curl "http://localhost:8125/spotfinder/dials.find_spots?file_name=/net/cci/dials/dectris/eiger16MNov2015/2015_11_10/insu6_1_master.h5&frame_number=44&stats=True"

Distl parameters
----------------

distl.image: path to file on disk.

Further parameters matching the name/value pairs used by distl are described here: http://cci.lbl.gov/labelit/html/spotfinder.html

Dials parameters
----------------

file_name: path to file on disk

frame_number: for multi-file data, such as Eiger data from Dectris detectors, this indicates which frame number to analyze.

stats: True or False. If True, reports detailed spotfinding statistics such as resolution estimates.

Further parameters matching the name/values pairs used by the dials spotfinder are described under :ref:`dials.find_spots <dials-find-spots>`.

