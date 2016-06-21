ega_client_galaxy_wrapper
=========================

Bio-molecular high throughput data is privacy sensitive and can not easily made accessible to the entire outside world. To manage access to long term-archival of such data the European Genome-phenome Archive (EGA) project was initiated to facilitate data access and management to funded projects after completion to enable continued access to these data. Strict protocols govern how information is managed, stored, transferred and distributed and each data provider is responsible for ensuring a Data Access Committee is in place to grant access to the data. Obtaining such data is achieved with the EGA download streamer.

License
-------
This embeds access to the EGA download streamer into Galaxy (16.07+). The EGA download streamer (https://www.ebi.ac.uk/ega/about/your_EGA_account/download_streaming_client) is free software released under the GNU GPL 3 license. By installating this tool you automatically agree with the terms of this license.


Installation
------------
The galaxy wrapper is available in the public tool sheds under the name *ega_download_streamer* and can be viewed at the following url: https://toolshed.g2.bx.psu.edu/view/yhoogstrate/ega_download_streamer

The galaxy web administrator can install this wrapper in the admin menu, under the tabs **Tools and Tool Shed**: *Search Tool Shed*. In here, the admin needs to go to the **Galaxy main tool shed** where the search term *ega_download_streamer* should show the tool installer maintained by *yhoogstrate*. By pressing on it, a new page will open. This page will indicate that it will automatically install the tool dependencies and shows this README, including the LICENSE, and a button that will install the entire package into the galaxy server.

Because currently galaxy does not (yet?) support a password input type, the tool has to be configured by setting up 1 generic account. This is done by adding lines, setting up the proper environment variables, to a destination (computer node) in galaxy's `./config/job_conf.xml`: 

``<destination id="local_dest" runner="local">``

``<env id="user">email</env>``

``<env id="pass">Secret12456...blabla</env>``

``</destination>``

The `runner="local"` should correspond to the runner that runs the ega_download_streamer (only for complex setups, local is default).


Maintainance
------------

The galaxy repository is maintained by Youri Hoogstrate and was written for an ELIXIR/EGA project.

Bugs and other contributions to the galaxy wrapper are more than welcome at:
https://github.com/ErasmusMC-Bioinformatics/ega_client_galaxy_wrapper
