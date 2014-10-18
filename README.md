climada_module_barisal_demo
===========================

Barisal, Bangladesh, demo module for tropical cyclone

this module contains the additional climada module to implement a demo tropical cyclone hazard event set for Barisal in Bangladesh.

please install climada core first, see https://github.com/davidnbresch/climada

In order to grant core climada access to additional modules, create a folder ‘modules’ in the core climada folder and copy/move any additional modules into climada/modules, without 'climada_module_' in the filename. 

E.g. if the addition module is named climada_module_MODULE_NAME, we should have
.../climada the core climada, with sub-folders as
.../climada/code
.../climada/data
.../climada/docs
and then
.../climada/modules/MODULE_NAME with contents such as 
.../climada/modules/MODULE_NAME/code
.../climada/modules/MODULE_NAME/data
.../climada/modules/MODULE_NAME/docs
this way, climada sources all modules' code upon startup

see climada/docs/climada_manual.pdf to get started

copyright (c) 2014, David N. Bresch, david.bresch@gmail.com all rights reserved.


see climada/docs/climada_manual.pdf to get started and

copyright (c) 2014, David N. Bresch, david.bresch@gmail.com, all rights reserved
