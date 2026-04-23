#!/bin/bash 

#run to produce shinyapp
R -e "shiny::runApp('dashboard/app.R', host='0.0.0.0', port=8080)"
