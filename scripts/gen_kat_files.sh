#!/bin/bash

echo 'Running Script for Level 1...' 
./build/apps/PQCgenKAT_sign_lvl1
mv PQCsignKAT_353_SQIsign_lvl1.req ./KAT/PQCsignKAT_353_SQIsign_lvl1.req
mv PQCsignKAT_353_SQIsign_lvl1.rsp ./KAT/PQCsignKAT_353_SQIsign_lvl1.rsp

echo 'Running Script for Level 3...' 
./build/apps/PQCgenKAT_sign_lvl3
mv PQCsignKAT_529_SQIsign_lvl3.req ./KAT/PQCsignKAT_529_SQIsign_lvl3.req
mv PQCsignKAT_529_SQIsign_lvl3.rsp ./KAT/PQCsignKAT_529_SQIsign_lvl3.rsp

echo 'Running Script for Level 5...' 
./build/apps/PQCgenKAT_sign_lvl5
mv PQCsignKAT_701_SQIsign_lvl5.req ./KAT/PQCsignKAT_701_SQIsign_lvl5.req
mv PQCsignKAT_701_SQIsign_lvl5.rsp ./KAT/PQCsignKAT_701_SQIsign_lvl5.rsp
