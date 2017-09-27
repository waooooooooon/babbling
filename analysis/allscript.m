cd /Users/takimoto/babbling/   % return to babbling

%%%%% nofeedback


cd nofeedback/create_firingdate/
script
clear all
 
cd ../pca_5000/
run('firing_script.m')
clear all
 
cd ../pca_outputneurons/
run('firing_script.m')
clear all
 
 
%%%%%%


cd /Users/takimoto/babbling/   % return to babbling


%%%%% onfeedback
cd onfeedback/create_firingdate/
run('script.m')
clear all
 
cd ../pca_5000/
run('firing_script.m')
clear all
 
cd ../pca_outputneurons/
run('firing_script.m')
clear all
 
%%%%%%

cd /Users/takimoto/babbling/   % return to babbling


%%%%%% analysis/dimention
cd analysis/dimention/
run('dimention.m')
clear all
 
%%%%%%


cd /Users/takimoto/babbling/   % return to babbling