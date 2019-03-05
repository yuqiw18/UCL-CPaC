There are three main files to run

<<SETUP>>

1.main_standalone
* This file includes all basic section tasks
* Place the 3.2gb flows.mat file under the /flow folder
* Type 'clear' in command window to clear the cache before running this file if you ran any task related to "myData_" image set before. The cache only works for "gjbLookAtTarget_" image set which is shared between file 1 and 2

2.main_standalone_advanced
* This file includes optional task [1]Better Path and [3]Slow Motion
* Place the 3.2gb flows.mat file under the /flow folder
* Type 'clear' in command window to clear the cache before running this file if you ran any task related to "myData_" image set before. The cache only works for "gjbLookAtTarget_" image set which is shared between file 1 and 2

3.main_standalone_advanced_own_data
* This file includes optional taks [1]Better Path and [2]Own Data
* Run GenerateFlowData.m and place the generated myflows.mat file under /flow folder
* Type 'clear' in command window to clear the cache before running this file if you ran any task related to "gjbLookAtTarget_" image set before. The cache only works for "myData_" image set which is not shared with any other file

<<INTERACTION>>
* The program will ask user to type a number to start from the specific frame
* Left click on a series of locations in the image to mark the path (at least 5 points)
* Press Enter to end drawing
* Output images are named based on the file you are running:
[1]output_basic_####.png
[2]output_adv_####.png
[3]output_adv_od_####.png
