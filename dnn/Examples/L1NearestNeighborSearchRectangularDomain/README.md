# Example of $L_1$ Nearest Neighbor Search in Rectangular Domain (L1NNS)
Sample execute file is available for $L_1$ Nearest Neighbor Search in Rectangular Domain (L1NNS).  
This field is for users who want to compile.   
This program is developed using Microsoft Visual Studio 2019.  
To compile and execute the program, you may need freeglut 3.0.0 MSVC Package in https://www.transmissionzero.co.uk/software/freeglut-devel/.   
Unzip freeglut-MSVC-3.0.0-2.mp.zip file. You may see three folders, **bin**, **include**, and **lib**.  
and put **include** and **lib** folder into the directory in which source codes are contained.  
Also, put freeglut.dll in *bin* folder into the above directory.   
Then do the following in Visual Studio.  
```
Project -> Properties -> C/C++ -> General -> Additional Include Directories -> add "./include;"  
Project -> Properties -> Linker -> General -> Additional Library Directories -> add "./lib;"  
Project -> Properties -> Linker -> Input -> Additional Dependencies -> add "freeglut.lib;"  
```

## Operations
It supports four operations using mouse and keyboard.  
|Mouse||
|:---|:---|
|[left mouse click]| inserting points, querying, or deleting objects |
|[drag and drop]| inserting rectangles |

|Keyboard||
|:---|:---|
|[esc]| Exit the program |
|[d]| Activate Deletion |
|[i]| Activate Insertion |
|[q]| Activate Query |
|[1]| From now on the query is the 1-nearest neighbor |
|[2]| From now on the query is the 2-nearest neighbor |
|[3]| From now on the query is the 3-nearest neighbor |
|[0]| From now on the query is the farthest neighbor |

### Query
	- Activate Query operation by 'q' key  
	- Set up query point  
	- If there is no data in rectangular domain, it does not work  
	- Highlight the query and its nearest data  
	- Also support 2-NNS, 3-NNS, and Farthest-NS  

### Insertion
	- Also initial state  
	- Activate Insertion operation by 'i' key  
	- Insertion of point and rectangle  
	- Overlap is impossible  

### Deletion
	- Activate Deletion operation by 'd' key  
	- Deletion of point and rectangle  
	- The nearest object from pointed location is deleted  
	- Pointed location is too far, there is no deletion  

### Exit
	- Activate Exit operation by esc key  
	- Exit the program