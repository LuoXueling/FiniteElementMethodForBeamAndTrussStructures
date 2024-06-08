# Finite Element Method for Beam and Truss Structures

This a coursework in 2021 when I was an undergraduate student at Shanghai Jiao Tong University. I implemented finite element method for 3D beam (2 nodes) and truss (2 nodes) structures, visualization of deformation and stress, and a simple UI (Chinese) using PyQt5.

## Requirement

I used Visual Studio for C++ finite element calculation, and Python 3.7 with `matplotlib`, `PyQt5`, and `PyQt5-tools` for visualizations and UI. 

## Examplary usage

1. Run `ui.py`
2. Click `New` (or run `newStructure.py`) to generate a new structure following guidance.
    * To generate a beam structure, for example:
        * When asked the type of elements: `beam`
        * When asked the coordinates of nodes:
            ```shell
            0,0.96,0
            1.44,0.96,0
            0,0,0
            1.44,0,0
            end
            ```
        * When asked the ID of nodes of elements:
            ```shell
            1,2
            1,3
            2,4
            end
            ```
        * When asked the ID of [support node](https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usb/default.htm?startat=pt06ch23s03alm08.html) (or additional node? I don't know its English name) of elements that defines the local y direction of sections together with the first node of the element (yz is the plane of the cross-sectional area): 
            ```shell
            0,0.96,1
            0,0.96,1
            1.44,0.96,1
            ```
        * When asked the constraints (When assigning constraints and if it's a 2D problem, all displacements and rotations in z direction should be constrained. ):
            ```shell
            3
            x,y,z,tx,ty,tz
            4 
            x,y,z,tx,ty,tz
            end
            ```
        * When asked about loads:
            ```shell
            1 3000,-3000,0,0,0,-720
            2 0,-3000,0,0,0,720
            end
            ```
        * When asked about properties of sections ( E,A,G,Ay,Az,Ix,Iy,Iz,K,a1,b1,c1,a2,b2,c2):
        ```shell
        300000000000,0.00068,1,1,1,0.00000065,0.00000065,0.00000065,0,0,0,0,0,0,0
        end
        ```
        * When asked the ID of section  of each element, use `all k` to assign the k-th section to all elements:
        ```shell
        all 1
        ```
    * To generate a truss structure, for example: 
        * When asked the type of elements: `bar`
        * Nodes:
            ```shell
            0,0,0
            4,0,0
            0,3,0
            4,3,0
            0,6,0
            4,6,0
            0,9,0
            4,9,0
            end
            ```
        * ID of nodes of each element:
            ```shell
            1,3
            1,4
            2,3
            2,4
            3,4
            3,5
            3,6
            4,5
            4,6
            5,6
            5,7
            5,8
            6,7
            6,8
            7,8
            end
            ```
        * Constraints (When assigning constraints and if it's a 2D problem, all z DOFs should be constrained. ):
            ```shell
            1
            x,y,z
            2 
            x,y,z
            3
            z 
            4 
            z 
            5 
            z 
            6 
            z 
            7 
            z 
            8 
            z 
            end
            ```
        * Loads:
            ```shell
            3 20000,0,0
            5 20000,0,0
            7 10000,0,0
            end
            ```
        * E and A of sections:
            ```shell
            200000000000.0 0.0002
            end
            ```
        * Assign sections to elements:
            ```shell
            all 1
            ```

3. Click `plot` or run `plotStucture.py` and select the input file to check the structure. Blue circles and line segments indicate constraints, and Red circles and line segments indicate loads. 

* The examplary beam structure

<image src=".assets/beam.png" width="50%"></image>

* The examplary truss structure

<image src=".assets/truss.png" width="50%"></image>

The figures are interactable in the UI. 

4. Click `compute` to call `csm.exe` compiled from the `csm` folder. `csm.exe` in the repo is in Chinese but I translate it into English in the source code `csm` without recompiling (since I don't have a windows PC currently and the project uses some libs like `windows.h`). 

5. Plot again to see the deformed configuration and stress

<image src=".assets/beam_res.png" width="50%"></image>

<image src=".assets/truss_res.png" width="50%"></image>

## More examples

* `bar_test_1.txt`:

<image src=".assets/bar_test_1.png" width="50%"></image>

* `beam_test_3.txt`:

<image src=".assets/beam_test_3.png" width="50%"></image>

* `beam_test_4.txt`:

<image src=".assets/beam_test_4.png" width="50%"></image>

* `beam_test_5.txt`:

<image src=".assets/beam_test_5.png" width="50%"></image>

* `final_test_beam.txt`:

<image src=".assets/final_test_beam.png" width="50%"></image>
