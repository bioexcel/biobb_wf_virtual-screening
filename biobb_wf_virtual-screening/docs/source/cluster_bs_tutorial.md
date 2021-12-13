# Protein-ligand Docking tutorial using BioExcel Building Blocks (biobb)
### -- *PDB Cluster90 Binding Site Version* --

***
This tutorial aims to illustrate the process of **protein-ligand docking**, step by step, using the **BioExcel Building Blocks library (biobb)**. The particular example used is the **Mitogen-activated protein kinase 14** (p38-α) protein (PDB code [3HEC](https://www.rcsb.org/structure/3HEC)), a well-known **Protein Kinase enzyme**, 
 in complex with the FDA-approved **Imatinib**, (PDB Ligand code [STI](https://www.rcsb.org/ligand/STI), DrugBank Ligand Code [DB00619](https://go.drugbank.com/drugs/DB00619)), a small molecule **kinase inhibitor** used to treat certain types of **cancer**. 
 
The tutorial will guide you through the process of identifying the **active site cavity** (pocket) without previous knowledge, and the final prediction of the **protein-ligand complex**. 

Please note that **docking algorithms**, and in particular, **AutoDock Vina** program used in this tutorial, are **non-deterministic**. That means that results obtained when running the workflow **could be diferent** from the ones we obtained during the writing of this tutorial (see [AutoDock Vina manual](http://vina.scripps.edu/manual.html)). We invite you to try the docking process several times to verify this behaviour. 
***

<div style="background:#b5e0dd; padding: 15px;"><strong>Important:</strong> it is recommended to execute this tutorial step by step (not as a single workflow execution, <strong><em>Run All</em></strong> mode), as it has interactive selections.</div>

## Settings

### Biobb modules used

 - [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
 - [biobb_structure_utils](https://github.com/bioexcel/biobb_structure_utils): Tools to modify or extract information from a PDB structure file.
 - [biobb_chemistry](https://github.com/bioexcel/biobb_chemistry): Tools to perform chemoinformatics processes.
 - [biobb_vs](https://github.com/bioexcel/biobb_vs): Tools to perform virtual screening studies.
 
### Auxiliar libraries used

 - [nb_conda_kernels](https://github.com/Anaconda-Platform/nb_conda_kernels): Enables a Jupyter Notebook or JupyterLab application in one conda environment to access kernels for Python, R, and other languages found in other environments.
 - [nglview](http://nglviewer.org/#nglview): Jupyter/IPython widget to interactively view molecular structures and trajectories in notebooks.
 - [ipywidgets](https://github.com/jupyter-widgets/ipywidgets): Interactive HTML widgets for Jupyter notebooks and the IPython kernel.

### Conda Installation and Launch

```console
git clone https://github.com/bioexcel/biobb_wf_virtual-screening.git
cd biobb_wf_virtual-screening
conda env create -f conda_env/environment.yml
conda activate biobb_VS_tutorial
jupyter-nbextension enable --py --user widgetsnbextension
jupyter-nbextension enable --py --user nglview
jupyter-notebook biobb_wf_virtual_screening/notebooks/clusterBindingSite/wf_vs_clusterBindingSite.ipynb
  ``` 

***
## Pipeline steps
 1. [Input Parameters](#input)
 2. [Fetching PDB Structure](#fetch)
 3. [Extract Protein Structure](#extractProtein)
 4. [Computing Protein Cavities (fpocket)](#fpocket)
 5. [Filtering Protein Cavities (fpocket output)](#fpocketFilter)
 6. [Extract Pocket Cavity ](#fpocketSelect)
 7. [Generating Cavity Box ](#cavityBox)
 8. [Downloading Small Molecule](#downloadSmallMolecule)
 9. [Converting Small Molecule](#sdf2pdb)
 10. [Preparing Small Molecule (ligand) for Docking](#ligand_pdb2pdbqt)
 11. [Preparing Target Protein for Docking](#protein_pdb2pdbqt)
 12. [Running the Docking](#docking)
 13. [Extract a Docking Pose](#extractPose)
 14. [Converting Ligand Pose to PDB format](#pdbqt2pdb)
 15. [Superposing Ligand Pose to the Target Protein Structure](#catPdb)
 16. [Comparing final result with experimental structure](#viewFinal)
 17. [Questions & Comments](#questions)
 
***
<img src="https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png" alt="Bioexcel2 logo"
	title="Bioexcel2 logo" width="400" />
***


<a id="input"></a>
## Input parameters
**Input parameters** needed:

 - **pdb_code**: PDB code of the experimental complex structure (if exists).<br>
In this particular example, the **p38α** structure in complex with the **Imatinib drug** was experimentally solved and deposited in the **PDB database** under the **3HEC** PDB code. The protein structure from this PDB file will be used as a **target protein** for the **docking process**, after stripping the **small molecule**. An **APO structure**, or any other structure from the **p38α** [cluster 100](https://www.rcsb.org/search?request=%7B%22query%22%3A%7B%22type%22%3A%22terminal%22%2C%22service%22%3A%22sequence%22%2C%22parameters%22%3A%7B%22target%22%3A%22pdb_protein_sequence%22%2C%22value%22%3A%22RPTFYRQELNKTIWEVPERYQNLSPVGSGAYGSVCAAFDTKTGLRVAVKKLSRPFQSIIHAKRTYRELRLLKHMKHENVIGLLDVFTPARSLEEFNDVYLVTHLMGADLNNIVKCQKLTDDHVQFLIYQILRGLKYIHSADIIHRDLKPSNLAVNEDCELKILDFGLARHTDDEMTGYVATRWYRAPEIMLNWMHYNQTVDIWSVGCIMAELLTGRTLFPGTDHIDQLKLILRLVGTPGAELLKKISSESARNYIQSLTQMPKMNFANVFIGANPLAVDLLEKMLVLDSDKRITAAQALAHAYFAQYHDPDDEPVADPYDQSFESRDLLIDEWKSLTYDEVISFVPPP%22%2C%22identity_cutoff%22%3A1%2C%22evalue_cutoff%22%3A0.1%7D%2C%22node_id%22%3A0%7D%2C%22return_type%22%3A%22polymer_entity%22%2C%22request_options%22%3A%7B%22pager%22%3A%7B%22start%22%3A0%2C%22rows%22%3A25%7D%2C%22scoring_strategy%22%3A%22combined%22%2C%22sort%22%3A%5B%7B%22sort_by%22%3A%22score%22%2C%22direction%22%3A%22desc%22%7D%5D%7D%2C%22request_info%22%3A%7B%22src%22%3A%22ui%22%2C%22query_id%22%3A%22bea5861f8b38a9e25a3e626b39d6bcbf%22%7D%7D) (sharing a 100% of sequence similarity with the **p38α** structure) could also be used as a **target protein**. This structure of the **protein-ligand complex** will be also used in the last step of the tutorial to check **how close** the resulting **docking pose** is from the known **experimental structure**. 
 -----
 - **ligandCode**: Ligand PDB code (3-letter code) for the small molecule (e.g. STI).<br>
In this particular example, the small molecule chosen for the tutorial is the FDA-approved drug **Imatinib** (PDB Code STI), a type of cancer growth blocker, used in [diferent types of leukemia](https://go.drugbank.com/drugs/DB00619).


```python
import nglview
import ipywidgets

pdb_code = "3HEC"         # P38 + Imatinib

ligand_code = "STI"       # Imatinib
```

<a id="fetch"></a>
***
## Fetching PDB structure
Downloading **PDB structure** with the **protein molecule** from the PDBe database.<br>
Alternatively, a **PDB file** can be used as starting structure. <br>
***
**Building Blocks** used:
 - [Pdb](https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb) from **biobb_io.api.pdb**
***


```python
from biobb_io.api.pdb import pdb

download_pdb = "download.pdb"
prop = {
  "pdb_code": pdb_code,
  "filter": ["ATOM", "HETATM"]
}

pdb(output_pdb_path=download_pdb,
    properties=prop)
```

<a id="vis3D"></a>
### Visualizing 3D structure
Visualizing the downloaded/given **PDB structure** using **NGL**.<br><br>
Note (and try to identify) the **Imatinib small molecule (STI)** and the **detergent (β-octyl glucoside) (BOG)** used in the experimental reservoir solution to obtain the crystal.


```python
view = nglview.show_structure_file(download_pdb, default=True)
view.center()
view._remote_call('setSize', target='Widget', args=['','600px'])

view
```

<img src='_static/cluster_bs/ngl1.png'></img>

<a id="extractProtein"></a>
***
## Extract Protein Structure
Extract **protein structure** from the **downloaded PDB file**. Removing **any extra molecule** (ligands, ions, water molecules). <br><br>
The **protein structure** will be used as a **target** in the **protein-ligand docking process**. 
***
**Building Blocks** used:
 - [extract_molecule](https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#module-utils.extract_molecule) from **biobb_structure_utils.utils.extract_molecule**
***


```python
from biobb_structure_utils.utils.extract_molecule import extract_molecule

pdb_protein = "pdb_protein.pdb"

extract_molecule(input_structure_path=download_pdb,
             output_molecule_path = pdb_protein)
```

<a id="vis3D"></a>
### Visualizing 3D structure
Visualizing the downloaded/given **PDB structure** using **NGL**.<br><br>
Note that the **small molecules** included in the original structure are now gone. The new structure only contains the **protein molecule**, which will be used as a **target** for the **protein-ligand docking**. 


```python
view = nglview.show_structure_file(pdb_protein, default=True)
view.center()
view._remote_call('setSize', target='Widget', args=['','600px'])

view
```

<img src='_static/cluster_bs/ngl2.png'></img>

<a id="bindingSite"></a>
***
## Computing Protein Cavities (Cluster90 Binding Site)
Computing the **protein cavities** (pockets) using information from the **PDB Cluster90**. The **PDB Cluster90** is a collection derived from the **PDB database**, containing structures having **less than 90% sequence identity** to each other. The **Cluster90 Binding Site** is using information from all the structures of the **Cluster90 collection** for a particular **input protein** to discover the possible **binding sites** from **small molecules** attached to **similar proteins**. <br>

These **cavities** will be then used in the **docking procedure** to try to find the **best region of the protein surface** where the small molecule can **bind**. <br><br>
Although in this particular example we already know the **binding site** region, as we started from a **protein-ligand complex** structure where the ligand was located in the same **binding site** as **Imatinib** is binding, this is not always the case. In the cases where we do not know these regions, the **Cluster90 binding site** utility will help us identifying the possible **binding sites** of our **target protein**.<br>
<br>

***
**Building Blocks** used:
 - [pdb_cluster_zip](https://biobb-io.readthedocs.io/en/latest/api.html#module-api.pdb_cluster_zip) from **biobb_io.api.pdb_cluster_zip**
 - [bindingsite](https://biobb-vs.readthedocs.io/en/latest/utils.html#module-utils.bindingsite) from **biobb_vs.utils.bindingsite**
***

<a id="cluster90"></a>
### Getting Cluster 90 collection 
Extracting the **Cluster90** collection from the input **PDB structure**. The collection will contain all the structures in the **PDB database** having **less than 90% sequence identity** to the input structure (our target protein).<br><br>
*Please note that depending on the size of the cluster, the execution can take a while (minutes).*


```python
from biobb_io.api.pdb_cluster_zip import pdb_cluster_zip

pdb_cluster = "pdb_cluster.zip"
prop = {
    "pdb_code": pdb_code,
    "filter": ["ATOM", "HETATM"],
    "cluster": 90
}

pdb_cluster_zip(output_pdb_zip_path = pdb_cluster,
            properties=prop)
```

<a id="bindingsite"></a>
### Extracting the Cluster90 Binding Site(s)  
Extracting the protein binding site(s) from the **Cluster90** collection. There are two ways in which the **Cluster90 Binding Site** can be used:

- **With information about the ligand**: In this particular example we already know the **ligand** we want to dock, and also that the **PDB database** contain an **experimental structure** where the **ligand** was solved in **complex** with our protein of interest (**p38-α Tyrosine kinase**). In that case, the **ligand id** can be used to guide the **binding site** tool. <br><br>
- **Without information about the ligand**: If the docking study starts from a **protein receptor** with no information about the **binding site** or any known **protein-ligand complex**, the **binding site** tool can not be guided. Still, the tool should be able to extract **binding site(s)** information from **similar proteins** solved with **small molecules** attached (if any).<br><br>

This example is using the **ligand id** to guide the **binding site tool**, but we invite you to try and explore the differences in the output removing this **input information** (ligand property in the next building block).

The **Cluster90 binding site** tool is internally running **sequence alignments** to structurally **superpose** the structures contained in the **Cluster90 collection** and extract the **residue numbering** corresponding to the **binding site** residues. Thus, **sequence alignment parameters** can be changed using building block input properties. In this example, we have chosen the well-known **blosum62** substitution matrix, with a **penalty** for opening a gap of -10.0, and a **penalty** to extend a gap of -0.5. The total number of **superimposed ligands** to be extracted from the cluster is limited at 15, and the **cut-off distance** around the ligand atoms to consider a residue to be part of a **binding site** is fixed at 5 Ångstroms.


```python
from biobb_vs.utils.bindingsite import bindingsite

output_bindingsite = "bindingsite.pdb"
prop = {
    "ligand": ligand_code,
    "matrix_name": "blosum62",
    "gap_open": -10.0,
    "gap_extend": -0.5,
    "max_num_ligands": 15,
    "radius": 5
}

bindingsite(input_pdb_path = pdb_protein,
            input_clusters_zip = pdb_cluster,
            output_pdb_path = output_bindingsite,
            properties=prop)
```

<a id="viewPockets"></a>
### Visualizing selected pockets (cavities)
Visualizing the selected **pockets** (cavities) from the generated list using **NGL viewer**.<br>

**Protein residues** forming the **cavity** are shown in **licorice** representation. **Pockets** are represented in a **greyish surface**. The **original ligand** (if exists) is shown in **green-colored ball and stick** representation. 


```python
view = nglview.show_structure_file(download_pdb, default=False)

# ligand
view[0].add_representation(repr_type='ball+stick', 
                          selection='STI',
                          aspect_ratio=4,
                          color='green')

view[0].add_representation(repr_type='cartoon', 
                        selection='not het',
                          opacity=.2,
                          color='#cccccc')

view.add_component(output_bindingsite, default=False)
view[1].add_representation(repr_type='surface', 
                           selection='*', 
                           opacity = .3,
                           radius='1.5',
                           lowResolution= True,
                           # 0: low resolution 
                           smooth=1,
                           useWorker= True,
                           wrap= True)
view[1].add_representation(repr_type='licorice', 
                        selection='*')

view[0].center()
view._remote_call('setSize', target='Widget', args=['','600px'])

view
```

<img src='_static/cluster_bs/ngl3.png'></img>

<a id="cavityBox"></a>
***
## Generating Cavity Box 
Generating a **box** surrounding the selected **protein cavity** (pocket), to be used in the **docking procedure**. The **box** is defining the region on the **surface** of the **protein target** where the **docking program** should explore a possible **ligand dock**.<br>
An offset of **12 Angstroms** is used to generate a **big enough box** to fit the **small molecule** and its possible rotations.<br>

***
**Building Blocks** used:
 - [box](https://biobb-vs.readthedocs.io/en/latest/utils.html#module-utils.box) from **biobb_vs.utils.box**
***


```python
from biobb_vs.utils.box import box

output_box = "box.pdb"
prop = {
    "offset": 12,
    "box_coordinates": True
}

box(input_pdb_path = output_bindingsite,
            output_pdb_path = output_box,
            properties=prop)
```

<a id="vis3D"></a>
### Visualizing binding site box in 3D structure
Visualizing the **protein structure**, the **selected cavity**, and the **generated box**, all together using **NGL** viewer. Using the **original structure** with the **small ligand** inside (Imatinib, [STI](https://www.rcsb.org/ligand/STI)), to check that the **selected cavity** is placed in the **same region** as the **original ligand**. 


```python
#view = nglview.show_structure_file(box, default=False)
view = nglview.NGLWidget()
#s = view.add_component(pdb_single_chain)
s = view.add_component(download_pdb)
b = view.add_component(output_box)
s = view.add_component(output_bindingsite)

atomPair = [
    [ "9999:Z.ZN1", "9999:Z.ZN2" ],
    [ "9999:Z.ZN2", "9999:Z.ZN4" ],
    [ "9999:Z.ZN4", "9999:Z.ZN3" ],
    [ "9999:Z.ZN3", "9999:Z.ZN1" ],
    
    [ "9999:Z.ZN5", "9999:Z.ZN6" ],
    [ "9999:Z.ZN6", "9999:Z.ZN8" ],
    [ "9999:Z.ZN8", "9999:Z.ZN7" ],
    [ "9999:Z.ZN7", "9999:Z.ZN5" ],
    
    [ "9999:Z.ZN1", "9999:Z.ZN5" ],
    [ "9999:Z.ZN2", "9999:Z.ZN6" ],
    [ "9999:Z.ZN3", "9999:Z.ZN7" ],
    [ "9999:Z.ZN4", "9999:Z.ZN8" ]
]

#view.shape.add_cylinder( [ 0, 2, 7 ], [ 10, 0, 9 ], [ 1, 0, 0 ], 0.1 )

# structure
s.add_representation(repr_type='cartoon', 
                        selection='not het',
                        color='#cccccc',
                       opacity=.2)
# ligands box
b.add_representation(repr_type='ball+stick', 
                        selection='9999',
                        color='pink',
                       aspectRatio = 10)
# lines box
b.add_representation(repr_type='distance', 
                        atomPair= atomPair,
                       labelColor= 'transparent',
                       color= 'black')

# output bindingsite
s.add_representation(repr_type='surface', 
                        selection='*', 
                        color='skyblue',
                        lowResolution= True,
                        # 0: low resolution 
                        smooth=1,
                        surfaceType= 'av', 
                        contour=True,
                        opacity=0.4,
                        useWorker= True,
                        wrap= True)


view.center()
view._remote_call('setSize', target='Widget', args=['','600px'])

view
```

<img src='_static/cluster_bs/ngl4.png'></img>

<a id="downloadSmallMolecule"></a>
***
## Downloading Small Molecule 
Downloading the desired **small molecule** to be used in the **docking procedure**. <br>
In this particular example, the small molecule of interest is the FDA-approved drug **Imatinib**, with PDB code **STI**.<br>

***
**Building Blocks** used:
 - [ideal_sdf](https://biobb-io.readthedocs.io/en/latest/api.html#module-api.ideal_sdf) from **biobb_io.api.ideal_sdf**
***


```python
from biobb_io.api.ideal_sdf import ideal_sdf

sdf_ideal = "ideal.sdf"
prop = {
  "ligand_code": ligand_code
}

ideal_sdf(output_sdf_path=sdf_ideal,
    properties=prop)
```

<a id="sdf2pdb"></a>
***
## Converting Small Molecule 
Converting the desired **small molecule** to be used in the **docking procedure**, from **SDF** format to **PDB** format using the **OpenBabel chemoinformatics** tool. <br>

***
**Building Blocks** used:
 - [babel_convert](https://biobb-chemistry.readthedocs.io/en/latest/babelm.html#module-babelm.babel_convert) from **biobb_chemistry.babelm.babel_convert**
***


```python
from biobb_chemistry.babelm.babel_convert import babel_convert

ligand = "ligand.pdb"
prop = {
    "input_format": "sdf",
    "output_format": "pdb",
    "obabel_path": "obabel"
}

babel_convert(input_path = sdf_ideal,
            output_path = ligand,
            properties=prop)
```

<a id="ligand_pdb2pdbqt"></a>
***
## Preparing Small Molecule (ligand) for Docking
Preparing the **small molecule** structure for the **docking procedure**. Converting the **PDB file** to a **PDBQT file** format (AutoDock PDBQT: Protein Data Bank, with Partial Charges (Q), & Atom Types (T)), needed by **AutoDock Vina**. <br><br>
The process adds **partial charges** and **atom types** to every atom. Besides, the **ligand flexibility** is also defined in the information contained in the file. The concept of **"torsion tree"** is used to represent the **rigid and rotatable** pieces of the **ligand**. A rigid piece (**"root"**) is defined, with zero or more rotatable pieces (**"branches"**), hanging from the root, and defining the **rotatable bonds**.<br><br>
More info about **PDBQT file format** can be found in the [AutoDock FAQ pages](http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file).

***
**Building Blocks** used:
 - [babel_convert](https://biobb-chemistry.readthedocs.io/en/latest/babelm.html#module-babelm.babel_convert) from **biobb_chemistry.babelm.babel_convert**
***


```python
from biobb_chemistry.babelm.babel_convert import babel_convert

prep_ligand = "prep_ligand.pdbqt"
prop = {
    "input_format": "pdb",
    "output_format": "pdbqt",
    "obabel_path": "obabel"
}

babel_convert(input_path = ligand,
            output_path = prep_ligand,
            properties=prop)
```

<a id="viewDrug"></a>
### Visualizing small molecule (drug)
Visualizing the desired **drug** to be docked to the **target protein**, using **NGL viewer**.<br>
- **Left panel**: **PDB-formatted** file, with all hydrogen atoms.
- **Right panel**: **PDBqt-formatted** file (AutoDock Vina-compatible), with **united atom model** (only polar hydrogens are placed in the structures to correctly type heavy atoms as hydrogen bond donors).



```python
from ipywidgets import HBox

v0 = nglview.show_structure_file(ligand)
v1 = nglview.show_structure_file(prep_ligand)

v0._set_size('500px', '')
v1._set_size('500px', '')

def on_change(change):
    v1._set_camera_orientation(change['new'])
    
v0.observe(on_change, ['_camera_orientation'])

HBox([v0, v1])
```

<img src='_static/cluster_bs/ngl5.png'></img>

<a id="protein_pdb2pdbqt"></a>
***
## Preparing Target Protein for Docking
Preparing the **target protein** structure for the **docking procedure**. Converting the **PDB file** to a **PDBqt file**, needed by **AutoDock Vina**. Similarly to the previous step, the process adds **partial charges** and **atom types** to every target protein atom. In this case, however, we are not taking into account **receptor flexibility**, although **Autodock Vina** allows some limited flexibility of selected **receptor side chains** [(see the documentation)](https://autodock-vina.readthedocs.io/en/latest/docking_flexible.html).<br>

***
**Building Blocks** used:
 - [str_check_add_hydrogens](https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#utils-str-check-add-hydrogens-module) from **biobb_structure_utils.utils.str_check_add_hydrogens**
***


```python
from biobb_structure_utils.utils.str_check_add_hydrogens import str_check_add_hydrogens

prep_receptor = "prep_receptor.pdbqt"
prop = {
    "charges": True,
    "mode": "auto"
}

str_check_add_hydrogens(input_structure_path = pdb_protein,
            output_structure_path = prep_receptor,
            properties=prop)
```

<a id="docking"></a>
***
## Running the Docking
Running the **docking process** with the prepared files:
- **ligand**
- **target protein**
- **binding site box**<br>

using **AutoDock Vina**. <br>

***
**Building Blocks** used:
 - [autodock_vina_run](https://biobb-vs.readthedocs.io/en/latest/vina.html#module-vina.autodock_vina_run) from **biobb_vs.vina.autodock_vina_run**
***


```python
from biobb_vs.vina.autodock_vina_run import autodock_vina_run

output_vina_pdbqt = "output_vina.pdbqt"
output_vina_log = "output_vina.log"

autodock_vina_run(input_ligand_pdbqt_path = prep_ligand,
             input_receptor_pdbqt_path = prep_receptor,
             input_box_path = output_box,
            output_pdbqt_path = output_vina_pdbqt,
             output_log_path = output_vina_log)
```

<a id="viewDocking"></a>
### Visualizing docking output poses
Visualizing the generated **docking poses** for the **ligand**, using **NGL viewer**. <br>
- **Left panel**: **Docking poses** displayed with atoms coloured by **partial charges** and **licorice** representation.
- **Right panel**: **Docking poses** displayed with atoms coloured by **element** and **ball-and-stick** representation.


```python
models = 'all'
#models = '/0 or /1 or /4'

v0 = nglview.show_structure_file(output_vina_pdbqt, default=False)
v0.add_representation(repr_type='licorice', 
                        selection=models,
                       colorScheme= 'partialCharge')
v0.center()
v1 = nglview.show_structure_file(output_vina_pdbqt, default=False)
v1.add_representation(repr_type='ball+stick', 
                        selection=models)
v1.center()

v0._set_size('500px', '')
v1._set_size('500px', '')

def on_change(change):
    v1._set_camera_orientation(change['new'])
    
v0.observe(on_change, ['_camera_orientation'])

HBox([v0, v1])
```

<img src='_static/cluster_bs/ngl6.png'></img>

<a id="selectPose"></a>
### Select Docking Pose
Select a specific **docking pose** from the output list for **visual inspection**.
<br>
Choose a **docking pose** from the **DropDown list**.


```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET = True)
structure = parser.get_structure("protein", output_vina_pdbqt)
models = []
for i, m in enumerate(structure):
    models.append(('model' + str(i), i))
    
mdsel = ipywidgets.Dropdown(
    options=models,
    description='Sel. model:',
    disabled=False,
)
display(mdsel)
```

<a id="extractPose"></a>
***
## Extract a Docking Pose
Extract a specific **docking pose** from the **docking** outputs. <br>

***
**Building Blocks** used:
 - [extract_model_pdbqt](https://biobb-vs.readthedocs.io/en/latest/utils.html#module-utils.extract_model_pdbqt) from **biobb_vs.utils.extract_model_pdbqt**
***


```python
from biobb_vs.utils.extract_model_pdbqt import extract_model_pdbqt

output_pdbqt_model = "output_model.pdbqt"
prop = {
    "model": 1
}

extract_model_pdbqt(input_pdbqt_path = output_vina_pdbqt,
             output_pdbqt_path = output_pdbqt_model,
            properties=prop)
```

<a id="pdbqt2pdb"></a>
***
## Converting Ligand Pose to PDB format
Converting **ligand pose** to **PDB format**. <br>

***
**Building Blocks** used:
 - [babel_convert](https://biobb-chemistry.readthedocs.io/en/latest/babelm.html#module-babelm.babel_convert) from **biobb_chemistry.babelm.babel_convert**
***


```python
from biobb_chemistry.babelm.babel_convert import babel_convert

output_pdb_model = "output_model.pdb"
prop = {
    "input_format": "pdbqt",
    "output_format": "pdb",
    "obabel_path": "obabel"
}

babel_convert(input_path = output_pdbqt_model,
             output_path = output_pdb_model,
            properties=prop)
```

<a id="catPdb"></a>
***
## Superposing Ligand Pose to the Target Protein Structure
Superposing **ligand pose** to the target **protein structure**, in order to see the **protein-ligand docking conformation**. <br><br>Building a new **PDB file** with both **target and ligand** (binding pose) structures. <br>

***
**Building Blocks** used:
 - [cat_pdb](https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#module-utils.cat_pdb) from **biobb_structure_utils.utils.cat_pdb**
***


```python
from biobb_structure_utils.utils.cat_pdb import cat_pdb

output_structure = "output_structure.pdb"

cat_pdb(input_structure1 = pdb_protein,
             input_structure2 = output_pdb_model,
             output_structure_path = output_structure)
```

<a id="viewFinal"></a>
### Comparing final result with experimental structure 
Visualizing and comparing the generated **protein-ligand** complex with the original **protein-ligand conformation** (downloaded from the PDB database), using **NGL viewer**. <br>
- **Licorice, element-colored** representation: **Experimental pose**.
- **Licorice, green-colored** representation: **Docking pose**.
<br>

Note that outputs from **AutoDock Vina** don't contain all the atoms, as the program works with a **united-atom representation** (i.e. only polar hydrogens).


```python
view = nglview.NGLWidget()

# v1 = Experimental Structure
v1 = view.add_component(download_pdb)

v1.clear()
v1.add_representation(repr_type='licorice', 
                     selection='STI',
                     radius=0.5)

# v2 = Docking result
v2 = view.add_component(output_structure)
v2.clear()
v2.add_representation(repr_type='cartoon', colorScheme = 'sstruc')
v2.add_representation(repr_type='licorice', radius=0.5, color= 'green', selection='UNL')

view._remote_call('setSize', target='Widget', args=['','600px'])
view

# align reference and output
code = """
var stage = this.stage;
var clist_len = stage.compList.length;
var i = 0;
var s = [];
for(i = 0; i <= clist_len; i++){
    if(stage.compList[i] != undefined && stage.compList[i].structure != undefined) {        
       s.push(stage.compList[i])
    }
}
NGL.superpose(s[0].structure, s[1].structure, true, ".CA")
s[ 0 ].updateRepresentations({ position: true })
s[ 0 ].autoView()
"""

view._execute_js_code(code)

view
```

<img src='_static/cluster_bs/ngl7.png'></img>

***
<a id="questions"></a>

## Questions & Comments

Questions, issues, suggestions and comments are really welcome!

* GitHub issues:
    * [https://github.com/bioexcel/biobb](https://github.com/bioexcel/biobb)

* BioExcel forum:
    * [https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library](https://ask.bioexcel.eu/c/BioExcel-Building-Blocks-library)

