{::nomarkdown}

<html>

<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>README.md</title>
  <link rel="stylesheet" href="https://stackedit.io/style.css" />
</head>

<body class="stackedit">
  <div class="stackedit__html"><h1 id="how-to-retrofitting-metal-organic-frameworks">How-To: Retrofitting Metal-Organic Frameworks</h1>
<h5 id="authors-christian-schneider-david-bodesheim-julian-keupp-rochus-schmid-gregor-kieslich">Authors: Christian Schneider, David Bodesheim, Julian Keupp, Rochus Schmid, Gregor Kieslich*</h5>
<p>* <a href="mailto:gregor.kieslich@tum.de">gregor.kieslich@tum.de</a></p>
<p>This is the README file to the python-based RetroFit script, which can evaluate the fit of a molecule (Cross Linker, CL) into a given MOF featuring open metal sites (OMS). This heuristic tool is aimed to guide experimental chemists in developing new materials by installing additional linkages in a given MOF. The use of RetroFit requires only little knowledge about computational chemistry and is explained in the following step by step.</p>
<h1 id="usage">Usage</h1>
<p>As explained in the main article and the supporting information, RetroFit is written in Python (available at <a href="http://www.python.org">http://www.python.org</a>) and uses three input data sets, i.e. the structure of the MOF, the structure of the CL and a model interaction potential (MIP) of a test system modelling the interaction between CL and OMS. In the following, the prerequisites for the Python environment, the three input data sets, the operation of the script and the output files are explained.</p>
<h2 id="prerequisites">Prerequisites</h2>
<p>The RetroFit code was tested on Python 2.7 and Python 3. Besides the standard Python packages, the Atomic Simulation Environment Package (ASE) must be installed. Refer to <a href="https://wiki.fysik.dtu.dk/ase/">https://wiki.fysik.dtu.dk/ase/</a> for further information.</p>
<h2 id="input">Input:</h2>
<h3 id="cl-structures">CL structures</h3>
<p>The structure of CL molecules to be used in RetroFit have to be provided as .xyz-files. DFT-optimized structures (e.g. using the Gaussian09 program package) can be converted to .xyz-files using Open Babel (<a href="https://openbabel.org/">https://openbabel.org/</a>). RetroFit reads the first four lines of the .xyz-file containing atomic coordinates. These have to be the ones of the bridging functional groups. E.g. for cis-1,2-Dicyanoethene:</p>
<ul>
<li>first line: N-atom of nitrile group in position 1</li>
<li>second line: C-atom of nitrile group in position 1</li>
<li>third line: N-atom of nitrile group in position 2</li>
<li>fourth line: C-atom of nitrile group in position 2</li>
</ul>
<p>If the order of atomic coordinates in the .xyz file is different from the description above, please change it manually.<br>
In order to provide the structure of the CLs for the RetroFit routine, all .xyz-files must be located in one folder and the location of the folder have to be provided in the config file (see below).</p>
<h3 id="energy-values">Energy values</h3>
<p>Through single point DFT calculations, the MIP describing the interaction between the CL and the OMS of the MOF can be accessed. The results of the DFT calculations are provided in the form of a text file that is used as input by the RetroFit code. Precisely, a file which has four coloums (distance R [angstrom] vs angle theta [deg] vs angle delta [deg] vs energy [Ha]) must be provided. The location of this file must be specified in the config file (see below).</p>
<h3 id="mof-structure">MOF structure</h3>
<p>The relative positions of the OMS within a MOF has to be provided for RetroFit. Therefor the R_MM distance and the angle enclosed by the vector of one metal center and its OMS (for a paddlewheel-based MOF, this corresponds to the metal-metal vector of the paddlewheel) and the vector between the two metal centers have to be manually extracted from the crystal structure of the MOF. R_MM and the angle have to be provided in the config file (see below).</p>
<h3 id="retrofit.config-file">RetroFit.config file</h3>
<p>A file named <code>RetroFit.config</code> has to be created which serves as the input for the RetroFit program an contains the information regarding the three input data sets described above plus a number of parameters for the operation of the RetroFit code, such as step sizes for the interpolation. All parameters of the config filed are specified in the following:</p>
<ul>
<li><strong>delta_steps</strong>:    Number of steps that the delta-testlist consists of. The higher, the more precise.</li>
<li><strong>a_steps</strong>:        Number of that the alpha-testlist consists of. The higher, the more precise. This value is only used if an energy-surface file is written out. If not, this value can be chosen arbitrarily but must be in a the config file.</li>
<li><strong>rnn_steps</strong>:      Number of that the rnn-testlist consists of. The higher, the more precise. This value is only used if an energy-surface file is written out. If not, this value can be chosen arbitrarily but must be in a the config file.</li>
<li><strong>guest_dir</strong>:      Path to directory with all the xyz files of guest-molecules.</li>
<li><strong>data_file</strong>:      Path to file of energy output of MIP screening.</li>
<li><strong>MMM_angle</strong>:      Angle between 3 metal atoms of 2 paddle-wheels [deg]. 2 of them being the open metal-site atoms.</li>
<li><strong>RMM</strong>:            distance between metal atoms of 2 open metal-sites [Å].</li>
</ul>
<p>Example of a <code>RetroFit.config</code> file:</p>
<pre><code>#
#
delta_steps     100
a_steps         100
rnn_steps       100
guest_dir       xyz-files
data_file       project/energies.out
MMM_angle       120.00
RMM             7.9997
</code></pre>
<h2 id="starting-the-script">Starting the script</h2>
<p>Open the comand window of your computer and navigate to the directory containing <code>RetroFit.py</code>. Write in command line:</p>
<pre class=" language-sh"><code class="prism  language-sh">$ python RetroFit.py PATH_TO_FOLDER_OF_CONFIG_FILE
</code></pre>
<p>Additional options are:</p>
<ul>
<li>
<pre><code> -sMOL: writes out an output-file which contains all screened molecules including the optimised values.
</code></pre>
</li>
<li>
<pre><code>  -sE: writes out and output-file for an Energysurface (RNN vs alpha). This can take a while.
</code></pre>
</li>
</ul>
<p>Command line input with all options:</p>
<pre class=" language-sh"><code class="prism  language-sh">$ python RetroFit.py PATH_TO_FOLDER_OF_CONFIG_FILE -sMOL -sE
</code></pre>
<h2 id="output">Output</h2>
<p>If the -sMOL flag is set, an outputfile called <code>MOL_LIST.out</code> is created in the directory where the .config file is located. It contains the information about the optimal values for R, theta and delta for a certain molecule geometry (defined by alpha and R_NN) and the respective Energy according to the input MIP (! Energy is converted internally from Ha to kcal/mol !). This file has the following format:</p>
<pre><code>name | alpha [deg] | R_NN [Å] | Ropt [Å] | theta [deg] | deltaopt [deg] |Emin [kcal/mol]
</code></pre>
<p>Note, Ropt refers to the optimal R_MD distance calculated by RetroFit.</p>
<p>The first 5 lines contain information about the input data and the parameters used as specified in <code>RetroFit.config</code>. These lines start with <code>###</code>. If the optimal values could not be calculated, Ropt, theta, deltaopt and Emin are replaced by <code>--</code>.</p>
<p>If the -sE flag is set, an outputfile called <code>E_SURF.out</code> is created in the directory where the .config file is located. It contains a the calculated energy values for sets of alpha and R_NN values and the respective optimal R, theta and delta values with the following format:</p>
<pre><code>alpha [deg] | R_NN [deg] | E_min [kcal/mol] | Ropt [Å] | theta [deg] | deltaopt [deg]
</code></pre>
<p>The first 4 lines contain information about the input data and the parameters used as specified in <code>RetroFit.config</code>. These lines start with <code>###</code>. If the optimal values could not be calculated, Ropt, theta, deltaopt and E_min are replaced by <code>--</code>.</p>
</div>
</body>

</html>

{:/}
