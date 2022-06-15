!> @mainpage KORC documentation
!! @author Leopoldo Carbajal @note KORC and its documentation are in
!! constant evolution. This makes it possible that there are issues
!! that need to be solved.  If you find an issue please report it
!! immediatly through the "issues" section of the github repository.
!!
!! @section s1 Introduction The Kinetic Orbit Runaway electrons Code
!! (KORC) is a full-orbit particle tracer that evolves relativistic
!! electrons in both analytical and grid-based pre-computed electric
!! and magnetic fields. KORC includes the effects of: the acceleration
!! of the electrons due to the electric field, synchrotron radiation
!! energy losses, and collisions with the background plasma containing
!! high-Z impurities.
!!
!! For details about the equations of motion of the relativistic
!! electrons followed by KOqRC we refer the user to <em>Carbajal et
!! al. Phys. Plasmas <b>24</b>, 042512 (2017)</em> and <em>Carbajal
!! and del-Castillo-Negrete, Nuclear Fusion, submitted
!! (2018)</em>. Also, for details about the KORC's synchrotron
!! radiation synthetic diagnostic we refer the user to <em>Carbajal
!! and del-Castillo-Negrete, Plasma Phys. Controll. Fusion <b>59</b>,
!! 124001 (2017)</em>.
!!
!! KORC is a modular Fortran 95 code that uses a hybrid MPI + open MP
!! parallelization paradigm to exploit multi-core nodes systems, such
!! as Cori and Edison NERSC systems (<a
!! href="https://www.nersc.gov">www.nersc.gov</a>).
!!
!!
!! @section s2 Installation @warning KORC has been installed and
!! tested in systems with OS 10.6 and higher, Ubuntu 14.04, and SuSe
!! 12. Though we have paid special attention to portability when
!! developing KORC, there is no guarantee that KORC can be compiled,
!! ran, and is accurate in other systems not listed above. As a rule
!! of thumb, we recommend to perform any convenient benchmark tests
!! when compiling and running KORC in a new system, this to make sure
!! that the external libraries and compilers do not modify the
!! simulation results.
!!
!! @subsection s1s2 Getting started To compile and run KORC you will
!! need the following: <ol> <li>The <a
!! href="https://support.hdfgroup.org/downloads/index.html">HDF5</a>
!! library for I/O.</li> <li>The <a
!! href="https://w3.pppl.gov/ntcc/PSPLINE/">PSPLINES</a> library for
!! interpolating the pre-computed electric and magnetic fields as well
!! aqs the plasma profiles.</li> <li>The GNU compilers suite, or</li>
!! <li>The INTEL compilers suite.</li> <li>The corresponding open MP
!! library</li> <li>The MPI or open MPI library</li> <li><a
!! href="https://git-scm.com/"> Git, for obtaining the latest version
!! of KORC. </a> </ol>
!!
!! The HDF5 and PSPLINES libraries need to be installed BEFORE
!! compiling KORC, or if these are already present in the system, they
!! need to be loaded to your environment and their installation paths
!! need to be added to the Makefile accordingly. We refer the user to
!! the documentation of HDF5 and PSPLINES for the specifics about
!! their installation.
!!
!! If you are installing the HDF5 and PSPLINES libraries in your
!! system, we recommend to perform a local installation. Then enter
!! the absolute path of the folder containing the <em>"lib"</em> and
!! <em>"bin"</em> folders of each library to the PSPLINE_INSTALL and
!! HDF5_INSTALL variables of the Makefile.
!!
!!
!! @note PSPLINES is known for having precision issues sometimes when
!! compiled using the INTEL compilers suite in Linux systems.  Please,
!! double-check that the interpolations are giving the correct
!! numbers.
!!
!! @subsection s2s2 Cloning and compiling KORC To obtain the latest
!! version of KORC you will need to clone the Github repository to the
!! system where you want to run KORC.  This can be done as follows:
!! <ol> <li> Make sure you have all the needed external libraries and
!! compilers in place.  <li>Using the terminal of your system,
!! <em>cd</em> to the directory where you want to save KORC.</li>
!! <li>Type the following in the terminal <em>"git
!! https://github.com/ORNL-Fusion/KORC.git"</em>. This will copy all
!! the KORC files from the Github repository to your local system.
!! <li><em>cd</em> to the KORC-FO folder.  <li>In the terminal type:
!! <em>"./compile.sh GNU"</em> if you are using GNU compilers, or
!! <em>"./compile.sh INTEL"</em> if you are using INTEL compilers.
!! </ol>
!!
!! @section s3 Running KORC If you are running KORC in a system that
!! uses a resource manager such as TORQUE or SLURM, you will need to
!! request the number of nodes and cores per node that you will need,
!! this along with any other resources such as wall time, memory per
!! node, etc. Decide the number of MPI processes (nmpi) and the number
!! of open MP threads per MPI process (nomp) that you will use in your
!! simulation.  Then, on the terminal or in your batch file type the
!! following:
!!
!! <em> mpirun -np nmpi -x OMP_NUM_THREADS=nomp ./bin/KORC
!! "full_path_to_input_file.korc" "full_path_to_the_outputs_folder"
!! </em>
!!
!! The above example assumes that you are using the open MPI wrapper
!! mpirun to execute the binary <em>"KORC"</em> generated after
!! compiling KORC. This wrapper can be replaced by any other wrapper
!! that you might want to use. Notice that when using a different MPI
!! wrapper the way to pass the number of MPI processes (nmpi) and the
!! number of open MP threads per MPI (nomp) to KORC might vary.
!!
!! @section s4 Generating the KORC's input file.
!! KORC uses...
